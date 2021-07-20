/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) by Kongsberg Oil & Gas Technologies.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg Oil & Gas Technologies
 *  about acquiring a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg Oil & Gas Technologies, Bygdoy Alle 5, 0257 Oslo, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

/*!
  \class SoVBO
  \brief The SoVBO class is used to handle OpenGL vertex buffer objects.

  It wraps the buffer handling, taking care of multi-context handling
  and allocation/deallocation of buffers. FIXME: more doc.

*/

#include "misc/SoVBO.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <Inventor/misc/SoContextHandler.h>
#include <Inventor/misc/SoGLDriverDatabase.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoShapeStyleElement.h>
#include <Inventor/C/tidbits.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/errors/SoDebugError.h>

#include "misc/SoVertexArrayIndexer.h"
#include "threads/threadsutilp.h"
#include "glue/glp.h"
#include "tidbitsp.h"
#include "share/gl/CoinGLPerformance.h"

static int vbo_vertex_count_min_limit = -1;
static int vbo_vertex_count_max_limit = -1;
static int vbo_render_as_vertex_arrays = -1;
static int vbo_enabled = -1;
static int vbo_debug = -1;

// VBO rendering seems to be faster than other rendering, even for
// large VBOs. Just set the default limit very high
static const int DEFAULT_MAX_LIMIT = 100000000;
static const int DEFAULT_MIN_LIMIT = 20;

static SbHash <SbBool, uint32_t> * vbo_isfast_hash;

/*!
  Constructor
*/
SoVBO::SoVBO(const GLenum target, const GLenum usage)
  : target(target),
    usage(usage),
    data(NULL),
    datasize(0),
    dataid(0),
    didalloc(FALSE),
    vbohash(5)
{
  SoContextHandler::addContextDestructionCallback(context_destruction_cb, this);
}


//
// Callback from SoGLCacheContextElement
//
void
SoVBO::vbo_delete(void * closure, uint32_t contextid)
{
  const cc_glglue * glue = cc_glglue_instance((int) contextid);
  GLuint id = (GLuint) ((uintptr_t) closure);
  cc_glglue_glDeleteBuffers(glue, 1, &id);
}

//
// Callback from SbHash
//
struct vbo_schedule :
  public SbHash <GLuint, uint32_t>::ApplyFunctor<void *>
{
  void operator()(uint32_t & key,
                GLuint & value,
                void * closure
                )
  {
    void * ptr = (void*) ((uintptr_t) value);
    SoGLCacheContextElement::scheduleDeleteCallback(key, SoVBO::vbo_delete, ptr);
  }
};

/*!
  Destructor
*/
SoVBO::~SoVBO()
{
  SoContextHandler::removeContextDestructionCallback(context_destruction_cb, this);
  // schedule delete for all allocated GL resources
  vbo_schedule functor;
  this->vbohash.apply(functor, static_cast<void *>(NULL));
  if (this->didalloc) {
    char * ptr = (char*) this->data;
    delete[] ptr;
  }
}

// atexit cleanup function
static void vbo_atexit_cleanup(void)
{
  delete vbo_isfast_hash;
  vbo_isfast_hash = NULL;
  vbo_vertex_count_min_limit = -1;
  vbo_vertex_count_max_limit = -1;
  vbo_render_as_vertex_arrays = -1;
  vbo_enabled = -1;
}

void
SoVBO::init(void)
{
  coin_glglue_add_instance_created_callback(context_created, NULL);

  vbo_isfast_hash = new SbHash <SbBool, uint32_t> (3);
  coin_atexit(vbo_atexit_cleanup, CC_ATEXIT_NORMAL);

  // use COIN_VBO_MAX_LIMIT to set the largest VBO we create
  if (vbo_vertex_count_max_limit < 0) {
    const char * env = coin_getenv("COIN_VBO_MAX_LIMIT");
    if (env) {
      vbo_vertex_count_max_limit = atoi(env);
    }
    else {
      vbo_vertex_count_max_limit = DEFAULT_MAX_LIMIT;
    }
  }

  // use COIN_VBO_MIN_LIMIT to set the smallest VBO we create
  if (vbo_vertex_count_min_limit < 0) {
    const char * env = coin_getenv("COIN_VBO_MIN_LIMIT");
    if (env) {
      vbo_vertex_count_min_limit = atoi(env);
    }
    else {
      vbo_vertex_count_min_limit = DEFAULT_MIN_LIMIT;
    }
  }

  // use COIN_VERTEX_ARRAYS to globally disable vertex array rendering
  if (vbo_render_as_vertex_arrays < 0) {
    const char * env = coin_getenv("COIN_VERTEX_ARRAYS");
    if (env) {
      vbo_render_as_vertex_arrays = atoi(env);
    }
    else {
      vbo_render_as_vertex_arrays = 1;
    }
  }

  // use COIN_VBO to globally disable VBOs when doing vertex array rendering
  if (vbo_enabled < 0) {
    const char * env = coin_getenv("COIN_VBO");
    if (env) {
      vbo_enabled = atoi(env);
    }
    else {
      vbo_enabled = 1;
    }
  }
  if (vbo_debug < 0) {
    const char * env = coin_getenv("COIN_DEBUG_VBO");
    if (env) {
      vbo_debug = atoi(env);
    }
    else {
      vbo_debug = 0;
    }
  }
}

/*!
  Used to allocate buffer data. The user is responsible for filling in
  the correct type of data in the buffer before the buffer is used.

  \sa setBufferData()
*/
void *
SoVBO::allocBufferData(intptr_t size, uint32_t dataid)
{
  // schedule delete for all allocated GL resources
  vbo_schedule functor;
  this->vbohash.apply(functor, static_cast<void *>(NULL));
  // clear hash table
  this->vbohash.clear();

  if (this->didalloc && this->datasize == size) {
    return (void*)this->data;
  }
  if (this->didalloc) {
    char * ptr = (char*) this->data;
    delete[] ptr;
  }

  char * ptr = new char[size];
  this->didalloc = TRUE;
  this->data = (const GLvoid*) ptr;
  this->datasize = size;
  this->dataid = dataid;
  return (void*) this->data;
}

/*!
  Sets the buffer data. \a dataid is a unique id used to identify
  the buffer data. In Coin it's possible to use the node id
  (SoNode::getNodeId()) to test if a buffer is valid for a node.
*/
void
SoVBO::setBufferData(const GLvoid * data, intptr_t size, uint32_t dataid)
{
  // schedule delete for all allocated GL resources
  vbo_schedule functor;
  this->vbohash.apply(functor, static_cast<void *>(NULL));
  // clear hash table
  this->vbohash.clear();

  // clean up old buffer (if any)
  if (this->didalloc) {
    char * ptr = (char*) this->data;
    delete[] ptr;
  }

  this->data = data;
  this->datasize = size;
  this->dataid = dataid;
  this->didalloc = FALSE;
}

/*!
  Returns the buffer data id.

  \sa setBufferData()
*/
uint32_t
SoVBO::getBufferDataId(void) const
{
  return this->dataid;
}

/*!
  Returns the data pointer and size.
*/
void
SoVBO::getBufferData(const GLvoid *& data, intptr_t & size)
{
  data = this->data;
  size = this->datasize;
}


/*!
  Binds the buffer for the context \a contextid.
*/
void
SoVBO::bindBuffer(uint32_t contextid)
{
  if ((this->data == NULL) ||
      (this->datasize == 0)) {
    assert(0 && "no data in buffer");
    return;
  }

  const cc_glglue * glue = cc_glglue_instance((int) contextid);

  GLuint buffer;
  if (!this->vbohash.get(contextid, buffer)) {
    // need to create a new buffer for this context
    cc_glglue_glGenBuffers(glue, 1, &buffer);
    cc_glglue_glBindBuffer(glue, this->target, buffer);
    cc_glglue_glBufferData(glue, this->target,
                           this->datasize,
                           this->data,
                           this->usage);
    this->vbohash.put(contextid, buffer);
  }
  else {
    // buffer already exists, bind it
    cc_glglue_glBindBuffer(glue, this->target, buffer);
  }

#if COIN_DEBUG
  if (vbo_debug) {
    if (this->target == GL_ELEMENT_ARRAY_BUFFER) {
      SoDebugError::postInfo("SoVBO::bindBuffer",
                             "Rendering using VBO. Index array size: %d",
                             this->datasize / sizeof(int32_t));
    }
    else {
      SoDebugError::postInfo("SoVBO::bindBuffer",
                             "Setting up buffer for rendering. Datasize: %d",
                             this->datasize);
    }
  }
#endif // COIN_DEBUG
}



//
// Callback from SoContextHandler
//
void
SoVBO::context_destruction_cb(uint32_t context, void * userdata)
{
  GLuint buffer;
  SoVBO * thisp = (SoVBO*) userdata;

  if (thisp->vbohash.get(context, buffer)) {
    const cc_glglue * glue = cc_glglue_instance((int) context);
    cc_glglue_glDeleteBuffers(glue, 1, &buffer);
    thisp->vbohash.remove(context);
  }
}


/*!
  Sets the global limits on the number of vertex data in a node before
  vertex buffer objects are considered to be used for rendering.
*/
void
SoVBO::setVertexCountLimits(const int minlimit, const int maxlimit)
{
  vbo_vertex_count_min_limit = minlimit;
  vbo_vertex_count_max_limit = maxlimit;
}

/*!
  Returns the vertex VBO minimum limit.

  \sa setVertexCountLimits()
 */
int
SoVBO::getVertexCountMinLimit(void)
{
  return vbo_vertex_count_min_limit;
}

/*!
  Returns the vertex VBO maximum limit.

  \sa setVertexCountLimits()
 */
int
SoVBO::getVertexCountMaxLimit(void)
{
  return vbo_vertex_count_max_limit;
}

SbBool
SoVBO::shouldCreateVBO(SoState * state, const uint32_t contextid, const int numdata)
{
  if (!vbo_enabled || !vbo_render_as_vertex_arrays) return FALSE;
  int minv = SoVBO::getVertexCountMinLimit();
  int maxv = SoVBO::getVertexCountMaxLimit();
  return
    (numdata >= minv) &&
    (numdata <= maxv) &&
    SoVBO::isVBOFast(contextid) &&
    !(SoShapeStyleElement::get(state)->getFlags() & SoShapeStyleElement::SHADOWMAP);

}

SbBool
SoVBO::shouldRenderAsVertexArrays(SoState * state,
                                  const uint32_t contextid,
                                  const int numdata)
{
  // FIXME: consider also using results from the performance tests

  // don't render as vertex arrays if there are very few elements to
  // be rendered. The VA setup overhead would make it slower than just
  // doing plain immediate mode rendering.
  return (numdata >= vbo_vertex_count_min_limit) && vbo_render_as_vertex_arrays;
}

SbBool
SoVBO::isVBOFast(const uint32_t contextid)
{
  SbBool result = TRUE;
  assert(vbo_isfast_hash != NULL);
  (void) vbo_isfast_hash->get(contextid, result);
  return result;
}

//
// callback from glglue (when a new glglue instance is created)
//
void
SoVBO::context_created(const uint32_t contextid, void * closure)
{
  SoVBO::testGLPerformance(contextid);
}

//
// test OpenGL performance for a context.
//
void
SoVBO::testGLPerformance(const uint32_t contextid)
{
  SbBool isfast;
  // did we alreay test this for this context?
  assert(vbo_isfast_hash != NULL);
  if (vbo_isfast_hash->get(contextid, isfast)) return;

  // Run time test disabled. Our old test seemed to be buggy, and
  // VBO should be fast on all platforms supporting it now. It was
  // really just one obscure laptop driver that caused problems for
  // us. This should be handled in the driver database anyway
  const cc_glglue * glue = cc_glglue_instance(contextid);
  if (SoGLDriverDatabase::isSupported(glue, SO_GL_VERTEX_BUFFER_OBJECT)) {
    vbo_isfast_hash->put(contextid, TRUE);
  }
  else {
    vbo_isfast_hash->put(contextid, FALSE);
  }
}
