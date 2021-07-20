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
  \class SoTextureCoordinateCube include/Inventor/nodes/SoTextureCoordinateCube.h
  \brief The SoTextureCoordinateCube class autogenerates cubemapped texture coordinated for shapes.
  \ingroup nodes

  The cube used for reference when mapping is the boundingbox for the shape.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinateCube {
    }
  \endcode

  \since Coin 2.3
*/
// FIXME: Add a better class description (20040123 handegar)

// *************************************************************************

#include <Inventor/nodes/SoTextureCoordinateCube.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG

#include <Inventor/C/glue/gl.h>
#include <Inventor/SbBox3f.h>
#include <Inventor/SoFullPath.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/caches/SoBoundingBoxCache.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoShape.h>
#include <Inventor/threads/SbStorage.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

typedef struct {
  SbVec3f origo;
  SbBox3f boundingbox;
  SoNode * currentshape;
  SoState * currentstate;
  SbVec4f texcoordreturn;
} so_texcoordcube_data;

static void
so_texcoordcube_construct_data(void * closure)
{
  so_texcoordcube_data * data = (so_texcoordcube_data *) closure;
  data->currentshape = NULL;
  data->currentstate = NULL;
  data->origo = SbVec3f(0,0,0);
}

static void
so_texcoordcube_destruct_data(void * closure)
{
}

// *************************************************************************

SO_NODE_SOURCE(SoTextureCoordinateCube);

// *************************************************************************

class SoTextureCoordinateCubeP {

public:
  SoTextureCoordinateCubeP(SoTextureCoordinateCube * texturenode)
    : master(texturenode) { }

  SbVec4f calculateTextureCoordinate(SbVec3f point, SbVec3f n);

  so_texcoordcube_data * so_texcoord_get_data() {
    so_texcoordcube_data * data = NULL;
    data = (so_texcoordcube_data *) this->so_texcoord_storage->get();
    assert(data && "Error retrieving thread data.");
    return data;
  }

  SbStorage * so_texcoord_storage;

private:
  SoTextureCoordinateCube * master;

};

static const SbVec4f & textureCoordinateCubeCallback(void * userdata, const SbVec3f & point, const SbVec3f & normal);

#define PRIVATE(p) (p->pimpl)
#define PUBLIC(p) (p->master)

// *************************************************************************

/*!
  Constructor.
*/
SoTextureCoordinateCube::SoTextureCoordinateCube(void)
{
  PRIVATE(this) = new SoTextureCoordinateCubeP(this);
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinateCube);

  pimpl->so_texcoord_storage = new SbStorage(sizeof(so_texcoordcube_data),
                                             so_texcoordcube_construct_data,
                                             so_texcoordcube_destruct_data);
}

/*!
  Destructor.
*/
SoTextureCoordinateCube::~SoTextureCoordinateCube()
{
  delete pimpl->so_texcoord_storage;
  delete PRIVATE(this);
}

// Documented in superclass.
void
SoTextureCoordinateCube::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinateCube, SO_FROM_COIN_2_3);

  SO_ENABLE(SoGLRenderAction, SoGLTextureCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoTextureCoordinateElement);
  SO_ENABLE(SoPickAction, SoTextureCoordinateElement);
}

const SbVec4f &
textureCoordinateCubeCallback(void * userdata,
                          const SbVec3f & point,
                          const SbVec3f & normal)
{

  SoTextureCoordinateCubeP * pimpl = (SoTextureCoordinateCubeP *) userdata;
  so_texcoordcube_data * data = pimpl->so_texcoord_get_data();

  SoState * state = data->currentstate;
  SoFullPath * path = (SoFullPath *) state->getAction()->getCurPath();
  SoNode * node = path->getTail();


  if (!node->isOfType(SoShape::getClassTypeId())) {
    // FIXME: A better way to handle this? (20040122 handegar)
    assert(FALSE && "TextureCoordinateCube callback called for a non-SoShape node.");
  }

  // Cast the node into a shape
  SoShape * shape = (SoShape *) node;

  if (shape != data->currentshape) {
    data->boundingbox.makeEmpty();
    const SoBoundingBoxCache * bboxcache = shape->getBoundingBoxCache();
    if (bboxcache && bboxcache->isValid(state)) {
      data->boundingbox = bboxcache->getProjectedBox();
      data->origo = data->boundingbox.getCenter();
    }
    else {
      shape->computeBBox(state->getAction(), data->boundingbox, data->origo);
      data->origo = data->boundingbox.getCenter();
    }
    data->currentshape = shape;

    // Expanding the bbox making it cube shaped
    float sx, sy, sz;
    data->boundingbox.getSize(sx, sy, sz);
    if (sy > sx) sx = sy;
    if (sz > sx) sx = sz;
    sx *= 0.5f;
    const SbVec3f c = data->origo;
    data->boundingbox.setBounds(c[0] - sx, c[1] - sx, c[2] - sx,
                                c[0] + sx, c[1] + sx, c[2] + sx);
  }

  const SbVec4f & ret = pimpl->calculateTextureCoordinate(point, normal);

  data->texcoordreturn = ret;
  return data->texcoordreturn;

}

SbVec4f
SoTextureCoordinateCubeP::calculateTextureCoordinate(SbVec3f point, SbVec3f n)
{

  so_texcoordcube_data * data = this->so_texcoord_get_data();

  double maxv = fabs(n[0]);
  int maxi = 0;

  if (fabs(n[1]) > maxv) { maxi = 1; maxv = fabs(n[1]); }
  if (fabs(n[2]) > maxv) { maxi = 2; }

  int i0 = (maxi + 1) % 3;
  int i1 = (maxi + 2) % 3;

  const SbVec3f bmax = data->boundingbox.getMax();
  const SbVec3f bmin = data->boundingbox.getMin();
  float d0 = bmax[i0] - bmin[i0];
  float d1 = bmax[i1] - bmin[i1];

  if (d0 == 0.0f) d0 = 1.0f;
  if (d1 == 0.0f) d1 = 1.0f;

  float s = (point[i0] - bmin[i0]) / d0;
  float t = (point[i1] - bmin[i1]) / d1;

  SbVec4f tc(s, t, 0.0f, 1.0f);
  switch (maxi) { // Flip textures according to projected cube-side
  case 0:
    tc[0] = 1.0f - t;
    tc[1] = s;
    break;
  case 1:
    tc[0] = t;
    tc[1] = 1.0f - s;
    break;
  }
  if (n[maxi] < 0.0f) {
    if (maxi == 1) {
      tc[1] = 1.0f - tc[1];
    }
    else {
      tc[0] = 1.0f - tc[0];
    }
  }

  return tc;
}


// Documented in superclass.
void
SoTextureCoordinateCube::doAction(SoAction * action)
{
  so_texcoordcube_data * data = PRIVATE(this)->so_texcoord_get_data();
  
  data->currentstate = action->getState();
  data->currentshape = NULL;

  int unit = SoTextureUnitElement::get(data->currentstate);
  if (unit == 0) {
    SoTextureCoordinateElement::setFunction(data->currentstate,
                                            this, textureCoordinateCubeCallback,
                                            PRIVATE(this));
  }
  else {
    SoMultiTextureCoordinateElement::setFunction(data->currentstate, this,
                                                 unit, textureCoordinateCubeCallback,
                                                 PRIVATE(this));
  }
}

// Documented in superclass.
void
SoTextureCoordinateCube::GLRender(SoGLRenderAction * action)
{
  so_texcoordcube_data * data = PRIVATE(this)->so_texcoord_get_data();

  data->currentstate = action->getState();
  data->currentshape = NULL;

  int unit = SoTextureUnitElement::get(data->currentstate);
  if (unit == 0) {
    SoTextureCoordinateElement::setFunction(data->currentstate,
                                            this, textureCoordinateCubeCallback,
                                            PRIVATE(this));
  }
  else {
    const cc_glglue * glue = cc_glglue_instance(SoGLCacheContextElement::get(action->getState()));
    int maxunits = cc_glglue_max_texture_units(glue);
    if (unit < maxunits) {        
      SoMultiTextureCoordinateElement::setFunction(data->currentstate, this,
                                                   unit, textureCoordinateCubeCallback,
                                                   PRIVATE(this));
    }
  }
}

// Documented in superclass.
void
SoTextureCoordinateCube::callback(SoCallbackAction * action)
{
  SoTextureCoordinateCube::doAction((SoAction *)action);
}

// Documented in superclass.
void
SoTextureCoordinateCube::pick(SoPickAction * action)
{
  SoTextureCoordinateCube::doAction((SoAction *)action);
}

#undef PRIVATE
#undef PUBLIC
