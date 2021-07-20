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
  \class SoContextHandler SoContextHandler.h Inventor/misc/SoContextHandler.h
  \brief The SoContextHandler class is for now to be treated as an internal class.
  \ingroup general

  \since Coin 2.0
*/

// FIXME: should be documented and be part of the Doxygen API doc,
// since it's a public class (and possibly useful
// externally). 20030225 mortene.
//
// UPDATE: there also a function in SoGLCacheContextElement which
// looks like it might do the same thing: scheduleDeleteCallback()..?
// Investigate.
//
// 20040723 mortene.

//
// UPDATE: No, scheduleDeleteCallback() will not do the same thing.
// It's only used for deleting SoGLDisplayLists.
//
// 20050209 pederb
//

// *************************************************************************

#include <Inventor/misc/SoContextHandler.h>

#include <stdlib.h>

#include <Inventor/errors/SoDebugError.h>
#include <Inventor/lists/SbList.h>

#include "misc/SbHash.h"
#include "threads/threadsutilp.h"
#include "tidbitsp.h"
#include "glue/glp.h"

// *************************************************************************

class socontexthandler_cbitem {
public:
  socontexthandler_cbitem(void) : func(NULL), closure(NULL), idx(0) { }
  
  int operator==(const socontexthandler_cbitem & theother) {
    return 
      this->func == theother.func &&
      this->closure == theother.closure;
  }
  
  operator unsigned long(void) const {
    unsigned long key = 0;
    // create an xor key
    const unsigned char * ptr = (const unsigned char *) this;
    
    // a bit hackish. Stop xor'ing at idx
    const unsigned char * stop = (const unsigned char*) &this->idx;
    
    const ptrdiff_t size = stop - ptr;
    
    for (int i = 0; i < size; i++) {
      int shift = (i%4) * 8;
      key ^= (ptr[i]<<shift);
    }
    return key;
  }
  
  SoContextHandler::ContextDestructionCB * func;
  void * closure;
  
  // this must be last!!
  int idx;
};

struct socontexthandler_sbhashcb :
  public SbHash <uint32_t, socontexthandler_cbitem>::ApplyFunctor<SbList <socontexthandler_cbitem> *>
{
  void operator()(socontexthandler_cbitem & key,
            uint32_t & obj, SbList <socontexthandler_cbitem> * list)
  {
    list->append(key);
  }
};

// "extern C" wrapper is needed with the OSF1/cxx compiler (probably a
// bug in the compiler, but it doesn't seem to hurt to do this
// anyway).
extern "C" { static int 
socontexthandler_qsortcb(const void * p0, const void * p1)
{
  socontexthandler_cbitem * i0 = (socontexthandler_cbitem *) p0; 
  socontexthandler_cbitem * i1 = (socontexthandler_cbitem *) p1; 

  return int(i0->idx) - int(i1->idx);
}
}

// *************************************************************************

static SbHash <uint32_t, socontexthandler_cbitem> * socontexthandler_hashlist;
static uint32_t socontexthandler_idx = 0;
static void * socontexthandler_mutex;

// *************************************************************************

static void
socontexthandler_cleanup(void)
{
#if COIN_DEBUG
  const int len = socontexthandler_hashlist ?
    socontexthandler_hashlist->getNumElements() : 0;
  if (len > 0) {
    // Can't use SoDebugError here, as SoError et al might have been
    // "cleaned up" already.
    (void)printf("Coin debug: socontexthandler_cleanup(): %d context-bound "
                 "resources not free'd before exit.\n", len);
  }
#endif // COIN_DEBUG  
  delete socontexthandler_hashlist;
  socontexthandler_hashlist = NULL;
  socontexthandler_idx = 0;
  CC_MUTEX_DESTRUCT(socontexthandler_mutex);
}

// *************************************************************************

/*!
  This method \e must be called by client code which destructs a
  context, to guarantee that there are no memory leaks upon context
  destruction.

  This will take care of correctly freeing context-bound resources,
  like OpenGL texture objects and display lists.

  Before calling this function, the context \e must be made current.

  Note that if you are using one of the standard GUI-binding libraries
  from Kongsberg Oil & Gas Technologies, this is taken care of
  automatically for contexts for canvases set up by SoQt, SoWin, etc.
*/
void
SoContextHandler::destructingContext(uint32_t contextid)
{
  CC_MUTEX_CONSTRUCT(socontexthandler_mutex);
  CC_MUTEX_LOCK(socontexthandler_mutex);
  if (socontexthandler_hashlist == NULL) { 
    CC_MUTEX_UNLOCK(socontexthandler_mutex);
    return; 
  }

  SbList <socontexthandler_cbitem> listcopy;
  socontexthandler_sbhashcb functor;
  socontexthandler_hashlist->apply(functor, &listcopy);
  CC_MUTEX_UNLOCK(socontexthandler_mutex);

  qsort((void*) listcopy.getArrayPtr(), 
        socontexthandler_hashlist->getNumElements(),
        sizeof(socontexthandler_cbitem), 
        socontexthandler_qsortcb);

  // process callbacks FILO-style so that callbacks registered first
  // are called last. HACK WARNING: SoGLCacheContextElement will add a
  // callback in initClass(). It's quite important that this callback
  // is called after all other callbacks (since the other callbacks
  // might schedule destruction of GL resources through the methods in
  // SoGLCacheContextElement). This criteria is met as it is now,
  // since it's the only callback added while initializing Coin
  // (SoDB::init()). 

  // FIXME: We should probably add a new method in
  // SoGLCacheContextElement which this class can call after all the
  // regular callbacks though. pederb, 2004-10-27
  for (int i = listcopy.getLength()-1; i >= 0; i--) {
    const socontexthandler_cbitem & item = listcopy[i];
    item.func(contextid, item.closure);
  }

  // tell glglue that this context is dead
  coin_glglue_destruct(contextid);
}

// *************************************************************************

/*!
  Add a callback which will be called every time a GL context is
  destructed. The callback should delete all GL resources tied to that
  context.
  
  All nodes/classes that allocate GL resources should set up a callback
  like this. Add the callback in the constructor of the node/class,
  and remove it in the destructor.

  \sa removeContextDestructionCallback()
*/
void
SoContextHandler::addContextDestructionCallback(ContextDestructionCB * func,
                                                void * closure)
{
  CC_MUTEX_CONSTRUCT(socontexthandler_mutex);
  CC_MUTEX_LOCK(socontexthandler_mutex);
  if (socontexthandler_hashlist == NULL) {
    socontexthandler_hashlist = new SbHash <uint32_t, socontexthandler_cbitem> (64);
    // make this callback trigger after the SoGLCacheContext cleanup function
    // by setting priority to -1
    coin_atexit((coin_atexit_f *)socontexthandler_cleanup, CC_ATEXIT_NORMAL_LOWPRIORITY);
  }
  socontexthandler_cbitem item;
  item.func = func;
  item.closure = closure;
  item.idx = socontexthandler_idx++;
  socontexthandler_hashlist->put(item, item.idx);
  CC_MUTEX_UNLOCK(socontexthandler_mutex);
}

/*!
  Remove a context destruction callback.
  
  \sa addContextDestructionCallback()
*/
void
SoContextHandler::removeContextDestructionCallback(ContextDestructionCB * func, void * closure)
{
  assert(socontexthandler_hashlist);

  socontexthandler_cbitem item;
  item.func = func;
  item.closure = closure;
  
  CC_MUTEX_LOCK(socontexthandler_mutex);
  int didremove = socontexthandler_hashlist->remove(item);
  assert(didremove);
  CC_MUTEX_UNLOCK(socontexthandler_mutex);
}

// *************************************************************************
