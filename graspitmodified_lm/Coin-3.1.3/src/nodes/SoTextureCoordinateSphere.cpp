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
  \class SoTextureCoordinateSphere include/Inventor/nodes/SoTextureCoordinateSphere.h
  \brief The SoTextureCoordinateSphere class autogenerates spheremapped texture coordinated for shapes.
  \ingroup nodes

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinateSphere {
    }
  \endcode

  \since Coin 2.3
*/
// FIXME: Add a better class description (20040123 handegar)

// *************************************************************************

#include <Inventor/nodes/SoTextureCoordinateSphere.h>

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
} so_texcoordsphere_data;

static void
so_texcoordsphere_construct_data(void * closure)
{
  so_texcoordsphere_data * data = (so_texcoordsphere_data *) closure;
  data->currentshape = NULL;
  data->currentstate = NULL;
  data->origo = SbVec3f(0,0,0);
}

static void
so_texcoordsphere_destruct_data(void * closure)
{
}

SO_NODE_SOURCE(SoTextureCoordinateSphere);

class SoTextureCoordinateSphereP {

public:
  SoTextureCoordinateSphereP(SoTextureCoordinateSphere * texturenode)
    : master(texturenode) { }

  SbVec4f calculateTextureCoordinate(SbVec3f point, SbVec3f n);

  so_texcoordsphere_data * so_texcoord_get_data() {
    so_texcoordsphere_data * data = NULL;
    data = (so_texcoordsphere_data *) this->so_texcoord_storage->get();
    assert(data && "Error retrieving thread data.");
    return data;
  }

  SbStorage * so_texcoord_storage;

private:
  SoTextureCoordinateSphere * master;
};


static const SbVec4f & textureCoordinateSphereCallback(void * userdata, const SbVec3f & point, const SbVec3f & normal);

#define PRIVATE(p) (p->pimpl)
#define PUBLIC(p) (p->master)


/*!
  Constructor.
*/
SoTextureCoordinateSphere::SoTextureCoordinateSphere(void)
{

  PRIVATE(this) = new SoTextureCoordinateSphereP(this);

  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinateSphere);

  pimpl->so_texcoord_storage = new SbStorage(sizeof(so_texcoordsphere_data),
                                             so_texcoordsphere_construct_data,
                                             so_texcoordsphere_destruct_data);
}

/*!
  Destructor.
*/
SoTextureCoordinateSphere::~SoTextureCoordinateSphere()
{
  delete pimpl->so_texcoord_storage;
  delete PRIVATE(this);
}

// Documented in superclass.
void
SoTextureCoordinateSphere::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinateSphere, SO_FROM_COIN_2_3);

  SO_ENABLE(SoGLRenderAction, SoGLTextureCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoTextureCoordinateElement);
  SO_ENABLE(SoPickAction, SoTextureCoordinateElement);

}

const SbVec4f &
textureCoordinateSphereCallback(void * userdata,
                          const SbVec3f & point,
                          const SbVec3f & normal)
{

  SoTextureCoordinateSphereP * pimpl = (SoTextureCoordinateSphereP *) userdata;
  so_texcoordsphere_data * data = pimpl->so_texcoord_get_data();

  SoState * state = data->currentstate;
  SoFullPath * path = (SoFullPath *) state->getAction()->getCurPath();
  SoNode * node = path->getTail();


  if (!node->isOfType(SoShape::getClassTypeId())) {
    // FIXME: A better way to handle this? (20040122 handegar)
    assert(FALSE && "TextureCoordinateSphere callback called for a non-SoShape node.");
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
  }

  const SbVec4f & ret = pimpl->calculateTextureCoordinate(point, normal);

  data->texcoordreturn = ret;
  return data->texcoordreturn;

}

SbVec4f
SoTextureCoordinateSphereP::calculateTextureCoordinate(SbVec3f point, SbVec3f n)
{

  // FIXME: This way of mapping will always lead to artifacts in the
  // change between 360 and 0 degrees around the Y-axis. This is
  // unavoidable as the callback cannot predict when the last vertex
  // will be received, and therefore be able to patch up the
  // transition. (20040127 handegar)

  SbVec4f tc((float) (atan2(point[0], point[2]) * (1.0/(2.0*M_PI)) + 0.5),
             (float) (atan2(point[1], sqrt(point[0]*point[0] + point[2]*point[2])) * (1.0/M_PI) + 0.5),
             0.0f, 1.0f);

  return tc;

}


// Documented in superclass.
void
SoTextureCoordinateSphere::doAction(SoAction * action)
{
  so_texcoordsphere_data * data = PRIVATE(this)->so_texcoord_get_data();
  
  data->currentstate = action->getState();
  data->currentshape = NULL;
  
  int unit = SoTextureUnitElement::get(data->currentstate);
  if (unit == 0) {
    SoTextureCoordinateElement::setFunction(data->currentstate,
                                            this, textureCoordinateSphereCallback,
                                            PRIVATE(this));
  }
  else {
    SoMultiTextureCoordinateElement::setFunction(data->currentstate, this,
                                                 unit, textureCoordinateSphereCallback,
                                                 PRIVATE(this));
  }
}

// Documented in superclass.
void
SoTextureCoordinateSphere::GLRender(SoGLRenderAction * action)
{
  so_texcoordsphere_data * data = PRIVATE(this)->so_texcoord_get_data();

  data->currentstate = action->getState();
  data->currentshape = NULL;

  int unit = SoTextureUnitElement::get(data->currentstate);
  if (unit == 0) {
    SoTextureCoordinateElement::setFunction(data->currentstate,
                                            this, textureCoordinateSphereCallback,
                                            PRIVATE(this));
  }
  else {
    const cc_glglue * glue = cc_glglue_instance(SoGLCacheContextElement::get(action->getState()));
    int maxunits = cc_glglue_max_texture_units(glue);
    if (unit < maxunits) {        
      SoMultiTextureCoordinateElement::setFunction(data->currentstate, this,
                                                   unit, textureCoordinateSphereCallback,
                                                   PRIVATE(this));
    }
  }

}

// Documented in superclass.
void
SoTextureCoordinateSphere::callback(SoCallbackAction * action)
{
  SoTextureCoordinateSphere::doAction((SoAction *)action);
}

// Documented in superclass.
void
SoTextureCoordinateSphere::pick(SoPickAction * action)
{
  SoTextureCoordinateSphere::doAction((SoAction *)action);
}

#undef PRIVATE
#undef PUBLIC
