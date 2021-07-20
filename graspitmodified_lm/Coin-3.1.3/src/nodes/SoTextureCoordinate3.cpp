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
  \class SoTextureCoordinate3 SoTextureCoordinate3.h Inventor/nodes/SoTextureCoordinate3.h
  \brief The SoTextureCoordinate3 class contains a set of coordinates for the mapping of 2D textures.
  \ingroup nodes

  When encountering nodes of this type during traversal, the
  coordinates it contains will be put on the state stack. Some shape
  nodes can then use these coordinates for explicit, detailed control
  of how 3D textures are mapped.

  (If 3D textures are used without any SoTextureCoordinate3 nodes in
  the scenegraph leading up to a shape node, the shape types have
  default fallbacks. So SoTextureCoordinate3 nodes are only necessary
  to use if you are not satisfied with the default mapping.)

  Note that an SoTextureCoordinate3 node will \e replace the
  coordinates already present in the state (if any).

  \COIN_CLASS_EXTENSION

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinate3 {
        point [  ]
    }
  \endcode

  \sa SoTextureCoordinate2
  \since Coin 2.0
  \since TGS Inventor 2.6
*/

// *************************************************************************

#include <Inventor/nodes/SoTextureCoordinate3.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/elements/SoGLVBOElement.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/C/glue/gl.h>

#include "nodes/SoSubNodeP.h"
#include "misc/SoVBO.h"

// *************************************************************************

/*!
  \var SoMFVec3f SoTextureCoordinate3::point

  The set of 3D texture coordinates. Default value of field is an
  empty set.

  Texture coordinates are usually specified in normalized coordinates,
  ie in the range [0, 1]. Coordinates outside the [0, 1] range can be
  used to repeat the texture across a surface.

  \sa SoTexture3::wrapR, SoTexure3::wrapS, SoTexture3::wrapT 
*/

// *************************************************************************

class SoTextureCoordinate3P {
 public:
  SoTextureCoordinate3P() : vbo(NULL) { }
  ~SoTextureCoordinate3P() { delete this->vbo; }
  SoVBO * vbo;
};

#define PRIVATE(obj) obj->pimpl

SO_NODE_SOURCE(SoTextureCoordinate3);

/*!
  Constructor.
*/
SoTextureCoordinate3::SoTextureCoordinate3(void)
{
  PRIVATE(this) = new SoTextureCoordinate3P;
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinate3);
  SO_NODE_ADD_FIELD(point, (NULL));
}

/*!
  Destructor.
*/
SoTextureCoordinate3::~SoTextureCoordinate3()
{
  delete PRIVATE(this);
}

// Documented in superclass.
void
SoTextureCoordinate3::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinate3, SO_FROM_INVENTOR_2_6|SO_FROM_COIN_2_0);

  SO_ENABLE(SoGLRenderAction, SoGLTextureCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoTextureCoordinateElement);
}

// Documented in superclass.
void
SoTextureCoordinate3::doAction(SoAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);
  
  if (unit == 0) {
    SoTextureCoordinateElement::set3(action->getState(), this,
                                     this->point.getNum(),
                                     this->point.getValues(0));
  }
  else {
    SoMultiTextureCoordinateElement::set3(action->getState(), this, unit,
                                          this->point.getNum(),
                                          this->point.getValues(0));
  }
}

// Documented in superclass.
void
SoTextureCoordinate3::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);

  if (unit == 0) {
    SoGLTextureCoordinateElement::setTexGen(action->getState(), this, NULL);
    SoTextureCoordinate3::doAction((SoAction *)action);
  }
  else {
    const cc_glglue * glue = cc_glglue_instance(SoGLCacheContextElement::get(state));
    int maxunits = cc_glglue_max_texture_units(glue);

    if (unit < maxunits) {
      SoGLMultiTextureCoordinateElement::setTexGen(action->getState(), this, unit, NULL);
      SoMultiTextureCoordinateElement::set3(action->getState(), this, unit,
                                            this->point.getNum(),
                                            this->point.getValues(0));
    }
  }

  SoBase::staticDataLock();
  const int num = this->point.getNum();
  SbBool setvbo = FALSE;
  if (SoGLVBOElement::shouldCreateVBO(state, num)) {
    setvbo = TRUE;
    SbBool dirty = FALSE;
    if (PRIVATE(this)->vbo == NULL) {
      PRIVATE(this)->vbo = new SoVBO(GL_ARRAY_BUFFER, GL_STATIC_DRAW); 
      dirty =  TRUE;
    }
    else if (PRIVATE(this)->vbo->getBufferDataId() != this->getNodeId()) {
      dirty = TRUE;
    }
    if (dirty) {
      PRIVATE(this)->vbo->setBufferData(this->point.getValues(0),
                                        num*sizeof(SbVec3f),
                                        this->getNodeId());
    }
  }
  else if (PRIVATE(this)->vbo && PRIVATE(this)->vbo->getBufferDataId()) {
    // clear buffers to deallocate VBO memory
    PRIVATE(this)->vbo->setBufferData(NULL, 0, 0);
  }
  SoBase::staticDataUnlock();
  SoGLVBOElement::setVertexVBO(state, setvbo ? PRIVATE(this)->vbo : NULL);
}

// Documented in superclass.
void
SoTextureCoordinate3::callback(SoCallbackAction * action)
{
  SoTextureCoordinate3::doAction((SoAction *)action);
}

// Documented in superclass.
void
SoTextureCoordinate3::pick(SoPickAction * action)
{
  SoTextureCoordinate3::doAction((SoAction *)action);
}

#undef PRIVATE
