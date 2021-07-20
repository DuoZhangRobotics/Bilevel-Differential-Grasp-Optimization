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
  \class SoTexture2Transform SoTexture2Transform.h Inventor/nodes/SoTexture2Transform.h
  \brief The SoTexture2Transform class is used to define 2D texture transformations.
  \ingroup nodes

  Textures applied to shapes in the scene can be transformed by
  "prefixing" in the state with instances of this node
  type. Translations, rotations and scaling in 2D can all be done.

  The default settings of this node's fields equals a "null
  transform", ie no transformation.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    Texture2Transform {
        translation 0 0
        rotation 0
        scaleFactor 1 1
        center 0 0
    }
  \endcode

  \sa SoTexture3Transform
*/

// *************************************************************************

#include <Inventor/nodes/SoTexture2Transform.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoGetMatrixAction.h>
#include <Inventor/elements/SoGLTextureMatrixElement.h>
#include <Inventor/elements/SoGLMultiTextureMatrixElement.h>
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoShapeStyleElement.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/C/glue/gl.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoSFVec2f SoTexture2Transform::translation

  Texture coordinate translation. Default value is [0, 0].
*/
/*!
  \var SoSFFloat SoTexture2Transform::rotation

  Texture coordinate rotation in radians (around z-axis, s is x-axis and t is
  y-axis).  Defaults to an identity rotation (ie zero rotation).
*/
/*!
  \var SoSFVec2f SoTexture2Transform::scaleFactor

  Texture coordinate scale factors. Default value is [1, 1].
*/
/*!
  \var SoSFVec2f SoTexture2Transform::center

  Center for scale and rotation. Default value is [0, 0].
*/

// *************************************************************************

SO_NODE_SOURCE(SoTexture2Transform);

/*!
  Constructor.
*/
SoTexture2Transform::SoTexture2Transform(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTexture2Transform);

  SO_NODE_ADD_FIELD(translation, (0.0f, 0.0f));
  SO_NODE_ADD_FIELD(rotation, (0.0f));
  SO_NODE_ADD_FIELD(scaleFactor, (1.0f, 1.0f));
  SO_NODE_ADD_FIELD(center, (0.0f, 0.0f));
}

/*!
  Destructor.
*/
SoTexture2Transform::~SoTexture2Transform()
{
}

// Documented in superclass.
void
SoTexture2Transform::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTexture2Transform, SO_FROM_INVENTOR_1|SoNode::VRML1);

  SO_ENABLE(SoGLRenderAction, SoGLTextureMatrixElement);
  SO_ENABLE(SoCallbackAction, SoTextureMatrixElement);
  SO_ENABLE(SoPickAction, SoTextureMatrixElement);

  SO_ENABLE(SoGLRenderAction, SoGLMultiTextureMatrixElement);
  SO_ENABLE(SoCallbackAction, SoMultiTextureMatrixElement);
  SO_ENABLE(SoPickAction, SoMultiTextureMatrixElement);
}


// Documented in superclass.
void
SoTexture2Transform::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();

  // don't modify the texture matrix while rendering the shadow map
  if (SoShapeStyleElement::get(state)->getFlags() & SoShapeStyleElement::SHADOWMAP) return;

  int unit = SoTextureUnitElement::get(state);

  if (unit == 0) {
    SoTexture2Transform::doAction(action);

  }
  else {
    const cc_glglue * glue =
      cc_glglue_instance(SoGLCacheContextElement::get(state));
    int maxunits = cc_glglue_max_texture_units(glue);

    if (unit < maxunits) {
      SbMatrix mat;
      this->makeMatrix(mat);
      SoMultiTextureMatrixElement::mult(state, this, unit, mat);
    }
    else {
      // we already warned in SoTextureUnit. I think it's best to just
      // ignore the texture here so that all textures for non-supported
      // units will be ignored. pederb, 2003-11-10
    }
  }
}

// Documented in superclass.
void
SoTexture2Transform::doAction(SoAction *action)
{
  SbMatrix mat;
  this->makeMatrix(mat);

  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);
  if (unit == 0) {
    SoTextureMatrixElement::mult(action->getState(), this,
                                 mat);
  }
  else {
    SoMultiTextureMatrixElement::mult(state, this, unit, mat);
  }
}

// Documented in superclass.
void
SoTexture2Transform::callback(SoCallbackAction *action)
{
  SoTexture2Transform::doAction(action);
}

// Documented in superclass.
void
SoTexture2Transform::getMatrix(SoGetMatrixAction * action)
{
  SbMatrix mat;
  this->makeMatrix(mat);
  action->getTextureMatrix().multLeft(mat);
  action->getTextureInverse().multRight(mat.inverse());
}

// Documented in superclass.
void
SoTexture2Transform::pick(SoPickAction * action)
{
  SoTexture2Transform::doAction(action);
}

//
// generate a matrix based on the fields
//
void
SoTexture2Transform::makeMatrix(SbMatrix & mat)
{
  SbMatrix tmp;
  SbVec2f c = this->center.isIgnored() ?
    SbVec2f(0.0f, 0.0f) :
    center.getValue();

  mat.makeIdentity();
  mat[3][0] = -c[0];
  mat[3][1] = -c[1];

  SbVec2f scale = this->scaleFactor.getValue();
  if (!this->scaleFactor.isIgnored() &&
      scale != SbVec2f(1.0f, 1.0f)) {
    tmp.makeIdentity();
    tmp[0][0] = scale[0];
    tmp[1][1] = scale[1];
    mat.multRight(tmp);
  }
  if (!this->rotation.isIgnored() && (this->rotation.getValue() != 0.0f)) {
    float cosa = (float)cos(this->rotation.getValue());
    float sina = (float)sin(this->rotation.getValue());
    tmp.makeIdentity();
    tmp[0][0] = cosa;
    tmp[1][0] = -sina;
    tmp[0][1] = sina;
    tmp[1][1] = cosa;
    mat.multRight(tmp);
  }
  if (!translation.isIgnored()) c+= this->translation.getValue();
  if (c != SbVec2f(0.0f, 0.0f)) {
    tmp.makeIdentity();
    tmp[3][0] = c[0];
    tmp[3][1] = c[1];
    mat.multRight(tmp);
  }
}
