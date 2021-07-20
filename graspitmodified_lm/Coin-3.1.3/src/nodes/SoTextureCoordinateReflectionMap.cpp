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
  \class SoTextureCoordinateReflectionMap SoTextureCoordinateReflectionMap.h Inventor/nodes/SoTextureCoordinateReflectionMap.h
  \brief The SoTextureCoordinateReflectionMap class generates 3D reflection texture coordinates.
  \ingroup nodes

  This node is usually used along with a SoCubeMapTexture node...
  
  FIXME: more doc.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinateReflectionMap {
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoTextureCoordinateReflectionMap.h>

#include <stdlib.h>
#include <float.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/system/gl.h>

#include "nodes/SoSubNodeP.h"
#include "tidbitsp.h"

// *************************************************************************

class SoTextureCoordinateReflectionMapP {
public:
  static SbVec4f * dummy_texcoords;
  static void cleanup_func(void);
};

SbVec4f * SoTextureCoordinateReflectionMapP::dummy_texcoords = NULL;

void
SoTextureCoordinateReflectionMapP::cleanup_func(void)
{
  delete SoTextureCoordinateReflectionMapP::dummy_texcoords;
}

// *************************************************************************

SO_NODE_SOURCE(SoTextureCoordinateReflectionMap);

/*!
  Constructor.
*/
SoTextureCoordinateReflectionMap::SoTextureCoordinateReflectionMap()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinateReflectionMap);
}

/*!
  Destructor.
*/
SoTextureCoordinateReflectionMap::~SoTextureCoordinateReflectionMap()
{
}

// doc in super
void
SoTextureCoordinateReflectionMap::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinateReflectionMap, SO_FROM_INVENTOR_1);

  SoTextureCoordinateReflectionMapP::dummy_texcoords = new SbVec4f(0.0f, 0.0f, 0.0f, 1.0f);
  coin_atexit((coin_atexit_f *)SoTextureCoordinateReflectionMapP::cleanup_func, CC_ATEXIT_NORMAL);
}

// generates texture coordinates for GLRender, callback and pick actions
const SbVec4f &
SoTextureCoordinateReflectionMap::generate(void *userdata,
                                         const SbVec3f & /* p */,
                                         const SbVec3f &n)
{
  SoState *state = (SoState*)userdata;
  SbVec3f wn; // normal in world (eye) coordinates
  SoModelMatrixElement::get(state).multDirMatrix(n, wn);
  SbVec3f u = n;

  u.normalize();
  wn.normalize();

  // reflection vector
  SbVec3f r = u - SbVec3f(2.0f*wn[0]*wn[0]*u[0],
                          2.0f*wn[1]*wn[1]*u[1],
                          2.0f*wn[2]*wn[2]*u[2]);
  r.normalize();

  (*SoTextureCoordinateReflectionMapP::dummy_texcoords)[0] = r[0];
  (*SoTextureCoordinateReflectionMapP::dummy_texcoords)[1] = r[1];
  (*SoTextureCoordinateReflectionMapP::dummy_texcoords)[2] = r[2];
  (*SoTextureCoordinateReflectionMapP::dummy_texcoords)[3] = 1.0f;
  return *SoTextureCoordinateReflectionMapP::dummy_texcoords;
}

// doc from parent
void
SoTextureCoordinateReflectionMap::doAction(SoAction * action)
{
  SoTextureCoordinateElement::setFunction(action->getState(), this,
                                          generate,
                                          action->getState());
}

// doc from parent
void
SoTextureCoordinateReflectionMap::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);
  if (unit == 0) {
    SoTextureCoordinateReflectionMap::doAction((SoAction *)action);
    SoGLTextureCoordinateElement::setTexGen(action->getState(),
                                            this, handleTexgen,
                                            action,
                                            generate,
                                            action->getState());
  }
  else {
    SoMultiTextureCoordinateElement::setFunction(action->getState(), this,
                                                 unit,
                                                 generate,
                                                 action->getState());
    SoGLMultiTextureCoordinateElement::setTexGen(action->getState(),
                                                 this, unit, handleTexgen, 
                                                 action,
                                                 generate,
                                                 action->getState());

  }
}

// doc from parent
void
SoTextureCoordinateReflectionMap::callback(SoCallbackAction * action)
{
  SoTextureCoordinateReflectionMap::doAction((SoAction *)action);
}

// doc from parent
void
SoTextureCoordinateReflectionMap::pick(SoPickAction * action)
{
  SoTextureCoordinateReflectionMap::doAction((SoAction *)action);
}

void
SoTextureCoordinateReflectionMap::handleTexgen(void * data)
{
  glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);
  glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);  
  glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);
}
