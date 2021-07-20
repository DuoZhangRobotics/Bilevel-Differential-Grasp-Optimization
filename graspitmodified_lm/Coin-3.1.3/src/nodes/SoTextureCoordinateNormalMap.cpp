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
  \class SoTextureCoordinateNormalMap SoTextureCoordinateNormalMap.h Inventor/nodes/SoTextureCoordinateNormalMap.h
  \brief The SoTextureCoordinateNormalMap class generates texture coordinates by projecting onto a surrounding texture.
  \ingroup nodes

  The texture specifying the enviroment will be mapped around the 
  scenegraph below this node using a sphere. The texture will be mapped
  onto the scenegraph taking camera position into account. This will
  lead to an object reflecting its enviroment.

  Here is a scenegraph example showing how enviroment mapping can be
  applied to an object:

  <code>
  #Inventor V2.1 ascii

  Separator {

    Texture2 {
      filename "ocean.jpg" # the enviroment, in this case ocean
    }
    TextureCoordinateNormalMap {}

    Cube {} # the enviromentally mapped object
  }
  </code>

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinateNormalMap {
    }
  \endcode
*/

// *************************************************************************

// FIXME: Can this somehow relate to 3D textures? (kintel 20020203)

#include <Inventor/nodes/SoTextureCoordinateNormalMap.h>

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
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/system/gl.h>

#include "tidbitsp.h"
#include "nodes/SoSubNodeP.h"

// *************************************************************************

class SoTextureCoordinateNormalMapP {
public:
  static SbVec4f * dummy_texcoords;
  static void cleanup_func(void);
};

SbVec4f * SoTextureCoordinateNormalMapP::dummy_texcoords = NULL;

void
SoTextureCoordinateNormalMapP::cleanup_func(void)
{
  delete SoTextureCoordinateNormalMapP::dummy_texcoords;
}

// *************************************************************************

SO_NODE_SOURCE(SoTextureCoordinateNormalMap);

/*!
  Constructor.
*/
SoTextureCoordinateNormalMap::SoTextureCoordinateNormalMap()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinateNormalMap);
}

/*!
  Destructor.
*/
SoTextureCoordinateNormalMap::~SoTextureCoordinateNormalMap()
{
}

// doc in super
void
SoTextureCoordinateNormalMap::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinateNormalMap, SO_FROM_INVENTOR_1);

  SoTextureCoordinateNormalMapP::dummy_texcoords = new SbVec4f(0.0f, 0.0f, 0.0f, 1.0f);
  coin_atexit((coin_atexit_f *)SoTextureCoordinateNormalMapP::cleanup_func, CC_ATEXIT_NORMAL);
}

// generates texture coordinates for GLRender, callback and pick actions
const SbVec4f &
SoTextureCoordinateNormalMap::generate(void *userdata,
                                         const SbVec3f & /* p */,
                                         const SbVec3f &n)
{
  return *SoTextureCoordinateNormalMapP::dummy_texcoords;
}

// doc from parent
void
SoTextureCoordinateNormalMap::doAction(SoAction * action)
{
  SoTextureCoordinateElement::setFunction(action->getState(), this,
                                          generate,
                                          action->getState());
}

// doc from parent
void
SoTextureCoordinateNormalMap::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);
  if (unit == 0) {
    SoTextureCoordinateNormalMap::doAction((SoAction *)action);
    SoGLTextureCoordinateElement::setTexGen(action->getState(),
                                            this, handleTexgen,
                                            NULL,
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
                                                 NULL,
                                                 generate,
                                                 action->getState());

  }
}

// doc from parent
void
SoTextureCoordinateNormalMap::callback(SoCallbackAction * action)
{
  SoTextureCoordinateNormalMap::doAction((SoAction *)action);
}

// doc from parent
void
SoTextureCoordinateNormalMap::pick(SoPickAction * action)
{
  SoTextureCoordinateNormalMap::doAction((SoAction *)action);
}

void
SoTextureCoordinateNormalMap::handleTexgen(void * /* data */)
{
  glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP);
  glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP);  
  glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP);
}
