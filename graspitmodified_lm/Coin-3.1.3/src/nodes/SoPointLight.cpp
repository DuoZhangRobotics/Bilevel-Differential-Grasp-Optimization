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
  \class SoPointLight SoPointLight.h Inventor/nodes/SoPointLight.h
  \brief The SoPointLight class is a node type for light sources.
  \ingroup nodes

  Pointlights emits light equally in all directions from a specified
  3D location.

  See also documentation of parent class for important information
  regarding light sources in general.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    PointLight {
        on TRUE
        intensity 1
        color 1 1 1
        location 0 0 1
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoPointLight.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/SbColor4f.h>
#include <Inventor/SbVec4f.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoEnvironmentElement.h>
#include <Inventor/elements/SoGLLightIdElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/elements/SoLightElement.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/system/gl.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoSFVec3f SoPointLight::location
  3D position of lightsource. Default value is <0, 0, 1>.
*/

// *************************************************************************

SO_NODE_SOURCE(SoPointLight);

// *************************************************************************

/*!
  Constructor.
*/
SoPointLight::SoPointLight(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoPointLight);

  SO_NODE_ADD_FIELD(location, (0.0f, 0.0f, 1.0f));
}

/*!
  Destructor.
*/
SoPointLight::~SoPointLight()
{
}

// Doc from superclass.
void
SoPointLight::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoPointLight, SO_FROM_INVENTOR_1|SoNode::VRML1);
}

// Doc from superclass.
void
SoPointLight::GLRender(SoGLRenderAction * action)
{
  if (!this->on.getValue()) return;

  int idx = SoGLLightIdElement::increment(action->getState());

  if (idx < 0) {
#if COIN_DEBUG
    SoDebugError::post("SoPointLight::GLRender()",
                       "Max # lights exceeded :(\n");
#endif // COIN_DEBUG
    return;
  }

  SoState * state = action->getState();

  SoLightElement::add(state, this, SoModelMatrixElement::get(state) *
                      SoViewingMatrixElement::get(state));
  
  GLenum light = (GLenum) (idx + GL_LIGHT0);

  SbVec3f attenuation = SoEnvironmentElement::getLightAttenuation(state);
  glLightf(light, GL_QUADRATIC_ATTENUATION, attenuation[0]);
  glLightf(light, GL_LINEAR_ATTENUATION, attenuation[1]);
  glLightf(light, GL_CONSTANT_ATTENUATION, attenuation[2]);

  SbColor4f lightcolor(0.0f, 0.0f, 0.0f, 1.0f);
  // disable ambient contribution from this light source
  glLightfv(light, GL_AMBIENT, lightcolor.getValue());

  lightcolor.setRGB(this->color.getValue());
  lightcolor *= this->intensity.getValue();

  glLightfv(light, GL_DIFFUSE, lightcolor.getValue());
  glLightfv(light, GL_SPECULAR, lightcolor.getValue());

  SbVec3f loc = this->location.getValue();

  // point (or spot) light when w = 1.0
  SbVec4f posvec(loc[0], loc[1], loc[2], 1.0f);
  glLightfv(light, GL_POSITION, posvec.getValue());

  // turning off spot light properties for ordinary lights
  glLightf(light, GL_SPOT_EXPONENT, 0.0);
  glLightf(light, GL_SPOT_CUTOFF, 180.0);
}
