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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#ifdef HAVE_VRML97

/*!
  \class SoVRMLDirectionalLight SoVRMLDirectionalLight.h Inventor/VRMLnodes/SoVRMLDirectionalLight.h
  \brief The SoVRMLDirectionalLight class is a node type for specifying directional light sources
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  DirectionalLight {
    exposedField SFFloat ambientIntensity  0        # [0,1]
    exposedField SFColor color             1 1 1    # [0,1]
    exposedField SFVec3f direction         0 0 -1   # (-inf,inf)
    exposedField SFFloat intensity         1        # [0,1]
    exposedField SFBool  on                TRUE
  }
  \endverbatim

  The DirectionalLight node defines a directional light source that
  illuminates along rays parallel to a given 3-dimensional vector. A
  description of the ambientIntensity, color, intensity, and on fields
  is in 4.6.6, Light sources
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.6>).

  The direction field specifies the direction vector of the
  illumination emanating from the light source in the local coordinate
  system. Light is emitted along parallel rays from an infinite
  distance away. A directional light source illuminates only the
  objects in its enclosing parent group.  The light may illuminate
  everything within this coordinate system, including all children and
  descendants of its parent group. The accumulated transformations of
  the parent nodes affect the light.  DirectionalLight nodes do not
  attenuate with distance. A precise description of VRML's lighting
  equations is contained in 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.14>).
  
*/

/*!
  \var SoSFVec3f SoVRMLDirectionalLight::direction
  The light direction.
*/

#include <Inventor/VRMLnodes/SoVRMLDirectionalLight.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/SbColor4f.h>
#include <Inventor/SbVec4f.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoEnvironmentElement.h>
#include <Inventor/elements/SoGLLightIdElement.h>
#include <Inventor/system/gl.h>
#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

#include "nodes/SoSubNodeP.h"

SO_NODE_SOURCE(SoVRMLDirectionalLight);

// Doc in parent
void
SoVRMLDirectionalLight::initClass(void) // static
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLDirectionalLight, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLDirectionalLight::SoVRMLDirectionalLight(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLDirectionalLight);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(direction, (0.0f, 0.0f, -1.0f));
}

/*!
  Destructor.
*/
SoVRMLDirectionalLight::~SoVRMLDirectionalLight()
{
}

// Doc in parent
void
SoVRMLDirectionalLight::GLRender(SoGLRenderAction * action)
{
  if (!this->on.getValue()) return;

  SoState * state = action->getState();
  int idx = SoGLLightIdElement::increment(state);

  if (idx < 0) {
#if COIN_DEBUG
    SoDebugError::postWarning("SoDirectionalLight::GLRender",
                              "Max # of OpenGL lights exceeded :(");
#endif // COIN_DEBUG
    return;
  }

  GLenum light = (GLenum) (idx + GL_LIGHT0);

  SbColor4f lightcolor(0.0f, 0.0f, 0.0f, 1.0f);
  lightcolor.setRGB(this->color.getValue());
  lightcolor *= this->ambientIntensity.getValue();
  glLightfv(light, GL_AMBIENT, lightcolor.getValue());

  lightcolor.setRGB(this->color.getValue());
  lightcolor *= this->intensity.getValue();

  glLightfv(light, GL_DIFFUSE, lightcolor.getValue());
  glLightfv(light, GL_SPECULAR, lightcolor.getValue());

  // GL directional light is specified towards light source
  SbVec3f dir = - this->direction.getValue();
  dir.normalize();

  // directional when w = 0.0
  SbVec4f dirvec(dir[0], dir[1], dir[2], 0.0f);
  glLightfv(light, GL_POSITION, dirvec.getValue());

  glLightf(light, GL_SPOT_EXPONENT, 0.0);
  glLightf(light, GL_SPOT_CUTOFF, 180.0);
  glLightf(light, GL_CONSTANT_ATTENUATION, 1);
  glLightf(light, GL_LINEAR_ATTENUATION, 0);
  glLightf(light, GL_QUADRATIC_ATTENUATION, 0);

}

#endif // HAVE_VRML97
