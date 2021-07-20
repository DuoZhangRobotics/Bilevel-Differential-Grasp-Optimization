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
  \class SoVRMLPointLight SoVRMLPointLight.h Inventor/VRMLnodes/SoVRMLPointLight.h
  \brief The SoVRMLPointLight class is used to represent a point light.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  PointLight {
    exposedField SFFloat ambientIntensity  0       # [0,1]
    exposedField SFVec3f attenuation       1 0 0   # [0, inf)
    exposedField SFColor color             1 1 1   # [0,1]
    exposedField SFFloat intensity         1       # [0,1]
    exposedField SFVec3f location          0 0 0   # (-inf, inf)
    exposedField SFBool  on                TRUE
    exposedField SFFloat radius            100     # [0, inf)
  }
  \endverbatim

  The PointLight node specifies a point light source at a 3D location
  in the local coordinate system. A point light source emits light
  equally in all directions; that is, it is
  omnidirectional. PointLight nodes are specified in the local
  coordinate system and are affected by ancestor transformations.
  Subclause 4.6.6, Light sources
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.6>),
  contains a detailed description of the ambientIntensity, color, and
  intensity fields.

  A PointLight node illuminates geometry within radius metres of its
  location. Both \e radius and \e location are affected by ancestors'
  transformations (scales affect \e radius and transformations affect
  \e location). 
  
  The \e radius field shall be greater than or equal to zero.
  PointLight node's illumination falls off with distance as specified
  by three attenuation coefficients. The attenuation factor is 1/
  max(attenuation[0] + attenuation[1]�r + attenuation[2]�r2, 1), where
  r is the distance from the light to the surface being illuminated.
  The default is no attenuation. An attenuation value of (0, 0, 0) is
  identical to (1, 0, 0). Attenuation values shall be greater than or
  equal to zero. A detailed description of VRML's lighting equations
  is contained in 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>).

*/

/*!
  \var SoSFVec3f SoVRMLPointLight::location
  Point light position.
*/

/*!
  \var SoSFFloat SoVRMLPointLight::radius
  Radius of point light.
*/

/*!
  \var SoSFVec3f SoVRMLPointLight::attenuation
  The point light attenuation.
*/

#include <Inventor/VRMLnodes/SoVRMLPointLight.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/SbColor4f.h>
#include <Inventor/SbVec4f.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLLightIdElement.h>
#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG
#include <Inventor/system/gl.h>

#include "nodes/SoSubNodeP.h"

SO_NODE_SOURCE(SoVRMLPointLight);

// doc in parent
void
SoVRMLPointLight::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLPointLight, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLPointLight::SoVRMLPointLight(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLPointLight);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(location, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(radius, (100.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(attenuation, (1.0f, 0.0f, 0.0f));
}

/*!
  Destructor.
*/
SoVRMLPointLight::~SoVRMLPointLight()
{
}

// Doc in parent
void
SoVRMLPointLight::GLRender(SoGLRenderAction * action)
{
  if (!this->on.getValue()) return;

  int idx = SoGLLightIdElement::increment(action->getState());

  if (idx < 0) {
#if COIN_DEBUG
    SoDebugError::post("SoVRMLPointLight::GLRender()",
                       "Max # lights exceeded :(");
#endif // COIN_DEBUG
    return;
  }

  GLenum light = (GLenum) (idx + GL_LIGHT0);

  SbVec3f att = this->attenuation.getValue();
  glLightf(light, GL_CONSTANT_ATTENUATION, att[0]);
  glLightf(light, GL_LINEAR_ATTENUATION, att[1]);
  glLightf(light, GL_QUADRATIC_ATTENUATION, att[2]);

  SbColor4f lightcolor(0.0f, 0.0f, 0.0f, 1.0f);
  lightcolor.setRGB(this->color.getValue());
  lightcolor *= this->ambientIntensity.getValue();
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

  // FIXME: consider radius
}

#endif // HAVE_VRML97
