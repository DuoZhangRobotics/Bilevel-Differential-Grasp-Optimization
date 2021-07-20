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
  \class SoVRMLSpotLight SoVRMLSpotLight.h Inventor/VRMLnodes/SoVRMLSpotLight.h
  \brief The SoVRMLSpotLight class defines a spot light source.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  SpotLight {
    exposedField SFFloat ambientIntensity  0         # [0,1]
    exposedField SFVec3f attenuation       1 0 0     # [0,inf)
    exposedField SFFloat beamWidth         1.570796  # (0,pi/2]
    exposedField SFColor color             1 1 1     # [0,1]
    exposedField SFFloat cutOffAngle       0.785398  # (0,pi/2]
    exposedField SFVec3f direction         0 0 -1    # (-inf, inf)
    exposedField SFFloat intensity         1         # [0,1]
    exposedField SFVec3f location          0 0 0     # (-inf, inf)
    exposedField SFBool  on                TRUE
    exposedField SFFloat radius            100       # [0, inf)
  }
  \endverbatim

  The SpotLight node defines a light source that emits light from a specific
  point along a specific direction vector and constrained within a solid angle.
  Spotlights may illuminate geometry nodes that respond to light sources and
  intersect the solid angle defined by the SpotLight. Spotlight nodes are
  specified in the local coordinate system and are affected by ancestors'
  transformations.

  A detailed description of ambientIntensity, color, intensity, and
  VRML's lighting equations is provided in 4.6.6, Light sources
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.6>).
  More information on lighting concepts can be found in 4.14, Lighting
  model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>),
  including a detailed description of the VRML lighting equations.

  The \e location field specifies a translation offset of the centre
  point of the light source from the light's local coordinate system
  origin.  This point is the apex of the solid angle which bounds
  light emission from the given light source.

  The \e direction field
  specifies the direction vector of the light's central axis defined
  in the local coordinate system.

  The \e on field specifies whether the
  light source emits light. If on is TRUE, the light source is
  emitting light and may illuminate geometry in the scene. If on is
  FALSE, the light source does not emit light and does not illuminate
  any geometry.

  The \e radius field specifies the radial extent of the
  solid angle and the maximum distance from location that may be
  illuminated by the light source. The light source does not emit
  light outside this radius.  The radius shall be greater than or
  equal to zero.

  Both \e radius and \e location are affected by ancestors'
  transformations (scales affect radius and transformations affect
  location).

  The \e cutOffAngle field specifies the outer bound of the
  solid angle.  The light source does not emit light outside of this
  solid angle.

  The \e beamWidth field specifies an inner solid angle in
  which the light source emits light at uniform full intensity. The
  light source's emission intensity drops off from the inner solid
  angle (beamWidth) to the outer solid angle (cutOffAngle) as
  described in the following equations:

  \verbatim
  angle = the angle between the Spotlight's direction vector
          and the vector from the Spotlight location to the point
          to be illuminated

  if (angle >= cutOffAngle):
    multiplier = 0
  else if (angle <= beamWidth):
    multiplier = 1
  else:
    multiplier = (angle - cutOffAngle) / (beamWidth - cutOffAngle)

  intensity(angle) = SpotLight.intensity � multiplier
  \endverbatim

  If the beamWidth
  is greater than the cutOffAngle, beamWidth is defined to be equal to
  the cutOffAngle and the light source emits full intensity within the
  entire solid angle defined by cutOffAngle.  Both beamWidth and
  cutOffAngle shall be greater than 0.0 and less than or equal to
  pi/2.

  Figure 6.16 depicts the beamWidth, cutOffAngle, direction, location,
  and radius fields of the SpotLight node.

  <center>
  <img src="http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/Images/spotlight.gif">
  Figure 6.16 -- SpotLight node
  </center>

  SpotLight illumination falls off with distance as specified by three
  attenuation coefficients. The attenuation factor is

  \verbatim
  1/max(attenuation[0] + attenuation[1]�r + attenuation[2]�r^2 , 1),
  \endverbatim

  where r is the distance from the light to the surface being
  illuminated. The default is no attenuation. An attenuation value of
  (0, 0, 0) is identical to (1, 0, 0). Attenuation values shall be
  greater than or equal to zero. A detailed description of VRML's
  lighting equations is contained in 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>).

*/

/*!
  \var SoSFVec3f SoVRMLSpotLight::location
  The light position. Default value is (0, 0, 0).
*/

/*!
  \var SoSFVec3f SoVRMLSpotLight::direction
  The light direction. Default value is (0, 0, 1).
*/

/*!
  \var SoSFFloat SoVRMLSpotLight::beamWidth
  The spot beam width. Default value is PI/2.
*/

/*!
  \var SoSFFloat SoVRMLSpotLight::cutOffAngle
  The spot light cut off angle. Default value is PI/4.
*/

/*!
  \var SoSFFloat SoVRMLSpotLight::radius
  The light radius. Light is not emitted past it. Default value is 100.
*/

/*!
  \var SoSFVec3f SoVRMLSpotLight::attenuation
  The attenuiation vector. Default value is (1, 0, 0).
*/

#include <Inventor/VRMLnodes/SoVRMLSpotLight.h>

#include <math.h>

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

SO_NODE_SOURCE(SoVRMLSpotLight);

// Doc in parent
void
SoVRMLSpotLight::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLSpotLight, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLSpotLight::SoVRMLSpotLight(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLSpotLight);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(location, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(direction,(0.0f, 0.0f, -1.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(beamWidth, (float(M_PI)/2.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(cutOffAngle, (float(M_PI)/4.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(radius, (100.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(attenuation, (1.0f, 0.0f, 0.0f));
}

/*!
  Destructor.
*/
SoVRMLSpotLight::~SoVRMLSpotLight()
{
}

// Doc in parent
void
SoVRMLSpotLight::GLRender(SoGLRenderAction * action)
{
  if (!this->on.getValue()) return;

  SoState * state = action->getState();
  int idx = SoGLLightIdElement::increment(state);

  if (idx < 0) {
#if COIN_DEBUG
    SoDebugError::post("SoSpotLight::GLRender()",
                       "Max # lights exceeded :(\n");
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
  glLightfv(light, GL_SPOT_DIRECTION, this->direction.getValue().getValue());

  float cutoff = SbClamp(this->cutOffAngle.getValue(), 0.0f, float(M_PI)*0.5f) * 180.0f / float(M_PI);
  glLightf(light, GL_SPOT_EXPONENT, 0.0f);
  glLightf(light, GL_SPOT_CUTOFF, cutoff);

  // FIXME: consider radius and beamWidth
}

#endif // HAVE_VRML97
