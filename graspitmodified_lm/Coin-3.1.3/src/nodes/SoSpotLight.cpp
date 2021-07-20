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
  \class SoSpotLight SoSpotLight.h Inventor/nodes/SoSpotLight.h
  \brief The SoSpotLight class is a node type for light sources with a cone shaped lightvolume.
  \ingroup nodes

  Spotlights are light sources with a position and a direction. They
  can be thought of as a pointlight with a lampshade.

  See also documentation of parent class for important information
  regarding light sources in general.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    SpotLight {
        on TRUE
        intensity 1
        color 1 1 1
        location 0 0 1
        direction 0 0 -1
        dropOffRate 0
        cutOffAngle 0.78539819
    }
  \endcode

  \sa SoSpotLight
*/

// *************************************************************************

#include <Inventor/nodes/SoSpotLight.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
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
  \var SoSFVec3f SoSpotLight::location

  3D position of light source. Default position is <0, 0, 1>.
*/
/*!
  \var SoSFVec3f SoSpotLight::direction

  Direction vector, where the light is pointing. Default is to point
  along the negative z-axis.
*/
/*!
  \var SoSFFloat SoSpotLight::dropOffRate

  The rate of intensity drop-off from the ray along the direction
  vector. Value must be between 0.0 (equal intensity for the whole
  cone of light), to 1.0 (a narrow intensity ray).

  Default value is 0.0.
*/
/*!
  \var SoSFFloat SoSpotLight::cutOffAngle

  The angle in radians from the direction vector where there will be
  no light outside (i.e. the angle of the "lampshade"). Default value
  is PI/4.0 (i.e. 45�). The value of this field will be clamped to
  [0.0, PI/2] before it is used.
*/


// *************************************************************************

SO_NODE_SOURCE(SoSpotLight);

// *************************************************************************

/*!
  Constructor.
*/
SoSpotLight::SoSpotLight(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoSpotLight);

  SO_NODE_ADD_FIELD(location, (0.0f, 0.0f, 1.0f));
  SO_NODE_ADD_FIELD(direction, (0.0f, 0.0f, -1.0f));
  SO_NODE_ADD_FIELD(dropOffRate, (0.0f));
  SO_NODE_ADD_FIELD(cutOffAngle, (float(M_PI)/4.0f));
}

/*!
  Destructor.
*/
SoSpotLight::~SoSpotLight()
{
}

// Doc in superclass.
void
SoSpotLight::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoSpotLight, SO_FROM_INVENTOR_1|SoNode::VRML1);
}

// Doc in superclass.
void
SoSpotLight::GLRender(SoGLRenderAction * action)
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
  glLightfv(light, GL_SPOT_DIRECTION, this->direction.getValue().getValue());

  float cutoff = this->cutOffAngle.getValue() * 180.0f / float(M_PI);
  float dropoff = SbClamp(this->dropOffRate.getValue(), 0.0f, 1.0f) * 128.0f;
  
#ifdef COIN_EXTRA_DEBUG // output a warning if the cutoff is invalid
                        // since we now clamp it (someone might have
                        // been setting it to 180.0, which would make
                        // this a PointLight)
  if (cutoff < 0.0f || cutoff > 90.0f) {
    SoDebugError::postWarning("SoSpotLight::GLRender",
                              "invalid cutOffAngle for SpotLight: %f, clamping to [0.0f, 90.0f]", cutoff);
  }
#endif // COIN_EXTRA_DEBUG

  cutoff = SbClamp(cutoff, 0.0f, 90.0f);

  glLightf(light, GL_SPOT_EXPONENT, dropoff);
  glLightf(light, GL_SPOT_CUTOFF, cutoff);
}
