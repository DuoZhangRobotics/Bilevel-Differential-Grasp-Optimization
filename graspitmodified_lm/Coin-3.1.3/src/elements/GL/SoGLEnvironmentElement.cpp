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
  \class SoGLEnvironmentElement Inventor/elements/SoGLEnvironmentElement.h
  \brief The SoGLEnvironmentElement class is for setting GL fog etc.
  \ingroup elements
*/


#include <Inventor/elements/SoGLEnvironmentElement.h>
#include <Inventor/elements/SoViewVolumeElement.h>
#include <Inventor/SbColor4f.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/system/gl.h>

#include <cassert>

SO_ELEMENT_SOURCE(SoGLEnvironmentElement);

/*!
  This static method initializes static data for the
  SoGLEnvironmentElement class.
*/

void
SoGLEnvironmentElement::initClass()
{
  SO_ELEMENT_INIT_CLASS(SoGLEnvironmentElement, inherited);
}

/*!
  The destructor.
*/

SoGLEnvironmentElement::~SoGLEnvironmentElement()
{
}

void
SoGLEnvironmentElement::init(SoState * state)
{
  inherited::init(state);
}

void
SoGLEnvironmentElement::pop(SoState * state,
                           const SoElement * prevTopElement)
{
  this->capture(state);
  this->updategl(state);
}

void
SoGLEnvironmentElement::setElt(SoState * const stateptr,
                               const float ambientIntensityarg,
                               const SbColor & ambientColorarg,
                               const SbVec3f & attenuationarg,
                               const int32_t fogTypearg,
                               const SbColor & fogColorarg,
                               const float fogVisibilityarg,
                               const float fogStartarg)
{
  inherited::setElt(stateptr, ambientIntensityarg, ambientColorarg,
                    attenuationarg, fogTypearg, fogColorarg, fogVisibilityarg,
                    fogStartarg);
  this->updategl(stateptr);
}


//! FIXME: write doc.

void
SoGLEnvironmentElement::updategl(SoState * const state)
{
  float ambient[4];
  ambient[0] = ambientColor[0] * ambientIntensity;
  ambient[1] = ambientColor[1] * ambientIntensity;
  ambient[2] = ambientColor[2] * ambientIntensity;
  ambient[3] = 1.0f;

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  if (fogType == (int)NONE) {
    glDisable(GL_FOG);
    return;
  }

  float nearval = 1.0f;
  float farval = 10.0f;

  const SbViewVolume &vv = SoViewVolumeElement::get(state);
  nearval = vv.getNearDist();
  farval = nearval + vv.getDepth();

  float dist = fogVisibility;
  if (dist > 0.0f) farval = dist;

  switch (fogType) {
  case HAZE:
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogf(GL_FOG_START, this->fogStart);
    glFogf(GL_FOG_END, farval);
    break;
  case FOG:
    glFogi(GL_FOG_MODE, GL_EXP);

    // formula used for finding density:
    //
    // 1/256 = e ^ -(density * farval)
    //
    // ln (1/256) = ln (e ^ -(density * farval))
    //
    // -5.545 = - density * farval
    //
    // density = 5.545 / farval;

    glFogf(GL_FOG_DENSITY, 5.545f / farval);
    break;
  case SMOKE:
    glFogi(GL_FOG_MODE, GL_EXP2);
    // formula used for finding density:
    //
    // 1/256 = e ^ (-(density * farval)^2)
    //
    // ln (1/256) = ln (e ^ (-(density * farval)^2))
    //
    // -5.545 = - (density * farval)^2
    //
    // sqrt(5.545) = density * farval
    //
    // density = 2.35 / farval
    //

    glFogf(GL_FOG_DENSITY, 2.35f / farval);
    break;
  default:
    assert(0 && "unknown fog type");
    break;
  }

  SbColor4f color(fogColor, 1.0f);
  glFogfv(GL_FOG_COLOR, color.getValue());
  glEnable(GL_FOG);
}
