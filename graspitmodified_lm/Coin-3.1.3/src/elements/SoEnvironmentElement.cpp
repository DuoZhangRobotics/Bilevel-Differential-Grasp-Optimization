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
  \class SoEnvironmentElement Inventor/elements/SoEnvironmentElement.h
  \brief The SoEnvironmentElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoEnvironmentElement.h>


#include <cassert>

#include "SbBasicP.h"

/*!
  \fn SoEnvironmentElement::FogType

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::ambientIntensity

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::ambientColor

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::attenuation

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::fogType

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::fogColor

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::fogVisibility

  FIXME: write doc.
*/

/*!
  \fn SoEnvironmentElement::fogStart

  FIXME: write doc.
*/

SO_ELEMENT_SOURCE(SoEnvironmentElement);

/*!
  This static method initializes static data for the
  SoEnvironmentElement class.
*/

void
SoEnvironmentElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoEnvironmentElement, inherited);
}

/*!
  The destructor.
*/

SoEnvironmentElement::~SoEnvironmentElement(void)
{
}

//! FIXME: write doc.

void
SoEnvironmentElement::set(SoState * const state,
                          SoNode * const node,
                          const float ambientIntensity,
                          const SbColor & ambientColor,
                          const SbVec3f & attenuation,
                          const int32_t fogType,
                          const SbColor & fogColor,
                          const float fogVisibility,
                          const float fogStart)
{
  SoEnvironmentElement * element =
    coin_safe_cast<SoEnvironmentElement *>
    (
     SoReplacedElement::getElement(state, classStackIndex, node)
     );
  if (element) {
    element->setElt(state, ambientIntensity, ambientColor, attenuation,
                    fogType, fogColor, fogVisibility, fogStart);
  }
}

//! FIXME: write doc.

void
SoEnvironmentElement::get(SoState * const state,
                          float & ambientIntensity,
                          SbColor & ambientColor,
                          SbVec3f & attenuation,
                          int32_t & fogType,
                          SbColor & fogColor,
                          float & fogVisibility,
                          float & fogStart)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
    SoElement::getConstElement(state, classStackIndex)
    );

  ambientIntensity = element->ambientIntensity;
  ambientColor = element->ambientColor;
  attenuation = element->attenuation;
  fogType = element->fogType;
  fogColor = element->fogColor;
  fogVisibility = element->fogVisibility;
  fogStart = element->fogStart;
}

//! FIXME: write doc.

void
SoEnvironmentElement::getDefault(float & ambientIntensity,
                                 SbColor & ambientColor,
                                 SbVec3f & attenuation,
                                 int32_t & fogType,
                                 SbColor & fogColor,
                                 float & fogVisibility,
                                 float & fogStart)
{
  ambientIntensity = 0.2f;
  ambientColor = SbColor(1.0f, 1.0f, 1.0f);
  attenuation = SbVec3f(0.0f, 0.0f, 1.0f);
  fogType = NONE;
  fogColor = SbColor(1.0f, 1.0f, 1.0f);
  fogVisibility = 0.0f;
  fogStart = 0.0f;
}

//! FIXME: write doc.

float
SoEnvironmentElement::getAmbientIntensity(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->ambientIntensity;
}

//! FIXME: write doc.

float
SoEnvironmentElement::getFogVisibility(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->fogVisibility;
}

//! FIXME: write doc.

const SbVec3f &
SoEnvironmentElement::getLightAttenuation(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->attenuation;
}

//! FIXME: write doc.

const SbColor &
SoEnvironmentElement::getAmbientColor(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->ambientColor;
}

//! FIXME: write doc.

const SbColor &
SoEnvironmentElement::getFogColor(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->fogColor;
}

//! FIXME: write doc.

int32_t
SoEnvironmentElement::getFogType(SoState * const state)
{
  const SoEnvironmentElement * element = coin_assert_cast<const SoEnvironmentElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->fogType;
}

//! FIXME: write doc.

void
SoEnvironmentElement::print(FILE * file) const
{
  fprintf(file, "SoEnvironmentElement[%p]\n", this);
}

//! FIXME: write doc.

void
SoEnvironmentElement::init(SoState * state)
{
  inherited::init(state);
  this->getDefault(this->ambientIntensity, this->ambientColor, this->attenuation,
                   fogType, fogColor, fogVisibility, fogStart);
}

//! FIXME: doc
void
SoEnvironmentElement::setElt(SoState * const,
                             const float ambientIntensityarg,
                             const SbColor & ambientColorarg,
                             const SbVec3f & attenuationarg,
                             const int32_t fogTypearg,
                             const SbColor & fogColorarg,
                             const float fogVisibilityarg,
                             const float fogStartarg)
{
  this->ambientIntensity = ambientIntensityarg;
  this->ambientColor = ambientColorarg;
  this->attenuation = attenuationarg;
  this->fogType = fogTypearg;
  this->fogColor = fogColorarg;
  this->fogVisibility = fogVisibilityarg;
  this->fogStart = fogStartarg;
}
