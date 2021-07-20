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
  \class SoSFColor SoSFColor.h Inventor/fields/SoSFColor.h
  \brief The SoSFColor class is a container for an SbColor value.
  \ingroup fields

  This field is used where nodes, engines or other field containers
  needs to store a single color value (i.e. Red + Green + Blue).

  Fields of this type stores their value to file as a "R G B" triple
  component, where each color component value is between 0.0 and 1.0.

  \sa SoMFColor
*/

// *************************************************************************

#include <Inventor/fields/SoSFColor.h>

#include <Inventor/SoInput.h>
#include <Inventor/SoOutput.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/errors/SoReadError.h>

#include "fields/SoSubFieldP.h"
#include "fields/shared.h"

// *************************************************************************

SO_SFIELD_SOURCE(SoSFColor, SbColor, const SbColor &);

// *************************************************************************

// Override from parent.
void
SoSFColor::initClass(void)
{
  SO_SFIELD_INTERNAL_INIT_CLASS(SoSFColor);
}

// No need to document readValue() and writeValue() here, as the
// necessary information is provided by the documentation of the
// parent classes.
#ifndef DOXYGEN_SKIP_THIS

SbBool
SoSFColor::readValue(SoInput * in)
{
  return
    in->read(this->value[0]) &&
    in->read(this->value[1]) &&
    in->read(this->value[2]);
}

void
SoSFColor::writeValue(SoOutput * out) const
{
  sosfvec3f_write_value(out, this->getValue());
}

#endif // DOXYGEN_SKIP_THIS


/*!
  Set color value from a vector with three elements. The three elements
  will be interpreted as red, green and blue, respectively.
*/
void
SoSFColor::setValue(const SbVec3f & vec)
{
  this->value.setValue(vec[0], vec[1], vec[2]);
  this->valueChanged();
}

/*!
  Set color value from \a red, \a green and \a blue. Value range
  for each component is between 0.0 and 1.0.
*/
void
SoSFColor::setValue(const float red, const float green, const float blue)
{
  this->value.setValue(red, green, blue);
  this->valueChanged();
}

/*!
  Set color value from a floating point number array with three
  elements. The three elements will be interpreted as red, green and
  blue, respectively.
*/
void
SoSFColor::setValue(const float rgb[3])
{
  this->value.setValue(rgb);
  this->valueChanged();
}

/*!
  Set color value from \a h, \a s and \a v, where \a is "hue", \a s
  is "saturation" and \a v is "value".
*/
void
SoSFColor::setHSVValue(const float h, const float s, const float v)
{
  this->value.setHSVValue(h, s, v);
  this->valueChanged();
}

/*!
  Set color value from a floating point number array with three
  elements. The three elements will be interpreted as hue,
  saturation and value, respectively.
*/
void
SoSFColor::setHSVValue(const float hsv[3])
{
  this->value.setHSVValue(hsv);
  this->valueChanged();
}

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoSFColor field;
  BOOST_CHECK_MESSAGE(SoSFColor::getClassTypeId() != SoType::badType(),
                      "SoSFColor class not initialized");
  BOOST_CHECK_MESSAGE(field.getTypeId() != SoType::badType(),
                      "missing class initialization");
}

#endif // COIN_TEST_SUITE
