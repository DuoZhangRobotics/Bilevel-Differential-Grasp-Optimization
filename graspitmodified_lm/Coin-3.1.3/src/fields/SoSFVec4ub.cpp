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
  \class SoSFVec4ub SoSFVec4ub.h Inventor/fields/SoSFVec4ub.h
  \brief The SoSFVec4ub class is a container for an SbVec4ub vector.
  \ingroup fields

  This field is used where nodes, engines or other field containers
  needs to store a single vector with four elements.

  \sa SbVec4ub, SoMFVec4ub
  \COIN_CLASS_EXTENSION
  \since Coin 2.5
*/

// *************************************************************************

#include <Inventor/fields/SoSFVec4ub.h>

#include <Inventor/SoInput.h>
#include <Inventor/SoOutput.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/errors/SoReadError.h>

#include "fields/SoSubFieldP.h"
#include "fields/shared.h"

// *************************************************************************

SO_SFIELD_SOURCE(SoSFVec4ub, SbVec4ub, SbVec4ub);

// *************************************************************************

// Override from parent class.
void
SoSFVec4ub::initClass(void)
{
  SO_SFIELD_INTERNAL_INIT_CLASS(SoSFVec4ub);
}

// *************************************************************************

// No need to document readValue() and writeValue() here, as the
// necessary information is provided by the documentation of the
// parent classes.
#ifndef DOXYGEN_SKIP_THIS

SbBool
SoSFVec4ub::readValue(SoInput * in)
{
  return
    in->readByte(this->value[0]) &&
    in->readByte(this->value[1]) &&
    in->readByte(this->value[2]) &&
    in->readByte(this->value[3]);
}

void
SoSFVec4ub::writeValue(SoOutput * out) const
{
  sosfvec4ub_write_value(out, this->getValue());
}

#endif // DOXYGEN_SKIP_THIS


/*!
  Set value of vector.
*/
void
SoSFVec4ub::setValue(uint8_t x, uint8_t y, uint8_t z, uint8_t w)
{
  this->setValue(SbVec4ub(x, y, z, w));
}

/*!
  Set value of vector.
*/
void
SoSFVec4ub::setValue(const uint8_t xyzw[4])
{
  this->setValue(SbVec4ub(xyzw));
}

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  BOOST_CHECK_MESSAGE(SoSFVec4ub::getClassTypeId() != SoType::badType(),
                      "SoSFVec4ub class not initialized");
  SoSFVec4ub field;
  BOOST_CHECK_MESSAGE(field.getTypeId() != SoType::badType(),
                      "SoSFVec4ub object wrongly initialized");
  // no default value initialization to test
  field.setValue(1, 2, 3, 4);
  BOOST_CHECK_EQUAL(field.getValue(), SbVec4ub(1, 2, 3, 4));
}

BOOST_AUTO_TEST_CASE(textinput)
{
  TestSuite::ResetReadErrorCount();
  SbBool ok;
  SoSFVec4ub field;
  ok = field.set("1 2 3 4");
  BOOST_CHECK_MESSAGE(ok == TRUE, "SoSFVec4ub read error");
  BOOST_CHECK_EQUAL(field.getValue(), SbVec4ub(1, 2, 3, 4));
  BOOST_CHECK_EQUAL(TestSuite::GetReadErrorCount(), 0);
  TestSuite::ResetReadErrorCount();
}

#endif // COIN_TEST_SUITE
