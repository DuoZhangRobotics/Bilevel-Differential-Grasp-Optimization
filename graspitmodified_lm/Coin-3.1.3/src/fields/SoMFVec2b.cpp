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
  \class SoMFVec2b SoMFVec2b.h Inventor/fields/SoMFVec2b.h
  \brief The SoMFVec2b class is a container for SbVec2b vectors.
  \ingroup fields

  This field is used where nodes, engines or other field containers
  needs to store an array of vectors with two elements.

  This field supports application data sharing through a
  setValuesPointer() method. See SoMField documentation for
  information on how to use this function.

  \sa SbVec2b, SoSFVec2b
  \COIN_CLASS_EXTENSION
  \since Coin 2.5
*/

// *************************************************************************

#include <Inventor/fields/SoMFVec2b.h>

#include <cassert>

#include <Inventor/SoInput.h>
#include <Inventor/errors/SoDebugError.h>

#include "fields/SoSubFieldP.h"
#include "fields/shared.h"

// *************************************************************************

SO_MFIELD_SOURCE(SoMFVec2b, SbVec2b, SbVec2b);

SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2b, SbVec2b, SbVec2b);
SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2b, SbVec2b, int8_t);

// *************************************************************************

// Override from parent class.
void
SoMFVec2b::initClass(void)
{
  SO_MFIELD_INTERNAL_INIT_CLASS(SoMFVec2b);
}

// *************************************************************************

// No need to document readValue() and writeValue() here, as the
// necessary information is provided by the documentation of the
// parent classes.
#ifndef DOXYGEN_SKIP_THIS

SbBool
SoMFVec2b::read1Value(SoInput * in, int idx)
{
  assert(idx < this->maxNum);
  return
    in->readByte(this->values[idx][0]) &&
    in->readByte(this->values[idx][1]);
}

void
SoMFVec2b::write1Value(SoOutput * out, int idx) const
{
  sosfvec2b_write_value(out, (*this)[idx]);
}

#endif // DOXYGEN_SKIP_THIS

// *************************************************************************

/*!
  Set \a num vector array elements from \a xy, starting at index
  \a start.
*/
void
SoMFVec2b::setValues(int start, int numarg, const int8_t xy[][2])
{
  if (start+numarg > this->maxNum) this->allocValues(start+numarg);
  else if (start+numarg > this->num) this->num = start+numarg;

  for(int i=0; i < numarg; i++) this->values[start+i] = SbVec2b(xy[i]);
  this->valueChanged();
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2b::set1Value(int idx, int8_t x, int8_t y)
{
  this->set1Value(idx, SbVec2b(x, y));
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2b::set1Value(int idx, const int8_t xy[2])
{
  this->set1Value(idx, SbVec2b(xy));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2b::setValue(int8_t x, int8_t y)
{
  this->setValue(SbVec2b(x, y));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2b::setValue(const int8_t xy[2])
{
  if (xy == NULL) this->setNum(0);
  else this->setValue(SbVec2b(xy));
}

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoMFVec2b field;
  BOOST_CHECK_MESSAGE(field.getTypeId() != SoType::badType(),
                      "missing class initialization");
  BOOST_CHECK_EQUAL(field.getNum(), 0);
}

#endif // COIN_TEST_SUITE
