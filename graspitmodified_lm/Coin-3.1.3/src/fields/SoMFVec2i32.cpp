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
  \class SoMFVec2i32 SoMFVec2i32.h Inventor/fields/SoMFVec2i32.h
  \brief The SoMFVec2i32 class is a container for SbVec2i32 vectors.
  \ingroup fields

  This field is used where nodes, engines or other field containers
  needs to store an array of vectors with two elements.

  This field supports application data sharing through a
  setValuesPointer() method. See SoMField documentation for
  information on how to use this function.

  \sa SbVec2i32, SoSFVec2i32
  \COIN_CLASS_EXTENSION
  \since Coin 2.5
*/

// *************************************************************************

#include <Inventor/fields/SoMFVec2i32.h>

#include <cassert>

#include <Inventor/SoInput.h>
#include <Inventor/errors/SoDebugError.h>

#include "fields/SoSubFieldP.h"
#include "fields/shared.h"

// *************************************************************************

SO_MFIELD_SOURCE(SoMFVec2i32, SbVec2i32, const SbVec2i32 &);

SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2i32, SbVec2i32, SbVec2i32);
SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2i32, SbVec2i32, int32_t);

// *************************************************************************

// Override from parent class.
void
SoMFVec2i32::initClass(void)
{
  SO_MFIELD_INTERNAL_INIT_CLASS(SoMFVec2i32);
}

// *************************************************************************

// No need to document readValue() and writeValue() here, as the
// necessary information is provided by the documentation of the
// parent classes.
#ifndef DOXYGEN_SKIP_THIS

SbBool
SoMFVec2i32::read1Value(SoInput * in, int idx)
{
  assert(idx < this->maxNum);
  return
    in->read(this->values[idx][0]) &&
    in->read(this->values[idx][1]);
}

void
SoMFVec2i32::write1Value(SoOutput * out, int idx) const
{
  sosfvec2i32_write_value(out, (*this)[idx]);
}

#endif // DOXYGEN_SKIP_THIS

// *************************************************************************

/*!
  Set \a num vector array elements from \a xy, starting at index
  \a start.
*/
void
SoMFVec2i32::setValues(int start, int numarg, const int32_t xy[][2])
{
  if (start+numarg > this->maxNum) this->allocValues(start+numarg);
  else if (start+numarg > this->num) this->num = start+numarg;

  for(int i=0; i < numarg; i++) this->values[start+i] = SbVec2i32(xy[i]);
  this->valueChanged();
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2i32::set1Value(int idx, int32_t x, int32_t y)
{
  this->set1Value(idx, SbVec2i32(x, y));
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2i32::set1Value(int idx, const int32_t xy[2])
{
  this->set1Value(idx, SbVec2i32(xy));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2i32::setValue(int32_t x, int32_t y)
{
  this->setValue(SbVec2i32(x, y));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2i32::setValue(const int32_t xy[2])
{
  if (xy == NULL) this->setNum(0);
  else this->setValue(SbVec2i32(xy));
}

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoMFVec2i32 field;
  BOOST_CHECK_MESSAGE(field.getTypeId() != SoType::badType(),
                      "missing class initialization");
  BOOST_CHECK_EQUAL(field.getNum(), 0);
}

#endif // COIN_TEST_SUITE
