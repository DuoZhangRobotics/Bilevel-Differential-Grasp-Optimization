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
  \class SoMFVec2d SoMFVec2d.h Inventor/fields/SoMFVec2d.h
  \brief The SoMFVec2d class is a container for SbVec2d vectors.
  \ingroup fields

  This field is used where nodes, engines or other field containers
  needs to store an array of vectors with two elements.

  This field supports application data sharing through a
  setValuesPointer() method. See SoMField documentation for
  information on how to use this function.

  \sa SbVec2d, SoSFVec2d
  \COIN_CLASS_EXTENSION
  \since Coin 2.5
*/

// *************************************************************************

#include <Inventor/fields/SoMFVec2d.h>

#include <cassert>

#include <Inventor/SoInput.h>
#include <Inventor/errors/SoDebugError.h>

#include "fields/SoSubFieldP.h"
#include "fields/shared.h"

// *************************************************************************

SO_MFIELD_SOURCE(SoMFVec2d, SbVec2d, const SbVec2d &);

SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2d, SbVec2d, SbVec2d);
SO_MFIELD_SETVALUESPOINTER_SOURCE(SoMFVec2d, SbVec2d, double);

// *************************************************************************

// Override from parent class.
void
SoMFVec2d::initClass(void)
{
  SO_MFIELD_INTERNAL_INIT_CLASS(SoMFVec2d);
}

// *************************************************************************

// No need to document readValue() and writeValue() here, as the
// necessary information is provided by the documentation of the
// parent classes.
#ifndef DOXYGEN_SKIP_THIS

SbBool
SoMFVec2d::read1Value(SoInput * in, int idx)
{
  assert(idx < this->maxNum);
  return
    in->read(this->values[idx][0]) &&
    in->read(this->values[idx][1]);
}

void
SoMFVec2d::write1Value(SoOutput * out, int idx) const
{
  sosfvec2d_write_value(out, (*this)[idx]);
}

#endif // DOXYGEN_SKIP_THIS

// *************************************************************************

/*!
  Set \a num vector array elements from \a xy, starting at index
  \a start.
*/
void
SoMFVec2d::setValues(int start, int numarg, const double xy[][2])
{
  if (start+numarg > this->maxNum) this->allocValues(start+numarg);
  else if (start+numarg > this->num) this->num = start+numarg;

  for(int i=0; i < numarg; i++) this->values[start+i] = SbVec2d(xy[i]);
  this->valueChanged();
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2d::set1Value(int idx, double x, double y)
{
  this->set1Value(idx, SbVec2d(x, y));
}

/*!
  Set the vector at \a idx.
*/
void
SoMFVec2d::set1Value(int idx, const double xy[2])
{
  this->set1Value(idx, SbVec2d(xy));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2d::setValue(double x, double y)
{
  this->setValue(SbVec2d(x, y));
}

/*!
  Set this field to contain a single vector with the given
  element values.
*/
void
SoMFVec2d::setValue(const double xy[2])
{
  if (xy == NULL) this->setNum(0);
  else this->setValue(SbVec2d(xy));
}

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoMFVec2d field;
  BOOST_CHECK_MESSAGE(field.getTypeId() != SoType::badType(),
                      "missing class initialization");
  BOOST_CHECK_EQUAL(field.getNum(), 0);
}

#endif // COIN_TEST_SUITE
