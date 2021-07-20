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
  \class SoTextureCoordinateElement Inventor/elements/SoTextureCoordinateElement.h
  \brief The SoTextureCoordinateElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.

  \COIN_CLASS_EXTENSION

  \since Coin 2.2
*/

#include "coindefs.h"
#include "SbBasicP.h"

#include <Inventor/elements/SoMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoGLVBOElement.h>
#include <Inventor/nodes/SoNode.h>
#include <cassert>

#define MAX_UNITS 16
#define PRIVATE(obj) obj->pimpl

class SoMultiTextureCoordinateElementP {
public:
  SoMultiTextureCoordinateElement::UnitData unitdata[MAX_UNITS];
};

SO_ELEMENT_CUSTOM_CONSTRUCTOR_SOURCE(SoMultiTextureCoordinateElement);

/*!
  This static method initializes static data for the
  SoMultiTextureCoordinateElement class.
*/

void
SoMultiTextureCoordinateElement::initClass()
{
  SO_ELEMENT_INIT_CLASS(SoMultiTextureCoordinateElement, inherited);
}


/*!
  The constructor.
*/
SoMultiTextureCoordinateElement::SoMultiTextureCoordinateElement(void)
{
  PRIVATE(this) = new SoMultiTextureCoordinateElementP;

  this->setTypeId(SoMultiTextureCoordinateElement::classTypeId);
  this->setStackIndex(SoMultiTextureCoordinateElement::classStackIndex);
}

/*!
  The destructor.
*/

SoMultiTextureCoordinateElement::~SoMultiTextureCoordinateElement()
{
  delete PRIVATE(this);
}

//! FIXME: write doc.

void
SoMultiTextureCoordinateElement::setDefault(SoState * const state,
                                            SoNode * const COIN_UNUSED_ARG(node),
                                            const int unit)
{
  if (state->isElementEnabled(SoGLVBOElement::getClassStackIndex())) {
    SoGLVBOElement::setTexCoordVBO(state, unit, NULL);
  }
  SoMultiTextureCoordinateElement * element =
    coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );

  assert(unit >= 0 && unit < MAX_UNITS);
  UnitData & ud = PRIVATE(element)->unitdata[unit];
  ud.nodeid = 0;
  ud.whatKind = SoTextureCoordinateElement::DEFAULT;
  ud.numCoords = 0;
}

//! FIXME: write doc.

void
SoMultiTextureCoordinateElement::setFunction(SoState * const state,
                                             SoNode * const node,
                                             const int unit,
                                             SoTextureCoordinateFunctionCB * const func,
                                             void * const userdata)
{
  if (state->isElementEnabled(SoGLVBOElement::getClassStackIndex())) {
    SoGLVBOElement::setTexCoordVBO(state, unit, NULL);
  }

  SoMultiTextureCoordinateElement * element =
    coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );

  assert(unit >= 0 && unit < MAX_UNITS);
  UnitData & ud = PRIVATE(element)->unitdata[unit];

  ud.nodeid = node->getNodeId();
  ud.funcCB = func;
  ud.funcCBData = userdata;
  ud.whatKind = SoTextureCoordinateElement::FUNCTION;
  ud.coords2 = NULL;
  ud.coords3 = NULL;
  ud.coords4 = NULL;
  ud.numCoords = 0;
}

//! FIXME: write doc.

void
SoMultiTextureCoordinateElement::set2(SoState * const state,
                                      SoNode * const node,
                                      const int unit,
                                      const int32_t numCoords,
                                      const SbVec2f * const coords)
{
  if (state->isElementEnabled(SoGLVBOElement::getClassStackIndex())) {
    SoGLVBOElement::setTexCoordVBO(state, unit, NULL);
  }
  SoMultiTextureCoordinateElement * element = coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );

  assert(unit >= 0 && unit < MAX_UNITS);
  UnitData & ud = PRIVATE(element)->unitdata[unit];

  ud.nodeid = node->getNodeId();
  ud.coordsDimension = 2;
  ud.numCoords = numCoords;
  ud.coords2 = coords;
  ud.coords3 = NULL;
  ud.coords4 = NULL;
  ud.whatKind = SoTextureCoordinateElement::EXPLICIT;
}

/*!
  FIXME: write doc.
*/
void
SoMultiTextureCoordinateElement::set3(SoState * const state,
                                      SoNode * const node,
                                      const int unit,
                                      const int32_t numCoords,
                                      const SbVec3f * const coords)
{
  if (state->isElementEnabled(SoGLVBOElement::getClassStackIndex())) {
    SoGLVBOElement::setTexCoordVBO(state, unit, NULL);
  }
  SoMultiTextureCoordinateElement * element =
    coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );

  assert(unit >= 0 && unit < MAX_UNITS);
  UnitData & ud = PRIVATE(element)->unitdata[unit];

  ud.nodeid = node->getNodeId();
  ud.coordsDimension = 3;
  ud.numCoords = numCoords;
  ud.coords2 = NULL;
  ud.coords3 = coords;
  ud.coords4 = NULL;
  ud.whatKind = SoTextureCoordinateElement::EXPLICIT;
}

//! FIXME: write doc.

void
SoMultiTextureCoordinateElement::set4(SoState * const state,
                                      SoNode * const node,
                                      const int unit,
                                      const int32_t numCoords,
                                      const SbVec4f * const coords)
{
  if (state->isElementEnabled(SoGLVBOElement::getClassStackIndex())) {
    SoGLVBOElement::setTexCoordVBO(state, unit, NULL);
  }
  SoMultiTextureCoordinateElement * element =
    coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );

  assert(unit >= 0 && unit < MAX_UNITS);
  UnitData & ud = PRIVATE(element)->unitdata[unit];

  ud.nodeid = node->getNodeId();
  ud.coordsDimension = 4;
  ud.numCoords = numCoords;
  ud.coords2 = NULL;
  ud.coords3 = NULL;
  ud.coords4 = coords;
  ud.whatKind = SoTextureCoordinateElement::EXPLICIT;
}

//! FIXME: write doc.

const SoMultiTextureCoordinateElement *
SoMultiTextureCoordinateElement::getInstance(SoState * const state)
{
  return coin_safe_cast<const SoMultiTextureCoordinateElement *>
    (
     getConstElement(state, classStackIndex)
     );
}

/*!
  This method returns texture coordinate for the given point and normal.
  The coordinate is returned as a 4D vector where the r and q coordinates
  may be set to 0 and 1 respecively depending on what texture coordinate
  dimension we're using.

  This method should only be used if the CoordType is FUNCTION.
*/

const SbVec4f &
SoMultiTextureCoordinateElement::get(const int unit,
                                     const SbVec3f & point,
                                     const SbVec3f & normal) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];

  assert((ud.whatKind == SoTextureCoordinateElement::FUNCTION ||
          ud.whatKind == SoTextureCoordinateElement::TEXGEN) && ud.funcCB);
  return (*(ud.funcCB))(ud.funcCBData, point, normal);
}

//! FIXME: write doc.

const SbVec2f &
SoMultiTextureCoordinateElement::get2(const int unit, const int index) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];

  assert(index >= 0 && index < ud.numCoords);
  assert(ud.whatKind == SoTextureCoordinateElement::EXPLICIT);
  if (ud.coordsDimension == 2) {
    return ud.coords2[index];
  }
  else {
    // need an instance we can write to
    SoMultiTextureCoordinateElement * elem = const_cast<SoMultiTextureCoordinateElement *>(this);

    if (ud.coordsDimension == 4) {
      float tmp = ud.coords4[index][3];
      float to2D = tmp == 0.0f ? 1.0f : 1.0f / tmp;

      elem->convert2.setValue(ud.coords4[index][0] * to2D,
                              ud.coords4[index][1] * to2D);
    }
    else { // coordsDimension == 3
      elem->convert2.setValue(ud.coords3[index][0],
                              ud.coords3[index][1]);
    }
    return this->convert2;
  }
}

/*!
  FIXME: write doc.

*/
const SbVec3f &
SoMultiTextureCoordinateElement::get3(const int unit, const int index) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];

  assert(index >= 0 && index < ud.numCoords);
  assert(ud.whatKind == SoTextureCoordinateElement::EXPLICIT);
  if (ud.coordsDimension == 3) {
    return ud.coords3[index];
  }
  else {
    // need an instance we can write to
    SoMultiTextureCoordinateElement * elem =
      const_cast<SoMultiTextureCoordinateElement *>(this);

    if (ud.coordsDimension==2) {
      elem->convert3.setValue(ud.coords2[index][0],
                              ud.coords2[index][1],
                              0.0f);
    }
    else { // this->coordsDimension==4
      ud.coords4[index].getReal(elem->convert3);
    }
    return this->convert3;
  }
}

//!  FIXME: write doc.

const SbVec4f &
SoMultiTextureCoordinateElement::get4(const int unit, const int index) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];

  assert(index >= 0 && index < ud.numCoords);
  assert(ud.whatKind == SoTextureCoordinateElement::EXPLICIT);
  if (ud.coordsDimension==4) {
    return ud.coords4[index];
  }
  else {
    // need an instance we can write to
    SoMultiTextureCoordinateElement * elem =
      const_cast<SoMultiTextureCoordinateElement *>(this);
    if (ud.coordsDimension == 2) {
      elem->convert4.setValue(ud.coords2[index][0],
                              ud.coords2[index][1],
                              0.0f,
                              1.0f);
    }
    else { // this->coordsDimension==3
      elem->convert4.setValue(ud.coords3[index][0],
                              ud.coords3[index][1],
                              ud.coords3[index][2],
                              1.0f);
    }
    return this->convert4;
  }
}

/*!
  This method is used by shapes.  Three return values are possible.

  DEFAULT means that the shapes should generate their own texture coordinates.

  EXPLICIT means that discrete texture coordinates are stored, and should be
  fetched with get2(), get3() or get4().

  FUNCTION means that get(point, normal) must be used to generate texture
  coordinates.
*/

SoTextureCoordinateElement::CoordType
SoMultiTextureCoordinateElement::getType(SoState * const state, const int unit)
{
  const SoMultiTextureCoordinateElement * element =
    coin_assert_cast<const SoMultiTextureCoordinateElement *>
    (getConstElement(state, classStackIndex));
  return element->getType(unit);
}

//! FIXME: write doc.

SoTextureCoordinateElement::CoordType
SoMultiTextureCoordinateElement::getType(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.whatKind;
}

//! FIXME: write doc.

void
SoMultiTextureCoordinateElement::init(SoState * state)
{
  inherited::init(state);
  for (int i = 0; i < MAX_UNITS; i++) {
    UnitData & ud = PRIVATE(this)->unitdata[i];
    ud.nodeid = 0;
    ud.whatKind = SoTextureCoordinateElement::DEFAULT;
    ud.funcCB = NULL;
    ud.funcCBData = NULL;
    ud.numCoords = 0;
    ud.coords2 = NULL;
    ud.coords3 = NULL;
    ud.coords4 = NULL;
    ud.coordsDimension = 2;
  }
}

//! FIXME: write doc.

//$ EXPORT INLINE
int32_t
SoMultiTextureCoordinateElement::getNum(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.numCoords;
}

//! FIXME: write doc. (for backwards compability. Use getDimension() instead).

//$ EXPORT INLINE
SbBool
SoMultiTextureCoordinateElement::is2D(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return (ud.coordsDimension==2);
}

/*!
  FIXME: write doc.
*/
int32_t
SoMultiTextureCoordinateElement::getDimension(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.coordsDimension;
}

/*!
  Returns a pointer to the 2D texture coordinate array. This method is not
  part of the OIV API.
*/
const SbVec2f *
SoMultiTextureCoordinateElement::getArrayPtr2(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.coords2;
}

/*!
  Returns a pointer to the 3D texture coordinate array.

*/
const SbVec3f *
SoMultiTextureCoordinateElement::getArrayPtr3(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.coords3;
}

/*!
  Returns a pointer to the 4D texture coordinate array. This method is not
  part of the OIV API.
*/
const SbVec4f *
SoMultiTextureCoordinateElement::getArrayPtr4(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  const UnitData & ud = PRIVATE(this)->unitdata[unit];
  return ud.coords4;
}

void
SoMultiTextureCoordinateElement::push(SoState * COIN_UNUSED_ARG(state))
{
  SoMultiTextureCoordinateElement * prev =
    coin_assert_cast<SoMultiTextureCoordinateElement *>
    (
     this->getNextInStack()
    );

  for (int i = 0; i < MAX_UNITS; i++) {
    PRIVATE(this)->unitdata[i] = PRIVATE(prev)->unitdata[i];
  }
}

SbBool
SoMultiTextureCoordinateElement::matches(const SoElement * elem) const
{
  const SoMultiTextureCoordinateElement * e =
    coin_assert_cast<const SoMultiTextureCoordinateElement *>(elem);
  for (int i = 0; i < MAX_UNITS; i++) {
    if (PRIVATE(e)->unitdata[i].nodeid != PRIVATE(this)->unitdata[i].nodeid) {
      return FALSE;
    }
  }
  return TRUE;
}

SoElement *
SoMultiTextureCoordinateElement::copyMatchInfo(void) const
{
  SoMultiTextureCoordinateElement * elem =
    static_cast<SoMultiTextureCoordinateElement *>(getTypeId().createInstance());
  for (int i = 0; i < MAX_UNITS; i++) {
    PRIVATE(elem)->unitdata[i].nodeid = PRIVATE(this)->unitdata[i].nodeid;
  }
  return elem;
}

/*!
  Returns the per-unit data for this element.
*/
SoMultiTextureCoordinateElement::UnitData &
SoMultiTextureCoordinateElement::getUnitData(const int unit)
{
  assert(unit >= 0 && unit < MAX_UNITS);
  return PRIVATE(this)->unitdata[unit];
}

const SoMultiTextureCoordinateElement::UnitData &
SoMultiTextureCoordinateElement::getUnitData(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  return PRIVATE(this)->unitdata[unit];
}



#undef MAX_UNITS
#undef PRIVATE
