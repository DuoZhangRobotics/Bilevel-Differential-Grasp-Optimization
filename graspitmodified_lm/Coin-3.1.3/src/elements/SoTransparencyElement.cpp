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
  \class SoTransparencyElement Inventor/elements/SoTransparencyElement.h
  \brief The SoTransparencyElement is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoTransparencyElement.h>

#include "SbBasicP.h"

#include <Inventor/elements/SoLazyElement.h>
#include <cassert>


SO_ELEMENT_SOURCE(SoTransparencyElement);

/*!
  This static method initializes static data for the
  SoTransparencyElement class.
*/

void
SoTransparencyElement::initClass()
{
  SO_ELEMENT_INIT_CLASS(SoTransparencyElement, inherited);
}

/*!
  The destructor.
*/

SoTransparencyElement::~SoTransparencyElement()
{
}

//! FIXME: write doc.

void
SoTransparencyElement::init(SoState * stateptr)
{
  this->state = stateptr;
}

//! FIXME: write doc.

void
SoTransparencyElement::set(SoState * const state, SoNode * const node,
                           const int32_t numvalues,
                           const float * const values)
{
  SoTransparencyElement * elem = const_cast<SoTransparencyElement *>
    (
     SoTransparencyElement::getInstance(state)
     );

  SoLazyElement::setTransparency(state, node, numvalues, values, &elem->colorpacker);
}

//! FIXME: write doc.

//$ EXPORT INLINE
int32_t
SoTransparencyElement::getNum(void) const
{
  SoLazyElement * lazy = SoLazyElement::getInstance(this->state);
  return lazy->getNumTransparencies();
}

//! FIXME: write doc.

float
SoTransparencyElement::get(const int index) const
{
  return SoLazyElement::getTransparency(this->state, index);
}

/*!
  Returns a pointer to the transparency values. This method is not part of the OIV API.
*/
const float *
SoTransparencyElement::getArrayPtr(void) const
{
  SoLazyElement * lazy = SoLazyElement::getInstance(this->state);
  return lazy->getTransparencyPointer();
}

//! FIXME: write doc.

const SoTransparencyElement *
SoTransparencyElement::getInstance(SoState *state)
{
  return coin_assert_cast<const SoTransparencyElement *>
    (
     state->getElementNoPush(classStackIndex)
     );
}
