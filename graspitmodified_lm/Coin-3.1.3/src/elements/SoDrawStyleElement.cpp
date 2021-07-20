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
  \class SoDrawStyleElement Inventor/elements/SoDrawStyleElement.h
  \brief The SoDrawStyleElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoDrawStyleElement.h>

#include <Inventor/elements/SoShapeStyleElement.h>


#include <cassert>

/*!
  \fn SoDrawStyleElement::Style

  FIXME: write doc.
*/

SO_ELEMENT_SOURCE(SoDrawStyleElement);

/*!
  This static method initializes static data for the
  SoDrawStyleElement class.
*/

void
SoDrawStyleElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoDrawStyleElement, inherited);
}

/*!
  The destructor.
*/

SoDrawStyleElement::~SoDrawStyleElement(void)
{
}

//! FIXME: write doc.

void
SoDrawStyleElement::init(SoState * state)
{
  inherited::init(state);
  this->data = getDefault();
}

//! FIXME: write doc.

void
SoDrawStyleElement::set(SoState * const state,
                        SoNode * const node,
                        const Style style)
{
  SoInt32Element::set(classStackIndex, state, node, static_cast<int32_t>(style));
  SoShapeStyleElement::setDrawStyle(state, static_cast<int32_t>(style));
}

//! FIXME: write doc.

//$ EXPORT INLINE
void
SoDrawStyleElement::set(SoState * const state, const Style style)
{
  set(state, NULL, style);
}

//! FIXME: write doc.

//$ EXPORT INLINE
SoDrawStyleElement::Style
SoDrawStyleElement::get(SoState * const state)
{
  return static_cast<Style>(inherited::get(classStackIndex, state));
}

//! FIXME: write doc.

//$ EXPORT INLINE
SoDrawStyleElement::Style
SoDrawStyleElement::getDefault()
{
  return FILLED;
}
