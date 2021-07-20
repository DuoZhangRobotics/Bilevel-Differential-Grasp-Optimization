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
  \class SoAnnoText3CharOrientElement Inventor/elements/SoAnnoText3CharOrientElement.h
  \brief The SoAnnoText3CharOrientElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoAnnoText3CharOrientElement.h>


#include <cassert>

SO_ELEMENT_SOURCE(SoAnnoText3CharOrientElement);

/*!
  This static method initializes static data for the
  SoAnnoText3CharOrientElement class.
*/

void
SoAnnoText3CharOrientElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoAnnoText3CharOrientElement, inherited);
}

/*!
  The destructor.
*/

SoAnnoText3CharOrientElement::~SoAnnoText3CharOrientElement(// virtual protected
    void)
{
}

//! FIXME: write doc.

void
SoAnnoText3CharOrientElement::init(SoState * state)
{
  inherited::init(state);
}

//! FIXME: write doc.

//$ EXPORT INLINE
void
SoAnnoText3CharOrientElement::set(SoState * const state, SbBool isOriented)
{
  inherited::set(classStackIndex, state, isOriented);
}

//! FIXME: write doc.

//$ EXPORT INLINE
SbBool
SoAnnoText3CharOrientElement::get(SoState * state)
{
  return static_cast<SbBool>(SoInt32Element::get(classStackIndex, state));
}

//! FIXME: write doc.

//$ EXPORT INLINE
SbBool
SoAnnoText3CharOrientElement::getDefault(void)
{
  return TRUE;
}
