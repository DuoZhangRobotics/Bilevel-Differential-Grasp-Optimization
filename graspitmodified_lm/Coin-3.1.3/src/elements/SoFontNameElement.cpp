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
  \class SoFontNameElement Inventor/elements/SoFontNameElement.h
  \brief The SoFontNameElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include "tidbitsp.h"
#include "SbBasicP.h"

#include <Inventor/elements/SoFontNameElement.h>

#include <cassert>



SbName * SoFontNameElement::defaultfontname = NULL;

/*!
  \fn SoFontNameElement::fontName

  FIXME: write doc.
*/

SO_ELEMENT_SOURCE(SoFontNameElement);

/*!
  This static method initializes static data for the SoFontNameElement
  class.
*/

void
SoFontNameElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoFontNameElement, inherited);

  SoFontNameElement::defaultfontname = new SbName("defaultFont");

  coin_atexit(reinterpret_cast<coin_atexit_f *>(SoFontNameElement::clean), CC_ATEXIT_NORMAL);
}

void
SoFontNameElement::clean(void)
{
  delete SoFontNameElement::defaultfontname;
}

/*!
  The destructor.
*/

SoFontNameElement::~SoFontNameElement()
{
}

//! FIXME: write doc.

void
SoFontNameElement::set(SoState * const state,
                       SoNode * const node,
                       const SbName fontName)
{
  SoFontNameElement * element =
    coin_safe_cast<SoFontNameElement *>
    (
     SoReplacedElement::getElement(state, classStackIndex, node)
     );

  if (element) {
    element->fontName = fontName;
  }
}

//! FIXME: write doc.

const SbName &
SoFontNameElement::get(SoState * const state)
{
  const SoFontNameElement * element = coin_assert_cast<const SoFontNameElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return element->fontName;
}

//! FIXME: write doc.

SbBool
SoFontNameElement::matches(const SoElement * element) const
{
  if (this == element)
    return TRUE;
  if (element->getTypeId() != SoFontNameElement::getClassTypeId())
    return FALSE;
  if (this->fontName != coin_assert_cast<const SoFontNameElement *>(element)->fontName)
    return FALSE;
  return TRUE;
}

//! FIXME: write doc.

SoElement *
SoFontNameElement::copyMatchInfo(void) const
{
  SoFontNameElement * element = static_cast<SoFontNameElement *>
    (
     SoFontNameElement::getClassTypeId().createInstance()
     );
  element->fontName = this->fontName;
  element->nodeId = this->nodeId;
  return element;
}

//! FIXME: write doc.

void
SoFontNameElement::print(FILE * file) const
{
  fprintf(file, "SoFontNameElement[%p]: font = %s\n", this,
           this->fontName.getString());
}

//! FIXME: write doc.

void
SoFontNameElement::init(SoState * state)
{
  inherited::init(state);
  this->fontName = *SoFontNameElement::defaultfontname;
}

//! FIXME: write doc.

SbName
SoFontNameElement::getDefault(void)
{
  return *SoFontNameElement::defaultfontname;
}
