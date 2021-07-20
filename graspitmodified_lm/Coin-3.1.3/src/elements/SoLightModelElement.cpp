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
  \class SoLightModelElement Inventor/elements/SoLightModelElement.h
  \brief The SoLightModelElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include "coindefs.h"
#include "SbBasicP.h"

#include <Inventor/elements/SoLightModelElement.h>
#include <Inventor/elements/SoLazyElement.h>
#include <cassert>

/*!
  \fn SoLightModelElement::Model

  FIXME: write doc.
*/

SO_ELEMENT_SOURCE(SoLightModelElement);

/*!
  This static method initializes static data for the
  SoLightModelElement class.
*/

void
SoLightModelElement::initClass()
{
  SO_ELEMENT_INIT_CLASS(SoLightModelElement, inherited);
}

/*!
  The destructor.
*/

SoLightModelElement::~SoLightModelElement()
{
}

//! FIXME: write doc.

void
SoLightModelElement::init(SoState * /* state */)
{
}

//! FIXME: write doc.

void
SoLightModelElement::set(SoState * const state, const Model model)
{
  SoLazyElement::setLightModel(state, static_cast<int32_t>(model));
}

//! FIXME: write doc.

void
SoLightModelElement::set(SoState * const state, SoNode * const COIN_UNUSED_ARG(node),
                         const Model model)
{
  SoLazyElement::setLightModel(state, static_cast<int32_t>(model));
}

//! FIXME: write doc.

SoLightModelElement::Model
SoLightModelElement::get(SoState * const state)
{
  return static_cast<SoLightModelElement::Model>(SoLazyElement::getLightModel(state));
}

//! FIXME: write doc.

SoLightModelElement::Model
SoLightModelElement::getDefault()
{
  return static_cast<SoLightModelElement::Model>(SoLazyElement::getDefaultLightModel());
}

//! FIXME: write doc

const SoLightModelElement *
SoLightModelElement::getInstance(SoState *state)
{
  //FIXME: Can this function or any of the similar functions ever
  //return NULL? BFG 20080916
  return coin_assert_cast<const SoLightModelElement *>
    (
     state->getElementNoPush(classStackIndex)
     );
}
