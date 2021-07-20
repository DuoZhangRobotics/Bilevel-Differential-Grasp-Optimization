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
  \class SoListenerOrientationElement Inventor/elements/SoListenerOrientationElement.h
  \brief The SoListenerOrientationElement holds the orientation of the current listener.
  \ingroup elements

  This orientation is set by SoListener nodes and SoCamera Nodes during
  audio rendering. When a SoListener is visited by the SoAudioRenderAction,
  it will add a new SoListenerOrientationElement to the state, holding it's
  orientation and with the setbylistener flag set. When a SoCamera is
  visited by SoAudioRenderAction, it will add a new
  SoListenerOrientationElement only if there are no previous elements with
  the setbylistener flag set.

  The SoListenerOrientationElement is used when the SoVRMLSound nodes render
  themselves.

  \COIN_CLASS_EXTENSION

  \since Coin 2.0
*/

#include <Inventor/elements/SoListenerOrientationElement.h>

#include "coindefs.h"
#include "SbBasicP.h"

#include <Inventor/nodes/SoNode.h>

/*!
  \fn SoListenerOrientationElement::orientation

  The orientation of the listener. Can be set by the SoListener class
  or the SoCamera class.
*/

SO_ELEMENT_SOURCE(SoListenerOrientationElement);

/*!
  This static method initializes static data for the
  SoListenerOrientationElement class.
*/

void
SoListenerOrientationElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoListenerOrientationElement, inherited);
}

/*!
  The destructor.
*/

SoListenerOrientationElement::~SoListenerOrientationElement(void)
{
}

/*!
  Initializes the element to it's default value. The default
  value for the orientation is (0.0f, 0.0f, 1.0f, 0.0f) and the
  default value for the setByListener flag is FALSE.
*/

void
SoListenerOrientationElement::init(SoState * state)
{
  inherited::init(state);
  this->orientation = SbRotation(0.0f, 0.0f, 1.0f, 0.0f);
  this->setbylistener = FALSE;
}

/*!
  Sets the current listener orientation, and indicates if it was set
  by a SoListener node or a SoCamera node.
*/

void
SoListenerOrientationElement::set(SoState * const state,
                                  SoNode * const COIN_UNUSED_ARG(node),
                                  const SbRotation & orientation,
                                  SbBool setbylistener)
{
  SoListenerOrientationElement *elem =
    coin_safe_cast<SoListenerOrientationElement *>
    (
     SoElement::getElement(state,classStackIndex)
     );
  if (elem) {
    elem->orientation = orientation;
    elem->setbylistener = setbylistener;
  }
}

//! Returns the current listener orientation

const SbRotation &
SoListenerOrientationElement::get(SoState * const state)
{
  const SoListenerOrientationElement * elem =
    coin_assert_cast<const SoListenerOrientationElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return elem->orientation;
}

/*!
  Returns TRUE if the orientation was set by a SoListener node,
  and FALSE if it was set by a SoCamera node
*/

SbBool
SoListenerOrientationElement::isSetByListener(SoState * const state)
{
  const SoListenerOrientationElement * elem =
    coin_assert_cast<const SoListenerOrientationElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
  return elem->setbylistener;
}

//! Prints the contents of the element (unimplemented)

void
SoListenerOrientationElement::print(FILE * /* file */) const
{
}
