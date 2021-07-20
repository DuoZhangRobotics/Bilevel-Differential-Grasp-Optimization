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
  \class SoAccumulatedElement SoAccumulatedElement.h Inventor/elements/SoAccumulatedElement.h
  \brief The SoAccumulatedElement class is an abstract class for storing accumulated state.
  \ingroup elements

  This is the superclass of elements where new element data \e
  accumulates with older data.

  The element stores node id values for all nodes accumulated during
  traversal for the current state. These id values are used to
  determine when to invalidate caches.

  \sa SoReplacedElement, SoFloatElement, SoInt32Element
*/

#include <Inventor/elements/SoAccumulatedElement.h>
#include <Inventor/nodes/SoNode.h>
#include <cassert>

#include "SbBasicP.h"

/*!
  \fn SoAccumulatedElement::nodeIds

  Stores the internal list of node id values for nodes accumulated on
  the stack for the element.
*/

SO_ELEMENT_ABSTRACT_SOURCE(SoAccumulatedElement);

// doc from parent
void
SoAccumulatedElement::initClass(void)
{
  SO_ELEMENT_INIT_ABSTRACT_CLASS(SoAccumulatedElement, inherited);
}

SoAccumulatedElement::~SoAccumulatedElement(void)
{
}

void
SoAccumulatedElement::init(SoState * state)
{
  inherited::init(state);
  // this is FALSE until node id's are copied
  this->recursecapture = FALSE;
}

void
SoAccumulatedElement::push(SoState * state)
{
  inherited::push(state);
  // this is FALSE until node id's are copied
  this->recursecapture = FALSE;
}

// Documented in superclass. Overridden to compare node ids.
SbBool
SoAccumulatedElement::matches(const SoElement * element) const
{
  const SoAccumulatedElement * elem = coin_assert_cast<const SoAccumulatedElement *>(element);
  return (elem->nodeIds == this->nodeIds);
}

/*!
  Empty the list of node ids.
*/
void
SoAccumulatedElement::clearNodeIds(void)
{
  this->nodeIds.truncate(0);
  // we do not depend on previous elements any more
  this->recursecapture = FALSE;
}

/*!
  Add the node id of \a node to the list of node ids.
*/
void
SoAccumulatedElement::addNodeId(const SoNode * const node)
{
  this->nodeIds.append(node->getNodeId());
}

/*!
  Empty the list of node ids, and add the id of \a node.
*/
void
SoAccumulatedElement::setNodeId(const SoNode * const node)
{
  this->clearNodeIds();
  this->addNodeId(node);
  // we do not depend on previous elements any more
  this->recursecapture = FALSE;
}

// Documented in superclass. Overridden to copy node ids.
SoElement *
SoAccumulatedElement::copyMatchInfo(void) const
{
  SoAccumulatedElement * element =
    static_cast<SoAccumulatedElement *>(this->getTypeId().createInstance());
  element->copyNodeIds(this);
  return element;
}

/*!
  Convenience method which copies the node ids from \a copyfrom to
  this element.

  \COIN_FUNCTION_EXTENSION
*/
void
SoAccumulatedElement::copyNodeIds(const SoAccumulatedElement * copyfrom)
{
  this->nodeIds = copyfrom->nodeIds;

  // this elements uses data from previous element in stack
  this->recursecapture = TRUE;
}

// Documented in superclass. Overridden to capture more elements.
void
SoAccumulatedElement::captureThis(SoState * state) const
{
  inherited::captureThis(state);

  // we need to recurse if element has copied data from previous
  // element in stack (or nextInStack as SGI was silly enough to call
  // it). This is because the depth of this element might not cause
  // cache to depend on this element, but the previous element(s)
  // might have a depth that will trigger a dependency.
  //                                              pederb, 2001-02-21
  if (this->recursecapture) {
    const SoAccumulatedElement * elem = coin_assert_cast<const SoAccumulatedElement *> (
      this->getNextInStack()
      );
    if (elem) elem->captureThis(state);
  }
}
