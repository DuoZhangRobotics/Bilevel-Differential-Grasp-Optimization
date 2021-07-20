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
  \class SoHandleEventAction SoHandleEventAction.h Inventor/actions/SoHandleEventAction.h
  \brief The SoHandleEventAction class distributes user events to the scene.
  \ingroup actions

  This is the action used by the GUI viewer classes to pass
  interaction events from the window system to the nodes in the scene
  graph.

  SoHandleEventAction also provides the functionality for tracking the
  object currently under the cursor, and functionality for "grabbing"
  the event focus.

  \sa SoEvent
*/

#include <Inventor/actions/SoHandleEventAction.h>

#include <Inventor/SbViewportRegion.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/elements/SoSwitchElement.h>
#include <Inventor/elements/SoViewVolumeElement.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/elements/SoWindowElement.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoInfo.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/misc/SoState.h>

#include "actions/SoSubActionP.h"

// *************************************************************************

// The private data for the SoHandleEventAction.

class SoHandleEventActionP {
public:
  SoHandleEventActionP(void) : owner(NULL) { }

  // Hidden private methods.

  void doPick(SoRayPickAction * ra);
  SoRayPickAction * getPickAction(void);

  // Hidden private variables.

  SbViewportRegion viewport;
  const SoEvent * event;
  SoNode * grabber;
  SoNode * pickroot;
  SbBool pickvalid;
  SbBool didpickall;
  SoRayPickAction * pickaction;

  SoHandleEventAction * owner;
};

#define PRIVATE(obj) ((obj)->pimpl)

// *************************************************************************

SO_ACTION_SOURCE(SoHandleEventAction);

// Overridden from parent class.
void
SoHandleEventAction::initClass(void)
{
  SO_ACTION_INTERNAL_INIT_CLASS(SoHandleEventAction, SoAction);

  SO_ENABLE(SoHandleEventAction, SoSwitchElement);
  SO_ENABLE(SoHandleEventAction, SoViewVolumeElement);
  SO_ENABLE(SoHandleEventAction, SoViewportRegionElement);
  SO_ENABLE(SoHandleEventAction, SoWindowElement);

}

/*!
  Constructor.

  SoHandleEventAction needs a \a viewportregion to pass on to the
  raypick action instance it uses for being able to track objects
  under the mouse cursor.
*/
SoHandleEventAction::SoHandleEventAction(const SbViewportRegion & viewportregion)
{
  PRIVATE(this)->owner = this;
  PRIVATE(this)->viewport = viewportregion;
  PRIVATE(this)->event = NULL;
  PRIVATE(this)->grabber = NULL;
  PRIVATE(this)->pickroot = NULL;
  PRIVATE(this)->pickvalid = FALSE;
  PRIVATE(this)->didpickall = FALSE;
  PRIVATE(this)->pickaction = NULL;

  SO_ACTION_CONSTRUCTOR(SoHandleEventAction);
}

/*!
  Destructor.
*/
SoHandleEventAction::~SoHandleEventAction()
{
  if (PRIVATE(this)->pickroot) PRIVATE(this)->pickroot->unref();
  delete PRIVATE(this)->pickaction;
}

/*!
  Set a new viewport region, replacing the one passed in the
  constructor.
*/
void
SoHandleEventAction::setViewportRegion(const SbViewportRegion & newregion)
{
  PRIVATE(this)->viewport = newregion;
  if (PRIVATE(this)->pickaction) PRIVATE(this)->pickaction->setViewportRegion(newregion);
}

/*!
  Returns the viewport region this action instance is using.

  Advanced Usage:

  You can also get the viewport region by accessing it through its element
  on the traversal state.  You do that the following way:

  \code
  #include <Inventor/elements/SoViewportRegionElement.h>

    SoState * state = action->getState();
    SbViewportRegion vp = SoViewportRegionElement::get(state);
  \endcode

  The reason for explaining this is that you can use this generic technique
  when you need access to state information when you can't seem to find the
  accessor function you need in the action implementation.  You can use
  it to for instance retrieve the view volume information, for which there
  are no accessor methods:

  \code
  #include <Inventor/elements/SoViewVolumeElement.h>

    SoState * state = action->getState();
    SbViewVolume vv = SoViewVolumeElement::get(state);
  \endcode

  When you do this on arbitrary action instances, you need to make sure
  that the given element is enabled for the action before you try to use it.
  The relevant functions for this are SoState::isElementEnabled() and
  SoElement::getClassStackIndex().
*/

const SbViewportRegion &
SoHandleEventAction::getViewportRegion(void) const
{
  return PRIVATE(this)->viewport;
}

/*!
  Set the event to distribute to the nodes of the scene.
*/
void
SoHandleEventAction::setEvent(const SoEvent * ev)
{
  PRIVATE(this)->event = ev;
}

/*!
  Returns the event this action is handling.
*/
const SoEvent *
SoHandleEventAction::getEvent(void) const
{
  return PRIVATE(this)->event;
}

/*!
  Marks the action instance as handled, hence terminates the action.

  The action is only marked as handled when a node in the graph
  "grabs" the event this action is carrying, so the handled flag will
  be \c FALSE after traversal if no nodes wanted the event.

  \sa isHandled()
*/
void
SoHandleEventAction::setHandled(void)
{
  this->setTerminated(TRUE);
}

/*!
  Returns whether or not the event has been handled by a node during
  scene graph traversal.

  \sa setHandled()
*/
SbBool
SoHandleEventAction::isHandled(void) const
{
  return this->hasTerminated();
}

/*!
  Set a \a node pointer which will get all future events handled by
  this action until releaseGrabber() is called.

  Note that since later SoHandleEventAction invokations are just applied
  directly on the grabber node, using SoHandleEventAction methods like
  getCurPath() will return bogus data.
*/
void
SoHandleEventAction::setGrabber(SoNode * node)
{
  // Check for inequality before executing code is not only good for
  // performance, but is also necessary to remove the potential for
  // infinite recursion. See comment in releaseGrabber().

  if (node != PRIVATE(this)->grabber) {
    this->releaseGrabber();
    PRIVATE(this)->grabber = node;
    if (node) node->grabEventsSetup();
  }
}

/*!
  Don't send the events to a "grabber" node anymore, use the default
  behavior of the action and pass them along to the scene graph again.

  \sa setGrabber()
*/
void
SoHandleEventAction::releaseGrabber(void)
{
  // Store old grabber node and set current node to NULL before
  // calling SoNode::grabEventsCleanup(), to avoid being vulnerable to
  // recursive calls from grabEventsCleanup() back to this method
  // (which happens from dragger classes).

  SoNode * old = PRIVATE(this)->grabber;
  PRIVATE(this)->grabber = NULL;
  if (old) old->grabEventsCleanup();
}

/*!
  Returns the grabber node, or \c NULL if no grabber is active.
*/
SoNode *
SoHandleEventAction::getGrabber(void) const
{
  return PRIVATE(this)->grabber;
}

/*!
  Sets the root \a node that is used for the pick action tracking the
  cursor.
*/
void
SoHandleEventAction::setPickRoot(SoNode * node)
{
  SoNode * oldroot = PRIVATE(this)->pickroot;
  PRIVATE(this)->pickroot = node;
  if (PRIVATE(this)->pickroot) PRIVATE(this)->pickroot->ref();
  if (oldroot) oldroot->unref();
  PRIVATE(this)->pickvalid = FALSE;
}

/*!
  Returns the root node that is used by nodes that is tracking the
  cursor.
*/
SoNode *
SoHandleEventAction::getPickRoot(void) const
{
  return PRIVATE(this)->pickroot;
}

/*!
  Sets the pick radius for cursor tracking.
*/
void
SoHandleEventAction::setPickRadius(const float radiusinpixels)
{
  PRIVATE(this)->getPickAction()->setRadius(radiusinpixels);
}

/*!
  Returns the SoPickedPoint information for the intersection point
  below the cursor.
*/
const SoPickedPoint *
SoHandleEventAction::getPickedPoint(void)
{
  SoRayPickAction * ra = PRIVATE(this)->getPickAction();
  if (!PRIVATE(this)->pickvalid || PRIVATE(this)->didpickall) {
    ra->setPickAll(FALSE);
    PRIVATE(this)->doPick(ra);
  }
  return ra->getPickedPoint();
}

/*!
  Returns a list of all intersection points below the mouse cursor.
*/
const SoPickedPointList &
SoHandleEventAction::getPickedPointList(void)
{
  SoRayPickAction * ra = PRIVATE(this)->getPickAction();
  if (!PRIVATE(this)->pickvalid || !PRIVATE(this)->didpickall) {
    ra->setPickAll(TRUE);
    PRIVATE(this)->doPick(ra);
  }
  return ra->getPickedPointList();
}

// Documented in superclass. Overridden to initialize local data
// members before executing the scene traversal.
void
SoHandleEventAction::beginTraversal(SoNode * node)
{
  assert(PRIVATE(this)->event);
  this->setPickRoot(node);

  this->getState()->push();
  SoViewportRegionElement::set(this->getState(), PRIVATE(this)->viewport);
  if (PRIVATE(this)->grabber) {
    this->traverse(PRIVATE(this)->grabber);
  }
  if (!this->isHandled()) {
    this->traverse(node);
  }
  this->getState()->pop();

  // clear the picked point list
  PRIVATE(this)->getPickAction()->reset();
  PRIVATE(this)->pickvalid = FALSE;
}

//////// Hidden private methods for //////////////////////////////////////
//////// SoHandleEventActionP (pimpl) ////////////////////////////////////

// Singleton pattern for the pick action instance.
SoRayPickAction *
SoHandleEventActionP::getPickAction(void)
{
  if (this->pickaction == NULL) {
    this->pickaction = new SoRayPickAction(this->viewport);
  }
  return this->pickaction;
}

void
SoHandleEventActionP::doPick(SoRayPickAction * ra)
{
  if (!this->event || !this->pickroot) return;

  SbBool didapply = FALSE;
  ra->setPoint(this->event->getPosition());
  if (this->owner->getWhatAppliedTo() == SoAction::PATH) {
    const SoPath * path = this->owner->getPathAppliedTo();
    if (path->getHead() == this->pickroot) {
      ra->apply(const_cast<SoPath *>(path));
      didapply = TRUE;
    }
    else { // make subpath if pickroot can be found in path
      int i, n = path->getLength();
      for (i = 1; i < n; i++) {
        if (path->getNode(i) == this->pickroot) break;
      }
      if (i < n) {
        SoPath * tmppath = path->copy(i);
        tmppath->ref();
        ra->apply(tmppath);
        tmppath->unref();
        didapply = TRUE;
      }
    }
  }
  if (!didapply) ra->apply(this->pickroot);
  this->didpickall = ra->isPickAll();
  this->pickvalid = TRUE;
}

#undef PRIVATE
