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
  \class SoEventCallback SoEventCallback.h Inventor/nodes/SoEventCallback.h
  \brief The SoEventCallback class provides functionality for catching events.
  \ingroup nodes

  Use SoEventCallback nodes in the scenegraph for catching user
  interaction events with the scenegraph's render canvas.


  This is how event handling works in Coin: when the user interacts
  with the render canvas, for instance by using the mouse pointer or
  by hitting the keyboard, the GUI interface toolkit (ie SoQt, SoWin,
  SoXt, Sc21 ...) will catch the event and translate it from a
  windowsystem-specific event to a generic Coin event. (For the types
  of generic Coin events, see the classes derived from SoEvent.)  This
  event will then be wrapped inside a SoHandleEventAction and applied
  to the scenegraph.  All this happens within the So[Qt|Xt|Win|...]
  toolkit.

  The SoHandleEventAction then traverses the scenegraph, delivering
  the event to any node type which "is interested" in it.  The
  SoEventCallback nodetype catches the action and forwards the event
  to a callback function set up by the application programmer.

  Be careful about which position in the scenegraph you insert
  SoEventCallback nodes if you are also using any of the built-in Coin
  library elements which are interested in user interaction events
  (like for instance the dragger and manipulator classes and the
  SoSelection nodes).  These Coin elements might catch the event for
  themselves, short-circuiting the SoHandleEventAction traversal so
  the event will never reach the SoEventCallback node(s) you insert.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    EventCallback {
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoEventCallback.h>

#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/events/SoEvent.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \typedef void SoEventCallbackCB(void * userdata, SoEventCallback * node)

  Callback functions for SoEventCallback::addEventCallback() must be
  of this type.  \a userdata is the last argument to
  SoEventCallback::addEventCallback(), and \a node is of course the
  SoEventCallback node in the scenegraph which caused the invocation
  of the callback.
*/

// *************************************************************************

SO_NODE_SOURCE(SoEventCallback);

// *************************************************************************

/*!
  Constructor.
*/
SoEventCallback::SoEventCallback(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoEventCallback);

  this->heaction = NULL;
  this->path = NULL;
}

/*!
  Destructor.
*/
SoEventCallback::~SoEventCallback()
{
  if (this->path) this->path->unref();
}

// Doc from superclass.
void
SoEventCallback::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoEventCallback, SO_FROM_INVENTOR_1);
}

// *************************************************************************

/*!
  Sets the path that must be picked before the registered callbacks
  are invoked. If \c NULL, callbacks will be invoked for every event
  that matches the callback event type.

  \sa getPath()
*/
void
SoEventCallback::setPath(SoPath * pathptr)
{
  if (this->path) {
    this->path->unref();
    this->path = NULL;
  }
  if (pathptr) {
#if COIN_DEBUG
    if (pathptr->getRefCount() == 0) {
      SoDebugError::postWarning("SoEventCallback::setPath",
                                "input path has reference count equal to zero");
    }
#endif // COIN_DEBUG
    this->path = pathptr->copy();
    this->path->ref();
  }
}

/*!
  Returns the path that must be picked before callbacks are invoked.

  \sa setPath()
*/
const SoPath *
SoEventCallback::getPath(void)
{
  return this->path;
}


/*!
  Set up a callback function \a f which will be invoked for the given
  \a eventtype. \a userdata will be given as the first argument to the
  function.
*/
void
SoEventCallback::addEventCallback(SoType eventtype, SoEventCallbackCB * f,
                                  void * userdata)
{
  struct CallbackInfo cb;
  cb.func = f;
  cb.eventtype = eventtype;
  cb.userdata = userdata;

  this->callbacks.append(cb);
}

/*!
  Unregister the given callback function \a f.
*/
void
SoEventCallback::removeEventCallback(SoType eventtype, SoEventCallbackCB * f,
                                     void * userdata)
{
  for (int i = 0; i < this->callbacks.getLength(); i++) {
    if (this->callbacks[i].func == f &&
        this->callbacks[i].eventtype == eventtype &&
        this->callbacks[i].userdata == userdata) {
      this->callbacks.remove(i);
      return;
    }
  }

#if COIN_DEBUG
  SoDebugError::postWarning("SoEventCallback::removeEventCallback",
                            "tried to remove non-existent callback function");
#endif // COIN_DEBUG
}


/*!
  Returns the SoHandleEventAction instance currently traversing the
  scene graph with the SoEvent-derived event object.
*/
SoHandleEventAction *
SoEventCallback::getAction(void) const
{
  return this->heaction;
}

/*!
  Returns a pointer to the event object which is currently being sent
  through the scenegraph.

  If your application code handles the event, you probably want to
  call SoEventCallback::setHandled() to notify the SoHandleEventAction
  that it should stop traversing the scenegraph with the event.
*/
const SoEvent *
SoEventCallback::getEvent(void) const
{
  return (this->heaction ? this->heaction->getEvent() : NULL);
}

/*!
  Returns the picked point for the current handle event traversal.

  This is obviously only related to events which can be considered
  "pick-style" events, like mousebutton presses.
*/
const SoPickedPoint *
SoEventCallback::getPickedPoint(void) const
{
  return this->heaction ? this->heaction->getPickedPoint() : NULL;
}


/*!
  Use this method from a callback function to notify the node that the
  event has been handled.

  The rest of the callbacks registered with the node will still be
  called, but further SoEventCallback nodes in the scene will not be
  notified about the event, neither will any other Coin elements in
  the scenegraph (like for instance SoDragger objects, SoSelection
  nodes or manipulators).

  Since callbacks registered within the same SoEventCallback node will
  still be invoked after the event has been handled, it is likely that
  you should use SoEventCallback::isHandled() to check for this
  condition from your callback functions.
*/
void
SoEventCallback::setHandled(void)
{
#if COIN_DEBUG
  if (!this->heaction) {
    SoDebugError::postWarning("SoEventCallback::setHandled",
                              "should only be called from event callbacks");
    return;
  }
#endif // COIN_DEBUG

  this->heaction->setHandled();
}

/*!
  Check whether or not the event from the SoHandleEventAction has been
  handled.
*/
SbBool
SoEventCallback::isHandled(void) const
{
#if COIN_DEBUG
  if (!this->heaction) {
    SoDebugError::postWarning("SoEventCallback::isHandled",
                              "should only be called from event callbacks");
    return TRUE;
  }
#endif // COIN_DEBUG

  return this->heaction->isHandled();
}


/*!
  Set up the node so all future events (until releaseEvents() is
  called) in Coin will be directly forwarded to this node.
*/
void
SoEventCallback::grabEvents(void)
{
#if COIN_DEBUG
  if (!this->heaction) {
    SoDebugError::postWarning("SoEventCallback::grabEvents",
                              "should only be called from event callbacks");
    return;
  }
#endif // COIN_DEBUG

  this->heaction->setGrabber(this);
}

/*!
  Do not grab event handling any more.

  \sa grabEvents()
*/
void
SoEventCallback::releaseEvents(void)
{
#if COIN_DEBUG
  if (!this->heaction) {
    SoDebugError::postWarning("SoEventCallback::releaseEvents",
                              "should only be called from event callbacks");
    return;
  }
#endif // COIN_DEBUG

  this->heaction->releaseGrabber();
}

/*!
  Invokes the registered callback functions.
 */
void
SoEventCallback::handleEvent(SoHandleEventAction * action)
{
  // check if correct path is picked
  if (this->path) {
    const SoPickedPoint * pp = action->getPickedPoint();
    if (!pp || !pp->getPath()->containsPath(this->path)) return;
  }

  // Make it possible to access the action object from the callback
  // code.
  SoHandleEventAction * tmphandle = this->heaction = action;

  SoType eventtype = this->heaction->getEvent()->getTypeId();

  // Invoke callbacks.
  for (int i = 0; i < this->callbacks.getLength(); i++) {
    if (eventtype.isDerivedFrom(this->callbacks[i].eventtype)) {
      SoEventCallbackCB * cb = this->callbacks[i].func;
      cb(this->callbacks[i].userdata, this);
    }
  }

  // Reset.
  //
  // FIXME: there is a fundamental flaw with this technique; if a new
  // SoHandleEventAction is applied to the scene graph while the
  // callbacks above are still processing, a new pointer will
  // overwrite the previous one (which is a problem that may cause
  // hard-to-find errors), and then it will be overwritten with a NULL
  // pointer below, which causes getAction(), getEvent() etc to return
  // NULL values.
  //
  // The fundamental design flaw here, as I see it, is that the
  // callback invocations should pass on the SoHandleEventAction
  // pointer, instead of temporarily storing it as a member of this
  // node. We can't really change this without breaking backward API
  // compatibility, though.
  //
  // I've added some detection code below to at least catch and warn
  // about the problem.
  //
  // 20041114 mortene.

  if (tmphandle != this->heaction) {
    SoDebugError::postWarning("SoEventCallback::handleEvent",
                              "detected additional SoHandleEventAction scene "
                              "graph traversal while still handling a previous "
                              "traversal -- this can cause hard-to-find bugs, "
                              "and should be avoided");
  }

  this->heaction = NULL;
}
