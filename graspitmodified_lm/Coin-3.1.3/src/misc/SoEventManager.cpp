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

#include <Inventor/SoEventManager.h>

#include <vector>
#include <algorithm>

#include <Inventor/SbViewportRegion.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/misc/SoState.h>

#include <Inventor/scxml/ScXML.h>
#include <Inventor/scxml/SoScXMLStateMachine.h>

#include "SbBasicP.h"

/*!
  \class SoEventManager SoEventManager.h Inventor/SoEventManager.h
  \brief The SoEventManager class provides event handling for a Coin3D viewer.

  It receives Coin events (typically translated from GUI events by the GUI 
  binding in use) and passes them to the state machine or directly to the
  scene graph depending on the navigation mode. SoEventManager implements the
  event interface for the Coin3D navigation system based on ScXML. 
  SoEventManager and SoRenderManager together supersede SoSceneManager, 
  which is now just a wrapper around these classes.
  The class is used by the QuarterWidget in Quarter (the Qt binding for Coin3D). 
*/

/*!
  \enum SoEventManager::NavigationState

  Sets how events are handled.
*/

/*!
  \var SoEventManager::NavigationState SoEventManager::NO_NAVIGATION

  Forwards the events only to the scene graph.
*/

/*!
  \var SoEventManager::NavigationState SoEventManager::JUST_NAVIGATION

  Forwards the events only to the state machines.
*/

/*!
  \var SoEventManager::NavigationState SoEventManager::MIXED_NAVIGATION

  Forwards the events to the scene graph first. If it does not get handled,
  the events get forwarded to the state machines.
*/

class SoEventManager::PImpl {
public:
  SoEventManager::NavigationState navigationstate;
  SoHandleEventAction * handleeventaction;
  SoSearchAction * searchaction;
  SbBool deletehandleeventaction;
  SoCamera * camera;
  SoNode * scene;

  std::vector<SoScXMLStateMachine *> statemachines;
}; // PImpl

#define PRIVATE(p) (p->pimpl)
#define PUBLIC(p) (p->publ)

SoEventManager::SoEventManager(void)
{
  PRIVATE(this)->navigationstate = SoEventManager::NO_NAVIGATION;
  PRIVATE(this)->handleeventaction = new SoHandleEventAction(SbViewportRegion(400, 400));
  PRIVATE(this)->deletehandleeventaction = TRUE;
  PRIVATE(this)->searchaction = new SoSearchAction;

  PRIVATE(this)->camera = NULL;
  PRIVATE(this)->scene = NULL;
}

/*!
  Destructor.

  The destructor destructs all subobjects it still has handles on when
  invoked.
*/
SoEventManager::~SoEventManager()
{
  this->setCamera(NULL);
  this->setSceneGraph(NULL);

  delete PRIVATE(this)->searchaction;

  if (PRIVATE(this)->deletehandleeventaction) {
    delete PRIVATE(this)->handleeventaction;
    PRIVATE(this)->handleeventaction = NULL;
  }
  for (int c = this->getNumSoScXMLStateMachines() - 1; c >= 0; --c) {
    SoScXMLStateMachine * sm = this->getSoScXMLStateMachine(c);
    this->removeSoScXMLStateMachine(sm);
    delete sm;
  }
}

/*!
  This method sets the navigation mode of the scene manager.  The navigation
  mode decides the system the user uses to navigate the 3D model.

  SoEventManager::NO_NAVIGATION is the default setting.

  \sa SoEventManager::NavigationState, getNavigationState
*/

void
SoEventManager::setNavigationState(NavigationState state)
{
  PRIVATE(this)->navigationstate = state;
}

/*!
  This method returns which state the camera navigation system is in.

  \sa SoEventManager::NavigationState, setNavigationState
*/
SoEventManager::NavigationState
SoEventManager::getNavigationState(void) const
{
  return PRIVATE(this)->navigationstate;
}

/*!
  Set the node which is top of the scene graph we're managing.  The \a
  sceneroot node reference count will be increased by 1, and any
  previously set scene graph top node will have it's reference count
  decreased by 1.

  \sa getSceneGraph
*/
void
SoEventManager::setSceneGraph(SoNode * const sceneroot)
{
  // Don't unref() until after we've set up the new root, in case the
  // old root == the new sceneroot. (Just to be that bit more robust.)
  SoNode * oldroot = PRIVATE(this)->scene;

  PRIVATE(this)->scene = sceneroot;

  if (PRIVATE(this)->scene)
    PRIVATE(this)->scene->ref();

  if (oldroot)
    oldroot->unref();

  for (int c = 0; c < this->getNumSoScXMLStateMachines(); ++c) {
    SoScXMLStateMachine * sm = this->getSoScXMLStateMachine(c);
    sm->setSceneGraphRoot(PRIVATE(this)->scene);
  }
}

/*!
  Returns pointer to root of scene graph.
*/
SoNode *
SoEventManager::getSceneGraph(void) const
{
  return PRIVATE(this)->scene;
}

/*!
  Sets the camera to be used.
*/
void
SoEventManager::setCamera(SoCamera * camera)
{
  SoCamera * oldcamera = PRIVATE(this)->camera;

  PRIVATE(this)->camera = camera;

  if (PRIVATE(this)->camera)
    PRIVATE(this)->camera->ref();

  if (oldcamera)
    oldcamera->unref();

  for (int i = this->getNumSoScXMLStateMachines() - 1; i >= 0; --i) {
    SoScXMLStateMachine * sm = this->getSoScXMLStateMachine(i);
    sm->setActiveCamera(PRIVATE(this)->camera);
  }
}

/*!
  Returns the current camera.
*/
SoCamera *
SoEventManager::getCamera(void) const
{
  return PRIVATE(this)->camera;
}


/*!
  Set size of rendering area for the viewport within the current
  window.
*/
void
SoEventManager::setSize(const SbVec2s & newsize)
{
  SbViewportRegion region = PRIVATE(this)->handleeventaction->getViewportRegion();
  SbVec2s origin = region.getViewportOriginPixels();
  region.setViewportPixels(origin, newsize);
  PRIVATE(this)->handleeventaction->setViewportRegion(region);
}

/*!
  Set \e only the origin of the viewport region within the rendering
  window.

  \sa setViewportRegion(), setWindowSize()
*/
void
SoEventManager::setOrigin(const SbVec2s & newOrigin)
{
  SbViewportRegion region = PRIVATE(this)->handleeventaction->getViewportRegion();
  SbVec2s size = region.getViewportSizePixels();
  region.setViewportPixels(newOrigin, size);
  PRIVATE(this)->handleeventaction->setViewportRegion(region);
}

/*!
  Sets the current viewport region. This will overwrite the default viewport region
  created in the constructor.

  \sa getViewportRegion()
*/
void
SoEventManager::setViewportRegion(const SbViewportRegion & newregion)
{
  PRIVATE(this)->handleeventaction->setViewportRegion(newregion);
}

/*!
  Returns the viewport region used by the event manager. 
*/
const SbViewportRegion &
SoEventManager::getViewportRegion(void) const
{
  return PRIVATE(this)->handleeventaction->getViewportRegion();
}

/*!
  Handles the event. Depending on the navigation state, this forwards the event
  to the state machines and/or the scene graph.
*/
SbBool
SoEventManager::processEvent(const SoEvent * const event)
{
  const SbViewportRegion & vp =
    PRIVATE(this)->handleeventaction->getViewportRegion();

  SbBool status = FALSE;

  int i = 0;
  switch (PRIVATE(this)->navigationstate) {
  case SoEventManager::NO_NAVIGATION:
    status = this->actuallyProcessEvent(event);
    break;
  case SoEventManager::JUST_NAVIGATION:
    for (i = this->getNumSoScXMLStateMachines() - 1; i >= 0; --i) {
      SoScXMLStateMachine * sm = this->getSoScXMLStateMachine(i);
      if (sm->isActive()) {
        sm->setViewportRegion(vp);
        if (sm->processSoEvent(event))
          status = TRUE;
      }
    }
    break;
  case SoEventManager::MIXED_NAVIGATION:
    if (this->actuallyProcessEvent(event)) {
      status = TRUE;
      break;
    }
    for (i = this->getNumSoScXMLStateMachines() - 1; i >= 0; --i) {
      SoScXMLStateMachine * sm = this->getSoScXMLStateMachine(i);
      if (sm->isActive()) {
        sm->setViewportRegion(vp);
        if (sm->processSoEvent(event))
          status = TRUE;
      }
    }
  }
  return status;
}

/*!
  Forwards the event to the scene graph.
*/
SbBool
SoEventManager::actuallyProcessEvent(const SoEvent * const event)
{
  assert(PRIVATE(this)->handleeventaction);

  SbBool handled = FALSE;
  if ( PRIVATE(this)->handleeventaction->getState() != NULL &&
       PRIVATE(this)->handleeventaction->getState()->getDepth() != 0 ) {
    // recursive invocation - action currently in use
#if COIN_DEBUG
    SoDebugError::post("SoEventManager::processEvent",
                       "Recursive invocation detected. Delay processing event "
                       "until the current event is finished processing.");
#endif // COIN_DEBUG
  } else if ( PRIVATE(this)->scene == NULL ) {
    // nothing
  } else {
    PRIVATE(this)->handleeventaction->setEvent(event);
    PRIVATE(this)->handleeventaction->apply(PRIVATE(this)->scene);
    handled = PRIVATE(this)->handleeventaction->isHandled();
  }

  return handled;
}

/*!
  Set the \a action to use for event handling. Overrides the default
  action made in the constructor.
 */
void
SoEventManager::setHandleEventAction(SoHandleEventAction * handleeventaction)
{
  assert(handleeventaction && "SoEventManager::setHandleEventAction, action == NULL");

  // remember old viewport region
  SbViewportRegion region = PRIVATE(this)->handleeventaction->getViewportRegion();
  if (PRIVATE(this)->deletehandleeventaction) {
    delete PRIVATE(this)->handleeventaction;
  }
  PRIVATE(this)->deletehandleeventaction = FALSE;
  PRIVATE(this)->handleeventaction = handleeventaction;
  PRIVATE(this)->handleeventaction->setViewportRegion(region);
}

/*!
  Returns pointer to event handler action.
 */
SoHandleEventAction *
SoEventManager::getHandleEventAction(void) const
{
  return PRIVATE(this)->handleeventaction;
}

/*!
  Returns the number of SCXML state machines registered on the SoEventManager
  object.
*/
int
SoEventManager::getNumSoScXMLStateMachines(void) const
{
  return static_cast<int>(PRIVATE(this)->statemachines.size());
}

/*!
  Returns the pointer to the Nth (\a idx) SCXML state machine registered
  on the SoSceneManager object.

  The \a idx argument must be a valid index, not outside the actual range
  of SCXML state machine object indices.
*/
SoScXMLStateMachine *
SoEventManager::getSoScXMLStateMachine(int idx) const
{
  assert(idx >= 0 &&
         idx < static_cast<int>(PRIVATE(this)->statemachines.size()));
  return PRIVATE(this)->statemachines.at(idx);
}

/*!
  Adds an SCXML state machine object to the SoEventManager's event pipeline.
  If it is not removed before SoSceneManager destruction, the SoSceneManager
  destructor will delete it.
*/
void
SoEventManager::addSoScXMLStateMachine(SoScXMLStateMachine * sm)
{
  PRIVATE(this)->statemachines.push_back(sm);
}

/*!
  Removes the SoSceneManager object reference to an SCXML state machine.
  It will just be removed, not destructed.

  If no reference to the given SCXML state machine exists, nothing will
  happen.
*/
void
SoEventManager::removeSoScXMLStateMachine(SoScXMLStateMachine * sm)
{
  std::vector<SoScXMLStateMachine *>::iterator findit =
    std::find(PRIVATE(this)->statemachines.begin(),
              PRIVATE(this)->statemachines.end(), sm);
  if (findit != PRIVATE(this)->statemachines.end()) {
    PRIVATE(this)->statemachines.erase(findit);
  }
}

#undef PRIVATE
#undef PUBLIC
