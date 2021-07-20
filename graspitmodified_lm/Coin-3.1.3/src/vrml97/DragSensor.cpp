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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#ifdef HAVE_VRML97

/*!
  \class SoVRMLDragSensor SoVRMLDragSensor.h Inventor/VRMLnodes/SoVRMLDragSensor.h
  \brief The SoVRMLDragSensor class is a superclass for VRML drag sensors.

  This class collects the two fields that are common for some of the
  sensor nodes, plus some common functions for these. Since this is an
  abstract "helper" class, it does not represent an actual node from
  the VRML97 specification, so don't use it as such.

  For more information, a detailed discussion of drag sensors is
  available in section 4.6.7.4 of the VRML97 specification:

  <http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.7.4>.
*/

#include <Inventor/VRMLnodes/SoVRMLDragSensor.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/elements/SoViewVolumeElement.h>
#include <Inventor/events/SoLocation2Event.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SbMatrix.h>

#include "nodes/SoSubNodeP.h"

/*!
  \var SoSFVec3f SoVRMLDragSensor::trackPoint_changed

  An event out that is generated when the trackpoint is changed.
*/

/*!
  \var SoSFBool SoVRMLDragSensor::autoOffset

  When \c TRUE, the offset field is updated when a dragging sequence
  is finished.  Default value is \c TRUE.
*/

/*!
  \fn SbBool SoVRMLDragSensor::dragStart(void)

  \COININTERNAL

  Called when dragger is selected (picked) by the user.
*/
/*!
  \fn void SoVRMLDragSensor::drag(void)

  \COININTERNAL

  Called when user drags the mouse after picking the dragger.
*/
/*!
  \fn void SoVRMLDragSensor::dragFinish(void)

  \COININTERNAL

  Called when mouse button is released after picking and interacting
  with the dragger.
*/


SO_NODE_ABSTRACT_SOURCE(SoVRMLDragSensor);

// Doc in parent
void
SoVRMLDragSensor::initClass(void)
{
  SO_NODE_INTERNAL_INIT_ABSTRACT_CLASS(SoVRMLDragSensor, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLDragSensor::SoVRMLDragSensor(void)
{
  SO_NODE_CONSTRUCTOR(SoVRMLDragSensor);

  SO_VRMLNODE_ADD_EVENT_OUT(trackPoint_changed);
  SO_VRMLNODE_ADD_EXPOSED_FIELD(autoOffset, (TRUE));
}

/*!
  Destructor.
*/
SoVRMLDragSensor::~SoVRMLDragSensor()
{
}

// Doc in parent
void
SoVRMLDragSensor::handleEvent(SoHandleEventAction * action)
{
  if (!this->enabled.getValue()) {
    if (this->isActive.getValue()) {
      this->isActive = FALSE;
      this->dragFinish();
    }
    inherited::handleEvent(action);
    return;
  }

  const SoEvent * event = action->getEvent();

  SbBool buttondown = SO_MOUSE_PRESS_EVENT(event, BUTTON1);
  SbBool buttonup = SO_MOUSE_RELEASE_EVENT(event, BUTTON1);
  SbBool mousemove = event->isOfType(SoLocation2Event::getClassTypeId());
  SbBool active = this->isActive.getValue();

  if ((!active && mousemove) ||
      (active && buttondown) ||
      (!active && buttonup)) {
    inherited::handleEvent(action);
    return;
  }

  SoState * state = action->getState();

  this->mousepos = event->getPosition();
  this->normpos = event->getNormalizedPosition(SoViewportRegionElement::get(state));
  this->viewvolume = SoViewVolumeElement::get(state);

  if (!active && buttondown) {
    SbBool hit = FALSE;
    const SoPickedPoint * pp = action->getPickedPoint();
    if (pp) {
      const SoFullPath * currpath = (const SoFullPath*) action->getCurPath();
      SoFullPath * parentpath = (SoFullPath*) currpath->copy(0, currpath->getLength()-1);
      SoNode * parentnode = parentpath->getTail();
      parentpath->ref();
      hit = pp->getPath()->containsPath(parentpath);
      parentpath->unref();

      if (hit) {
        // set object space hit point
        this->obj2world = pp->getObjectToWorld(parentnode);
        this->world2obj = this->obj2world.inverse();
        this->world2obj.multVecMatrix(pp->getPoint(), this->hitpt);
        if (this->dragStart()) {
          this->isActive = TRUE;
        }
      }
    }
  }
  else if (active && buttonup) {
    this->dragFinish();
    this->isActive = FALSE;
  }
  else if (active && mousemove) {
    this->drag();
  }
  inherited::handleEvent(action);
}

/*!
  Returns the interaction starting point.
*/
const SbVec3f &
SoVRMLDragSensor::getLocalStartingPoint(void) const
{
  return this->hitpt;
}

/*!
  Returns the matrix that transforms to the world coordinate system.
*/
const SbMatrix &
SoVRMLDragSensor::getLocalToWorldMatrix(void) const
{
  return this->obj2world;
}

/*!
  Returns the matrix that transforms from the world coordinate system.
*/
const SbMatrix &
SoVRMLDragSensor::getWorldToLocalMatrix(void) const
{
  return this->world2obj;
}

/*!
  Returns the current view volume.
*/
const SbViewVolume &
SoVRMLDragSensor::getViewVolume(void) const
{
  return this->viewvolume;
}

/*!
  Returns the current normalized pointer position.
*/
const SbVec2f &
SoVRMLDragSensor::getNormalizedLocaterPosition(void) const
{
  return this->normpos;
}

#endif // HAVE_VRML97
