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
  \class SoVRMLVisibilitySensor SoVRMLVisibilitySensor.h Inventor/VRMLnodes/SoVRMLVisibilitySensor.h
  \brief The SoVRMLVisibilitySensor class will generate events based on visibility.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  VisibilitySensor {
    exposedField SFVec3f center   0 0 0      # (-,)
    exposedField SFBool  enabled  TRUE
    exposedField SFVec3f size     0 0 0      # [0,)
    eventOut     SFTime  enterTime
    eventOut     SFTime  exitTime
    eventOut     SFBool  isActive
  }
  \endverbatim

  The VisibilitySensor node detects visibility changes of a
  rectangular box as the user navigates the world. VisibilitySensor is
  typically used to detect when the user can see a specific object or
  region in the scene in order to activate or deactivate some
  behaviour or animation. The purpose is often to attract the
  attention of the user or to improve performance.

  The \e enabled field enables and disables the VisibilitySensor node.
  If enabled is set to FALSE, the VisibilitySensor node does not send
  events. If enabled is TRUE, the VisibilitySensor node detects
  changes to the visibility status of the box specified and sends
  events through the isActive eventOut. A TRUE event is output to
  isActive when any portion of the box impacts the rendered view. A
  FALSE event is sent when the box has no effect on the view. Browsers
  shall guarantee that, if isActive is FALSE, the box has absolutely
  no effect on the rendered view. Browsers may err liberally when
  isActive is TRUE. For example, the box may affect the rendering.

  The exposed fields \e center and \e size specify the object space
  location of the box centre and the extents of the box (i.e., width,
  height, and depth). The VisibilitySensor node's box is affected by
  hierarchical transformations of its parents. The components of the
  size field shall be greater than or equal to zero.

  The \e enterTime event is generated whenever the isActive TRUE event
  is generated, and exitTime events are generated whenever isActive
  FALSE events are generated. A VisibilitySensor read from a VRML file
  shall generate isActive TRUE and enterTime events if the sensor is
  enabled and the visibility box is visible. A VisibilitySensor
  inserted into the transformation hierarchy shall generate isActive
  TRUE and enterTime events if the sensor is enabled and the
  visibility box is visible. A VisibilitySensor removed from the
  transformation hierarchy shall generate isActive FALSE and exitTime
  events if the sensor is enabled and the visibility box is visible.

  Each VisibilitySensor node behaves independently of all other
  VisibilitySensor nodes. Every enabled VisibilitySensor node that is
  affected by the user's movement receives and sends events, possibly
  resulting in multiple VisibilitySensor nodes receiving and sending
  events simultaneously. Unlike TouchSensor nodes, there is no notion
  of a VisibilitySensor node lower in the scene graph "grabbing"
  events. Multiply instanced VisibilitySensor nodes (i.e., DEF/USE)
  use the union of all the boxes defined by their instances. An
  instanced VisibilitySensor node shall detect visibility changes for
  all instances of the box and send events appropriately.
  
*/

/*!
  \var SoSFVec3f SoVRMLVisibilitySensor::center
  Visibility area center. Default value is (0, 0, 0).
*/

/*!
  \var SoSFVec3f SoVRMLVisibilitySensor::size
  Visibility area size. Default value is (0, 0, 0).
*/

/*!
  \var SoSFBool SoVRMLVisibilitySensor::enabled
  Enable/disable sensor. Default value is TRUE.
*/

/*!
  \var SoSFTime SoVRMLVisibilitySensor::enterTime
  An event out that is triggered when the region becomes visible.
*/

/*!
  \var SoSFTime SoVRMLVisibilitySensor::exitTime
  An event out that is triggered when the region becomes not visible.
*/

/*!
  \var SoSFBool SoVRMLVisibilitySensor::isActive
  An event out that is generated when the visibility state changes.
*/

#include <Inventor/VRMLnodes/SoVRMLVisibilitySensor.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoCullElement.h>
#include <Inventor/SoDB.h>
#include <Inventor/SbBox3f.h>
#include <Inventor/fields/SoSFTime.h>

#include "nodes/SoSubNodeP.h"

SO_NODE_SOURCE(SoVRMLVisibilitySensor);

//
// returns the current time. First tries the realTime field, then
// SbTime::getTimeOfDay() if field is not found.
//
static SbTime
visibilitysensor_get_current_time(void)
{
  SoField * realtime = SoDB::getGlobalField("realTime");
  if (realtime && realtime->isOfType(SoSFTime::getClassTypeId())) {
    return ((SoSFTime*)realtime)->getValue();
  }
  return SbTime::getTimeOfDay();
}

// Doc in parent
void
SoVRMLVisibilitySensor::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLVisibilitySensor, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLVisibilitySensor::SoVRMLVisibilitySensor(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLVisibilitySensor);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(center, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(size, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(enabled, (TRUE));

  SO_VRMLNODE_ADD_EVENT_OUT(enterTime);
  SO_VRMLNODE_ADD_EVENT_OUT(exitTime);
  SO_VRMLNODE_ADD_EVENT_OUT(isActive);
  this->isActive = FALSE;
}

/*!
  Destructor.
*/
SoVRMLVisibilitySensor::~SoVRMLVisibilitySensor()
{
}

// Doc in parent
void
SoVRMLVisibilitySensor::GLRender(SoGLRenderAction * action)
{
  SbVec3f c = this->center.getValue();
  SbVec3f s = this->size.getValue();

  SbBool wasvisible = this->isActive.getValue();
  SbBool visible = FALSE;

  if (s != SbVec3f(0.0f, 0.0f, 0.0f)) {
    SbBox3f box(c[0]-s[0], c[1]-s[1], c[2]-s[2],
                c[0]+s[0], c[1]+s[1], c[2]+s[2]);
    if (!SoCullElement::cullTest(action->getState(), box, TRUE)) {
      // FIXME: the SoCullElement cull test only tests if box is outside
      // one of the planes, and the box might not be culled even if it's
      // not visible for some cases. pederb, 2002-05-16
      visible = TRUE;
    }
  }
  if (visible && !wasvisible) {
    this->enterTime = visibilitysensor_get_current_time();
    this->isActive = TRUE;
  }
  else if (!visible && wasvisible) {
    this->exitTime = visibilitysensor_get_current_time();
    this->isActive = FALSE;
  }
}

#endif // HAVE_VRML97
