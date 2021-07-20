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
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_DRAGGERS

/*!
  \class SoScale1Dragger SoScale1Dragger.h Inventor/draggers/SoScale1Dragger.h
  \brief The SoScale1Dragger class provides a mechanism for the end-user to scale in one dimension.
  \ingroup draggers

  \DRAGGER_DEFAULT_SCREENSHOT

  <center>
  <img src="http://doc.coin3d.org/images/Coin/draggers/scale1.png">
  </center>

  Use this dragger to allow the end-user of your application to scale
  along the X-axis. (Use a transformation node in front of the dragger
  to position it and re-orient it to scale along any vector.)

  \sa SoScale2Dragger, SoScaleUniformDragger, SoScale2UniformDragger
*/

#include <Inventor/draggers/SoScale1Dragger.h>

#include <cstring>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/SbRotation.h>
#include <Inventor/projectors/SbLineProjector.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <data/draggerDefaults/scale1Dragger.h>

#include "nodekits/SoSubKitP.h"
#include "coindefs.h"

/*!
  \var SoSFVec3f SoScale1Dragger::scaleFactor

  Continuously updated to contain the current vector of scaling along
  the X, Y and Z axes.

  For the SoScale1Dragger, only the X component is used, the Y and Z
  components will always be equal to 1 (ie no scaling).
*/

/*!
  \var SoFieldSensor * SoScale1Dragger::fieldSensor
  \COININTERNAL
*/
/*!
  \var SbLineProjector * SoScale1Dragger::lineProj
  \COININTERNAL
*/

#define THISP(d) static_cast<SoScale1Dragger *>(d)

class SoScale1DraggerP {
public:
};

SO_KIT_SOURCE(SoScale1Dragger);


// doc in superclass
void
SoScale1Dragger::initClass(void)
{
  SO_KIT_INTERNAL_INIT_CLASS(SoScale1Dragger, SO_FROM_INVENTOR_1);
}

// FIXME: document which parts need to be present in the geometry
// scenegraph, and what role they play in the dragger. 20010913 mortene.
/*!
  \DRAGGER_CONSTRUCTOR

  \NODEKIT_PRE_DIAGRAM

  \verbatim
  CLASS SoScale1Dragger
  -->"this"
        "callbackList"
        "topSeparator"
           "motionMatrix"
           "geomSeparator"
  -->         "scalerSwitch"
  -->            "scaler"
  -->            "scalerActive"
  -->         "feedbackSwitch"
  -->            "feedback"
  -->            "feedbackActive"
  \endverbatim

  \NODEKIT_POST_DIAGRAM


  \NODEKIT_PRE_TABLE

  \verbatim
  CLASS SoScale1Dragger
  PVT   "this",  SoScale1Dragger  ---
        "callbackList",  SoNodeKitListPart [ SoCallback, SoEventCallback ]
  PVT   "topSeparator",  SoSeparator  ---
  PVT   "motionMatrix",  SoMatrixTransform  ---
  PVT   "geomSeparator",  SoSeparator  ---
  PVT   "scalerSwitch",  SoSwitch  ---
        "scaler",  SoSeparator  ---
        "scalerActive",  SoSeparator  ---
  PVT   "feedbackSwitch",  SoSwitch  ---
        "feedback",  SoSeparator  ---
        "feedbackActive",  SoSeparator  ---
  \endverbatim

  \NODEKIT_POST_TABLE
*/
SoScale1Dragger::SoScale1Dragger(void)
{
  SO_KIT_INTERNAL_CONSTRUCTOR(SoScale1Dragger);

  SO_KIT_ADD_CATALOG_ENTRY(scalerSwitch, SoSwitch, TRUE, geomSeparator, feedbackSwitch, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(scaler, SoSeparator, TRUE, scalerSwitch, scalerActive, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(scalerActive, SoSeparator, TRUE, scalerSwitch, "", TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(feedbackSwitch, SoSwitch, TRUE, geomSeparator, "", FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(feedback, SoSeparator, TRUE, feedbackSwitch, feedbackActive, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(feedbackActive, SoSeparator, TRUE, feedbackSwitch, "", TRUE);

  if (SO_KIT_IS_FIRST_INSTANCE()) {
    SoInteractionKit::readDefaultParts("scale1Dragger.iv",
                                       SCALE1DRAGGER_draggergeometry,
                                       static_cast<int>(strlen(SCALE1DRAGGER_draggergeometry)));
  }

  SO_KIT_ADD_FIELD(scaleFactor, (1.0f, 1.0f, 1.0f));

  SO_KIT_INIT_INSTANCE();

  // initialize default parts
  this->setPartAsDefault("scaler", "scale1Scaler");
  this->setPartAsDefault("scalerActive", "scale1ScalerActive");
  this->setPartAsDefault("feedback", "scale1Feedback");
  this->setPartAsDefault("feedbackActive", "scale1FeedbackActive");

  // initialize swich values
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);

  // setup projector
  this->lineProj = new SbLineProjector();
  this->addStartCallback(SoScale1Dragger::startCB);
  this->addMotionCallback(SoScale1Dragger::motionCB);
  this->addFinishCallback(SoScale1Dragger::finishCB);

  this->addValueChangedCallback(SoScale1Dragger::valueChangedCB);

  this->fieldSensor = new SoFieldSensor(SoScale1Dragger::fieldSensorCB, this);
  this->fieldSensor->setPriority(0);

  this->setUpConnections(TRUE, TRUE);
}

/*!
  Protected destructor.

  (Dragger classes are derived from SoBase, so they are reference
  counted and automatically destroyed when their reference count goes
  to 0.)
 */
SoScale1Dragger::~SoScale1Dragger()
{
  delete this->lineProj;
  delete this->fieldSensor;
}

// Doc in superclass.
SbBool
SoScale1Dragger::setUpConnections(SbBool onoff, SbBool doitalways)
{
  if (!doitalways && this->connectionsSetUp == onoff) return onoff;

  SbBool oldval = this->connectionsSetUp;

  if (onoff) {
    inherited::setUpConnections(onoff, doitalways);

    SoScale1Dragger::fieldSensorCB(this, NULL);

    if (this->fieldSensor->getAttachedField() != &this->scaleFactor) {
      this->fieldSensor->attach(&this->scaleFactor);
    }
  }
  else {
    if (this->fieldSensor->getAttachedField() != NULL) {
      this->fieldSensor->detach();
    }
    inherited::setUpConnections(onoff, doitalways);
  }
  this->connectionsSetUp = onoff;
  return oldval;
}

/*! \COININTERNAL */
void
SoScale1Dragger::fieldSensorCB(void * d, SoSensor *)
{
  assert(d);
  SoScale1Dragger * thisp = THISP(d);
  SbMatrix matrix = thisp->getMotionMatrix();

  SbVec3f t, s;
  SbRotation r, so;

  matrix.getTransform(t, r, s, so);
  s = thisp->scaleFactor.getValue();
  matrix.setTransform(t, r, s, so);
  thisp->setMotionMatrix(matrix);
}

/*! \COININTERNAL */
void
SoScale1Dragger::valueChangedCB(void * COIN_UNUSED_ARG(f), SoDragger * d)
{
  SoScale1Dragger * thisp = THISP(d);
  SbMatrix matrix = thisp->getMotionMatrix();

  SbVec3f trans, scale;
  SbRotation rot, scaleOrient;
  matrix.getTransform(trans, rot, scale, scaleOrient);
  thisp->fieldSensor->detach();
  if (thisp->scaleFactor.getValue() != scale)
    thisp->scaleFactor = scale;
  thisp->fieldSensor->attach(&thisp->scaleFactor);
}

/*! \COININTERNAL */
void
SoScale1Dragger::startCB(void * COIN_UNUSED_ARG(f), SoDragger * d)
{
  SoScale1Dragger * thisp = THISP(d);
  thisp->dragStart();
}

/*! \COININTERNAL */
void
SoScale1Dragger::motionCB(void * COIN_UNUSED_ARG(f), SoDragger * d)
{
  SoScale1Dragger * thisp = THISP(d);
  thisp->drag();
}

/*! \COININTERNAL */
void
SoScale1Dragger::finishCB(void * COIN_UNUSED_ARG(f), SoDragger * d)
{
  SoScale1Dragger * thisp = THISP(d);
  thisp->dragFinish();
}

/*! \COININTERNAL
  Called when dragger is selected (picked) by the user.
*/
void
SoScale1Dragger::dragStart(void)
{
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 1);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 1);

  SbVec3f hitPt = this->getLocalStartingPoint();
  this->lineProj->setLine(SbLine(hitPt, hitPt + SbVec3f(1.0f, 0.0f, 0.0f)));
}

/*! \COININTERNAL
  Called when user drags the mouse after picking the dragger.
*/
void
SoScale1Dragger::drag(void)
{
  this->lineProj->setViewVolume(this->getViewVolume());
  this->lineProj->setWorkingSpace(this->getLocalToWorldMatrix());

  SbVec3f projPt = lineProj->project(this->getNormalizedLocaterPosition());
  SbVec3f startPt = this->getLocalStartingPoint();

  float motion = projPt[0];
  if (startPt[0] != 0.0f)
    motion /= startPt[0];
  else
    motion = 0.0f;

  if (motion < 0.0f) motion = 0.0f;

  this->setMotionMatrix(this->appendScale(this->getStartMotionMatrix(),
                                          SbVec3f(static_cast<float>(fabs(motion)), 1.0f, 1.0f),
                                          SbVec3f(0.0f, 0.0f, 0.0f)));
}

/*! \COININTERNAL
  Called when mouse button is released after picking and interacting
  with the dragger.
*/
void
SoScale1Dragger::dragFinish(void)
{
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
}

#undef THISP
#endif // HAVE_DRAGGERS
