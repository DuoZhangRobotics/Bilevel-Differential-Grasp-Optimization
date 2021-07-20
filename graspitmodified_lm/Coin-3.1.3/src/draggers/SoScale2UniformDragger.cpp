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
  \class SoScale2UniformDragger SoScale2UniformDragger.h Inventor/draggers/SoScale2UniformDragger.h
  \brief The SoScale2UniformDragger class provides a mechanism for the end-user to scale in two dimensions.
  \ingroup draggers

  \DRAGGER_DEFAULT_SCREENSHOT

  <center>
  <img src="http://doc.coin3d.org/images/Coin/draggers/scale2uniform.png">
  </center>

  Use this dragger to allow the end-user of your application to scale
  along the X-axis and the Y-axis. (Use a transformation node in front
  of the dragger to position it and re-orient it to scale in any
  plane.)

  Scaling with this dragger can only be done in a uniform manner, ie
  the X component of the SoScale2UniformDragger::scaleFactor will
  always equal the Y component.

  For non-uniform scaling operations in 2 dimensions, use the
  SoScale2Dragger.

  \sa SoScaleUniformDragger
*/

#include <Inventor/draggers/SoScale2UniformDragger.h>

#include <cstring>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include <Inventor/SbRotation.h>
#include <Inventor/projectors/SbLineProjector.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <data/draggerDefaults/scale2UniformDragger.h>

#include "nodekits/SoSubKitP.h"

/*!
  \var SoSFVec3f SoScale2UniformDragger::scaleFactor

  Continuously updated to contain the current vector of scaling along
  the X, Y and Z axes.

  For the SoScale2UniformDragger, only the X and Y components are
  used, the Z component will always be equal to 1 (ie no scaling).
*/

/*!
  \var SoFieldSensor * SoScale2UniformDragger::fieldSensor
  \COININTERNAL
*/
/*!
  \var SbLineProjector * SoScale2UniformDragger::lineProj
  \COININTERNAL
*/

#define THISP(d) static_cast<SoScale2UniformDragger *>(d)

class SoScale2UniformDraggerP {
public:
};

SO_KIT_SOURCE(SoScale2UniformDragger);


// doc in superclass
void
SoScale2UniformDragger::initClass(void)
{
  SO_KIT_INTERNAL_INIT_CLASS(SoScale2UniformDragger, SO_FROM_INVENTOR_1);
}

// FIXME: document which parts need to be present in the geometry
// scenegraph, and what role they play in the dragger. 20010913 mortene.
/*!
  \DRAGGER_CONSTRUCTOR

  \NODEKIT_PRE_DIAGRAM

  \verbatim
  CLASS SoScale2UniformDragger
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
  CLASS SoScale2UniformDragger
  PVT   "this",  SoScale2UniformDragger  --- 
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
SoScale2UniformDragger::SoScale2UniformDragger(void)
{
  SO_KIT_INTERNAL_CONSTRUCTOR(SoScale2UniformDragger);

  SO_KIT_ADD_CATALOG_ENTRY(scalerSwitch, SoSwitch, TRUE, geomSeparator, feedbackSwitch, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(scaler, SoSeparator, TRUE, scalerSwitch, scalerActive, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(scalerActive, SoSeparator, TRUE, scalerSwitch, "", TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(feedbackSwitch, SoSwitch, TRUE, geomSeparator, "", FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(feedback, SoSeparator, TRUE, feedbackSwitch, feedbackActive, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(feedbackActive, SoSeparator, TRUE, feedbackSwitch, "", TRUE);

  if (SO_KIT_IS_FIRST_INSTANCE()) {
    SoInteractionKit::readDefaultParts("scale2UniformDragger.iv",
                                       SCALE2UNIFORMDRAGGER_draggergeometry,
                                       static_cast<int>(strlen(SCALE2UNIFORMDRAGGER_draggergeometry)));
  }

  SO_KIT_ADD_FIELD(scaleFactor, (1.0f, 1.0f, 1.0f));
  SO_KIT_INIT_INSTANCE();

  // initialize default parts
  this->setPartAsDefault("scaler", "scale2UniformScaler");
  this->setPartAsDefault("scalerActive", "scale2UniformScalerActive");
  this->setPartAsDefault("feedback", "scale2UniformFeedback");
  this->setPartAsDefault("feedbackActive", "scale2UniformFeedbackActive");

  // initialize swich values
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);

  // setup projector
  this->lineProj = new SbLineProjector();

  this->addStartCallback(SoScale2UniformDragger::startCB);
  this->addMotionCallback(SoScale2UniformDragger::motionCB);
  this->addFinishCallback(SoScale2UniformDragger::finishCB);

  this->addValueChangedCallback(SoScale2UniformDragger::valueChangedCB);

  this->fieldSensor = new SoFieldSensor(SoScale2UniformDragger::fieldSensorCB, this);
  this->fieldSensor->setPriority(0);

  this->setUpConnections(TRUE, TRUE);
}

/*!
  Protected destructor.

  (Dragger classes are derived from SoBase, so they are reference
  counted and automatically destroyed when their reference count goes
  to 0.)
 */
SoScale2UniformDragger::~SoScale2UniformDragger()
{
  delete this->lineProj;
  delete this->fieldSensor;
}

// Doc in superclass.
SbBool
SoScale2UniformDragger::setUpConnections(SbBool onoff, SbBool doitalways)
{
  if (!doitalways && this->connectionsSetUp == onoff) return onoff;

  SbBool oldval = this->connectionsSetUp;

  if (onoff) {
    inherited::setUpConnections(onoff, doitalways);

    SoScale2UniformDragger::fieldSensorCB(this, NULL);

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
SoScale2UniformDragger::fieldSensorCB(void * d, SoSensor *)
{
  assert(d);
  SoScale2UniformDragger * thisp = THISP(d);
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
SoScale2UniformDragger::valueChangedCB(void *, SoDragger * d)
{
  SoScale2UniformDragger * thisp = THISP(d);
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
SoScale2UniformDragger::startCB(void *, SoDragger * d)
{
  SoScale2UniformDragger * thisp = THISP(d);
  thisp->dragStart();
}

/*! \COININTERNAL */
void
SoScale2UniformDragger::motionCB(void *, SoDragger * d)
{
  SoScale2UniformDragger * thisp = THISP(d);
  thisp->drag();
}

/*! \COININTERNAL */
void
SoScale2UniformDragger::finishCB(void *, SoDragger * d)
{
  SoScale2UniformDragger * thisp = THISP(d);
  thisp->dragFinish();
}

/*! \COININTERNAL
  Called when dragger is selected (picked) by the user.
*/
void
SoScale2UniformDragger::dragStart(void)
{
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 1);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 1);

  SbVec3f startPt = this->getLocalStartingPoint();
  startPt[2] = 0.0f;
  this->lineProj->setLine(SbLine(SbVec3f(0.0f, 0.0f, 0.0f), startPt));
}

/*! \COININTERNAL
  Called when user drags the mouse after picking the dragger.
*/
void
SoScale2UniformDragger::drag(void)
{
  this->lineProj->setViewVolume(this->getViewVolume());
  this->lineProj->setWorkingSpace(this->getLocalToWorldMatrix());
  SbVec3f startPt = this->getLocalStartingPoint();
  SbVec3f projPt = lineProj->project(this->getNormalizedLocaterPosition());
  startPt[2] = 0.0f;
  projPt[2] = 0.0f;

  float orglen = startPt.length();
  float currlen = projPt.length();

  float scale = 0.0f;
  if (orglen > 0.0f) scale = currlen / orglen;

  if (scale > 0.0f && startPt.dot(projPt) < 0.0f) scale = 0.0f;

  this->setMotionMatrix(this->appendScale(this->getStartMotionMatrix(),
                                          SbVec3f(scale, scale, 1.0f),
                                          SbVec3f(0.0f, 0.0f, 0.0f)));
}

/*! \COININTERNAL
  Called when mouse button is released after picking and interacting
  with the dragger.
*/
void
SoScale2UniformDragger::dragFinish(void)
{
  SoSwitch *sw;
  sw = SO_GET_ANY_PART(this, "scalerSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
  sw = SO_GET_ANY_PART(this, "feedbackSwitch", SoSwitch);
  SoInteractionKit::setSwitchValue(sw, 0);
}

#undef THISP
#endif // HAVE_DRAGGERS
