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
  \class SoPointLightDragger SoPointLightDragger.h Inventor/draggers/SoPointLightDragger.h
  \brief The SoPointLightDragger class provides interactive geometry for manipulating a point light source.
  \ingroup draggers

  \DRAGGER_DEFAULT_SCREENSHOT

  <center>
  <img src="http://doc.coin3d.org/images/Coin/draggers/pointlight.png">
  </center>

  This dragger is well suited to use for setting up the fields of a
  SoPointLight node, as it provides geometry for the end-user to
  translate a point in 3D space.

  The Coin library also includes a manipulator class,
  SoPointLightManip, which wraps the functionality provided by this
  class inside the necessary mechanisms for connecting it to
  SoPointLight node instances in a scenegraph.

  \sa SoPointLightManip
*/

#include <Inventor/draggers/SoPointLightDragger.h>

#include <cstring>

#include <Inventor/draggers/SoDragPointDragger.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <data/draggerDefaults/pointLightDragger.h>

#include "nodekits/SoSubKitP.h"
#include "SbBasicP.h"

/*!
  \var SoSFVec3f SoPointLightDragger::translation

  This field is continuously updated to contain the translation of the
  pointlight dragger. The application programmer will typically
  connect this to the SoPointLight::location field of a SoPointLight
  node (unless using the SoPointLightManip class, where this is taken
  care of automatically).

  It may also of course be connected to any other location /
  translation field controlling the position of scenegraph geometry,
  it does not have to part of a SoPointLight node specifically.
*/

class SoPointLightDraggerP {
public:
};

SO_KIT_SOURCE(SoPointLightDragger);

// doc in superclass
void
SoPointLightDragger::initClass(void)
{
  SO_KIT_INTERNAL_INIT_CLASS(SoPointLightDragger, SO_FROM_INVENTOR_1);
}

// FIXME: document which parts need to be present in the geometry
// scenegraph, and what role they play in the dragger. 20010913 mortene.
/*!
  \DRAGGER_CONSTRUCTOR

  \NODEKIT_PRE_DIAGRAM

  \verbatim
  CLASS SoPointLightDragger
  -->"this"
        "callbackList"
        "topSeparator"
           "motionMatrix"
  -->      "material"
  -->      "translator"
           "geomSeparator"
  \endverbatim

  \NODEKIT_POST_DIAGRAM


  \NODEKIT_PRE_TABLE

  \verbatim
  CLASS SoPointLightDragger
  PVT   "this",  SoPointLightDragger  --- 
        "callbackList",  SoNodeKitListPart [ SoCallback, SoEventCallback ] 
  PVT   "topSeparator",  SoSeparator  --- 
  PVT   "motionMatrix",  SoMatrixTransform  --- 
        "material",  SoMaterial  --- 
        "translator",  SoDragPointDragger  --- 
  PVT   "geomSeparator",  SoSeparator  --- 
  \endverbatim

  \NODEKIT_POST_TABLE
*/
SoPointLightDragger::SoPointLightDragger(void)
{
  SO_KIT_INTERNAL_CONSTRUCTOR(SoPointLightDragger);

  SO_KIT_ADD_CATALOG_ENTRY(material, SoMaterial, TRUE, topSeparator, translator, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(translator, SoDragPointDragger, TRUE, topSeparator, geomSeparator, TRUE);

  if (SO_KIT_IS_FIRST_INSTANCE()) {
    SoInteractionKit::readDefaultParts("pointLightDragger.iv",
                                       POINTLIGHTDRAGGER_draggergeometry,
                                       static_cast<int>(strlen(POINTLIGHTDRAGGER_draggergeometry)));
  }

  SO_KIT_ADD_FIELD(translation, (0.0f, 0.0f, 0.0f));
  SO_KIT_INIT_INSTANCE();

  SoDragger *pdragger = SO_GET_ANY_PART(this, "translator", SoDragPointDragger);
  assert(pdragger);

  this->setPartAsDefault("material", "pointLightOverallMaterial");

  this->addValueChangedCallback(SoPointLightDragger::valueChangedCB);
  this->fieldSensor = new SoFieldSensor(SoPointLightDragger::fieldSensorCB, this);
  this->fieldSensor->setPriority(0);
  this->setUpConnections(TRUE, TRUE);
}

/*!
  Protected destructor.

  (Dragger classes are derived from SoBase, so they are reference
  counted and automatically destroyed when their reference count goes
  to 0.)
 */
SoPointLightDragger::~SoPointLightDragger()
{
  delete this->fieldSensor;
}

// Doc in superclass.
SbBool
SoPointLightDragger::setUpConnections(SbBool onoff, SbBool doitalways)
{
  if (!doitalways && this->connectionsSetUp == onoff) return onoff;

  if (onoff) {
    inherited::setUpConnections(onoff, doitalways);
    SoDragger * child = coin_assert_cast<SoDragger *>(this->getAnyPart("translator", FALSE));
    child->setPartAsDefault("yzTranslator.translator",
                            "pointLightTranslatorPlaneTranslator");
    child->setPartAsDefault("xzTranslator.translator",
                            "pointLightTranslatorPlaneTranslator");
    child->setPartAsDefault("xyTranslator.translator",
                            "pointLightTranslatorPlaneTranslator");

    child->setPartAsDefault("yzTranslator.translatorActive",
                            "pointLightTranslatorPlaneTranslatorActive");
    child->setPartAsDefault("xzTranslator.translatorActive",
                            "pointLightTranslatorPlaneTranslatorActive");
    child->setPartAsDefault("xyTranslator.translatorActive",
                            "pointLightTranslatorPlaneTranslatorActive");

    child->setPartAsDefault("xTranslator.translator",
                            "pointLightTranslatorLineTranslator");
    child->setPartAsDefault("yTranslator.translator",
                            "pointLightTranslatorLineTranslator");
    child->setPartAsDefault("zTranslator.translator",
                            "pointLightTranslatorLineTranslator");

    child->setPartAsDefault("xTranslator.translatorActive",
                            "pointLightTranslatorLineTranslatorActive");
    child->setPartAsDefault("yTranslator.translatorActive",
                            "pointLightTranslatorLineTranslatorActive");
    child->setPartAsDefault("zTranslator.translatorActive",
                            "pointLightTranslatorLineTranslatorActive");
    this->registerChildDragger(child);
    if (this->fieldSensor->getAttachedField() != &this->translation) {
      this->fieldSensor->attach(&this->translation);
    }
  }
  else {
    SoDragger * child = coin_assert_cast<SoDragger *>(this->getAnyPart("translator", FALSE));
    this->unregisterChildDragger(child);
    if (this->fieldSensor->getAttachedField() != NULL) {
      this->fieldSensor->detach();
    }
    inherited::setUpConnections(onoff, doitalways);
  }
  return !(this->connectionsSetUp = onoff);
}

// Doc in superclass.
void
SoPointLightDragger::setDefaultOnNonWritingFields(void)
{
  this->translator.setDefault(TRUE);

  inherited::setDefaultOnNonWritingFields();
}

/*! \COININTERNAL */
void
SoPointLightDragger::fieldSensorCB(void *d, SoSensor *)
{
  SoPointLightDragger * thisp = static_cast<SoPointLightDragger *>(d);
  SbMatrix matrix = thisp->getMotionMatrix();
  thisp->workFieldsIntoTransform(matrix);
  thisp->setMotionMatrix(matrix);
}

/*! \COININTERNAL */
void
SoPointLightDragger::valueChangedCB(void *, SoDragger * d)
{
  SoPointLightDragger * thisp = static_cast<SoPointLightDragger *>(d);

  SbMatrix matrix = thisp->getMotionMatrix();
  SbVec3f t;
  t[0] = matrix[3][0];
  t[1] = matrix[3][1];
  t[2] = matrix[3][2];

  thisp->fieldSensor->detach();
  if (thisp->translation.getValue() != t) {
    thisp->translation = t;
  }
  thisp->fieldSensor->attach(&thisp->translation);
}

#endif // HAVE_DRAGGERS
