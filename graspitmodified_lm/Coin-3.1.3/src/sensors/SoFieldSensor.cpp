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
  \class SoFieldSensor SoFieldSensor.h Inventor/sensors/SoFieldSensor.h
  \brief The SoFieldSensor class detects changes to a field.
  \ingroup sensors

  Attach a field to a sensor of this type to put it under
  surveillance, so you can act upon changes to the field.

  An SoFieldSensor can also act for delete-callback purposes alone and
  does not need a regular notification-based callback.
*/

#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/fields/SoField.h>


/*!
  Constructor.
 */
SoFieldSensor::SoFieldSensor(void)
{
  this->convict = NULL;
}

/*!
  Constructor taking as parameters the sensor callback function and
  the userdata which will be passed the callback.

  \sa setFunction(), setData()
 */
SoFieldSensor::SoFieldSensor(SoSensorCB * func, void * data)
  : inherited(func, data)
{
  this->convict = NULL;
}

/*!
  Destructor.
*/
SoFieldSensor::~SoFieldSensor(void)
{
  if (this->convict) this->detach();
}

/*!
  Attach sensor to a field. Whenever the field's value changes, the
  sensor will be triggered and call the callback function.

  A field sensor can be attached to only a single field at a
  time. When this method is invoked multiple times, each subsequent
  call will replace the field the sensor is monitoring with the new \a
  field.

  When the given field is deleted, the sensor will automatically be
  detached.

  \sa detach()
 */
void
SoFieldSensor::attach(SoField * field)
{
  if (this->convict) this->detach();
  this->convict = field;
  field->addAuditor(this, SoNotRec::SENSOR);
  field->evaluate();
}

/*!
  Detach sensor from field. As long as an SoFieldSensor is detached,
  it will never call its callback function.

  \sa attach()
 */
void
SoFieldSensor::detach(void)
{
  if (this->convict) {
    this->convict->removeAuditor(this, SoNotRec::SENSOR);
    this->convict = NULL;
    if (this->isScheduled()) this->unschedule();
  }
}

/*!
  Returns a pointer to the field connected to the sensor.

  \sa attach(), detach()
 */
SoField *
SoFieldSensor::getAttachedField(void) const
{
  return this->convict;
}

// Doc from superclass.
void
SoFieldSensor::trigger(void)
{
  this->convict->evaluate();
  inherited::trigger();
}

// Doc from superclass.
void
SoFieldSensor::notify(SoNotList * l)
{
  // Overridden to only propagate if the field that caused the
  // notification is the one this sensor is attached to.
  if (l->getLastField() == this->convict) {
    inherited::notify(l);
  }
}

// Doc from superclass.
void
SoFieldSensor::dyingReference(void)
{
  SoFieldContainer * dyingcontainer = this->getAttachedField()->getContainer();

  this->invokeDeleteCallback();

  if (this->getAttachedField() != NULL &&
      this->getAttachedField()->getContainer() == dyingcontainer) {
    // sensor is attached, and to the same field-container
    this->detach();
    // FIXME: we could maybe do an exception for the globalfield-container,
    // and for loose fields I assume have NULL for getContainer() - those cases
    // should be checked at the field-pointer level instead.
  }
}
