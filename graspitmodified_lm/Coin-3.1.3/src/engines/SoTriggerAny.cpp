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
  \class SoTriggerAny SoTriggerAny.h Inventor/engines/SoTriggerAny.h
  \brief The SoTriggerAny class is a fan-in engine for triggers.
  \ingroup engines

  When any one of the input triggers are "pulsed", any field connected
  as a slave to the engine output will be notified.
*/

#include <Inventor/engines/SoTriggerAny.h>
#include <Inventor/lists/SoEngineOutputList.h>

#include "engines/SoSubEngineP.h"

/*!
  \var SoSFTrigger SoTriggerAny::input0
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input1
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input2
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input3
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input4
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input5
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input6
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input7
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input8
  Input trigger.
*/
/*!
  \var SoSFTrigger SoTriggerAny::input9
  Input trigger.
*/

/*!
  \var SoEngineOutput SoTriggerAny::output

  (SoSFTrigger) Connect to the output with the field(s) you want
  notified upon any input trigger "pulses".
*/


SO_ENGINE_SOURCE(SoTriggerAny);

// Documented in superclass.
void
SoTriggerAny::initClass(void)
{
  SO_ENGINE_INTERNAL_INIT_CLASS(SoTriggerAny);
}

/*!
  Default constructor.
*/
SoTriggerAny::SoTriggerAny(void)
{
  SO_ENGINE_INTERNAL_CONSTRUCTOR(SoTriggerAny);

  SO_ENGINE_ADD_INPUT(input0, ());
  SO_ENGINE_ADD_INPUT(input1, ());
  SO_ENGINE_ADD_INPUT(input2, ());
  SO_ENGINE_ADD_INPUT(input3, ());
  SO_ENGINE_ADD_INPUT(input4, ());
  SO_ENGINE_ADD_INPUT(input5, ());
  SO_ENGINE_ADD_INPUT(input6, ());
  SO_ENGINE_ADD_INPUT(input7, ());
  SO_ENGINE_ADD_INPUT(input8, ());
  SO_ENGINE_ADD_INPUT(input9, ());

  SO_ENGINE_ADD_OUTPUT(output, SoSFTrigger);
}

/*!
  Destructor is protected because explicit destruction of engines is
  not allowed.
*/
SoTriggerAny::~SoTriggerAny()
{
}

// Documented in superclass.
void
SoTriggerAny::evaluate(void)
{
  SO_ENGINE_OUTPUT(output, SoSFTrigger, setValue());
}
