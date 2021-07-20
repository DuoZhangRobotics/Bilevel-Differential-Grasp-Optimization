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
  \class SoComposeRotation SoCompose.h Inventor/engines/SoCompose.h
  \brief The SoComposeRotation class is used to compose rotations from angle and axis.
  \ingroup engines

  Simple usage example:

  \code
  #Inventor V2.1 ascii
  
  Separator {
     Transform {
        rotation =
        ComposeRotation { axis 0 1 0  angle =
           ElapsedTime { }.timeOut
        }.rotation
     }
     Cube { }
  }
  \endcode
*/

#include <Inventor/engines/SoComposeRotation.h>
#include <Inventor/lists/SoEngineOutputList.h>
#include <Inventor/fields/SoMFRotation.h>

#include "engines/SoSubEngineP.h"

// *************************************************************************

/*!
  \var SoMFVec3f SoComposeRotation::axis
  Set of axis vectors for the output rotations. Default value is (0.0f, 0.0f, 1.0f).
*/
/*!
  \var SoMFFloat SoComposeRotation::angle
  Set of scalar rotation values for the output rotations. Default value is 0.0.
*/
/*!
  \var SoEngineOutput SoComposeRotation::rotation

  (SoMFRotation) Rotations generated from the angle and axis input
  fields.
*/

// *************************************************************************

SO_ENGINE_SOURCE(SoComposeRotation);

// *************************************************************************

SoComposeRotation::SoComposeRotation()
{
  SO_ENGINE_INTERNAL_CONSTRUCTOR(SoComposeRotation);

  SO_ENGINE_ADD_INPUT(axis,(0.0f,0.0f,1.0f));
  SO_ENGINE_ADD_INPUT(angle,(0.0f));

  SO_ENGINE_ADD_OUTPUT(rotation,SoMFRotation);
}

// Documented in superclass.
void
SoComposeRotation::initClass()
{
  SO_ENGINE_INTERNAL_INIT_CLASS(SoComposeRotation);
}

//
// private members
//
SoComposeRotation::~SoComposeRotation()
{
}

// *************************************************************************

// Documented in superclass.
void
SoComposeRotation::evaluate()
{
  int numAxis=axis.getNum();
  int numAngle=angle.getNum();

  int numOut=numAxis>numAngle?numAxis:numAngle;

  SO_ENGINE_OUTPUT(rotation,SoMFRotation,setNum(numOut));

  int i;

  float angleVal;
  for (i=0;i<numOut;i++) {
    const SbVec3f axisVal=i<numAxis?axis[i]:axis[numAxis-1];
    angleVal=i<numAngle?angle[i]:angle[numAngle-1];

    SO_ENGINE_OUTPUT(rotation,SoMFRotation,set1Value(i,axisVal,angleVal));
  }
}

// *************************************************************************
