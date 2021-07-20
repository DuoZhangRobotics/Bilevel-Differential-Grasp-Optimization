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
  \class SoComposeVec3f SoCompose.h Inventor/engines/SoCompose.h
  \brief The SoComposeVec3f class is used to compose 3D vectors from floats.
  \ingroup engines
*/

#include <Inventor/engines/SoComposeVec3f.h>
#include <Inventor/lists/SoEngineOutputList.h>

#include "engines/SoSubEngineP.h"

SO_ENGINE_SOURCE(SoComposeVec3f);

/*!
  \var SoMFFloat SoComposeVec3f::x
  First coordinates of the output vectors.
*/
/*!
  \var SoMFFloat SoComposeVec3f::y
  Second coordinates of the output vectors.
*/
/*!
  \var SoMFFloat SoComposeVec3f::z
  Third coordinates of the output vectors.
*/
/*!
  \var SoEngineOutput SoComposeVec3f::vector
  (SoMFVec3f) 3D vectors.
*/

#ifndef DOXYGEN_SKIP_THIS // No need to document these.

SoComposeVec3f::SoComposeVec3f()
{
  SO_ENGINE_INTERNAL_CONSTRUCTOR(SoComposeVec3f);

  SO_ENGINE_ADD_INPUT(x,(0.0f));
  SO_ENGINE_ADD_INPUT(y,(0.0f));
  SO_ENGINE_ADD_INPUT(z,(0.0f));

  SO_ENGINE_ADD_OUTPUT(vector,SoMFVec3f);
}

// Documented in superclass.
void
SoComposeVec3f::initClass()
{
  SO_ENGINE_INTERNAL_INIT_CLASS(SoComposeVec3f);
}

//
// private members
//
SoComposeVec3f::~SoComposeVec3f()
{
}

// Documented in superclass.
void
SoComposeVec3f::evaluate()
{
  int numX=x.getNum();
  int numY=y.getNum();
  int numZ=z.getNum();

  int numOut=numX>numY?numX:numY;
  numOut=numZ>numOut?numZ:numOut;

  SO_ENGINE_OUTPUT(vector,SoMFVec3f,setNum(numOut));

  int i;
  float xVal,yVal,zVal;
  for (i=0;i<numOut;i++) {
    xVal=i<numX?x[i]:x[numX-1];
    yVal=i<numY?y[i]:y[numY-1];
    zVal=i<numZ?z[i]:z[numZ-1];

    SO_ENGINE_OUTPUT(vector,SoMFVec3f,set1Value(i,xVal,yVal,zVal));
  }
}

#endif // !DOXYGEN_SKIP_THIS
