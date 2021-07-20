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
  \class SoDecomposeVec4f SoCompose.h Inventor/engines/SoCompose.h
  \brief The SoDecomposeVec4f class is used to decompose 4D vectors into four floats.
  \ingroup engines
*/

#include <Inventor/engines/SoDecomposeVec4f.h>
#include <Inventor/lists/SoEngineOutputList.h>

#include "engines/SoSubEngineP.h"

SO_ENGINE_SOURCE(SoDecomposeVec4f);

/*!
  \var SoMFVec4f SoDecomposeVec4f::vector
  Set of input vectors to be decomposed into their coordinate
  elements.
*/
/*!
  \var SoEngineOutput SoDecomposeVec4f::x
  (SoMFFloat) First coordinates of the input vectors.
*/
/*!
  \var SoEngineOutput SoDecomposeVec4f::y
  (SoMFFloat) Second coordinates of the input vectors.
*/
/*!
  \var SoEngineOutput SoDecomposeVec4f::z
  (SoMFFloat) Third coordinates of the input vectors.
*/
/*!
  \var SoEngineOutput SoDecomposeVec4f::w
  (SoMFFloat) Fourth coordinates of the input vectors.
*/


#ifndef DOXYGEN_SKIP_THIS // No need to document these.

SoDecomposeVec4f::SoDecomposeVec4f()
{
  SO_ENGINE_INTERNAL_CONSTRUCTOR(SoDecomposeVec4f);

  SO_ENGINE_ADD_INPUT(vector,(0,0,0,0));

  SO_ENGINE_ADD_OUTPUT(x,SoMFFloat);
  SO_ENGINE_ADD_OUTPUT(y,SoMFFloat);
  SO_ENGINE_ADD_OUTPUT(z,SoMFFloat);
  SO_ENGINE_ADD_OUTPUT(w,SoMFFloat);
}

// Documented in superclass.
void
SoDecomposeVec4f::initClass()
{
  SO_ENGINE_INTERNAL_INIT_CLASS(SoDecomposeVec4f);
}

//
// private members
//
SoDecomposeVec4f::~SoDecomposeVec4f()
{
}

// Documented in superclass.
void
SoDecomposeVec4f::evaluate()
{
  int num = this->vector.getNum();

  SO_ENGINE_OUTPUT(x,SoMFFloat,setNum(num));
  SO_ENGINE_OUTPUT(y,SoMFFloat,setNum(num));
  SO_ENGINE_OUTPUT(z,SoMFFloat,setNum(num));
  SO_ENGINE_OUTPUT(w,SoMFFloat,setNum(num));

  int i;
  for (i = 0; i < num; i++) {
    SO_ENGINE_OUTPUT(x,SoMFFloat,set1Value(i,vector[i][0]));
    SO_ENGINE_OUTPUT(y,SoMFFloat,set1Value(i,vector[i][1]));
    SO_ENGINE_OUTPUT(z,SoMFFloat,set1Value(i,vector[i][2]));
    SO_ENGINE_OUTPUT(w,SoMFFloat,set1Value(i,vector[i][3]));
  }
}

#endif // !DOXYGEN_SKIP_THIS
