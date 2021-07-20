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
  \class SoDecomposeMatrix SoCompose.h Inventor/engines/SoCompose.h
  \brief The SoDecomposeMatrix class is used to decompose a matrix into simple transformations.
  \ingroup engines
*/

#include <Inventor/engines/SoDecomposeMatrix.h>
#include <Inventor/lists/SoEngineOutputList.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoMFRotation.h>

#include "engines/SoSubEngineP.h"

SO_ENGINE_SOURCE(SoDecomposeMatrix);

/*!
  \var SoMFMatrix SoDecomposeMatrix::matrix
  Set of transformation matrices to decompose into their
  translation/rotation/scale parts.
*/
/*!
  \var SoMFVec3f SoDecomposeMatrix::center
  Center points of transform matrices.
*/
/*!
  \var SoEngineOutput SoDecomposeMatrix::translation
  (SoMFVec3f) Translation parts of input matrices.
*/
/*!
  \var SoEngineOutput SoDecomposeMatrix::rotation
  (SoMFRotation) Rotation parts of input matrices.
*/
/*!
  \var SoEngineOutput SoDecomposeMatrix::scaleFactor
  (SoMFVec3f) Scale vectors of input matrices.
*/
/*!
  \var SoEngineOutput SoDecomposeMatrix::scaleOrientation
  (SoMFRotation) Scale orientation values of the input matrices.
*/


#ifndef DOXYGEN_SKIP_THIS // No need to document these.

// Default constructor.
SoDecomposeMatrix::SoDecomposeMatrix()
{
  SO_ENGINE_INTERNAL_CONSTRUCTOR(SoDecomposeMatrix);

  SO_ENGINE_ADD_INPUT(matrix,(SbMatrix()));
  SO_ENGINE_ADD_INPUT(center,(SbVec3f()));

  SO_ENGINE_ADD_OUTPUT(translation,SoMFVec3f);
  SO_ENGINE_ADD_OUTPUT(rotation,SoMFRotation);
  SO_ENGINE_ADD_OUTPUT(scaleFactor,SoMFVec3f);
  SO_ENGINE_ADD_OUTPUT(scaleOrientation,SoMFRotation);
}

// Documented in superclass.
void
SoDecomposeMatrix::initClass()
{
  SO_ENGINE_INTERNAL_INIT_CLASS(SoDecomposeMatrix);
}

//
// private members
//
SoDecomposeMatrix::~SoDecomposeMatrix()
{
}

// Documented in superclass.
void
SoDecomposeMatrix::evaluate()
{
  int num = this->matrix.getNum();

  SO_ENGINE_OUTPUT(translation,SoMFVec3f,setNum(num));
  SO_ENGINE_OUTPUT(rotation,SoMFRotation,setNum(num));
  SO_ENGINE_OUTPUT(scaleFactor,SoMFVec3f,setNum(num));
  SO_ENGINE_OUTPUT(scaleOrientation,SoMFRotation,setNum(num));

  int i;
  SbVec3f translationVal,scaleFactorVal;
  SbRotation rotationVal,scaleOrientationVal;
  for (i = 0; i < num; i++) {
    SbVec3f c = (i < center.getNum()) ? center[i] : SbVec3f(0.0f, 0.0f, 0.0f);
    this->matrix[i].getTransform(translationVal,rotationVal,scaleFactorVal,
                                 scaleOrientationVal, c);
    SO_ENGINE_OUTPUT(translation,SoMFVec3f,set1Value(i,translationVal));
    SO_ENGINE_OUTPUT(rotation,SoMFRotation,set1Value(i,rotationVal));
    SO_ENGINE_OUTPUT(scaleFactor,SoMFVec3f,set1Value(i,scaleFactorVal));
    SO_ENGINE_OUTPUT(scaleOrientation,SoMFRotation,
                     set1Value(i,scaleOrientationVal));
  }
}

#endif // !DOXYGEN_SKIP_THIS
