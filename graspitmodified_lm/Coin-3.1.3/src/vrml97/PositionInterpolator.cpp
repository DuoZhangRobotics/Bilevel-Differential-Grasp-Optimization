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
  \class SoVRMLPositionInterpolator SoVRMLPositionInterpolator.h Inventor/VRMLnodes/SoVRMLPositionInterpolator.h
  \brief The SoVRMLPositionInterpolator class is used to interpolate 3D points.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  PositionInterpolator {
    eventIn      SFFloat set_fraction        # (-,)
    exposedField MFFloat key           []    # (-,)
    exposedField MFVec3f keyValue      []    # (-,)
    eventOut     SFVec3f value_changed
  }
  \endverbatim
  
  The PositionInterpolator node linearly interpolates among a list of
  3D vectors. The keyValue field shall contain exactly as many values
  as in the key field.  4.6.8, Interpolator nodes
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.8>),
  contains a more detailed discussion of interpolators.

*/

/*!
  \var SoMFVec3f SoVRMLPositionInterpolator::keyValue
  The keyValue vector.
*/

/*!
  \var SoEngineOutput SoVRMLPositionInterpolator::value_changed
  The eventOut which is sent every time the interpolator has calculated a new value.
*/

#include <Inventor/VRMLnodes/SoVRMLPositionInterpolator.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>

#include "engines/SoSubNodeEngineP.h"

SO_NODEENGINE_SOURCE(SoVRMLPositionInterpolator);

// Doc in parent
void
SoVRMLPositionInterpolator::initClass(void) // static
{
  SO_NODEENGINE_INTERNAL_INIT_CLASS(SoVRMLPositionInterpolator);
}

/*!
  Constructor.
*/
SoVRMLPositionInterpolator::SoVRMLPositionInterpolator(void)
{
  SO_NODEENGINE_INTERNAL_CONSTRUCTOR(SoVRMLPositionInterpolator);

  SO_VRMLNODE_ADD_EMPTY_EXPOSED_MFIELD(keyValue);
  SO_NODEENGINE_ADD_OUTPUT(value_changed, SoSFVec3f);
}

/*!
  Destructor.
*/
SoVRMLPositionInterpolator::~SoVRMLPositionInterpolator()
{
}

// Doc in parent
void
SoVRMLPositionInterpolator::evaluate(void)
{
  float interp;
  int idx = this->getKeyValueIndex(interp);
  if (idx < 0) return;

  const SbVec3f * v = this->keyValue.getValues(0);

  SbVec3f v0 = v[idx];
  if (interp > 0.0f) {
    SbVec3f v1 = v[idx+1];
    v0 = v0 + (v1-v0)*interp;
  }
  SO_ENGINE_OUTPUT(value_changed, SoSFVec3f, setValue(v0));
}

#endif // HAVE_VRML97
