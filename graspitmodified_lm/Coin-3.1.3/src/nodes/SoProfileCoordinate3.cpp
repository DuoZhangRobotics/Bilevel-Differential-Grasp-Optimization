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
  \class SoProfileCoordinate3 SoProfileCoordinate3.h Inventor/nodes/SoProfileCoordinate3.h
  \brief The SoProfileCoordinate3 class is a node specifying a set of 3D coordinates for profiles.
  \ingroup nodes

  Use nodes of this type to provide coordinates to profiles.

  The third element of the coordinate vectors is used for
  normalization. A node of this type where all the normalization
  values are equal to 1.0 is the equivalent of setting up an
  SoProfileCoordinate2 node.

  Note that an SoProfileCoordinate3 node will \e replace the profile
  coordinates already present in the state (if any).

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    ProfileCoordinate3 {
        point 0 0 1
    }
  \endcode

  \sa SoProfile, SoProfileCoordinate2
*/

// *************************************************************************

#include <Inventor/nodes/SoProfileCoordinate3.h>

#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoGetPrimitiveCountAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/elements/SoProfileCoordinateElement.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoMFVec3f SoProfileCoordinate3::point

  Pool of coordinate points for the traversal state.
*/

// *************************************************************************

SO_NODE_SOURCE(SoProfileCoordinate3);

/*!
  Constructor.
*/
SoProfileCoordinate3::SoProfileCoordinate3(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoProfileCoordinate3);

  SO_NODE_ADD_FIELD(point, (0.0f, 0.0f, 1.0f));
}

/*!
  Destructor.
*/
SoProfileCoordinate3::~SoProfileCoordinate3()
{
}

// doc in super
void
SoProfileCoordinate3::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoProfileCoordinate3, SO_FROM_INVENTOR_1);

  SO_ENABLE(SoGLRenderAction, SoProfileCoordinateElement);
  SO_ENABLE(SoPickAction, SoProfileCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoProfileCoordinateElement);
  SO_ENABLE(SoGetBoundingBoxAction, SoProfileCoordinateElement);
  SO_ENABLE(SoGetPrimitiveCountAction, SoProfileCoordinateElement);
}

void
SoProfileCoordinate3::getBoundingBox(SoGetBoundingBoxAction * action)
{
  SoProfileCoordinate3::doAction(action);
}

void
SoProfileCoordinate3::doAction(SoAction * action)
{
  SoProfileCoordinateElement::set3(action->getState(), this,
                                   this->point.getNum(),
                                   this->point.getValues(0));
}

void
SoProfileCoordinate3::GLRender(SoGLRenderAction * action)
{
  SoProfileCoordinate3::doAction(action);
}

void
SoProfileCoordinate3::callback(SoCallbackAction * action)
{
  SoProfileCoordinate3::doAction(action);
}

void
SoProfileCoordinate3::pick(SoPickAction * action)
{
  SoProfileCoordinate3::doAction(action);
}

void
SoProfileCoordinate3::getPrimitiveCount(SoGetPrimitiveCountAction * action)
{
  SoProfileCoordinate3::doAction(action);
}
