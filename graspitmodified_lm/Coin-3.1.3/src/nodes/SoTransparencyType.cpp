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
  \class SoTransparencyType SoTransparencyType.h Inventor/nodes/SoTransparencyType.h
  \brief The SoTransparencyType class is a node for setting the transparency type for shapes.
  \ingroup nodes

  In earlier versions of Coin/Open Inventor it was only possible to
  set the transparency mode globally for an entire scene graph, which
  could be inconvenient if different transparency types was wanted for
  different shapes.

  Here is a screenshot of the different transparency modes used in a
  single scene.
 
  <center>
  <img src="http://doc.coin3d.org/images/Coin/nodes/transparencytype.png">
  </center>

  \COIN_CLASS_EXTENSION

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TransparencyType {
        value SCREEN_DOOR
    }
  \endcode

  \sa SoGLRenderAction::TransparencyType
  \since Coin 2.0
*/

// *************************************************************************

#include <Inventor/nodes/SoTransparencyType.h>

#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetPrimitiveCountAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/elements/SoOverrideElement.h>
#include <Inventor/elements/SoShapeStyleElement.h>
#include <Inventor/elements/SoLazyElement.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \enum SoTransparencyType::Type
  Enumeration of available transparency types. See documentation in
  SoGLRenderAction for a description of the different types.
*/

/*!
  \var SoSFEnum SoTransparencyType::value

  The transparency type to use for subsequent shape nodes in the scene
  graph.
*/


// *************************************************************************

SO_NODE_SOURCE(SoTransparencyType);

/*!
  Constructor.
*/
SoTransparencyType::SoTransparencyType(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTransparencyType);

  SO_NODE_ADD_FIELD(value, (SCREEN_DOOR));

  SO_NODE_DEFINE_ENUM_VALUE(Type, SCREEN_DOOR);
  SO_NODE_DEFINE_ENUM_VALUE(Type, ADD);
  SO_NODE_DEFINE_ENUM_VALUE(Type, DELAYED_ADD);
  SO_NODE_DEFINE_ENUM_VALUE(Type, BLEND);
  SO_NODE_DEFINE_ENUM_VALUE(Type, DELAYED_BLEND);
  SO_NODE_DEFINE_ENUM_VALUE(Type, SORTED_OBJECT_ADD);
  SO_NODE_DEFINE_ENUM_VALUE(Type, SORTED_OBJECT_BLEND);
  SO_NODE_DEFINE_ENUM_VALUE(Type, SORTED_OBJECT_SORTED_TRIANGLE_ADD);
  SO_NODE_DEFINE_ENUM_VALUE(Type, SORTED_OBJECT_SORTED_TRIANGLE_BLEND);
  SO_NODE_DEFINE_ENUM_VALUE(Type, NONE);

  SO_NODE_SET_SF_ENUM_TYPE(value, Type);
}


/*!
  Destructor.
*/
SoTransparencyType::~SoTransparencyType()
{
}

// Doc from superclass.
void
SoTransparencyType::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTransparencyType, SO_FROM_INVENTOR_1);
}

// Doc from superclass.
void
SoTransparencyType::GLRender(SoGLRenderAction * action)
{
  SoTransparencyType::doAction(action);
}

// Doc from superclass.
void
SoTransparencyType::doAction(SoAction * action)
{
  if (!this->value.isIgnored()
      && !SoOverrideElement::getTransparencyTypeOverride(action->getState())) {
    SoShapeStyleElement::setTransparencyType(action->getState(),
                                             this->value.getValue());
    SoLazyElement::setTransparencyType(action->getState(),
                                       this->value.getValue());
    if (this->isOverride()) {
      SoOverrideElement::setTransparencyTypeOverride(action->getState(), this, TRUE);
    }
  }
}

// Doc from superclass.
void
SoTransparencyType::callback(SoCallbackAction * action)
{
  SoTransparencyType::doAction((SoAction *)action);
}
