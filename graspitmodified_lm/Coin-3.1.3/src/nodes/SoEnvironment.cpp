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
  \class SoEnvironment SoEnvironment.h Inventor/nodes/SoEnvironment.h
  \brief The SoEnvironment class is a node for specifying global rendering parameters.
  \ingroup nodes

  This node type provides the application programmer with the ability
  to set global parameters influencing lighting and fog.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    Environment {
        ambientIntensity 0.2
        ambientColor 1 1 1
        attenuation 0 0 1
        fogType NONE
        fogColor 1 1 1
        fogVisibility 0
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoEnvironment.h>

#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLEnvironmentElement.h>
#include <Inventor/elements/SoLightAttenuationElement.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \enum SoEnvironment::FogType
  Enumeration of available types of fog.
*/
/*!
  \var SoEnvironment::FogType SoEnvironment::NONE

  No fog. Visibility will be equal for all objects, independent of
  distance to camera.
*/
/*!
  \var SoEnvironment::FogType SoEnvironment::HAZE

  Fog where visibility will decrease linearly with distance to camera.
*/
/*!
  \var SoEnvironment::FogType SoEnvironment::FOG

  Fog where visibility will decrease exponentially with distance to
  camera.
*/
/*!
  \var SoEnvironment::FogType SoEnvironment::SMOKE

  Fog where visibility will decrease exponentially with the square of
  the distance to camera (simulating really thick fog).
*/


/*!
  \var SoSFFloat SoEnvironment::ambientIntensity

  A global ambient value for the light intensity for the complete
  scene. This will provide some light even when there are no light
  sources defined for the scene.

  Valid values is from 0.0 (no ambient light) to 1.0 (full ambient
  light intensity). Default value is 0.2.
*/
/*!
  \var SoSFColor SoEnvironment::ambientColor

  The color of the global ambient light. Defaults to full intensity
  white, ie <1.0, 1.0, 1.0>.

  Each component of this RGB vector is multiplied with the scalar
  value of SoEnvironment::ambientIntensity.
*/
/*!
  \var SoSFVec3f SoEnvironment::attenuation

  Squared, linear, and constant light attenuation coefficients (in
  that order). Please note that this is the reverse order of how
  attenuation is usually specified.
*/

/*!
  \var SoSFEnum SoEnvironment::fogType

  The fog model. See SoEnvironment::FoType.
*/
/*!
  \var SoSFColor SoEnvironment::fogColor

  Color of fog. Defaults to full intensity white.
*/
/*!
  \var SoSFFloat SoEnvironment::fogVisibility

  The "cut-off" distance from the camera where objects will be totally
  obscured by fog. If set to 0.0, the far plane distance will be used
  instead.

  Default value is 0.0.
*/

// *************************************************************************

SO_NODE_SOURCE(SoEnvironment);

/*!
  Constructor.
*/
SoEnvironment::SoEnvironment()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoEnvironment);

  SO_NODE_ADD_FIELD(ambientIntensity, (0.2f));
  SO_NODE_ADD_FIELD(ambientColor, (1.0f, 1.0f, 1.0f));
  SO_NODE_ADD_FIELD(attenuation, (0.0f, 0.0f, 1.0f));
  SO_NODE_ADD_FIELD(fogType, (SoEnvironment::NONE));
  SO_NODE_ADD_FIELD(fogColor, (1.0f, 1.0f, 1.0f));
  SO_NODE_ADD_FIELD(fogVisibility, (0.0f));

  SO_NODE_DEFINE_ENUM_VALUE(FogType, NONE);
  SO_NODE_DEFINE_ENUM_VALUE(FogType, HAZE);
  SO_NODE_DEFINE_ENUM_VALUE(FogType, FOG);
  SO_NODE_DEFINE_ENUM_VALUE(FogType, SMOKE);
  SO_NODE_SET_SF_ENUM_TYPE(fogType, FogType);
}

/*!
  Destructor.
*/
SoEnvironment::~SoEnvironment()
{
}

// Doc from superclass.
void
SoEnvironment::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoEnvironment, SO_FROM_INVENTOR_2_0);

  SO_ENABLE(SoGLRenderAction, SoGLEnvironmentElement);
  SO_ENABLE(SoCallbackAction, SoEnvironmentElement);
  SO_ENABLE(SoGLRenderAction, SoLightAttenuationElement);
  SO_ENABLE(SoCallbackAction, SoLightAttenuationElement);
}

// Doc from superclass.
void
SoEnvironment::GLRender(SoGLRenderAction * action)
{
  SoLightAttenuationElement::set(action->getState(), this,
                                 this->attenuation.getValue());
  SoEnvironmentElement::set(action->getState(),
                            this,
                            this->ambientIntensity.getValue(),
                            this->ambientColor.getValue(),
                            this->attenuation.getValue(),
                            (int32_t)fogType.getValue(),
                            this->fogColor.getValue(),
                            this->fogVisibility.getValue());
}

// Doc from superclass.
void
SoEnvironment::callback(SoCallbackAction *action)
{
  SoLightAttenuationElement::set(action->getState(), this,
                                 this->attenuation.getValue());
  SoEnvironmentElement::set(action->getState(),
                            this,
                            this->ambientIntensity.getValue(),
                            this->ambientColor.getValue(),
                            this->attenuation.getValue(),
                            (int32_t)fogType.getValue(),
                            this->fogColor.getValue(),
                            this->fogVisibility.getValue());
}
