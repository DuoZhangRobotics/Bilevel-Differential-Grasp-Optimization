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
  \class SoShadowSpotLight SoShadowSpotLight.h Inventor/annex/FXViz/nodes/SoShadowSpotLight.h
  \brief The SoShadowSpotLight class is a node for setting up a spot light which casts shadows.

  This node can be used instead of a normal SpotLight if you need to
  improve the performance by supplying a simplified scene graph to be
  used when rendering the shadow map(s). For instance, the shadow map
  scene graph doesn't need any textures or materials, and any
  non-casters can also be excluded from this scene graph. It's more
  optimal to use this node than to use the SoShadowStyle node to
  control this, at the cost of some extra application complexity.

  It's especially useful if you have a scene with few shadow caster
  nodes and lots of shadow receiver nodes.

  Currently, this node must be placed somewhere in the SoShadowGroup
  subgraph to cast shadows.

  \ingroup nodes

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    ShadowSpotLight {
      shadowMapScene NULL
      nearDistance -1
      farDistance -1
    }
  \endcode

  Here is the the example from SoShadowGroup, modified to use SoShadowSpotLight
  instead of a normal SoSpotLight. Notice that only the sphere casts shadows.

  \code

  DirectionalLight { direction 0 0 -1 intensity 0.2 }

  ShadowGroup {
    quality 1 # to get per pixel lighting

    ShadowSpotLight {
      location -8 -8 8.0
      direction 1 1 -1
      cutOffAngle 0.35
      dropOffRate 0.7

      shadowMapScene
      DEF sphere Separator {
          Complexity { value 1.0 }
          Material { diffuseColor 1 1 0 specularColor 1 1 1 shininess 0.9 }
          Shuttle { translation0 -3 1 0 translation1 3 -5 0 speed 0.25 on TRUE }
          Translation { translation -5 0 2 }
          Sphere { radius 2.0 }
        }
    }
    # need to insert the sphere in the regular scene graph as well
    USE sphere

    Separator {
      Material { diffuseColor 1 0 0 specularColor 1 1 1 shininess 0.9 }
      Shuttle { translation0 0 -5 0 translation1 0 5 0 speed 0.15 on TRUE }
      Translation { translation 0 0 -3 }
      Cube { depth 1.8 }
    }
    Separator {
      Material { diffuseColor 0 1 0 specularColor 1 1 1 shininess 0.9 }
      Shuttle { translation0 -5 0 0 translation1 5 0 0 speed 0.3 on TRUE }
      Translation { translation 0 0 -3 }
      Cube { }
    }

    Coordinate3 { point [ -10 -10 -3, 10 -10 -3, 10 10 -3, -10 10 -3 ] }
    Material { specularColor 1 1 1 shininess 0.9 }

    Complexity { textureQuality 0.1 }
    Texture2 { image 2 2 3 0xffffff 0x225588 0x225588 0xffffff }
    Texture2Transform { scaleFactor 4 4 }
    FaceSet { numVertices 4 }
  }

  \endcode

  \since Coin 3.0
*/

/*!
  \var SoSFNode SoShadowGroup::shadowMapScene
  
  The shadow map scene graph. If this is NULL (the default), the node
  will behave as a normal SoSpotLight node.  

*/

/*!
  \var SoSFFloat SoShadowGroup::nearDistance

  Can be used to set a fixed near distance for this spot light. The value in this
  field will be used if it's set to > 0.0. Default value is -1.0.
*/

/*!
  \var SoSFFloat SoShadowGroup::farDistance

  Can be used to set a fixed far distance for this spot light. The value in this
  field will be used if it's set to > 0.0. Default value is -1.0.
*/



// *************************************************************************

#include <Inventor/annex/FXViz/nodes/SoShadowSpotLight.h>

#include <stdio.h>
#include <Inventor/actions/SoGLRenderAction.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************


SO_NODE_SOURCE(SoShadowSpotLight);

/*!
  Constructor.
*/
SoShadowSpotLight::SoShadowSpotLight(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShadowSpotLight);
  SO_NODE_ADD_FIELD(shadowMapScene, (NULL));
  SO_NODE_ADD_FIELD(nearDistance, (-1.0f));
  SO_NODE_ADD_FIELD(farDistance, (-1.0f));
}

/*!
  Destructor.
*/
SoShadowSpotLight::~SoShadowSpotLight()
{
}

// Doc from superclass.
void
SoShadowSpotLight::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShadowSpotLight, SO_FROM_COIN_3_0);
}

// Doc from superclass.
void
SoShadowSpotLight::GLRender(SoGLRenderAction * action)
{
  inherited::GLRender(action);
}



#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoShadowSpotLight * node = new SoShadowSpotLight;
  assert(node);
  node->ref();
  BOOST_CHECK_MESSAGE(node->getTypeId() != SoType::badType(),
                      "missing class initialization");
  node->unref();
}

#endif // COIN_TEST_SUITE
