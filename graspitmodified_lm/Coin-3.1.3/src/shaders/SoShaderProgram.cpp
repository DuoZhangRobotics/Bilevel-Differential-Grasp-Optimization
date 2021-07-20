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

  \class SoShaderProgram SoShaderProgram.h Inventor/nodes/SoShaderProgram.h

  \brief The SoShaderProgram class is used to specify a set of
  vertex/geometry/fragment objects.

  This node can store one of each of SoVertexShader, SoGeometryShader
  and SoFragmentShader in its shaderObject field. Coin will load all
  shader objects specified there, and attach all objects into a
  program before binding it as the current shader program.

  \ingroup shaders

  A typical scene graph with shaders will look something like this:

  \code

  Separator {
    ShaderProgram {
      shaderObject [
        VertexShader {
          sourceProgram "myvertexshader.glsl"
          parameter [
            ShaderParameter1f { name "myvertexparam" value 1.0 }
          ]
        }
        FragmentShader {
          sourceProgram "myfragmentshader.glsl"
          parameter [
            ShaderParameter1f { name "myfragmentparam" value 2.0 }
          ]
        }
      ]
    }
    Cube { }
  }

  \endcode

  This will render the Cube with the vertex and fragment shaders
  specified in myvertexshader.glsl and myfragmentshader.glsl. Coin
  also supports ARB shaders and Cg shaders (if the Cg library is
  installed). However, we recommend using GLSL since we will focus
  mostly on support this shader language.

  Coin defines some named parameters that can be added by the
  application programmer, and which will be automatically updated by
  Coin while traversing the scene graph.

  \li coin_texunit[n]_model - Set to 0 when texturing is disabled, and
  to SoTextureImageElement::Model if there's a current texture on the
  state for unit \a n.

  \li coin_light_model - Set to 1 for PHONG, 0 for BASE_COLOR lighting.

  Example scene graph that renders per-fragment OpenGL Phong lighting
  for one light source. The shaders assume the first light source is a
  directional light. This is the case if you open the file in a
  standard examiner viewer.

  The iv-file:
  \code
  Separator {
    ShaderProgram {
      shaderObject [
        VertexShader {
          sourceProgram "perpixel_vertex.glsl"
        }
        FragmentShader {
          sourceProgram "perpixel_fragment.glsl"
        }
      ]
    }
    Complexity { value 1.0 }
    Material { diffuseColor 1 0 0 specularColor 1 1 1 shininess 0.9 }
    Sphere { }

    Translation { translation 3 0 0 }
    Material { diffuseColor 0 1 0 specularColor 1 1 1 shininess 0.9 }
    Cone { }

    Translation { translation 3 0 0 }
    Material { diffuseColor 0.8 0.4 0.1 specularColor 1 1 1 shininess 0.9 }
    Cylinder { }
  }
  \endcode

  The vertex shader (perpixel_vertex.glsl):
  \code
  varying vec3 ecPosition3;
  varying vec3 fragmentNormal;

  void main(void)
  {
    vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
    ecPosition3 = ecPosition.xyz / ecPosition.w;
    fragmentNormal = normalize(gl_NormalMatrix * gl_Normal);

    gl_Position = ftransform();
    gl_FrontColor = gl_Color;
  }
  \endcode

  The fragment shader (perpixel_vertex.glsl):
  \code
  varying vec3 ecPosition3;
  varying vec3 fragmentNormal;

  void DirectionalLight(in int i,
                        in vec3 normal,
                        inout vec4 ambient,
                        inout vec4 diffuse,
                        inout vec4 specular)
  {
    float nDotVP; // normal . light direction
    float nDotHV; // normal . light half vector
    float pf;     // power factor

    nDotVP = max(0.0, dot(normal, normalize(vec3(gl_LightSource[i].position))));
    nDotHV = max(0.0, dot(normal, vec3(gl_LightSource[i].halfVector)));

    if (nDotVP == 0.0)
      pf = 0.0;
    else
      pf = pow(nDotHV, gl_FrontMaterial.shininess);

    ambient += gl_LightSource[i].ambient;
    diffuse += gl_LightSource[i].diffuse * nDotVP;
    specular += gl_LightSource[i].specular * pf;
  }

  void main(void)
  {
    vec3 eye = -normalize(ecPosition3);
    vec4 ambient = vec4(0.0);
    vec4 diffuse = vec4(0.0);
    vec4 specular = vec4(0.0);
    vec3 color;

    DirectionalLight(0, normalize(fragmentNormal), ambient, diffuse, specular);

    color =
      gl_FrontLightModelProduct.sceneColor.rgb +
      ambient.rgb * gl_FrontMaterial.ambient.rgb +
      diffuse.rgb * gl_Color.rgb +
      specular.rgb * gl_FrontMaterial.specular.rgb;

    gl_FragColor = vec4(color, gl_Color.a);
  }
  \endcode

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    ShaderProgram {
      shaderObject []
    }
  \endcode

  \sa SoShaderObject
  \sa SoShaderProgram
  \since Coin 2.5
*/

/*!
  \var SoMFNode SoShaderProgram::shaderObject

  The shader objects.

*/

#include <Inventor/nodes/SoShaderProgram.h>

#include <assert.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLShaderProgramElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/nodes/SoShaderObject.h>
#include <Inventor/sensors/SoNodeSensor.h>

#include "nodes/SoSubNodeP.h"
#include "shaders/SoGLShaderProgram.h"

// *************************************************************************

class SoShaderProgramP
{
public:
  SoShaderProgramP(SoShaderProgram * ownerptr);
  ~SoShaderProgramP();

  void GLRender(SoGLRenderAction * action);

  SoShaderProgramEnableCB * enablecb;
  void * enablecbclosure;

private:
  SoShaderProgram * owner;
  SoGLShaderProgram glShaderProgram;

  static void sensorCB(void * data, SoSensor *);
  SoNodeSensor * sensor;
};

#define PRIVATE(p) ((p)->pimpl)
#define PUBLIC(p) ((p)->owner)

// *************************************************************************

SO_NODE_SOURCE(SoShaderProgram);

// *************************************************************************

// doc from parent
void
SoShaderProgram::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderProgram,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);

  SO_ENABLE(SoGLRenderAction, SoGLShaderProgramElement);
}

/*!
  Constructor.
*/
SoShaderProgram::SoShaderProgram(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderProgram);

  SO_NODE_ADD_FIELD(shaderObject, (NULL));
  this->shaderObject.setNum(0);
  this->shaderObject.setDefault(TRUE);


  PRIVATE(this) = new SoShaderProgramP(this);
  PRIVATE(this)->enablecb = NULL;
  PRIVATE(this)->enablecbclosure = NULL;
}

/*!
  Destructor.
*/
SoShaderProgram::~SoShaderProgram()
{
  delete PRIVATE(this);
}

// doc from parent
void
SoShaderProgram::GLRender(SoGLRenderAction * action)
{
  PRIVATE(this)->GLRender(action);
}

// doc from parent
void
SoShaderProgram::search(SoSearchAction * action)
{
  // Include this node in the search.
  SoNode::search(action);
  if (action->isFound()) return;

  // we really can't do this since this node hasn't got an SoChildList
  // instance
#if 0 // disabled, not possible to search under this node
  int numindices;
  const int * indices;
  if (action->getPathCode(numindices, indices) == SoAction::IN_PATH) {
    // FIXME: not implemented -- 20050129 martin
  }
  else { // traverse all shader objects
    int num = this->shaderObject.getNum();
    for (int i=0; i<num; i++) {
      SoNode * node = this->shaderObject[i];
      action->pushCurPath(i, node);
      SoNodeProfiling profiling;
      profiling.preTraversal(action);
      node->search(action);
      profiling.postTraversal(action);
      action->popCurPath();
      if (action->isFound()) return;
    }
  }
#endif // disabled
}

/*!
  Adds a callback which is called every time this program is enabled/disabled.
*/
void
SoShaderProgram::setEnableCallback(SoShaderProgramEnableCB * cb,
                                   void * closure)
{
  PRIVATE(this)->enablecb = cb;
  PRIVATE(this)->enablecbclosure = closure;
}

// *************************************************************************

SoShaderProgramP::SoShaderProgramP(SoShaderProgram * ownerptr)
{
  PUBLIC(this) = ownerptr;
  this->sensor = new SoNodeSensor(SoShaderProgramP::sensorCB, this);
  this->sensor->attach(ownerptr);
}

SoShaderProgramP::~SoShaderProgramP()
{
  delete this->sensor;
}

void
SoShaderProgramP::GLRender(SoGLRenderAction * action)
{
  SoState *state = action->getState();

  // FIXME: (Martin 2004-09-21) find an alternative to invalidating the cache
  SoCacheElement::invalidate(state);

  this->glShaderProgram.removeShaderObjects();
  this->glShaderProgram.setEnableCallback(this->enablecb,
                                          this->enablecbclosure);

  SoGLShaderProgramElement::set(state, PUBLIC(this), &this->glShaderProgram);

  int i, cnt = PUBLIC(this)->shaderObject.getNum();

  // load shader objects
  for (i = 0; i <cnt; i++) {
    SoNode * node = PUBLIC(this)->shaderObject[i];
    if (node->isOfType(SoShaderObject::getClassTypeId())) {
      ((SoShaderObject *)node)->GLRender(action);
    }
  }

  // enable shader after all shader objects have been loaded
  SoGLShaderProgramElement::enable(state, TRUE);

  // update parameters after all shader objects have been added and enabled

  for (i = 0; i <cnt; i++) {
    SoNode * node = PUBLIC(this)->shaderObject[i];
    if (node->isOfType(SoShaderObject::getClassTypeId())) {
      ((SoShaderObject *)node)->updateParameters(state);
    }
  }
}

void
SoShaderProgramP::sensorCB(void * data, SoSensor *)
{
  // nothing to do now
}

#undef PRIVATE
#undef PUBLIC

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  SoShaderProgram * node = new SoShaderProgram;
  assert(node);
  node->ref();
  BOOST_CHECK_MESSAGE(node->getTypeId() != SoType::badType(),
                      "missing class initialization");
  node->unref();
}

#endif // COIN_TEST_SUITE
