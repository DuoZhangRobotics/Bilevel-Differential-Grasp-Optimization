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
  \class SoShapeHints SoShapeHints.h Inventor/nodes/SoShapeHints.h
  \brief The SoShapeHints class is a node containing hints about how to render geometry.
  \ingroup nodes

  The SoShapeHints node is used to set up clues to the rendering
  subsystem about how particular aspects of the subsequent geometry in
  the scene graph should be drawn.

  The default settings of the rendering system is tuned towards
  programmer convenience. For instance, the default is to render both
  sides of all polygons to make sure we avoid any "holes" in the
  geometry if the vertex ordering should happen to differ for
  different polygons.

  If the programmer gives up some of this convenience and uses
  SoShapeHints nodes to more closely specify information about the
  scene graph geometry, the clues given by the SoShapeHints node(s)
  will then be used by the rendering subsystem to avoid certain costly
  operations. Significant gains in rendering speed could be seen as a
  result.

  Be aware that the way backface culling and two-sided lighting is
  decided by the rendering system is not at all intuitive.  Here are
  the common rules of how primitive shapes will render themselves with
  regard to how the SoShapeHints::vertexOrdering and
  SoShapeHints::shapeType fields are set:

  <ul>

  <li>vertexOrdering == CLOCKWISE or COUNTERCLOCKWISE, shapeType ==
      SOLID: causes primitives to be backface culled and rendered with
      one-sided lighting.

  <li>vertexOrdering == CLOCKWISE or COUNTERCLOCKWISE, shapeType ==
      UNKNOWN_SHAPE_TYPE: primitives are \e not backface culled, and
      they are rendered with two-sided lighting.

  <li>vertexOrdering == UNKNOWN_ORDERING, any shapeType: primitives
      are \e not backface culled, and they are rendered with one-sided
      lighting. The OpenGL vertex ordering will be set to counter clockwise
      ordering.
  </ul>

  The UNKNOWN_ORDERING enum has a special and non-intuitive meaning.
  The ordering is really set to counter clockwise -- in OpenGL and
  when generating normals. However, if you want to render your
  geometry with one-sided lighting and backface culling disabled, you
  have to use this enum value, and your polygons need to be in counter
  clockwise ordering.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    ShapeHints {
        vertexOrdering UNKNOWN_ORDERING
        shapeType UNKNOWN_SHAPE_TYPE
        faceType CONVEX
        creaseAngle 0
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/actions/SoCallbackAction.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/elements/SoCreaseAngleElement.h>
#include <Inventor/elements/SoGLShapeHintsElement.h>
#include <Inventor/elements/SoOverrideElement.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \enum SoShapeHints::VertexOrdering

  Enumeration of available ways to specify ordering of vertices for a
  polygon.
*/
/*!
  \var SoShapeHints::VertexOrdering SoShapeHints::UNKNOWN_ORDERING

  Ordering not known, render both sides of the polygon.
*/
/*!
  \var SoShapeHints::VertexOrdering SoShapeHints::CLOCKWISE

  Vertices are specified in a clockwise order.
*/
/*!
  \var SoShapeHints::VertexOrdering SoShapeHints::COUNTERCLOCKWISE

  Vertices are specified in a counter-clockwise order.
*/

/*!
  \enum SoShapeHints::ShapeType
  Enumeration of different shape types.
*/
/*!
  \var SoShapeHints::ShapeType SoShapeHints::UNKNOWN_SHAPE_TYPE
  Nothing known about the shape, be conservative when rendering.
*/
/*!
  \var SoShapeHints::ShapeType SoShapeHints::SOLID

  The subsequent shapes in the graph are all known to be completely
  "closed", solid 3D shapes. Backface culling will be done if
  vertexOrdering is known.
*/

/*!
  \enum SoShapeHints::FaceType
  Enumeration of polygon face types.
*/
/*!
  \var SoShapeHints::FaceType SoShapeHints::UNKNOWN_FACE_TYPE

  Signifies: nothing is known about subsequent polygon data, be
  conservative when rendering.

  All polygons in the scene will be analyzed to see if they needs to
  be tessellated (broken up) into triangles before passed on to the
  underlying immediate mode rendering system.

  The OpenGL rendering system handles most complex polygon types, but
  not all: it can for instance have problems with many-sided, concave
  polygons (concave polygons are "hollow", that is: rounded
  inwards). Coin's internal tessellator will most often handle the
  cases that OpenGL fails on.

  So if you are seeing weird artifacts in how complex polygons are
  rendered, try to change the SoShapeHints::faceType field to this
  value and see if they are then rendered correctly.

  Beware that turning on this functionality might have the effect of
  making the rendering performance worse. If it has a noticable effect
  on your particular scenegraph, we advise that you investigate
  whether you could change how the polygons are generated for Coin
  rendering and then avoid using this flag.
*/
/*!
  \var SoShapeHints::FaceType SoShapeHints::CONVEX

  Subsequent faces are all convex, so turn off the check for and
  tessellation of inconvex faces.

  Subsequent polygons from faceset-type nodes (like SoFaceSet and
  SoIndexedFaceSet) will be sent unmodified to OpenGL, thereby
  assuming that the polygons are in a form handled by OpenGL.
*/


/*!
  \var SoSFEnum SoShapeHints::vertexOrdering

  Specifies how vertices are ordered for polygon faces.

  Set this field to SoShapeHints::CLOCKWISE or
  SoShapeHints::COUNTERCLOCKWISE if possible to turn on backface
  culling and thereby optimize rendering.

  Default value is SoShapeHints::UNKNOWN_ORDERING.
*/
/*!
  \var SoSFEnum SoShapeHints::shapeType

  Hint about whether or not shapes are known to be "closed".  Default
  value is SoShapeHints::UNKNOWN_SHAPE_TYPE.
*/
/*!
  \var SoSFEnum SoShapeHints::faceType

  Hint about whether or not polygon faces are known to be convex.
  Default value is SoShapeHints::CONVEX.
*/
/*!
  \var SoSFFloat SoShapeHints::creaseAngle

  When normals are automatically generated by Coin (i.e. SoNormal
  nodes are not used), this is the smallest angle between two faces
  where we would still calculate normals to do flatshading.

  If the angle between the normals of two neighboring faces is less
  than the value of this field, the faces will be smoothshaded around
  their common edge.
  
  The angle is specified in radians, and the default value is 0.0, 
  meaning no smoothing will be done by default.
*/

// *************************************************************************

SO_NODE_SOURCE(SoShapeHints);

/*!
  Constructor.
*/
SoShapeHints::SoShapeHints(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShapeHints);

  SO_NODE_ADD_FIELD(vertexOrdering, (UNKNOWN_ORDERING));
  SO_NODE_ADD_FIELD(shapeType, (UNKNOWN_SHAPE_TYPE));
  SO_NODE_ADD_FIELD(faceType, (CONVEX));
  SO_NODE_ADD_FIELD(creaseAngle, (0.0f));

  SO_NODE_DEFINE_ENUM_VALUE(VertexOrdering, UNKNOWN_ORDERING);
  SO_NODE_DEFINE_ENUM_VALUE(VertexOrdering, CLOCKWISE);
  SO_NODE_DEFINE_ENUM_VALUE(VertexOrdering, COUNTERCLOCKWISE);

  SO_NODE_DEFINE_ENUM_VALUE(ShapeType, UNKNOWN_SHAPE_TYPE);
  SO_NODE_DEFINE_ENUM_VALUE(ShapeType, SOLID);

  SO_NODE_DEFINE_ENUM_VALUE(FaceType, UNKNOWN_FACE_TYPE);
  SO_NODE_DEFINE_ENUM_VALUE(FaceType, CONVEX);

  SO_NODE_SET_SF_ENUM_TYPE(vertexOrdering, VertexOrdering);
  SO_NODE_SET_SF_ENUM_TYPE(shapeType, ShapeType);
  SO_NODE_SET_SF_ENUM_TYPE(faceType, FaceType);

}

/*!
  Destructor.
*/
SoShapeHints::~SoShapeHints()
{
}

// doc in super
void
SoShapeHints::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShapeHints, SO_FROM_INVENTOR_2_0|SoNode::VRML1);

  SO_ENABLE(SoCallbackAction, SoCreaseAngleElement);
  SO_ENABLE(SoCallbackAction, SoShapeHintsElement);
  SO_ENABLE(SoGLRenderAction, SoCreaseAngleElement);
  SO_ENABLE(SoGLRenderAction, SoGLShapeHintsElement);
  SO_ENABLE(SoGetBoundingBoxAction, SoCreaseAngleElement);
  SO_ENABLE(SoGetBoundingBoxAction, SoShapeHintsElement);
  SO_ENABLE(SoPickAction, SoCreaseAngleElement);
  SO_ENABLE(SoPickAction, SoShapeHintsElement);
}

void
SoShapeHints::doAction(SoAction * action)
{
  SoState * state = action->getState();

  uint32_t flags = SoOverrideElement::getFlags(state);
#define TEST_OVERRIDE(bit) ((SoOverrideElement::bit & flags) != 0)

  // store current values, in case some are overridden or ignored
  SoShapeHintsElement::VertexOrdering vo;
  SoShapeHintsElement::ShapeType st;
  SoShapeHintsElement::FaceType ft;
  SoShapeHintsElement::get(state, vo, st, ft);

  if (!this->vertexOrdering.isIgnored() && !TEST_OVERRIDE(SHAPE_HINTS)) {
    vo = (SoShapeHintsElement::VertexOrdering) this->vertexOrdering.getValue();
    if (this->isOverride()) {
      SoOverrideElement::setShapeHintsOverride(state, this, TRUE);
    }
  }
  if (!this->shapeType.isIgnored() && !TEST_OVERRIDE(SHAPE_HINTS)) {
    st = (SoShapeHintsElement::ShapeType) this->shapeType.getValue();
    if (this->isOverride()) {
      SoOverrideElement::setShapeHintsOverride(state, this, TRUE);
    }
  }
  if (!this->faceType.isIgnored() && !TEST_OVERRIDE(SHAPE_HINTS)) {
    ft = (SoShapeHintsElement::FaceType) this->faceType.getValue();
    if (this->isOverride()) {
      SoOverrideElement::setShapeHintsOverride(state, this, TRUE);
    }
  }
  SoShapeHintsElement::set(action->getState(), this,
                           vo, st, ft);

  if (!this->creaseAngle.isIgnored() && !TEST_OVERRIDE(CREASE_ANGLE)) {
    float ca = this->creaseAngle.getValue();
    // Fix to handle VRML1 ShapeHints nodes correctly. The default
    // creaseAngle value for VRML1 is 0.5, while it's 0.0 for
    // Inventor 2.1
    if (this->creaseAngle.isDefault() &&
        (this->getNodeType() == SoNode::VRML1) &&
        (ca == 0.0f)) {
      ca = 0.5f;
    }
    SoCreaseAngleElement::set(state, this, ca);
    if (this->isOverride()) {
      SoOverrideElement::setCreaseAngleOverride(state, this, TRUE);
    }
  }
#undef TEST_OVERRIDE
}

void
SoShapeHints::GLRender(SoGLRenderAction * action)
{
  SoShapeHints::doAction(action);
}

void
SoShapeHints::callback(SoCallbackAction * action)
{
  SoShapeHints::doAction(action);
}

void
SoShapeHints::pick(SoPickAction * action)
{
  SoShapeHints::doAction(action);
}

void
SoShapeHints::getBoundingBox(SoGetBoundingBoxAction * action)
{
  SoShapeHints::doAction(action);
}
