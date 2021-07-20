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

// Note: the class documentation for the basic primitive shapes
// SoSphere, SoCylinder, SoCone and SoCube have many common, or at
// least close to common, paragraphs. If you make any changes, check
// those other shapes too, to see if your updates / fixes should be
// migrated. <mortene@sim.no>.
/*!
  \class SoCone SoCone.h Inventor/nodes/SoCone.h
  \brief The SoCone class is for rendering cone shapes.
  \ingroup nodes

  Insert a cone shape into the scenegraph. The cone is rendered with
  the current material, texture and drawstyle settings (if any,
  otherwise the default settings are used).

  The SoCone node class is provided as a convenient abstraction for
  the application programmer to use "complex" shapes of this type
  without having to do the tessellation to polygons and other
  low-level programming herself.

  A cone is visualized by the underlying rendering system by first
  tessellating the conceptual cone into a set of polygons. To control
  the trade-off between an as much as possible correct visual
  appearance of the cone versus fast rendering, use an SoComplexity
  node to influence the number of polygons generated from the
  tessellation process. (The higher the complexity value, the more
  polygons will be generated, the more \e rounded the sides of the
  cone will look.) Set the SoComplexity::value field to what you
  believe would be a good trade-off between correctness and speed for
  your particular application.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    Cone {
        bottomRadius 1
        height 2
        parts (SIDES | BOTTOM)
    }
  \endcode

  \sa SoCylinder, SoSphere, SoCube
*/

#include <Inventor/nodes/SoCone.h>

#include <assert.h>
#include <math.h>

#include <Inventor/SbLine.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetPrimitiveCountAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/bundles/SoMaterialBundle.h>
#include <Inventor/details/SoConeDetail.h>
#include <Inventor/elements/SoMaterialBindingElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoGLTexture3EnabledElement.h>
#include <Inventor/elements/SoTextureCoordinateElement.h>
#include <Inventor/misc/SoState.h>

#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

#include "nodes/SoSubNodeP.h"
#include "misc/SoGL.h"
#include "misc/SoPick.h"
#include "misc/SoGenerate.h"

/*!
  \enum SoCone::Part
  Enumerates the various parts of the cone, for setting inclusion or
  exclusion from the shape.
*/

/*!
  \var SoSFBitMask SoCone::parts
  The parts to use for the cone shape. Defaults to SoCone::ALL.
*/
/*!
  \var SoSFFloat SoCone::bottomRadius
  Radius of the cone's bottom disc. Default value is 1.0.
*/
/*!
  \var SoSFFloat SoCone::height
  Height of cone. Default value is 2.0.
*/

#define CONE_SIDE_NUMTRIS 40.0f

// *************************************************************************

SO_NODE_SOURCE(SoCone);

/*!
  Constructor.
*/
SoCone::SoCone(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoCone);

  SO_NODE_ADD_FIELD(bottomRadius, (1.0f));
  SO_NODE_ADD_FIELD(height, (2.0f));
  SO_NODE_ADD_FIELD(parts, (ALL));

  SO_NODE_DEFINE_ENUM_VALUE(Part, SIDES);
  SO_NODE_DEFINE_ENUM_VALUE(Part, BOTTOM);
  SO_NODE_DEFINE_ENUM_VALUE(Part, ALL);
  SO_NODE_SET_SF_ENUM_TYPE(parts, Part);
}

/*!
  Destructor.
*/
SoCone::~SoCone()
{
}


// Doc from parent.
void
SoCone::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoCone, SO_FROM_INVENTOR_1|SoNode::VRML1);
}

// Doc from parent.
void
SoCone::computeBBox(SoAction * action, SbBox3f & box, SbVec3f & center)
{
  float r = this->bottomRadius.getValue();
  float h = this->height.getValue();

  // Allow negative values.
  if (h < 0.0f) h = -h;
  if (r < 0.0f) r = -r;

  float half_height = h/2.0f;

  // The SIDES are present, so just find the middle point and enclose
  // everything.
  if (this->parts.getValue() & SoCone::SIDES) {
    center.setValue(0.0f, 0.0f, 0.0f);
    box.setBounds(SbVec3f(-r, -half_height, -r), SbVec3f(r, half_height, r));
  }
  // ..no SIDES, but we've still got the bottom (NB: OIV misses this case).
  else if (this->parts.getValue() & SoCone::BOTTOM) {
    center.setValue(0.0f, -half_height, 0.0f);
    box.setBounds(SbVec3f(-r, -half_height, -r), SbVec3f(r, -half_height, r));
  }
  // ..no parts present. My confidence is shot -- I feel very small.
  else {
    center.setValue(0.0f, 0.0f, 0.0f);
    box.setBounds(SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(0.0f, 0.0f, 0.0f));
  }
}

// Doc from parent.
void
SoCone::GLRender(SoGLRenderAction * action)
{
  if (!shouldGLRender(action)) return;

  SoState * state = action->getState();

  SbBool doTextures = FALSE;
  SbBool do3DTextures = FALSE;
  if (SoGLTextureEnabledElement::get(state)) doTextures = TRUE;
  else if (SoGLTexture3EnabledElement::get(state)) do3DTextures = TRUE;

  SoCone::Part p = (SoCone::Part) this->parts.getValue();

  SoMaterialBundle mb(action);
  SbBool sendNormals = !mb.isColorOnly() || 
    (SoTextureCoordinateElement::getType(state) == SoTextureCoordinateElement::FUNCTION);

  unsigned int flags = 0;
  if (doTextures) flags |= SOGL_NEED_TEXCOORDS;
  else if (do3DTextures) flags |= SOGL_NEED_3DTEXCOORDS;
  if (sendNormals) flags |= SOGL_NEED_NORMALS;
  if (p & SoCone::SIDES) flags |= SOGL_RENDER_SIDE;
  if (p & SoCone::BOTTOM) flags |= SOGL_RENDER_BOTTOM;

  SoMaterialBindingElement::Binding bind =
    SoMaterialBindingElement::get(state);
  if (bind == SoMaterialBindingElement::PER_PART ||
      bind == SoMaterialBindingElement::PER_PART_INDEXED)
    flags |= SOGL_MATERIAL_PER_PART;

  mb.sendFirst();

  float complexity = this->getComplexityValue(action);

  sogl_render_cone(this->bottomRadius.getValue(),
                   this->height.getValue(),
                   (int)(CONE_SIDE_NUMTRIS * complexity),
                   &mb,
                   flags, state);

#if COIN_DEBUG && 0 // debug
  SoDebugError::postInfo("SoCone::GLRender", "end");
#endif // debug
}

/*!
  Add a \a part to the cone.

  \sa removePart(), hasPart()
*/
void
SoCone::addPart(SoCone::Part part)
{
  if (this->hasPart(part)) {
#if COIN_DEBUG
    SoDebugError::postWarning("SoCone::addPart", "part already set");
#endif // COIN_DEBUG
    return;
  }

  this->parts.setValue(this->parts.getValue() | part);
}

/*!
  Remove a \a part from the cone.

  \sa addPart(), hasPart()
*/
void
SoCone::removePart(SoCone::Part part)
{
  if (!this->hasPart(part)) {
#if COIN_DEBUG
    SoDebugError::postWarning("SoCone::removePart", "part was not set");
#endif // COIN_DEBUG
    return;
  }

  this->parts.setValue(this->parts.getValue() & ~part);
}

/*!
  Returns \c TRUE if rendering of the given \a part is currently
  turned on.

  \sa addPart(), removePart()
*/
SbBool
SoCone::hasPart(SoCone::Part part) const
{
  return (this->parts.getValue() & part) ? TRUE : FALSE;
}

// Doc from parent.
void
SoCone::rayPick(SoRayPickAction * action)
{
  if (!this->shouldRayPick(action)) return;

  SoCone::Part p = (SoCone::Part) this->parts.getValue();
  unsigned int flags = 0;
  if (p & SoCone::SIDES) flags |= SOPICK_SIDES;
  if (p & SoCone::BOTTOM) flags |= SOPICK_BOTTOM;

  SoMaterialBindingElement::Binding bind =
    SoMaterialBindingElement::get(action->getState());
  if (bind == SoMaterialBindingElement::PER_PART ||
      bind == SoMaterialBindingElement::PER_PART_INDEXED)
    flags |= SOPICK_MATERIAL_PER_PART;
  
  sopick_pick_cone(this->bottomRadius.getValue(),
                   this->height.getValue(),
                   flags, this, action);
}

// Doc from parent.
void
SoCone::getPrimitiveCount(SoGetPrimitiveCountAction * action)
{
  if (!this->shouldPrimitiveCount(action)) return;

  float complexity = this->getComplexityValue(action);
  int numtris = (int)(complexity * CONE_SIDE_NUMTRIS);

  if (this->parts.getValue() & SoCone::BOTTOM) {
    action->addNumTriangles(numtris-2);
  }
  if (this->parts.getValue() & SoCone::SIDES) {
    action->addNumTriangles(numtris);
  }
}

// Doc from parent.
void
SoCone::generatePrimitives(SoAction * action)
{
  SoCone::Part p = (SoCone::Part) this->parts.getValue();
  unsigned int flags = 0;
  if (p & SoCone::SIDES) flags |= SOGEN_GENERATE_SIDE;
  if (p & SoCone::BOTTOM) flags |= SOGEN_GENERATE_BOTTOM;

  SoMaterialBindingElement::Binding bind =
    SoMaterialBindingElement::get(action->getState());
  if (bind == SoMaterialBindingElement::PER_PART ||
      bind == SoMaterialBindingElement::PER_PART_INDEXED)
    flags |= SOGL_MATERIAL_PER_PART;

  float complexity = this->getComplexityValue(action);

  sogen_generate_cone(this->bottomRadius.getValue(),
                      this->height.getValue(),
                      (int)(CONE_SIDE_NUMTRIS * complexity),
                      flags,
                      this,
                      action);
}

#undef CONE_SIDE_NUMTRIS
