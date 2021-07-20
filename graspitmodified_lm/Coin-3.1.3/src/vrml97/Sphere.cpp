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
  \class SoVRMLSphere SoVRMLSphere.h Inventor/VRMLnodes/SoVRMLSphere.h
  \brief The SoVRMLSphere class is used to represent a spherical 3D object.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  Sphere {
    field SFFloat radius  1    # (0, inf)
  }
  \endverbatim
  
  The Sphere node specifies a sphere centred at (0, 0, 0) in the local
  coordinate system. The radius field specifies the radius of the
  sphere and shall be greater than zero. Figure 6.15 depicts the
  fields of the Sphere node.

  <center>
  <img src="http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/Images/sphere.gif">
  Figure 6.15 -- Sphere node
  </center>

  When a texture is applied to a sphere, the texture covers the entire
  surface, wrapping counterclockwise from the back of the sphere
  (i.e., longitudinal arc intersecting the -Z-axis) when viewed from
  the top of the sphere. The texture has a seam at the back where the
  X=0 plane intersects the sphere and Z values are
  negative. TextureTransform affects the texture coordinates of the
  Sphere.  The Sphere node's geometry requires outside faces
  only. When viewed from the inside the results are undefined.  

*/

/*!
  \var SoSFFloat SoVRMLSphere::radius
  Sphere radius. Default value is 1.0.
*/

#include <Inventor/VRMLnodes/SoVRMLSphere.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/bundles/SoMaterialBundle.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoLazyElement.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/actions/SoGetPrimitiveCountAction.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/elements/SoGLShapeHintsElement.h>
#include <Inventor/elements/SoTextureCoordinateElement.h>

#include "nodes/SoSubNodeP.h"
#include "misc/SoGenerate.h"
#include "misc/SoPick.h"
#include "misc/SoGL.h"

#define SPHERE_NUM_SLICES 30.0f
#define SPHERE_NUM_STACKS 30.0f

SO_NODE_SOURCE(SoVRMLSphere);

// Doc in parent
void
SoVRMLSphere::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLSphere, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLSphere::SoVRMLSphere(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLSphere);

  SO_VRMLNODE_ADD_FIELD(radius, (1.0f));
}

/*!
  Destructor.
*/
SoVRMLSphere::~SoVRMLSphere()
{
}

// Doc in parent
void
SoVRMLSphere::GLRender(SoGLRenderAction * action)
{
  if (!shouldGLRender(action)) return;

  SoState * state = action->getState();

  SoMaterialBundle mb(action);
  mb.sendFirst();

  SbBool doTextures = SoGLTextureEnabledElement::get(state);

  SbBool sendNormals = !mb.isColorOnly() ||
    (SoTextureCoordinateElement::getType(state) == SoTextureCoordinateElement::FUNCTION);
  
  float complexity = this->getComplexityValue(action);

  unsigned int flags = 0;
  if (sendNormals) flags |= SOGL_NEED_NORMALS;
  if (doTextures) flags |= SOGL_NEED_TEXCOORDS;

  // enable back face culling
  SoGLShapeHintsElement::forceSend(state, TRUE, TRUE);

  sogl_render_sphere(this->radius.getValue(),
                     (int)(SPHERE_NUM_SLICES * complexity),
                     (int)(SPHERE_NUM_STACKS * complexity),
                     &mb,
                     flags, state);
}

// Doc in parent
void
SoVRMLSphere::rayPick(SoRayPickAction * action)
{
  if (!shouldRayPick(action)) return;

  sopick_pick_sphere(this->radius.getValue(),
                     action);
}

// Doc in parent
void
SoVRMLSphere::getPrimitiveCount(SoGetPrimitiveCountAction * action)
{
  if (!this->shouldPrimitiveCount(action)) return;

  float complexity = this->getComplexityValue(action);
  action->addNumTriangles((int)(complexity*2.0f*SPHERE_NUM_SLICES*(SPHERE_NUM_STACKS-1)));
}

void
SoVRMLSphere::generatePrimitives(SoAction * action)
{
  float complexity = this->getComplexityValue(action);

  sogen_generate_sphere(this->radius.getValue(),
                        (int)(SPHERE_NUM_SLICES * complexity),
                        (int)(SPHERE_NUM_STACKS * complexity),
                        this,
                        action);
}

// Doc in parent
void
SoVRMLSphere::computeBBox(SoAction * action,
                          SbBox3f & box,
                          SbVec3f & center)
{
  float r = this->radius.getValue();

  // Allow negative values.
  if (r < 0.0f) r = -r;

  box.setBounds(SbVec3f(-r, -r, -r), SbVec3f(r, r, r));
  center.setValue(0.0f, 0.0f, 0.0f);
}

#undef SPHERE_NUM_SLICES
#undef SPHERE_NUM_STACKS

#endif // HAVE_VRML97
