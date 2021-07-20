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

#include <Inventor/VRMLnodes/SoVRML.h>
#include <Inventor/VRMLnodes/SoVRMLNodes.h>

void
so_vrml_init(void)
{
  SoVRMLGeometry::initClass();
  SoVRMLVertexPoint::initClass();
  SoVRMLVertexShape::initClass();
  SoVRMLIndexedShape::initClass();

  SoVRMLParent::initClass();
  SoVRMLGroup::initClass();

  SoVRMLTexture::initClass();

  SoVRMLInterpolator::initClass();

  SoVRMLLight::initClass();
  
  SoVRMLSensor::initClass();
  SoVRMLDragSensor::initClass();

  SoVRMLAnchor::initClass();
  SoVRMLAppearance::initClass();
  SoVRMLAudioClip::initClass();
  SoVRMLBackground::initClass();
  SoVRMLBillboard::initClass();
  SoVRMLBox::initClass();
  SoVRMLCollision::initClass();
  SoVRMLColor::initClass();
  SoVRMLColorInterpolator::initClass();
  SoVRMLCone::initClass();
  SoVRMLCoordinate::initClass();
  SoVRMLCoordinateInterpolator::initClass();
  SoVRMLCylinder::initClass();
  SoVRMLCylinderSensor::initClass();
  SoVRMLDirectionalLight::initClass();
  SoVRMLElevationGrid::initClass();
  SoVRMLExtrusion::initClass();
  SoVRMLFog::initClass();
  SoVRMLFontStyle::initClass();
  SoVRMLImageTexture::initClass();
  SoVRMLIndexedFaceSet::initClass();

  SoVRMLVertexLine::initClass();
  SoVRMLIndexedLine::initClass();
  SoVRMLIndexedLineSet::initClass();
  SoVRMLInline::initClass();
  SoVRMLLOD::initClass();
  SoVRMLShape::initClass();
  SoVRMLMaterial::initClass();
  SoVRMLMovieTexture::initClass();
  SoVRMLNavigationInfo::initClass();
  SoVRMLNormal::initClass();
  SoVRMLNormalInterpolator::initClass();
  SoVRMLOrientationInterpolator::initClass();
  SoVRMLPixelTexture::initClass();
  SoVRMLPlaneSensor::initClass();
  SoVRMLPointLight::initClass();
  SoVRMLPointSet::initClass();
  SoVRMLPositionInterpolator::initClass();
  SoVRMLProximitySensor::initClass();
  SoVRMLScalarInterpolator::initClass();
  SoVRMLScript::initClass();
  SoVRMLSound::initClass();
  SoVRMLSphere::initClass();
  SoVRMLSphereSensor::initClass();
  SoVRMLSpotLight::initClass();
  SoVRMLSwitch::initClass();
  SoVRMLText::initClass();
  SoVRMLTextureCoordinate::initClass();
  SoVRMLTextureTransform::initClass();
  SoVRMLTimeSensor::initClass();
  SoVRMLTouchSensor::initClass();
  SoVRMLTransform::initClass();
  SoVRMLViewpoint::initClass();
  SoVRMLVisibilitySensor::initClass();
  SoVRMLWorldInfo::initClass();
}

#endif // HAVE_VRML97
