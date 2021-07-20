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
  \class SoOrthographicCamera SoOrthographicCamera.h Inventor/nodes/SoOrthographicCamera.h
  \brief The SoOrthographicCamera class defines a camera node with orthographic rendering.
  \ingroup nodes

  Orthographic rendering will not give a particularly realistic
  impression of the scene, but non-realistic rendering is for various
  reasons widely used in applications for e.g. Computer Aided Design.

  Also, a simple technique for drawing overlay / HUD style graphics
  (often appearing to be in 2D) can be implemented by setting up a
  "sub scene graph" with an SoOrthographicCamera and the geometry.  As
  a simple demonstration of this concept, load this file into e.g. the
  ExaminerViewer:

  \verbatim
  #Inventor V2.1 ascii
  
  Separator {
     PerspectiveCamera { position 0 0 5 }
     Cone { }
  
     Separator {
        OrthographicCamera {
           position -0.75 -0.75 2
        }
        BaseColor { rgb 1 1 0 }
        Sphere { radius 0.2 }
     }
  }
  \endverbatim

  You will likely encounter Z-buffer issues with this technique which
  makes the overlay / HUD graphics end up interspersed with the "real"
  geometry. If so, this can be solved by e.g. inserting an SoCallback
  node in the sub-scene, where you let the callback disable the depth
  buffer with glDisable(GL_DEPTH_TEST).

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    OrthographicCamera {
        viewportMapping ADJUST_CAMERA
        position 0 0 1
        orientation 0 0 1  0
        nearDistance 1
        farDistance 10
        aspectRatio 1
        focalDistance 5
        height 2
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoOrthographicCamera.h>

#include <Inventor/SbSphere.h>
#include <Inventor/errors/SoDebugError.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoSFFloat SoOrthographicCamera::height
  Height of viewport in world-space scale. Defaults to 2.0 units.
*/

// *************************************************************************

SO_NODE_SOURCE(SoOrthographicCamera);

/*!
  Constructor.
*/
SoOrthographicCamera::SoOrthographicCamera()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoOrthographicCamera);

  SO_NODE_ADD_FIELD(height, (2.0f));
}

/*!
  Destructor.
*/
SoOrthographicCamera::~SoOrthographicCamera()
{
}

// Doc in superclass.
void
SoOrthographicCamera::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoOrthographicCamera, SO_FROM_INVENTOR_1|SoNode::VRML1);
}

/*!
  Scale SoOrthographicCamera::height by multiplying with \a
  scalefactor.
*/
void
SoOrthographicCamera::scaleHeight(float scalefactor)
{
  this->height.setValue(this->height.getValue() * scalefactor);
}

// Doc in superclass.
SbViewVolume
SoOrthographicCamera::getViewVolume(float useaspectratio) const
{
  SbViewVolume volume;

  if (useaspectratio == 0.0f) useaspectratio = this->aspectRatio.getValue();
  float halfheight = this->height.getValue() * 0.5f;
  float halfwidth = halfheight * useaspectratio;
  volume.ortho(-halfwidth, halfwidth,
               -halfheight, halfheight,
               this->nearDistance.getValue(), this->farDistance.getValue());
  
  volume.rotateCamera(this->orientation.getValue());
  volume.translateCamera(this->position.getValue());
  return volume;
}

// Doc in superclass.
void
SoOrthographicCamera::viewBoundingBox(const SbBox3f & box,
                                      float aspect, float slack)
{
#if COIN_DEBUG
  if (box.isEmpty()) {
    SoDebugError::postWarning("SoOrthographicCamera::viewBoundingBox",
                              "bounding box empty");
    return;
  }
#endif // COIN_DEBUG

  // First, we want to move the camera in such a way that it is
  // pointing straight at the center of the scene bounding box -- but
  // without modifiying the rotation value (so we can't use
  // SoCamera::pointAt()).
  SbVec3f cameradirection;
  this->orientation.getValue().multVec(SbVec3f(0, 0, -1), cameradirection);
  this->position.setValue(box.getCenter() + -cameradirection);


  // Get the radius of the bounding sphere.
  SbSphere bs;
  bs.circumscribe(box);
  float radius = bs.getRadius();


  // Make sure that everything will still be inside the viewing volume
  // even if the aspect ratio "favorizes" width over height.
  if (aspect < 1.0f)
    this->height = 2 * radius / aspect;
  else
    this->height = 2 * radius;


  // Move the camera to the edge of the bounding sphere, while still
  // pointing at the scene.
  SbVec3f direction = this->position.getValue() - box.getCenter();
  (void) direction.normalize(); // we know this is not a null vector
  this->position.setValue(box.getCenter() + direction * radius);


  // Set up the clipping planes according to the slack value (a value
  // of 1.0 will yield clipping planes that are tangent to the
  // bounding sphere of the scene).
  float distance_to_midpoint =
    (this->position.getValue() - box.getCenter()).length();
  this->nearDistance = distance_to_midpoint - radius * slack;
  this->farDistance = distance_to_midpoint + radius * slack;


  // The focal distance is simply the distance from the camera to the
  // scene midpoint. This field is not used in rendering, its just
  // provided to make it easier for the user to do calculations based
  // on the distance between the camera and the scene.
  this->focalDistance = distance_to_midpoint;
}
