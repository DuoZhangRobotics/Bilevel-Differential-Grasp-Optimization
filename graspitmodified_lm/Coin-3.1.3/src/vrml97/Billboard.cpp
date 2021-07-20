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
  \class SoVRMLBillboard SoVRMLBillboard.h Inventor/VRMLnodes/SoVRMLBillboard.h
  \brief The SoVRMLBillboard class is used for rotating geometry towards the viewpoint.
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  Billboard {
    eventIn      MFNode   addChildren
    eventIn      MFNode   removeChildren
    exposedField SFVec3f  axisOfRotation 0 1 0     # (-inf, inf)
    exposedField MFNode   children       []
    field        SFVec3f  bboxCenter     0 0 0     # (-inf, inf)
    field        SFVec3f  bboxSize       -1 -1 -1  # (0, inf) or -1,-1,-1
  }
  \endverbatim

  The Billboard node is a grouping node which modifies its coordinate
  system so that the Billboard node's local Z-axis turns to point at
  the viewer.  The Billboard node has children which may be other
  children nodes.  The axisOfRotation field specifies which axis to
  use to perform the rotation. This axis is defined in the local
  coordinate system.  When the axisOfRotation field is not (0, 0, 0),
  the following steps describe how to rotate the billboard to face the
  viewer:

  \li a. Compute the vector from the Billboard node's origin to the
  viewer's position. This vector is called the billboard-to-viewer
  vector.

  \li b. Compute the plane defined by the axisOfRotation and the
  billboard-to-viewer vector.

  \li c. Rotate the local Z-axis of the billboard into the plane from b., pivoting
  around the axisOfRotation.

  When the axisOfRotation field is set to (0, 0, 0), the special case
  of viewer-alignment is indicated. In this case, the object rotates
  to keep the billboard's local Y-axis parallel with the Y-axis of the
  viewer.  This special case is distinguished by setting the
  axisOfRotation to (0, 0, 0). The following steps describe how to
  align the billboard's Y-axis to the Y-axis of the viewer:

  \li d. Compute the billboard-to-viewer vector.

  \li e. Rotate the Z-axis of the billboard to be collinear with the
  billboard-to-viewer vector and pointing towards the viewer's
  position.

  \li f. Rotate the Y-axis of the billboard to be parallel and
  oriented in the same direction as the Y-axis of the viewer.

  If the axisOfRotation and the billboard-to-viewer line are
  coincident, the plane cannot be established and the resulting
  rotation of the billboard is undefined. For example, if the
  axisOfRotation is set to (0,1,0) (Y-axis) and the viewer flies over
  the billboard and peers directly down the Y-axis, the results are
  undefined.  Multiple instances of Billboard nodes (DEF/USE) operate
  as expected: each instance rotates in its unique coordinate system
  to face the viewer.  Subclause 4.6.5, Grouping and children nodes
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.5>),
  provides a description of the children, addChildren, and
  removeChildren fields and eventIns.  The bboxCenter and bboxSize
  fields specify a bounding box that encloses the Billboard node's
  children. This is a hint that may be used for optimization
  purposes. The results are undefined if the specified bounding box is
  smaller than the actual bounding box of the children at any time. A
  default bboxSize value, (-1, -1, -1), implies that the bounding box
  is not specified and if needed shall be calculated by the browser. A
  description of the bboxCenter and bboxSize fields is contained in
  4.6.4, Bounding boxes
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.4>),

  \ENDWEB3D

  The following example VRML scene is a simple illustration of how
  SoVRMLBillboard::axisOfRotation constrains rotation around the given
  vector:

  \verbatim
  #VRML V2.0 utf8
  
  Transform {
     translation -2 0 0
     children [
        Billboard {
          children [ Box { } ]
        }
     ]
  }
  
  Transform {
     translation 2 0 0
     children [
        Billboard {
          axisOfRotation 0 1 0
          children [ Box { } ]
        }
     ]
  }
  
  Transform {
     translation 0 -2 0
     children [ Box { size 10 0.1 10 } ]
  }
  \endverbatim

*/

// *************************************************************************

/*!
  \var SoSFVec3f SoVRMLBillboard::axisOfRotation

  The axis of rotation for the geometry.
*/

/*!
  \var SoSFVec3f SoVRMLBillboard::bboxCenter
  The bounding box center hint. Default value is (0, 0, 0).
*/

/*!
  \var SoSFVec3f SoVRMLBillboard::bboxSize
  The bounding box size hint. Default value is (-1, -1, -1).
*/

// *************************************************************************

#include <Inventor/VRMLnodes/SoVRMLBillboard.h>

#include <math.h>
#include <float.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/elements/SoViewVolumeElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/SbPlane.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/SbRotation.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/system/gl.h>

#include "nodes/SoSubNodeP.h"
#include "misc/SoGL.h"
#include "glue/glp.h"
#include "profiler/SoNodeProfiling.h"

// *************************************************************************

SO_NODE_SOURCE(SoVRMLBillboard);

// *************************************************************************

// Doc in parent
void
SoVRMLBillboard::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLBillboard, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLBillboard::SoVRMLBillboard(void)
{
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLBillboard);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(axisOfRotation, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_FIELD(bboxCenter, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_FIELD(bboxSize, (-1.0f, -1.0f, -1.0f));
}

/*!
  Constructor. \a numchildren is the expected number of children.
*/
SoVRMLBillboard::SoVRMLBillboard(int numchildren)
  : inherited(numchildren)
{
}

/*!
  Destructor.
*/
SoVRMLBillboard::~SoVRMLBillboard()
{
}

// *************************************************************************

// Doc in parent
void
SoVRMLBillboard::doAction(SoAction * action)
{
  SoState * state = action->getState();
  state->push();
  this->performRotation(state);
  SoGroup::doAction(action);
  state->pop();
}

// Doc in parent
void
SoVRMLBillboard::callback(SoCallbackAction * action)
{
  SoVRMLBillboard::doAction((SoAction*) action);
}

// Doc in parent
void
SoVRMLBillboard::GLRender(SoGLRenderAction * action)
{
  switch (action->getCurPathCode()) {
  case SoAction::NO_PATH:
  case SoAction::BELOW_PATH:
    this->GLRenderBelowPath(action);
    break;
  case SoAction::OFF_PATH:
    // do nothing. Separator will reset state.
    break;
  case SoAction::IN_PATH:
    this->GLRenderInPath(action);
    break;
  }
}

// Doc in parent
void
SoVRMLBillboard::getBoundingBox(SoGetBoundingBoxAction * action)
{
  SoState * state = action->getState();
  state->push();
  this->performRotation(state);
  SoGroup::getBoundingBox(action);
  state->pop();
}

// Doc in parent
void
SoVRMLBillboard::getMatrix(SoGetMatrixAction * action)
{
  // FIXME: hmmm, what should I do here. pederb, 2001-10-30
}

// Doc in parent
void
SoVRMLBillboard::pick(SoPickAction * action)
{
  SoVRMLBillboard::doAction((SoAction*) action);
}

// Doc in parent
void
SoVRMLBillboard::search(SoSearchAction * action)
{
  SoNode::search(action);
  if (action->isFound()) return;
  SoGroup::doAction(action);
}

// Doc in parent
void
SoVRMLBillboard::GLRenderBelowPath(SoGLRenderAction * action)
{
  SoState * state = action->getState();

  // never cache this node
  SoCacheElement::invalidate(state);

  state->push();
  this->performRotation(state);

  int n = this->getChildren()->getLength();
  SoNode ** childarray = (SoNode**) this->getChildren()->getArrayPtr();

  action->pushCurPath();
  for (int i = 0; i < n && !action->hasTerminated(); i++) {
    action->popPushCurPath(i, childarray[i]);
    if (action->abortNow()) {
      // only cache if we do a full traversal
      break;
    }
    SoNodeProfiling profiling;
    profiling.preTraversal(action);
    childarray[i]->GLRenderBelowPath(action);
    profiling.postTraversal(action);

#if COIN_DEBUG
    // The GL error test is default disabled for this optimized
    // path.  If you get a GL error reporting an error in the
    // Separator node, enable this code by setting the environment
    // variable COIN_GLERROR_DEBUGGING to "1" to see exactly which
    // node caused the error.
    static SbBool chkglerr = sogl_glerror_debugging();
    if (chkglerr) {
      cc_string str;
      cc_string_construct(&str);
      const unsigned int errs = coin_catch_gl_errors(&str);
      if (errs > 0) {
        SoDebugError::post("SoVRMLBillboard::GLRenderBelowPath",
                           "GL error: '%s', nodetype: %s",
                           cc_string_get_text(&str),
                           (*this->getChildren())[i]->getTypeId().getName().getString());
      }
      cc_string_clean(&str);
    }
#endif // COIN_DEBUG
  }
  action->popCurPath();
  state->pop();
}

// Doc in parent
void
SoVRMLBillboard::GLRenderInPath(SoGLRenderAction * action )
{
  int numindices;
  const int * indices;

  SoAction::PathCode pathcode = action->getPathCode(numindices, indices);

  if (pathcode == SoAction::IN_PATH) {
    SoState * state = action->getState();
    SoCacheElement::invalidate(state);
    SoNode ** childarray = (SoNode**) this->getChildren()->getArrayPtr();
    state->push();
    this->performRotation(state);

    int childidx = 0;
    for (int i = 0; i < numindices; i++) {
      for (; childidx < indices[i] && !action->hasTerminated(); childidx++) {
        SoNode * offpath = childarray[childidx];
        if (offpath->affectsState()) {
          action->pushCurPath(childidx, offpath);
          if (!action->abortNow()) {
            SoNodeProfiling profiling;
            profiling.preTraversal(action);
            offpath->GLRenderOffPath(action);
            profiling.postTraversal(action);
          }
          action->popCurPath(pathcode);
        }
      }
      SoNode * inpath = childarray[childidx];
      action->pushCurPath(childidx, inpath);
      if (!action->abortNow()) {
        SoNodeProfiling profiling;
        profiling.preTraversal(action);
        inpath->GLRenderInPath(action);
        profiling.postTraversal(action);
      }
      action->popCurPath(pathcode);
      childidx++;
    }
    state->pop();
  }
  else {
    // we got to the end of the path
    assert(action->getCurPathCode() == SoAction::BELOW_PATH);
    this->GLRenderBelowPath(action);
  }
}

// Doc in parent
void
SoVRMLBillboard::GLRenderOffPath(SoGLRenderAction * action)
{
  // do nothing
}

// Doc in parent
void
SoVRMLBillboard::notify(SoNotList * list)
{
  inherited::notify(list);
}

//
// private method that appends the needed rotation to the state
//
void
SoVRMLBillboard::performRotation(SoState * state) const
{
  SbVec3f rotaxis = this->axisOfRotation.getValue();

  SbMatrix imm = SoModelMatrixElement::get(state).inverse();
  const SbViewVolume & vv = SoViewVolumeElement::get(state);
  SbVec3f up, look, right;
  
  imm.multDirMatrix(vv.getViewUp(), up);
  imm.multDirMatrix(-vv.getProjectionDirection(), look);

  if (rotaxis == SbVec3f(0.0f, 0.0f, 0.0f)) {
    // always orient the billboard towards the viewer
    right = up.cross(look);
    up = look.cross(right);
  } else { 
    // The VRML97 spec calls for rotating the local z-axis of the
    // billboard to face the viewer, pivoting around the axis of
    // rotation. If the axis of rotation is the z axis, this angle
    // will be zero, and no rotation can happen. We don't actually
    // bother to compute this angle at all, but set up = rotaxis and
    // use cross products from there to construct the rotation matrix.
    // This is more numerically stable, more general, and the code is
    // much nicer, but it also means that we must check specifically
    // for this case.
    if (rotaxis == SbVec3f(0.0f, 0.0f, 1.0f)) { return; }
    
    up = rotaxis;
    right = up.cross(look);
    look = right.cross(up);
  }

  // construct the rotation matrix with the vectors defining the
  // desired orientation
  SbMatrix matrix = SbMatrix::identity();

  right.normalize();
  up.normalize();
  look.normalize();
  
  matrix[0][0] = right[0];
  matrix[0][1] = right[1];
  matrix[0][2] = right[2];
  
  matrix[1][0] = up[0];
  matrix[1][1] = up[1];
  matrix[1][2] = up[2];
  
  matrix[2][0] = look[0];
  matrix[2][1] = look[1];
  matrix[2][2] = look[2];
  
  // append the desired rotation to the state
  SoModelMatrixElement::rotateBy(state, (SoNode*) this, SbRotation(matrix));
}

#endif // HAVE_VRML97
