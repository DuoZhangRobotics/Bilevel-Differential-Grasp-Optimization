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
  \class SoGetMatrixAction SoGetMatrixAction.h Inventor/actions/SoGetMatrixAction.h
  \brief The SoGetMatrixAction class is an action for accumulating the transformation matrix of a subgraph.
  \ingroup actions

  This action makes it easy to calculate and convert to and from the
  global coordinate system of your scene and local coordinates of
  parts in a hierarchical model.

  As opposed to most other action types, the SoGetMatrixAction does
  not traverse children of the node it is applied to -- just the node
  itself. When applied to paths, it stops at the last node and does
  not continue further with the children of the tail node.

  Typical usage when querying for world space position, orientation
  and/or scaling would be as follows:

  \code

  // First get hold of an SoPath through the scenegraph down to the
  // node ("mynode") you want to query about its current world space
  // transformation(s).

  SoSearchAction * searchaction = new SoSearchAction;
  searchaction->setNode(mynode);
  searchaction->apply(myscenegraphroot);

  SoPath * path = searchaction->getPath();
  assert(path != NULL);

  // Then apply the SoGetMatrixAction to get the full transformation
  // matrix from world space.

  const SbViewportRegion vpr = myviewer->getViewportRegion();
  SoGetMatrixAction * getmatrixaction = new SoGetMatrixAction(vpr);
  getmatrixaction->apply(path);

  SbMatrix transformation = getmatrixaction->getMatrix();

  // And if you want to access the individual transformation
  // components of the matrix:

  SbVec3f translation;
  SbRotation rotation;
  SbVec3f scalevector;
  SbRotation scaleorientation;

  transformation.getTransform(translation, rotation, scalevector, scaleorientation);
  
  \endcode
*/


// Implementation note: nodes with special behavior on this action
// must set both the matrix and the inverse matrix explicitly, as no
// tracking is done to see when the matrices have been modified.


#include <Inventor/actions/SoGetMatrixAction.h>

#include <cassert>

#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/nodes/SoNode.h>

#include "actions/SoSubActionP.h"

class SoGetMatrixActionP {
public:
};

SO_ACTION_SOURCE(SoGetMatrixAction);


// Overridden from parent.
void
SoGetMatrixAction::initClass(void)
{
  SO_ACTION_INTERNAL_INIT_CLASS(SoGetMatrixAction, SoAction);

  SO_ENABLE(SoGetMatrixAction, SoViewportRegionElement);
}


/*!
  Constructor.

  The \a region viewport specification is not used by this action, and
  is passed along in case it is needed by any nodes.
*/
SoGetMatrixAction::SoGetMatrixAction(const SbViewportRegion & region)
  : viewportregion(region)
{
  SO_ACTION_CONSTRUCTOR(SoGetMatrixAction);
}

/*!
  Destructor.
*/
SoGetMatrixAction::~SoGetMatrixAction()
{
}

/*!
  Set the viewport \a region.

  \sa SoGetMatrixAction::SoGetMatrixAction()
*/
void
SoGetMatrixAction::setViewportRegion(const SbViewportRegion & region)
{
  this->viewportregion = region;
}

/*!
  Returns the viewport region for the action instance.
*/
const SbViewportRegion &
SoGetMatrixAction::getViewportRegion(void) const
{
  return this->viewportregion;
}

/*!
  Returns the accumulated transformation matrix.

  Note: don't modify the returned matrix. This should only be done if
  you are implementing your own transformation type node
  extensions. This advice is also valid for the other matrix access
  methods documented below.
*/
SbMatrix &
SoGetMatrixAction::getMatrix(void)
{
  return this->matrix;
}

/*!
  Returns the inverse of the accumulated transformation matrix.
*/
SbMatrix &
SoGetMatrixAction::getInverse(void)
{
  return this->invmatrix;
}

/*!
  Returns the accumulated texture matrix.
*/
SbMatrix &
SoGetMatrixAction::getTextureMatrix(void)
{
  return this->texmatrix;
}

/*!
  Returns the inverse of the accumulated texture matrix.
*/
SbMatrix &
SoGetMatrixAction::getTextureInverse(void)
{
  return this->invtexmatrix;
}

// Documented in superclass. Overridden from parent class to
// initialize the matrices before traversal starts.
void
SoGetMatrixAction::beginTraversal(SoNode * node)
{
  assert(this->traversalMethods);

  SoViewportRegionElement::set(state, this->viewportregion);
  this->matrix.makeIdentity();
  this->invmatrix.makeIdentity();
  this->texmatrix.makeIdentity();
  this->invtexmatrix.makeIdentity();

  this->traverse(node);
}
