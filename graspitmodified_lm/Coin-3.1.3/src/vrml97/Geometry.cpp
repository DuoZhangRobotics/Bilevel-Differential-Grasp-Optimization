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
  \class SoVRMLGeometry SoVRMLGeometry.h Inventor/VRMLnodes/SoVRMLGeometry.h
  \brief The SoVRMLGeometry class is a superclass for VRML shapes.
*/

#include <Inventor/VRMLnodes/SoVRMLGeometry.h>

#include <Inventor/VRMLnodes/SoVRMLParent.h>
#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/elements/SoGLShapeHintsElement.h>
#include <Inventor/elements/SoOverrideElement.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/actions/SoSearchAction.h>
#ifdef HAVE_THREADS
#include <Inventor/threads/SbMutex.h>
#endif // HAVE_THREADS

#include "nodes/SoSubNodeP.h"

// *************************************************************************

class SoVRMLGeometryP {
public:  
  SoChildList * childlist;
  SbBool childlistvalid;

#ifdef COIN_THREADSAFE
  SbMutex childlistmutex;
  void lockChildList(void) { this->childlistmutex.lock(); }
  void unlockChildList(void) { this->childlistmutex.unlock(); }
#else // !COIN_THREADSAFE
  void lockChildList(void) { }
  void unlockChildList(void) { }
#endif // !COIN_THREADSAFE
};

#define PRIVATE(thisp) ((thisp)->pimpl)

// *************************************************************************

SO_NODE_ABSTRACT_SOURCE(SoVRMLGeometry);

// *************************************************************************

// Doc in parent
void
SoVRMLGeometry::initClass(void)
{
  SO_NODE_INTERNAL_INIT_ABSTRACT_CLASS(SoVRMLGeometry, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLGeometry::SoVRMLGeometry(void)
{
  PRIVATE(this) = new SoVRMLGeometryP;
  // supply a NULL-pointer as parent, since notifications will be 
  // handled by the fields that actually contain the node(s)
  PRIVATE(this)->childlist = new SoChildList(NULL);
  PRIVATE(this)->childlistvalid = FALSE;
  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLGeometry);
}

/*!
  Destructor.
*/
SoVRMLGeometry::~SoVRMLGeometry()
{
  delete PRIVATE(this)->childlist;
  delete PRIVATE(this);
}

// Doc in parent
SbBool
SoVRMLGeometry::shouldGLRender(SoGLRenderAction * action)
{
  return inherited::shouldGLRender(action);
}
 
/*!
  Convenience method that updates the shape hints element.
*/
void
SoVRMLGeometry::setupShapeHints(SoState * state, const SbBool ccw, const SbBool solid)
{
  if (!(SoOverrideElement::getFlags(state) & SoOverrideElement::SHAPE_HINTS)) {
    SoShapeHintsElement::set(state, this, 
                             ccw ? SoShapeHintsElement::COUNTERCLOCKWISE : SoShapeHintsElement::CLOCKWISE,
                             solid ? SoShapeHintsElement::SOLID : SoShapeHintsElement::UNKNOWN_SHAPE_TYPE,
                             SoShapeHintsElement::FACE_TYPE_AS_IS);
  }
}

// Doc in parent
SoChildList *
SoVRMLGeometry::getChildren(void) const
{
  if (!PRIVATE(this)->childlistvalid) {
    // this is not 100% thread safe. The assumption is that no nodes
    // will be added or removed while a scene graph is being
    // traversed. For Coin, this is an ok assumption.
    PRIVATE(this)->lockChildList();
    // test again after we've locked
    if (!PRIVATE(this)->childlistvalid) {
      SoVRMLGeometry * thisp = (SoVRMLGeometry*) this;
      SoVRMLParent::updateChildList(thisp, *(PRIVATE(thisp)->childlist));
      PRIVATE(thisp)->childlistvalid = TRUE;
    }
    PRIVATE(this)->unlockChildList();
  }
  return PRIVATE(this)->childlist;
}

// Doc in parent
void
SoVRMLGeometry::search(SoSearchAction * action)
{
  SoNode::search(action);
  if (action->isFound() || this->getChildren() == NULL) return;

  int numindices;
  const int * indices;
  if (action->getPathCode(numindices, indices) == SoAction::IN_PATH) {
    this->getChildren()->traverseInPath(action, numindices, indices);
  }
  else {
    this->getChildren()->traverse(action); // traverse all children
  }
}

// Doc in parent
void
SoVRMLGeometry::notify(SoNotList * list)
{
  SoField * f = list->getLastField();
  if (f && f->getTypeId() == SoSFNode::getClassTypeId()) {
    PRIVATE(this)->childlistvalid = FALSE;
  }
  inherited::notify(list);
}

// Doc in parent
void
SoVRMLGeometry::copyContents(const SoFieldContainer * from,
                             SbBool copyConn)
{
  inherited::copyContents(from, copyConn);
  PRIVATE(this)->childlistvalid = FALSE;
  PRIVATE(this)->childlist->truncate(0);
}

#undef PRIVATE

#endif // HAVE_VRML97
