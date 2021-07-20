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
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_VRML97

/*!
  \class SoVRMLAppearance SoVRMLAppearance.h Inventor/VRMLnodes/SoVRMLAppearance.h
  \brief The SoVRMLAppearance class specifies visual properties for shapes.
  \ingroup VRMLnodes
  
  \WEB3DCOPYRIGHT

  \verbatim

  Appearance {
    exposedField SFNode material          NULL
    exposedField SFNode texture           NULL
    exposedField SFNode textureTransform  NULL
  }
  \endverbatim

  The Appearance node specifies the visual properties of geometry. The
  value for each of the fields in this node may be NULL. However, if
  the field is non-NULL, it shall contain one node of the appropriate
  type.  The material field, if specified, shall contain a VRMLMaterial
  node. If the material field is NULL or unspecified, lighting is off
  (all lights are ignored during rendering of the object that
  references this Appearance) and the unlit object colour is (1, 1,
  1). Details of the VRML lighting model are in 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>).

  The texture field, if specified, shall contain one of the various
  types of texture nodes (VRMLImageTexture, VRMLMovieTexture, or
  VRMLPixelTexture).  If the texture node is NULL or the texture field
  is unspecified, the object that references this Appearance is not
  textured.  The textureTransform field, if specified, shall contain a
  VRMLTextureTransform node. If the textureTransform is NULL or
  unspecified, the textureTransform field has no effect.

*/

/*!
  \var SoSFNode SoVRMLAppearance::material
  Can contain an SoVRMLMaterial node. Is NULL by default.
*/

/*!
  \var SoSFNode SoVRMLAppearance::texture
  Can contain a texture node. Is NULL by default.
*/

/*!
  \var SoSFNode SoVRMLAppearance::textureTransform
  Can contain an SoVRMLTextureTransform node. Is NULL by default.
*/

#include <Inventor/VRMLnodes/SoVRMLAppearance.h>

#include <stddef.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/VRMLnodes/SoVRMLParent.h>
#include <Inventor/VRMLnodes/SoVRMLMaterial.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/actions/SoActions.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/elements/SoLazyElement.h>
#include <Inventor/elements/SoTextureQualityElement.h>
#include <Inventor/elements/SoTextureImageElement.h>

#ifdef HAVE_THREADS
#include <Inventor/threads/SbMutex.h>
#endif // HAVE_THREADS

#include "nodes/SoSubNodeP.h"
#include "profiler/SoNodeProfiling.h"

// *************************************************************************

class SoVRMLAppearanceP {
public:
  SoChildList * childlist;
  SbBool childlistvalid;

#ifdef COIN_THREADSAFE
  SbMutex mutex;
  void lock(void) { this->mutex.lock(); }
  void unlock(void) { this->mutex.unlock(); }
#else // !COIN_THREADSAFE
  void lock(void) {  }
  void unlock(void) { }
#endif // !COIN_THREADSAFE
  uint32_t fakecolor;
};

#define PRIVATE(thisp) ((thisp)->pimpl)

// *************************************************************************

SO_NODE_SOURCE(SoVRMLAppearance);

// *************************************************************************

// doc in parent
void
SoVRMLAppearance::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLAppearance, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLAppearance::SoVRMLAppearance(void)
{
  PRIVATE(this) = new SoVRMLAppearanceP;
  // supply a NULL-pointer as parent, since notifications will be 
  // handled by the fields that actually contain the node(s)
  PRIVATE(this)->childlist = new SoChildList(NULL);
  PRIVATE(this)->childlistvalid = FALSE;

  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLAppearance);
  
  SO_VRMLNODE_ADD_EXPOSED_FIELD(material, (NULL));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(texture, (NULL));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(textureTransform, (NULL));
}

/*!
  Destructor.
*/
SoVRMLAppearance::~SoVRMLAppearance()
{
  delete PRIVATE(this)->childlist;
  delete PRIVATE(this);
}

// doc in parent
void
SoVRMLAppearance::doAction(SoAction * action)
{
  int numindices;
  const int * indices;
  if (action->getPathCode(numindices, indices) == SoAction::IN_PATH) {
    this->getChildren()->traverseInPath(action, numindices, indices);
  }
  else {
    this->getChildren()->traverse(action); // traverse all children
  }
}

// doc in parent
void
SoVRMLAppearance::callback(SoCallbackAction * action)
{
  SoVRMLAppearance::doAction(action);
}

// doc in parent
void
SoVRMLAppearance::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();
  int numindices;
  const int * indices;
  SoAction::PathCode pathcode = action->getPathCode(numindices, indices);

  SoNode ** childarray = (SoNode**) this->getChildren()->getArrayPtr();

  if (pathcode == SoAction::IN_PATH) {
    int lastchild = indices[numindices - 1];
    for (int i = 0; i <= lastchild && !action->hasTerminated(); i++) {
      SoNode * child = childarray[i];
      action->pushCurPath(i, child);
      if (action->getCurPathCode() != SoAction::OFF_PATH ||
          child->affectsState()) {
        if (!action->abortNow()) {
          SoNodeProfiling profiling;
          profiling.preTraversal(action);
          child->GLRender(action);
          profiling.postTraversal(action);
        }
        else {
          SoCacheElement::invalidate(state);
        }
      }
      action->popCurPath(pathcode);
    }
  }
  else {
    action->pushCurPath();
    int n = this->getChildren()->getLength();
    for (int i = 0; i < n && !action->hasTerminated(); i++) {
      action->popPushCurPath(i, childarray[i]);
      if (action->abortNow()) {
        // only cache if we do a full traversal
        SoCacheElement::invalidate(state);
        break;
      }
      SoNodeProfiling profiling;
      profiling.preTraversal(action);
      childarray[i]->GLRender(action);
      profiling.postTraversal(action);
    }
    action->popCurPath();
  }

  // workaround for weird (IMO) VRML97 texture/material handling. For 
  // RGB[A] textures, the texture color should replace the diffuse color.
  SbVec2s size;
  int nc;
  (void) SoTextureImageElement::getImage(state, size, nc);
  
  if (this->texture.getValue() && 
      SoTextureQualityElement::get(state) > 0.0f && 
      size != SbVec2s(0,0) && 
      nc >= 3) {

    float t = SoLazyElement::getTransparency(state, 0);
    uint32_t alpha = (uint32_t) ((1.0f - t) * 255.0f);
    
    // lock just in case two threads get here at the same time
    PRIVATE(this)->lock();
    PRIVATE(this)->fakecolor = 0xffffff00 | alpha;
    PRIVATE(this)->unlock();
    SoLazyElement::setPacked(state, this, 1, &PRIVATE(this)->fakecolor, alpha != 255);
  }
}

// doc in parent
void
SoVRMLAppearance::search(SoSearchAction * action)
{
  SoNode::search(action);
  if (action->isFound()) return;
  SoVRMLAppearance::doAction(action);
}

// doc in parent
SoChildList *
SoVRMLAppearance::getChildren(void) const
{
  if (!PRIVATE(this)->childlistvalid) {
    // this is not 100% thread safe. The assumption is that no nodes
    // will be added or removed while a scene graph is being
    // traversed. For Coin, this is an ok assumption.
    PRIVATE(this)->lock();
    // test again after we've locked
    if (!PRIVATE(this)->childlistvalid) {
      SoVRMLAppearance * thisp = (SoVRMLAppearance*) this;
      SoVRMLParent::updateChildList(thisp, *(PRIVATE(thisp)->childlist));
      PRIVATE(thisp)->childlistvalid = TRUE;
    }
    PRIVATE(this)->unlock();
  }
  return PRIVATE(this)->childlist;
}

// doc in parent
void
SoVRMLAppearance::notify(SoNotList * list)
{
  SoField * f = list->getLastField();
  if (f && f->getTypeId() == SoSFNode::getClassTypeId()) {
    PRIVATE(this)->childlistvalid = FALSE;
  }
  inherited::notify(list);
}

// doc in parent
void
SoVRMLAppearance::copyContents(const SoFieldContainer * from,
                               SbBool copyConn)
{
  inherited::copyContents(from, copyConn);
  PRIVATE(this)->childlistvalid = FALSE;
  PRIVATE(this)->childlist->truncate(0);
}

#undef PRIVATE

#endif // HAVE_VRML97
