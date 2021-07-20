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
  \class SoChildList SoChildList.h Inventor/misc/SoChildList.h
  \brief The SoChildList class is a container for node children.
  \ingroup general

  This class does automatic notification on the parent nodes upon
  adding or removing children.

  Methods for action traversal of the children are also provided.
*/

#include <Inventor/misc/SoChildList.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/SbName.h>

#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG



/*!
  Default constructor, sets parent container and initializes a minimal
  list.
*/
SoChildList::SoChildList(SoNode * const parentptr)
  : SoNodeList()
{
  this->parent = parentptr;
}

/*!
  Constructor with hint about list size.

  \sa SoNodeList::SoNodeList(const int)
*/
SoChildList::SoChildList(SoNode * const parentptr, const int size)
  : SoNodeList(size)
{
  this->parent = parentptr;
}

/*!
  Copy constructor.

  \sa SoNodeList::SoNodeList(const SoNodeList &)
*/
SoChildList::SoChildList(SoNode * const parentptr, const SoChildList & cl)
  : SoNodeList()
{
  this->parent = parentptr;
  this->copy(cl);
}

/*!
  Destructor.
*/
SoChildList::~SoChildList()
{
  this->truncate(0);
}

/*!
  Append a new \a node instance as a child of our parent container.

  Automatically notifies parent node and any SoPath instances auditing
  paths with nodes from this list.
*/
void
SoChildList::append(SoNode * const node)
{
  if (this->parent) {
    node->addAuditor(this->parent, SoNotRec::PARENT);
  }
  SoNodeList::append(node);

  if (this->parent) {
    this->parent->startNotify();
  }
  // Doesn't need to notify SoPath auditors, as adding a new node at
  // _the end_ won't affect any path "passing through" this childlist.
}

/*!
  Insert a new \a node instance as a child of our parent container at
  position \a addbefore.

  Automatically notifies parent node and any SoPath instances auditing
  paths with nodes from this list.
*/
void
SoChildList::insert(SoNode * const node, const int addbefore)
{
  assert(addbefore <= this->getLength());
  if (this->parent) {
    node->addAuditor(this->parent, SoNotRec::PARENT);
  }
  SoNodeList::insert(node, addbefore);

  // FIXME: shouldn't we move this startNotify() call to the end of
  // the function?  pederb, 2002-10-02
  if (this->parent) {
    this->parent->startNotify();
    for (int i=0; i < this->auditors.getLength(); i++) {
      this->auditors[i]->insertIndex(this->parent, addbefore);
    }
  }
}

/*!
  Remove the child node pointer at \a index.

  Automatically notifies parent node and any SoPath instances auditing
  paths with nodes from this list.
*/
void
SoChildList::remove(const int index)
{
  assert(index >= 0 && index < this->getLength());
  if (this->parent) {
    SoNodeList::operator[](index)->removeAuditor(this->parent, SoNotRec::PARENT);
  }
  // FIXME: we experienced memory corruption if the
  // SoNodeList::remove(index) statement was placed here (before
  // updating paths). It seems to be working ok now, but we should
  // figure out exactly why we can't remove the node before updating
  // the paths.  pederb, 2002-10-02
  if (this->parent) {
    for (int i=0; i < this->auditors.getLength(); i++) {
      this->auditors[i]->removeIndex(this->parent, index);
    }
  }
  SoNodeList::remove(index);
  if (this->parent) {
    this->parent->startNotify();
  }
}

// Documented in superclass. Overridden to handle notification.
void
SoChildList::truncate(const int length)
{
  const int n = this->getLength();
  assert(length >= 0 && length <= n);

  if (length != n) {
    if (this->parent) {
      for (int i = length; i < n; i++) {
        SoNodeList::operator[](i)->removeAuditor(this->parent, SoNotRec::PARENT);
      }
    }
    SoNodeList::truncate(length);

    // FIXME: shouldn't we move this startNotify() call to the end of
    // the function?  pederb, 2002-10-02
    if (this->parent) {
      this->parent->startNotify();
      for (int k=0; k < this->auditors.getLength(); k++) {
        for (int j=length; j < n; j++) {
          this->auditors[k]->removeIndex(this->parent, j);
        }
      }
    }
  }
}

/*!
  Copy contents of \a cl into this list.
*/
void
SoChildList::copy(const SoChildList & cl)
{
  // Overridden from superclass to handle notification.

  if (this == &cl) return;

  // Call truncate() explicitly here to get the path notification.
  this->truncate(0);
  SoBaseList::copy(cl);

  // it's important to add parent as auditor for all nodes (this is
  // usually done in SoChildList::append/insert)
  if (this->parent) {
    for (int i = 0; i < this->getLength(); i++) {
      (*this)[i]->addAuditor(this->parent, SoNotRec::PARENT);
    }
    this->parent->startNotify();
  }
}

/*!
  Index operator to set element at \a index. Does \e not expand array
  bounds if \a index is outside the list.
*/
void
SoChildList::set(const int index, SoNode * const node)
{
  // Overridden from superclass to handle notification.

#if COIN_DEBUG && 0 // debug
  SoDebugError::postInfo("SoChildList::set",
                         "(%p) index=%d, node=%p, oldnode=%p",
                         this, index, node, (*this)[index]);
#endif // debug

  assert(index >= 0 && index < this->getLength());
  if (this->parent) {
    SoNodeList::operator[](index)->removeAuditor(this->parent, SoNotRec::PARENT);
  }
  SoBaseList::set(index, (SoBase *)node);
  if (this->parent) {
    node->addAuditor(this->parent, SoNotRec::PARENT);
  }
  // FIXME: shouldn't we move this startNotify() call to the end of
  // the function?  pederb, 2002-10-02
  if (this->parent) {
    this->parent->startNotify();
    for (int i=0; i < this->auditors.getLength(); i++)
      this->auditors[i]->replaceIndex(this->parent, index, node);
  }
}

/*!
  Optimized IN_PATH traversal method.
  
  This method is an extension versus the Open Inventor API.
*/
void
SoChildList::traverseInPath(SoAction * const action,
                            const int numindices,
                            const int * indices)
{
  assert(action->getCurPathCode() == SoAction::IN_PATH);

  // only traverse nodes in path list, and nodes off path that
  // affects state.
  int childidx = 0;

  for (int i = 0; i < numindices && !action->hasTerminated(); i++) {
    int stop = indices[i];
    for (; childidx < stop && !action->hasTerminated(); childidx++) {
      // we are off path. Check if node affects state before traversing
      SoNode * node = (*this)[childidx];
      if (node->affectsState()) {
        action->pushCurPath(childidx, node);
        action->traverse(node);
        action->popCurPath(SoAction::IN_PATH);
      }
    }

    if (!action->hasTerminated()) {
      // here we are in path. Always traverse
      SoNode * node = (*this)[childidx];
      action->pushCurPath(childidx, node);
      action->traverse(node);
      action->popCurPath(SoAction::IN_PATH);
      childidx++;
    }
  }
}

/*!
  Traverse child nodes in the list from index \a first up to and
  including index \a last, or until the SoAction::hasTerminated() flag
  of \a action has been set.
*/
void
SoChildList::traverse(SoAction * const action, const int first, const int last)
{
  int i;
  SoNode * node = NULL;

  assert((first >= 0) && (first < this->getLength()) && "index out of bounds");
  assert((last >= 0) && (last < this->getLength()) && "index out of bounds");
  assert((last >= first) && "erroneous indices");

#if COIN_DEBUG
  // Calculate a checksum over the children node pointers, to later
  // catch attempts at changing the scene graph layout mid-traversal
  // with an assert. (chksum reversed to initial value and controlled
  // at the bottom end of this function.)
  //
  // Note: we might find this to be overly strict, because there are
  // cases where this will stop an unharmful attempt at changing the
  // current group node's set of children. But that's only if the
  // application programmer _really_, _really_ know what he is doing,
  // and it's still a slippery slope.. so "better safe than sorry" and
  // all that.
  //
  // mortene.
  uintptr_t chksum = 0xdeadbeef;
  for (i = first; i <= last; i++) { chksum ^= (uintptr_t)(*this)[i]; }
  SbBool changedetected = FALSE;
#endif // COIN_DEBUG

  SoAction::PathCode pathcode = action->getCurPathCode();

  switch (pathcode) {
  case SoAction::NO_PATH:
  case SoAction::BELOW_PATH:
    // always traverse all nodes.
    action->pushCurPath();
    for (i = first; (i <= last) && !action->hasTerminated(); i++) {
#if COIN_DEBUG
      if (i >= this->getLength()) {
        changedetected = TRUE;
        break;
      }
#endif // COIN_DEBUG
      node = (*this)[i];
      action->popPushCurPath(i, node);
      action->traverse(node);
    }
    action->popCurPath();
    break;
  case SoAction::OFF_PATH:
    for (i = first; (i <= last) && !action->hasTerminated(); i++) {      
#if COIN_DEBUG
      if (i >= this->getLength()) {
        changedetected = TRUE;
        break;
      }
#endif // COIN_DEBUG
      node = (*this)[i];
      // only traverse nodes that affects state
      if (node->affectsState()) {
        action->pushCurPath(i, node);
        action->traverse(node);
        action->popCurPath(pathcode);
      }
    }
    break;
  case SoAction::IN_PATH:
    for (i = first; (i <= last) && !action->hasTerminated(); i++) {
#if COIN_DEBUG
      if (i >= this->getLength()) {
        changedetected = TRUE;
        break;
      }
#endif // COIN_DEBUG
      node = (*this)[i];
      action->pushCurPath(i, node);
      // if we're OFF_PATH after pushing, we only traverse if the node
      // affects the state.
      if ((action->getCurPathCode() != SoAction::OFF_PATH) ||
          node->affectsState()) {
        action->traverse(node);
      }
      action->popCurPath(pathcode);
    }
    break;
  default:
    assert(0 && "unknown path code.");
    break;
  }

#if COIN_DEBUG
  if (!changedetected) {
    for (i = last; i >= first; i--) { chksum ^= (uintptr_t)(*this)[i]; }
    if (chksum != 0xdeadbeef) changedetected = TRUE;
  }
  if (changedetected) {
    SoDebugError::postWarning("SoChildList::traverse",
                              "Detected modification of scene graph layout "
                              "during action traversal. This is considered to "
                              "be hazardous and error prone, and we "
                              "strongly advice you to change your code "
                              "and/or reorganize your scene graph so that "
                              "this is not necessary.");
  }
#endif // COIN_DEBUG
}

/*!
  Traverse all nodes in the list, invoking their methods for the given
  \a action.
*/
void
SoChildList::traverse(SoAction * const action)
{
  if (this->getLength() == 0) return;
  this->traverse(action, 0, this->getLength() - 1);
}

/*!
  Traverse the node at \a index (and possibly its children, if its a
  group node), applying the nodes' method for the given \a action.
*/
void
SoChildList::traverse(SoAction * const action, const int index)
{
  assert((index >= 0) && (index < this->getLength()) && "index out of bounds");
  this->traverse(action, index, index);
}

/*!
  Traverse the \a node (and possibly its children, if its a group
  node), applying the nodes' method for the given \a action.
*/
void
SoChildList::traverse(SoAction * const action, SoNode * node)
{
  int idx = this->find(node);
  assert(idx != -1);
  this->traverse(action, idx);
}

/*!
  Notify \a path whenever this list of node children changes.
*/
void
SoChildList::addPathAuditor(SoPath * const path)
{
#if COIN_DEBUG && 0 // debug
  SoDebugError::postInfo("SoChildList::addPathAuditor",
                         "add SoPath auditor %p to list %p", path, this);
#endif // debug

  this->auditors.append(path);
}

/*!
  Remove \a path as an auditor for our list of node children.
*/
void
SoChildList::removePathAuditor(SoPath * const path)
{
#if COIN_DEBUG && 0 // debug
  SoDebugError::postInfo("SoChildList::removePathAuditor",
                         "remove SoPath auditor %p from list %p", path, this);
#endif // debug

  const int index = this->auditors.find(path);
#if COIN_DEBUG
  if (index == -1) {
    SoDebugError::post("SoChildList::removePathAuditor",
                       "no SoPath %p is auditing list %p! (of parent %p (%s))",
                       path,
                       this,
                       this->parent,
                       this->parent ? this->parent->getTypeId().getName().getString() : "<no type>");
    return;
  }
#endif // COIN_DEBUG
  this->auditors.remove(index);
}
