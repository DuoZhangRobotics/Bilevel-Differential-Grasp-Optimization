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
  \class SoActionMethodList SoActionMethodList.h Inventor/lists/SoActionMethodList.h
  \brief The SoActionMethodList class contains function pointers for action methods.
  \ingroup actions

  An SoActionMethodList contains one function pointer per node
  type. Each action contains an SoActioMethodList to know which
  functions to call during scene graph traversal.
*/

#include <Inventor/lists/SoActionMethodList.h>
#include <Inventor/lists/SoTypeList.h>
#include <Inventor/lists/SbList.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/nodes/SoNode.h>
#include <assert.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#ifdef COIN_THREADSAFE
#include <Inventor/threads/SbMutex.h>
#endif // COIN_THREADSAFE

#ifndef DOXYGEN_SKIP_THIS

class SoActionMethodListP {
public:
  SoActionMethodList * parent;
  int setupnumtypes;
  SbList <SoType> addedtypes;
  SbList <SoActionMethod> addedmethods;

#ifdef COIN_THREADSAFE
  SbMutex mutex;
#endif // COIN_THREADSAFE
  void lock(void) {
#ifdef COIN_THREADSAFE
    this->mutex.lock();
#endif
  }
  void unlock(void) {
#ifdef COIN_THREADSAFE
    this->mutex.unlock();
#endif
  }
};

#endif // DOXYGEN_SKIP_THIS

#define PRIVATE(obj) ((obj)->pimpl)

/*!
  The constructor.  The \a parentlist argument is the parent action's
  action method list.  It can be \c NULL for action method lists that
  are not based on inheriting from a parent action.
*/
SoActionMethodList::SoActionMethodList(SoActionMethodList * const parentlist)
{
  PRIVATE(this) = new SoActionMethodListP;
  PRIVATE(this)->parent = parentlist;
  PRIVATE(this)->setupnumtypes = 0;
}

/*!
  Destructor.
*/
SoActionMethodList::~SoActionMethodList()
{
  delete PRIVATE(this);
}

// Documented in superclass. Overridden from parent to cast from \c
// void pointer.
SoActionMethod &
SoActionMethodList::operator[](const int index)
{
  return (SoActionMethod&)SbPList::operator[](index);
}

/*!
  Add a function pointer to a node type's action method.
*/
void
SoActionMethodList::addMethod(const SoType node, const SoActionMethod method)
{
  assert(node != SoType::badType());
  PRIVATE(this)->lock();
  PRIVATE(this)->addedtypes.append(node);
  PRIVATE(this)->addedmethods.append(method);
  PRIVATE(this)->setupnumtypes = 0; // force a new setUp
  PRIVATE(this)->unlock();
}

// dummy method used for detecting unset action methods
static void unsetActionMethod(SoAction *, SoNode *)
{
}

/*!
  This method must be called as the last initialization step before
  using the list. It fills in \c NULL entries with the parent's
  method.
*/
void
SoActionMethodList::setUp(void)
{
  PRIVATE(this)->lock();
  if (PRIVATE(this)->setupnumtypes != SoType::getNumTypes()) {
    int i, n;

    this->truncate(0); // clear action method list

    // first set all methods that have been set directly through SO_ACTION_ADD_METHOD()
    n = PRIVATE(this)->addedtypes.getLength();
    for (i = 0; i < n; i++) {
      (*this)[SoNode::getActionMethodIndex(PRIVATE(this)->addedtypes[i])] = PRIVATE(this)->addedmethods[i];
    }
    
    // make sure SoNode's action method is set to avoid a NULL action method
    i = SoNode::getActionMethodIndex(SoNode::getClassTypeId());
    if ((*this)[i] == NULL) {
      if (PRIVATE(this)->parent == NULL) {
        (*this)[i] = SoAction::nullAction;
      }
      else {
        // set to a dummy method to detect unset methods in the final pass
        (*this)[i] = unsetActionMethod;
      }
    }

    // for node types with no action method, inherit from parent nodetype(s)
    SoTypeList allnodes;
    SoType::getAllDerivedFrom(SoNode::getClassTypeId(), allnodes);
    n = allnodes.getLength();

    for (i = 0; i < n; i++) {
      SoType type = allnodes[i];
      int idx = SoNode::getActionMethodIndex(type);
      SoActionMethod m = (*this)[idx];
      if (m == NULL) {
        do {
          type = type.getParent();
          m = (*this)[SoNode::getActionMethodIndex(type)];
        } while (m == NULL);
        (*this)[idx] = m;
      }
    }

    // inherit unset methods from parent action
    if (PRIVATE(this)->parent != NULL) {
      PRIVATE(this)->parent->setUp();
      n = this->getLength();
      for (i = 0; i < n; i++) {
        if ((*this)[i] == unsetActionMethod) {
          (*this)[i] = (*PRIVATE(this)->parent)[i];
        }
      }
    }
    // used to detect when a new node has been added
    PRIVATE(this)->setupnumtypes = SoType::getNumTypes();
  }
  PRIVATE(this)->unlock();
}

#undef PRIVATE
