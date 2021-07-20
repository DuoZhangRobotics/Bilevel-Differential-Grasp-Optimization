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
  \class SoTextureCoordinateBundle include/Inventor/bundles/SoTextureCoordinateBundle.h
  \brief The SoTextureCoordinateBundle class simplifies texture coordinate handling.
  \ingroup bundles

  It is unlikely that application programmers should need to know how
  to use this class, as it is mostly intended for internal use.
*/
// FIXME: document class better.


#include <Inventor/bundles/SoTextureCoordinateBundle.h>

#include <Inventor/misc/SoState.h>
#include <Inventor/elements/SoTextureImageElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoGLTexture3EnabledElement.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoMultiTextureEnabledElement.h>
#include <Inventor/elements/SoBumpMapElement.h>
#include <Inventor/elements/SoBumpMapCoordinateElement.h>

#include <Inventor/nodes/SoShape.h>
#include <Inventor/nodes/SoVertexProperty.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/caches/SoBoundingBoxCache.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Inventor/system/gl.h>

#include <cassert>
#include "SbBasicP.h"
#include "coindefs.h"

#define FLAG_FUNCTION           0x01
#define FLAG_NEEDCOORDS         0x02
#define FLAG_DEFAULT            0x04
#define FLAG_DIDPUSH            0x08
#define FLAG_3DTEXTURES         0x10
#define FLAG_DIDINITDEFAULT     0x20
#define FLAG_NEEDINDICES        0x40

/*!
  Constructor with \a action being the action applied to the node.
  The \a forRendering parameter must be \e TRUE if the bundle is to
  be used for sending texture coordinates to GL during rendering.
  The \a setUpDefault must be \e TRUE if default texture coordinates
  should be generated.
*/
SoTextureCoordinateBundle::
SoTextureCoordinateBundle(SoAction * const action,
                          const SbBool forRendering,
                          const SbBool setUpDefault)
  : SoBundle(action)
{
  this->flags = 0;
  //
  // return immediately if there is no texture
  //
  SbBool needinit =
    SoTextureEnabledElement::get(this->state) ||
    SoTexture3EnabledElement::get(this->state);

  SbBool glrender = forRendering || action->isOfType(SoGLRenderAction::getClassTypeId());

  const SbBool * multienabled = NULL;
  int multimax = 0;
  if (glrender) {
    multienabled =
      SoMultiTextureEnabledElement::getEnabledUnits(this->state, multimax);
  }
  SbBool bumpenabled = glrender && (SoBumpMapElement::get(this->state) != NULL);

  if (!needinit && !multienabled && !bumpenabled) return;

  // It is safe to assume that shapenode is of type SoShape, so we
  // cast to SoShape before doing any operations on the node.
  this->shapenode = coin_assert_cast<SoShape *>(action->getCurPathTail());

  this->coordElt = SoTextureCoordinateElement::getInstance(this->state);
  switch (this->coordElt->getType()) {
  case SoTextureCoordinateElement::DEFAULT:
    this->initDefault(action, forRendering);
    break;
  case SoTextureCoordinateElement::EXPLICIT:
    if (this->coordElt->getNum() > 0) {
      this->flags |= FLAG_NEEDCOORDS;
    }
    else {
      this->initDefault(action, forRendering);
    }
    break;
  case SoTextureCoordinateElement::FUNCTION:
    this->flags |= FLAG_FUNCTION;
    this->flags |= FLAG_NEEDCOORDS; // not automatically generated
    break;

  case SoTextureCoordinateElement::TEXGEN:
    // texcoord won't be needed. This will only happen during SoGLRenderAction,
    // when GL generates texture coorinates. Therefore, we will not set
    // the FLAG_NEEDCOORDS here.
    this->flags |= FLAG_FUNCTION;
    if (!forRendering) {
      this->flags |= FLAG_NEEDCOORDS;
    }
    break;
  default:
    assert(0 && "unknown CoordType");
    break;
  }

  this->glElt = NULL;
  if (glrender) {
    SbBool needindices = !this->isFunction();
    if (multienabled) {
      const SoMultiTextureCoordinateElement * melem =
        SoMultiTextureCoordinateElement::getInstance(state);
      for (int i = 1; i <= multimax; i++) {
        if (multienabled[i]) {
          if ((melem->getType(i) == SoTextureCoordinateElement::DEFAULT) ||
              ((melem->getType(state, i) == SoTextureCoordinateElement::EXPLICIT) &&
               (melem->getNum(i) == 0))) {
            this->initDefaultMulti(action, i);
          }
          else if (!needindices &&
                   (melem->getType(state, i) == SoTextureCoordinateElement::EXPLICIT)) {
            needindices = TRUE;
          }
        }
      }
    }
    if (!needindices && this->isFunction()) {
      // check if bump mapping needs texture coordinate indices
      if (bumpenabled &&
          (SoBumpMapCoordinateElement::getInstance(state)->getNum())) {
        needindices = TRUE;
      }
    }
    if (needindices && this->isFunction()) this->flags |= FLAG_NEEDINDICES;
    this->glElt = static_cast<const SoGLTextureCoordinateElement *>(this->coordElt);
    this->glElt->initMulti(action->getState());
  }
  if ((this->flags & FLAG_DEFAULT) && !setUpDefault) {
    // FIXME: I couldn't be bothered to support this yet. It is for picking
    // optimization only, I think. pederb, 20000218
  }
}

/*!
  Destructor.
*/
SoTextureCoordinateBundle::~SoTextureCoordinateBundle()
{
  if (this->flags & FLAG_DIDPUSH) this->state->pop();
}

/*!
  Returns \e TRUE if texture coordinates is needed during rendering.
*/
SbBool
SoTextureCoordinateBundle::needCoordinates() const
{
  return (this->flags & FLAG_NEEDCOORDS) != 0;
}

/*!
  Returns \e TRUE if a texture coordinate function should be used.
*/
SbBool
SoTextureCoordinateBundle::isFunction() const
{
  return (this->flags & FLAG_FUNCTION) != 0;
}

/*!

  Returns \e TRUE if isFunction() is \e TRUE, but the texture
  coordinate indices are needed either by bump mapping or by one of
  the other texture units.

  \since Coin 2.2
*/
SbBool
SoTextureCoordinateBundle::needIndices(void) const
{
  return (this->flags & FLAG_NEEDINDICES) != 0;
}


/*!
  Returns the texture coordinates based on \a point and \a normal.
  Should only be used if SoTextureCoordinateBundle::isFunction() is \a TRUE.
*/
const SbVec4f &
SoTextureCoordinateBundle::get(const SbVec3f &point, const SbVec3f &normal)
{
  assert(this->coordElt != NULL && (this->flags & FLAG_FUNCTION));
  if (this->flags & FLAG_DEFAULT) {
    SbVec3f pt;
    if (this->flags & FLAG_3DTEXTURES) {
      pt = point - this->defaultorigo;
      this->dummyInstance[2] = pt[2]/this->defaultsize[2];
    }
    else {
      pt.setValue(point[this->defaultdim0]-this->defaultorigo[0],
                  point[this->defaultdim1]-this->defaultorigo[1],
                  0.0f);
    }
    this->dummyInstance[0] = pt[0]/this->defaultsize[0];
    this->dummyInstance[1] = pt[1]/this->defaultsize[1];
    return this->dummyInstance;
  }
  else {
    return coordElt->get(point, normal);
  }
}

/*!
  Returns the texture coordinates at index \a index.
  Should only be used if SoTextureCoordinateBundle::isFunction() is \a FALSE.
*/
const SbVec4f &
SoTextureCoordinateBundle::get(const int index)
{
  assert(coordElt && !(this->flags & FLAG_FUNCTION));
  return coordElt->get4(index);
}


/*!
  \fn void SoTextureCoordinateBundle::send(const int index) const
  Send texture coordinates to GL. Should only be used if
  SoTextureCoordinateBundle::isFunction() is \a FALSE.
*/

/*!
  \fn void SoTextureCoordinateBundle::send(const int index, const SbVec3f &point, const SbVec3f &normal) const
  Convenience function that will make it transparent to the rendering
  code if ordinary texture coordinates or function texture coordinates
  are used.
*/


//
// initialize default texture coordinates for unit 0
//
void
SoTextureCoordinateBundle::initDefault(SoAction * const action,
                                       const SbBool COIN_UNUSED_ARG(forRendering))
{
  this->flags |= FLAG_NEEDCOORDS;
  this->flags |= FLAG_DEFAULT;
  this->flags |= FLAG_FUNCTION;

  if (!(this->flags & FLAG_DIDPUSH)) {
    this->state->push();
    this->flags |= FLAG_DIDPUSH;
  }
  // have coordelt generate the default texture coordinates using a
  // callback to this instance.
  SoTextureCoordinateElement::setFunction(this->state, this->shapenode,
                                          SoTextureCoordinateBundle::defaultCB,
                                          this);
  this->coordElt = SoTextureCoordinateElement::getInstance(this->state);

  if (!(this->flags & FLAG_DIDINITDEFAULT)) {
    this->initDefaultCallback(action);
  }
}

//
// initialize default texture coordinates for unit > 0
//
void
SoTextureCoordinateBundle::initDefaultMulti(SoAction * action, const int unit)
{
  this->flags |= FLAG_NEEDCOORDS;

  if (!(this->flags & FLAG_DIDPUSH)) {
    this->state->push();
    this->flags |= FLAG_DIDPUSH;
  }
  SoMultiTextureCoordinateElement::setFunction(this->state, this->shapenode, unit,
                                               SoTextureCoordinateBundle::defaultCBMulti,
                                               this);
  if (!(this->flags & FLAG_DIDINITDEFAULT)) {
    this->initDefaultCallback(action);
  }
}

//
// callback for default texture coordinates (for texture unit 0)
//
const SbVec4f &
SoTextureCoordinateBundle::defaultCB(void * userdata,
                                     const SbVec3f & point,
                                     const SbVec3f & normal)
{
  return static_cast<SoTextureCoordinateBundle *>(userdata)->get(point, normal);
}

//
// callback for default texture coordinates (for texture units > 0)
//
const SbVec4f &
SoTextureCoordinateBundle:: defaultCBMulti(void * userdata,
                                           const SbVec3f & point,
                                           const SbVec3f & COIN_UNUSED_ARG(normal))
{
  SoTextureCoordinateBundle * thisp = static_cast<SoTextureCoordinateBundle *>(userdata);

  SbVec3f pt;
  if (thisp->flags & FLAG_3DTEXTURES) {
    pt = point - thisp->defaultorigo;
    thisp->dummyInstance[2] = pt[2]/thisp->defaultsize[2];
  }
  else {
    pt.setValue(point[thisp->defaultdim0]-thisp->defaultorigo[0],
                point[thisp->defaultdim1]-thisp->defaultorigo[1],
                0.0f);
  }
  thisp->dummyInstance[0] = pt[0]/thisp->defaultsize[0];
  thisp->dummyInstance[1] = pt[1]/thisp->defaultsize[1];
  return thisp->dummyInstance;
}

//
// Set up stuff needed for default texture coordinate mapping callback
//
void
SoTextureCoordinateBundle::initDefaultCallback(SoAction * action)
{
  this->flags |= FLAG_DIDINITDEFAULT;
  //
  // calculate needed stuff for default mapping
  //
  SbBox3f box;
  SbVec3f center;
  // this could be very slow, but if you're looking for speed, default
  // texture coordinate mapping shouldn't be used. We might optimize this
  // by using a SoTextureCoordinateCache soon though. pederb, 20000218

  SoShape * shape = coin_assert_cast<SoShape *>(this->shapenode);
  const SoBoundingBoxCache * bboxcache = shape->getBoundingBoxCache();
  if (bboxcache && bboxcache->isValid(action->getState())) {
    box = bboxcache->getProjectedBox();
    if (bboxcache->isCenterSet()) center = bboxcache->getCenter();
    else center = box.getCenter();
  }
  else {
    shape->computeBBox(action, box, center);
  }

  // just use som default values if the shape bbox is empty
  SbVec3f size(1.0f, 1.0f, 1.0f);
  SbVec3f origo(0.f, 0.0f, 0.0f);
  if (!box.isEmpty()) {
    box.getSize(size[0], size[1], size[2]);
    origo = box.getMin();
  }

  // Map S,T,R to X,Y,Z for 3D texturing
  if (SoTexture3EnabledElement::get(this->state)) {
    this->flags |= FLAG_3DTEXTURES;
    this->defaultdim0 = 0;
    this->defaultdim1 = 1;

    this->defaultorigo[2] = origo[2];
    this->defaultsize[2] = size[2];
  }
  else { // 2D textures
    this->defaultsize[2] = 1.0f;
    this->flags &= ~FLAG_3DTEXTURES;
    // find the two biggest dimensions
    int smallest = 0;
    float smallval = size[0];
    if (size[1] < smallval) {
      smallest = 1;
      smallval = size[1];
    }
    if (size[2] < smallval) {
      smallest = 2;
    }

    this->defaultdim0 = (smallest + 1) % 3;
    this->defaultdim1 = (smallest + 2) % 3;

    if (size[this->defaultdim0] == size[this->defaultdim1]) {
      // FIXME: this is probably an OIV bug. The OIV man pages are not
      // clear on this point (surprise), but the VRML specification states
      // that if the two dimensions are equal, the ordering X>Y>Z should
      // be used.
#if 0 // the correct way to do it
      if (this->defaultdim0 > this->defaultdim1) {
        SbSwap(this->defaultdim0, this->defaultdim1);
      }
#else // the OIV way to do it.
      if (this->defaultdim0 < this->defaultdim1) {
        SbSwap(this->defaultdim0, this->defaultdim1);
      }
#endif // OIV compatibility fix
    }
    else if (size[this->defaultdim0] < size[this->defaultdim1]) {
      SbSwap(this->defaultdim0, this->defaultdim1);
    }
  }

  this->defaultorigo[0] = origo[this->defaultdim0];
  this->defaultorigo[1] = origo[this->defaultdim1];
  this->defaultsize[0] = size[this->defaultdim0];
  this->defaultsize[1] = size[this->defaultdim1];

  // if bbox is empty in one dimension we just want to set it to
  // 0.0. The size should be set to 1.0 to avoid division by zero.
  for (int i = 0; i < 3; i++) {
    if (this->defaultsize[i] <= 0.0f) {
      this->defaultsize[i] = 1.0f;
    }
  }

  this->dummyInstance[2] = 0.0f;
  this->dummyInstance[3] = 1.0f;

  assert(this->defaultsize[0] > 0.0f);
  assert(this->defaultsize[1] > 0.0f);
}

#undef FLAG_FUNCTION
#undef FLAG_NEEDCOORDS
#undef FLAG_DEFAULT
#undef FLAG_DIDPUSH
#undef FLAG_3DTEXTURES
#undef FLAG_DIDINITDEFAULT
#undef FLAG_NEEDINDICES
