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
  \class SoShapeStyleElement Inventor/elements/SoShapeStyleElement.h
  \brief The SoShapeStyleElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoShapeStyleElement.h>

#include "SbBasicP.h"

#include <Inventor/elements/SoLazyElement.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoDrawStyleElement.h>
#include <Inventor/elements/SoComplexityTypeElement.h>

#include <coindefs.h> // COIN_OBSOLETED()
#include <cassert>


#define DELAYRENDER_MASK \
  (SoShapeStyleElement::BBOXCMPLX| \
   SoShapeStyleElement::INVISIBLE| \
   SoShapeStyleElement::ABORTCB|   \
   SoShapeStyleElement::BIGIMAGE|  \
   SoShapeStyleElement::BUMPMAP|   \
   SoShapeStyleElement::VERTEXARRAY)

#define TRANSPTYPE_MASK 0x0001f

SO_ELEMENT_SOURCE(SoShapeStyleElement);

/*!
  This static method initializes static data for the SoShapeStyleElement class.
*/

void
SoShapeStyleElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoShapeStyleElement, inherited);
}

/*!
  The destructor.
*/

SoShapeStyleElement::~SoShapeStyleElement()
{
}

//! FIXME: write doc.

void
SoShapeStyleElement::init(SoState * state)
{
  inherited::init(state);
  this->flags = LIGHTING;
}

//! FIXME: write doc.

void
SoShapeStyleElement::push(SoState * COIN_UNUSED_ARG(state))
{
  SoShapeStyleElement * prev = coin_assert_cast<SoShapeStyleElement *>(this->getNextInStack());
  this->flags = prev->flags;
}

//! FIXME: write doc.

void
SoShapeStyleElement::pop(SoState * state, const SoElement * prevTopElement)
{
  inherited::pop(state, prevTopElement);
}

//! FIXME: write doc.

SbBool
SoShapeStyleElement::matches(const SoElement * element) const
{
  const SoShapeStyleElement * elem =
    coin_assert_cast<const SoShapeStyleElement *>(element);
  return this->flags == elem->flags;
}

//! FIXME: write doc.

SoElement *
SoShapeStyleElement::copyMatchInfo(void) const
{
  SoShapeStyleElement * elem =
    static_cast<SoShapeStyleElement *>(this->getTypeId().createInstance());
  elem->flags = this->flags;
  return elem;
}

//! FIXME: write doc.

const SoShapeStyleElement *
SoShapeStyleElement::get(SoState * const state)
{
  return coin_assert_cast<const SoShapeStyleElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
}

//! FIXME: write doc.

void
SoShapeStyleElement::setDrawStyle(SoState * const state,
                                  const int32_t value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value == static_cast<int32_t>(SoDrawStyleElement::INVISIBLE)) {
    elem->flags |= INVISIBLE;
  }
  else {
    elem->flags &= ~INVISIBLE;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setComplexityType(SoState * const state,
                                       const int32_t value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value == static_cast<int32_t>(SoComplexityTypeElement::BOUNDING_BOX)) {
    elem->flags |= BBOXCMPLX;
  }
  else {
    elem->flags &= ~BBOXCMPLX;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setTransparencyType(SoState * const state,
                                         const int32_t value)
{
  SoShapeStyleElement * elem = getElement(state);

  elem->flags &= ~TRANSPTYPE_MASK;
  assert(value <= TRANSPTYPE_MASK);
  elem->flags |= (value & TRANSPTYPE_MASK);

  if ((value == int(SoGLRenderAction::SORTED_OBJECT_SORTED_TRIANGLE_BLEND)) ||
      (value == int(SoGLRenderAction::SORTED_OBJECT_SORTED_TRIANGLE_ADD))) {
    elem->flags |= TRANSP_SORTED_TRIANGLES;
  }
  else {
    elem->flags &= ~TRANSP_SORTED_TRIANGLES;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setTextureEnabled(SoState * const state,
                                       const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= TEXENABLED;
  }
  else {
    elem->flags &= ~TEXENABLED;
  }
}

/*!
  FIXME: write doc.

  \COIN_FUNCTION_EXTENSION

  \since Coin 2.0
*/
void
SoShapeStyleElement::setTexture3Enabled(SoState * const state,
                                       const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= TEX3ENABLED;
  }
  else {
    elem->flags &= ~TEX3ENABLED;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setTextureFunction(SoState * const state,
                                        const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= TEXFUNC;
  }
  else {
    elem->flags &= ~TEXFUNC;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setLightModel(SoState * const state,
                                   const int32_t value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value != static_cast<int32_t>(SoLazyElement::BASE_COLOR)) {
    elem->flags |= LIGHTING;
  }
  else {
    elem->flags &= ~LIGHTING;
  }
}

//! FIXME: write doc.

void
SoShapeStyleElement::setOverrides(SoState * const state,
                                  const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= OVERRIDE;
  }
  else {
    elem->flags &= ~OVERRIDE;
  }
}

//! FIXME: write doc.

SbBool
SoShapeStyleElement::isScreenDoor(SoState * const state)
{
  const SoShapeStyleElement * elem = getConstElement(state);
  return ((elem->flags & TRANSPTYPE_MASK) == SoGLRenderAction::SCREEN_DOOR);
}

/*!
  Returns the current transparency type.

  \COIN_FUNCTION_EXTENSION

  \since Coin 2.0
*/
int
SoShapeStyleElement::getTransparencyType(SoState * const state)
{
  const SoShapeStyleElement * elem = getConstElement(state);
  return static_cast<int>(elem->flags & TRANSPTYPE_MASK);
}

/*!
  FIXME: write doc.
*/

SbBool
SoShapeStyleElement::mightNotRender() const
{
  if ((this->flags & DELAYRENDER_MASK) != 0) return TRUE;
  return FALSE;
}

/*!
  FIXME: write doc.
*/

SbBool
SoShapeStyleElement::needNormals() const
{
  return (this->flags & LIGHTING) != 0;
}

/*!
  FIXME: write doc.
*/

SbBool
SoShapeStyleElement::needTexCoords(void) const
{
  return (this->flags&(TEXENABLED|TEX3ENABLED)) != 0;
}

/*!
  Not implemented in Coin. It is used by SoVertexProperty in SGI OIV.
*/
int
SoShapeStyleElement::getRenderCaseMask(void) const
{
  COIN_OBSOLETED();
  return 0;
}

/*!
  Returns if texture function is currently enabled.
*/
SbBool
SoShapeStyleElement::isTextureFunction(void) const
{
  return (this->flags&TEXFUNC) != 0;
}

/*!
  Returns the current modifiable instance (might cause a push())
*/
SoShapeStyleElement *
SoShapeStyleElement::getElement(SoState * const state)
{
  return coin_assert_cast<SoShapeStyleElement *>
    (
     SoElement::getElement(state, classStackIndex)
     );
}
/*!
  Returns the current read-only instance.
*/
const SoShapeStyleElement *
SoShapeStyleElement::getConstElement(SoState * const state)
{
  return coin_assert_cast<const SoShapeStyleElement *>
    (
     SoElement::getConstElement(state, classStackIndex)
     );
}

/*!
  Sets bumpmap enabled.

  \since Coin 2.4
*/
void
SoShapeStyleElement::setBumpmapEnabled(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= BUMPMAP;
  }
  else {
    elem->flags &= ~BUMPMAP;
  }
}

/*!
  Sets bigimage enabled.

  \since Coin 2.4
*/
void
SoShapeStyleElement::setBigImageEnabled(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= BIGIMAGE;
  }
  else {
    elem->flags &= ~BIGIMAGE;
  }
}

/*!
  Sets if vertex array rendering might be used.

  \since Coin 2.4
*/
void
SoShapeStyleElement::setVertexArrayRendering(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= VERTEXARRAY;
  }
  else {
    elem->flags &= ~VERTEXARRAY;
  }
}

/*!
  Sets material transparency.

  \since Coin 2.4
*/
void
SoShapeStyleElement::setTransparentMaterial(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= TRANSP_MATERIAL;
  }
  else {
    elem->flags &= ~TRANSP_MATERIAL;
  }
}

/*!
  Sets texture transparency.

  \since Coin 2.4
*/
void
SoShapeStyleElement::setTransparentTexture(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= TRANSP_TEXTURE;
  }
  else {
    elem->flags &= ~TRANSP_TEXTURE;
  }
}

/*!
  Sets whether we are rendering to a shadow (depth) map or not.

  \since Coin 2.5
*/
void
SoShapeStyleElement::setShadowMapRendering(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= SHADOWMAP;
  }
  else {
    elem->flags &= ~SHADOWMAP;
  }
}

/*!
  Sets whether we are rendering with shadows or not.

  \since Coin 2.5
*/
void
SoShapeStyleElement::setShadowsRendering(SoState * state, const SbBool value)
{
  SoShapeStyleElement * elem = getElement(state);
  if (value) {
    elem->flags |= SHADOWS;
  }
  else {
    elem->flags &= ~SHADOWS;
  }
}

/*!
  Returns the state flags. Used internally to optimize rendering.

  \ since Coin 2.4
*/
unsigned int
SoShapeStyleElement::getFlags(void) const
{
  return this->flags;
}


#undef DELAYRENDER_MASK
#undef TRANSPTYPE_MASK
