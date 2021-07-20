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
  \class SoVRMLMaterial SoVRMLMaterial.h Inventor/VRMLnodes/SoVRMLMaterial.h
  \brief The SoVRMLMaterial class is used to assign a material to geometry.
  \ingroup VRMLnodes
  
  \WEB3DCOPYRIGHT

  \verbatim
  Material {
    exposedField SFFloat ambientIntensity  0.2         # [0,1]
    exposedField SFColor diffuseColor      0.8 0.8 0.8 # [0,1]
    exposedField SFColor emissiveColor     0 0 0       # [0,1]
    exposedField SFFloat shininess         0.2         # [0,1]
    exposedField SFColor specularColor     0 0 0       # [0,1]
    exposedField SFFloat transparency      0           # [0,1]
  }
  \endverbatim

  The Material node specifies surface material properties for
  associated geometry nodes and is used by the VRML lighting equations
  during rendering.  Subclause 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>),
  contains a detailed description of the VRML lighting model
  equations.  All of the fields in the Material node range from 0.0 to
  1.0.  The fields in the Material node determine how light reflects
  off an object to create colour:

  - The ambientIntensity field specifies how much ambient light from
    light sources this surface shall reflect. Ambient light is
    omnidirectional and depends only on the number of light sources, not
    their positions with respect to the surface. Ambient colour is
    calculated as ambientIntensity � diffuseColor.
 
  - The diffuseColor field reflects all VRML light sources depending
    on the angle of the surface with respect to the light source. The more
    directly the surface faces the light, the more diffuse light reflects.
 
  - The emissiveColor field models "glowing" objects.
    This can be useful for displaying pre-lit models (where the light energy
    of the room is computed explicitly), or for displaying scientific data.
 
  - The specularColor and shininess fields determine the
    specular highlights (e.g., the shiny spots on an apple). When the
    angle from the light to the surface is close to the angle from the surface
    to the viewer, the specularColor is added to the diffuse and ambient
    colour calculations. Lower shininess values produce soft glows, while higher
    values result in sharper, smaller highlights.
 
  - The transparency field specifies how "clear" an object
    is, with 1.0 being completely transparent, and 0.0 completely opaque.  

*/

/*!
  \var SoSFColor SoVRMLMaterial::diffuseColor
  The diffuse color component. Default value is (0.8, 0.8, 0.8).
*/

/*!
  \var SoSFFloat SoVRMLMaterial::ambientIntensity
  The ambient intensity. Default value is 0.2.
*/

/*!
  \var SoSFColor SoVRMLMaterial::specularColor
  The specular color component. Default value is (0, 0, 0).
*/

/*!
  \var SoSFColor SoVRMLMaterial::emissiveColor
  The emissive color component.  Default value is (0, 0, 0).
*/

/*!
  \var SoSFFloat SoVRMLMaterial::shininess
  The shininess value. A number between 0 and 1. Default value is 0.2.
*/

/*!
  \var SoSFFloat SoVRMLMaterial::transparency
  The material transparency. Valid range is from 0.0 (completely opaque) to 1.0 (completely transparent).
  Default value is 0.0.
*/


#include <Inventor/VRMLnodes/SoVRMLMaterial.h>

#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/elements/SoOverrideElement.h>
#include <Inventor/elements/SoShapeStyleElement.h>
#include <Inventor/elements/SoGLLazyElement.h>
#include <Inventor/errors/SoDebugError.h>
#ifdef COIN_THREADSAFE
#include <Inventor/threads/SbStorage.h>
#endif // COIN_THREADSAFE

#include "nodes/SoSubNodeP.h"

#ifndef DOXYGEN_SKIP_THIS

class SoVRMLMaterialP {
public:
  
  SoVRMLMaterialP() 
#ifdef COIN_THREADSAFE
    : colorpacker_storage(sizeof(void*), alloc_colorpacker, free_colorpacker)
#endif // COIN_THREADSAFE
  {}

  SbColor tmpambient;
  float tmptransparency;

#ifdef COIN_THREADSAFE
  SbStorage colorpacker_storage;
#else // COIN_THREADSAFE
  SoColorPacker single_colorpacker;
#endif // COIN_THREADSAFE
  
  SoColorPacker * getColorPacker(void) {
#ifdef COIN_THREADSAFE
    SoColorPacker ** cptr = (SoColorPacker**) this->colorpacker_storage.get();
    return * cptr;
#else // COIN_THREADSAFE
    return &this->single_colorpacker;
#endif // COIN_THREADSAFE
  }

#ifdef COIN_THREADSAFE
private:
  static void alloc_colorpacker(void * data) {
    SoColorPacker ** cptr = (SoColorPacker**) data;
    *cptr = new SoColorPacker;
  }
  static void free_colorpacker(void * data) {
    SoColorPacker ** cptr = (SoColorPacker**) data;
    delete *cptr;
  }
#endif // COIN_THREADSAFE
};

#endif // DOXYGEN_SKIP_THIS

#define PRIVATE(obj) ((obj)->pimpl)

SO_NODE_SOURCE(SoVRMLMaterial);

// Doc in parent
void
SoVRMLMaterial::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLMaterial, SO_VRML97_NODE_TYPE);
}

/*!
  Constructor.
*/
SoVRMLMaterial::SoVRMLMaterial(void)
{
  PRIVATE(this) = new SoVRMLMaterialP;

  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLMaterial);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(diffuseColor, (0.8f, 0.8f, 0.8f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(ambientIntensity, (0.2f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(specularColor, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(emissiveColor, (0.0f, 0.0f, 0.0f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(shininess, (0.2f));
  SO_VRMLNODE_ADD_EXPOSED_FIELD(transparency, (0.0f));
}

/*!
  Destructor.
*/
SoVRMLMaterial::~SoVRMLMaterial()
{
  delete PRIVATE(this);
}

// Doc in parent
void
SoVRMLMaterial::doAction(SoAction * action)
{
  SoState * state = action->getState();

  uint32_t bitmask = 0;
  uint32_t flags = SoOverrideElement::getFlags(state);

#define TEST_OVERRIDE(bit) ((SoOverrideElement::bit & flags) != 0)

  if (!this->diffuseColor.isIgnored() &&
      !TEST_OVERRIDE(AMBIENT_COLOR)) {
    PRIVATE(this)->tmpambient = this->diffuseColor.getValue();
    if (!this->ambientIntensity.isIgnored())
      PRIVATE(this)->tmpambient *= this->ambientIntensity.getValue();
    bitmask |= SoLazyElement::AMBIENT_MASK;
    if (this->isOverride()) {
      SoOverrideElement::setAmbientColorOverride(state, this, TRUE);
    }
  }
  if (!this->diffuseColor.isIgnored() &&
      !TEST_OVERRIDE(DIFFUSE_COLOR)) {
    // Note: the override flag bit values for diffuseColor and
    // transparency are equal (done like that to match SGI/TGS
    // Inventor behavior), so overriding one will also override the
    // other.
    bitmask |= SoLazyElement::DIFFUSE_MASK;
    if (this->isOverride()) {
      SoOverrideElement::setDiffuseColorOverride(state, this, TRUE);
    }
  }
  if (!this->emissiveColor.isIgnored() &&
      !TEST_OVERRIDE(EMISSIVE_COLOR)) {

    bitmask |= SoLazyElement::EMISSIVE_MASK;
    if (this->isOverride()) {
      SoOverrideElement::setEmissiveColorOverride(state, this, TRUE);
    }

  }
  if (!this->specularColor.isIgnored() &&
      !TEST_OVERRIDE(SPECULAR_COLOR)) {
    bitmask |= SoLazyElement::SPECULAR_MASK;
    if (this->isOverride()) {
      SoOverrideElement::setSpecularColorOverride(state, this, TRUE);
    }
  }
  if (!this->shininess.isIgnored() &&
      !TEST_OVERRIDE(SHININESS)) {
    bitmask |= SoLazyElement::SHININESS_MASK;
    if (this->isOverride()) {
      SoOverrideElement::setShininessOverride(state, this, TRUE);
    }
  }
  if (!this->transparency.isIgnored() &&
      !TEST_OVERRIDE(TRANSPARENCY)) {
    PRIVATE(this)->tmptransparency = this->transparency.getValue();
    bitmask |= SoLazyElement::TRANSPARENCY_MASK;
    // Note: the override flag bit values for diffuseColor and
    // transparency are equal (done like that to match SGI/TGS
    // Inventor behavior), so overriding one will also override the
    // other.
    if (this->isOverride()) {
      SoOverrideElement::setTransparencyOverride(state, this, TRUE);
    }
  }
#undef TEST_OVERRIDE

  if (bitmask) {
#if COIN_DEBUG
    if (bitmask & SoLazyElement::SHININESS_MASK) {
      static int didwarn = 0;
      if (!didwarn && (this->shininess.getValue() < 0.0f || this->shininess.getValue() > 1.0f)) {
        SoDebugError::postWarning("SoMaterial::GLRender",
                                  "Shininess out of range [0-1]. "
                                  "The shininess value will be clamped. "
                                  "This warning will be printed only once, but there might be more errors. "
                                  "You should check and fix your code and/or VRML exporter.");
        didwarn = 1;
      }
    }
#endif // COIN_DEBUG
    
    SoLazyElement::setMaterials(state, this, bitmask,
                                PRIVATE(this)->getColorPacker(),
                                &this->diffuseColor.getValue(), 1,
                                &PRIVATE(this)->tmptransparency, 1,
                                PRIVATE(this)->tmpambient,
                                this->emissiveColor.getValue(),
                                this->specularColor.getValue(),
                                SbClamp(this->shininess.getValue(), 0.0f, 1.0f),
                                PRIVATE(this)->tmptransparency > 0.0f);
  }
}

// Doc in parent
void
SoVRMLMaterial::GLRender(SoGLRenderAction * action)
{
  SoVRMLMaterial::doAction(action);
}

// Doc in parent
void
SoVRMLMaterial::callback(SoCallbackAction * action)
{
  SoVRMLMaterial::doAction(action);
}

#undef PRIVATE

#endif // HAVE_VRML97
