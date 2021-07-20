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
  \class SoVRMLPixelTexture SoVRMLPixelTexture.h Inventor/VRMLnodes/SoVRMLPixelTexture.h
  \brief The SoVRMLPixelTexture class is used for mapping a texture image onto geometry..
  \ingroup VRMLnodes

  \WEB3DCOPYRIGHT

  \verbatim
  PixelTexture {
    exposedField SFImage  image      0 0 0    # see SoSFImage
    field        SFBool   repeatS    TRUE
    field        SFBool   repeatT    TRUE
  }
  \endverbatim

  The PixelTexture node defines a 2D image-based texture map as an
  explicit array of pixel values (image field) and parameters
  controlling tiling repetition of the texture onto geometry.  Texture
  maps are defined in a 2D coordinate system (s, t) that ranges from
  0.0 to 1.0 in both directions. The bottom edge of the pixel image
  corresponds to the S-axis of the texture map, and left edge of the
  pixel image corresponds to the T-axis of the texture map. The
  lower-left pixel of the pixel image corresponds to s=0.0, t=0.0, and
  the top-right pixel of the image corresponds to s = 1.0, t = 1.0.
  See 4.6.11, Texture maps
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.6.11>),
  for a general description of texture
  maps. Figure 6.13 depicts an example PixelTexture.

  <center>
  <img src="http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/Images/PixelTexture.gif">
  Figure 6.13 -- PixelTexture node
  </center>

  See 4.14, Lighting model
  (<http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/part1/concepts.html#4.14>),
  for a description of how the texture values interact with the
  appearance of the geometry.  SoSFImage, describes the
  specification of an image.  The repeatS and repeatT fields specify
  how the texture wraps in the S and T directions. If repeatS is TRUE
  (the default), the texture map is repeated outside the 0-to-1
  texture coordinate range in the S direction so that it fills the
  shape. If repeatS is FALSE, the texture coordinates are clamped in
  the S direction to lie within the 0.0 to 1.0 range. The repeatT
  field is analogous to the repeatS field.

*/

/*!
  SoSFImage SoVRMLPixelTexture::image
  The image data.
*/

#include <Inventor/VRMLnodes/SoVRMLPixelTexture.h>

#include <assert.h>
#include <stddef.h>

#include <Inventor/VRMLnodes/SoVRMLMacros.h>
#include <Inventor/SoInput.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoGLTexture3EnabledElement.h>
#include <Inventor/elements/SoGLTextureImageElement.h>
#include <Inventor/elements/SoTextureQualityElement.h>
#include <Inventor/elements/SoTextureOverrideElement.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/errors/SoReadError.h>
#include <Inventor/misc/SoGLBigImage.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/SbImage.h>
#ifdef HAVE_THREADS
#include <Inventor/threads/SbMutex.h>
#endif // HAVE_THREADS

#include "nodes/SoSubNodeP.h"
#include "elements/SoTextureScalePolicyElement.h"

// *************************************************************************

class SoVRMLPixelTextureP {
public:
  SoGLImage * glimage;
  SbBool glimagevalid;
  int readstatus;

#ifdef COIN_THREADSAFE
  SbMutex glimagemutex;
  void lock_glimage(void) { this->glimagemutex.lock(); }
  void unlock_glimage(void) { this->glimagemutex.unlock(); }
#else // !COIN_THREADSAFE
  void lock_glimage(void) { }
  void unlock_glimage(void) { }
#endif // !COIN_THREADSAFE
};

// *************************************************************************

SO_NODE_SOURCE(SoVRMLPixelTexture);

// *************************************************************************

// Doc in parent
void
SoVRMLPixelTexture::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoVRMLPixelTexture, SO_VRML97_NODE_TYPE);

  SoType type = SoVRMLPixelTexture::getClassTypeId();
  SoRayPickAction::addMethod(type, SoNode::rayPickS);
}

#define PRIVATE(obj) ((obj)->pimpl)

/*!
  Constructor.
*/
SoVRMLPixelTexture::SoVRMLPixelTexture(void)
{
  PRIVATE(this) = new SoVRMLPixelTextureP;

  SO_VRMLNODE_INTERNAL_CONSTRUCTOR(SoVRMLPixelTexture);

  SO_VRMLNODE_ADD_EXPOSED_FIELD(image, (SbVec2s(0,0), 0, NULL));

  PRIVATE(this)->glimage = NULL;
  PRIVATE(this)->glimagevalid = FALSE;
  PRIVATE(this)->readstatus = 1;
}

/*!
  Destructor.
*/
SoVRMLPixelTexture::~SoVRMLPixelTexture()
{
  if (PRIVATE(this)->glimage) PRIVATE(this)->glimage->unref(NULL);
  delete PRIVATE(this);
}

static SoGLImage::Wrap
pixeltexture_translate_wrap(const SbBool repeat)
{
  if (repeat) return SoGLImage::REPEAT;
  return SoGLImage::CLAMP_TO_EDGE;
}

// Doc in parent
void
SoVRMLPixelTexture::doAction(SoAction * action)
{
  SoState * state = action->getState();

  if (SoTextureOverrideElement::getImageOverride(state))
    return;

  SoTexture3EnabledElement::set(state, this, FALSE);
  int nc;
  SbVec2s size;
  const unsigned char * bytes = this->image.getValue(size, nc);

  if (size == SbVec2s(0, 0)) {
    SoTextureEnabledElement::set(state, this, FALSE);    
  }
  else {
    SoTextureImageElement::set(state, this,
                               size, nc, bytes,
                               pixeltexture_translate_wrap(this->repeatS.getValue()),
                               pixeltexture_translate_wrap(this->repeatT.getValue()),
                               SoTextureImageElement::MODULATE,
                               SbColor(1.0f, 1.0f, 1.0f));
    SoTextureEnabledElement::set(state, this, TRUE);
  }
  if (this->isOverride()) {
    SoTextureOverrideElement::setImageOverride(state, TRUE);
  }
}

void 
SoVRMLPixelTexture::rayPick(SoRayPickAction * action)
{
  SoVRMLPixelTexture::doAction(action);
}

// Doc in parent
void
SoVRMLPixelTexture::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();

  if (SoTextureOverrideElement::getImageOverride(state))
    return;

  SoGLTexture3EnabledElement::set(state, this, FALSE);

  float quality = SoTextureQualityElement::get(state);

  PRIVATE(this)->lock_glimage();

  if (!PRIVATE(this)->glimagevalid) {
    int nc;
    SbVec2s size;
    const unsigned char * bytes =
      this->image.getValue(size, nc);
    SoTextureScalePolicyElement::Policy scalepolicy =
      SoTextureScalePolicyElement::get(state);      
    SbBool needbig = (scalepolicy == SoTextureScalePolicyElement::FRACTURE);

    if (needbig &&
        (PRIVATE(this)->glimage == NULL ||
         PRIVATE(this)->glimage->getTypeId() != SoGLBigImage::getClassTypeId())) {
      if (PRIVATE(this)->glimage) PRIVATE(this)->glimage->unref(state);
      PRIVATE(this)->glimage = new SoGLBigImage();
    }
    else if (!needbig &&
             (PRIVATE(this)->glimage == NULL ||
              PRIVATE(this)->glimage->getTypeId() != SoGLImage::getClassTypeId())) {
      if (PRIVATE(this)->glimage) PRIVATE(this)->glimage->unref(state);
      PRIVATE(this)->glimage = new SoGLImage();
    }

    if (scalepolicy == SoTextureScalePolicyElement::SCALE_DOWN) {
      PRIVATE(this)->glimage->setFlags(PRIVATE(this)->glimage->getFlags()|SoGLImage::SCALE_DOWN);
    }

    if (bytes && size != SbVec2s(0,0)) {
      PRIVATE(this)->glimage->setData(bytes, size, nc,
                             pixeltexture_translate_wrap(this->repeatS.getValue()),
                             pixeltexture_translate_wrap(this->repeatT.getValue()),
                             quality);
      PRIVATE(this)->glimagevalid = TRUE;
      // don't cache while creating a texture object
      SoCacheElement::setInvalid(TRUE);
      if (state->isCacheOpen()) {
        SoCacheElement::invalidate(state);
      }
    }
  }

  PRIVATE(this)->unlock_glimage();

  SoGLTextureImageElement::set(state, this,
                               PRIVATE(this)->glimagevalid ? PRIVATE(this)->glimage : NULL,
                               SoTextureImageElement::MODULATE,
                               SbColor(1.0f, 1.0f, 1.0f));

  SoGLTextureEnabledElement::set(state,
                                 this, PRIVATE(this)->glimagevalid &&
                                 quality > 0.0f);

  if (this->isOverride()) {
    SoTextureOverrideElement::setImageOverride(state, TRUE);
  }
}

// doc in parent
void
SoVRMLPixelTexture::callback(SoCallbackAction * action)
{
  SoVRMLPixelTexture::doAction(action);
}

// doc in parent
SbBool
SoVRMLPixelTexture::readInstance(SoInput * in,
                                 unsigned short flags)
{
  PRIVATE(this)->glimagevalid = FALSE;
  return inherited::readInstance(in, flags);
}

/*!
  Overloaded to detect when fields change.
*/
void
SoVRMLPixelTexture::notify(SoNotList * list)
{
  PRIVATE(this)->glimagevalid = FALSE;
  SoNode::notify(list);
}

#undef PRIVATE

#endif // HAVE_VRML97
