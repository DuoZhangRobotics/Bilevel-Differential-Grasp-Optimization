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
  \class SoTextureCoordinate2 SoTextureCoordinate2.h Inventor/nodes/SoTextureCoordinate2.h
  \brief The SoTextureCoordinate2 class contains a set of coordinates for the mapping of 2D textures.
  \ingroup nodes

  When encountering nodes of this type during traversal, the
  coordinates it contains will be put on the state stack. Some shape
  nodes (for instance SoIndexedFaceSet, among many others) can then
  use these coordinates for explicit, detailed control of how textures
  are mapped to it's surfaces.

  (If texturemapping is used without any SoTextureCoordinate2 nodes in
  the scenegraph leading up to a shape node, all shape types have
  default fallbacks. So SoTextureCoordinate2 nodes are only necessary
  to use if you are not satisfied with the default mapping.)

  Note that an SoTextureCoordinate2 node will \e replace the
  coordinates already present in the state (if any).

  Here's a very simple example (in Inventor scenegraph file format --
  mapping it to sourcecode is straightforward) that shows how to set
  up two quadratic polygons, one mapped 1:1 to the texture, the other
  using only the upper left quarter of the texture:

\code

Separator {
   Texture2 {
      image 6 8 3
      0x00ff0000 0x00ff0000 0x000000ff 0x000000ff 0x00ff00ff 0x00ff00ff
      0x00ff0000 0x00ff0000 0x000000ff 0x000000ff 0x00ff00ff 0x00ff00ff
      0x00ff0000 0x00ff0000 0x000000ff 0x000000ff 0x00ff00ff 0x00ff00ff
      0x0000ff00 0x0000ff00 0x0000ffff 0x0000ffff 0x0000ff00 0x0000ff00
      0x0000ff00 0x0000ff00 0x0000ffff 0x0000ffff 0x0000ff00 0x0000ff00
      0x00ffff00 0x00ffff00 0x000000ff 0x000000ff 0x00ffffff 0x00ffffff
      0x00ffff00 0x00ffff00 0x000000ff 0x000000ff 0x00ffffff 0x00ffffff
      0x00ffff00 0x00ffff00 0x000000ff 0x000000ff 0x00ffffff 0x00ffffff
   }

   Coordinate3 { point [ -1 -1 0, 1 -1 0, 1 1 0, -1 1 0 ] }

   # "1:1 mapping" to actual texture appearance. (Note that Y goes
   # from bottom to top, versus the common way of specifying bitmap
   # data from top to bottom.)
   TextureCoordinate2 { point [ 0 1, 1 1, 1 0, 0 0 ] }

   IndexedFaceSet {
      coordIndex [ 0, 1, 2, 3, -1 ]
      textureCoordIndex [ 0, 1, 2, 3, -1 ]
   }

   Translation { translation +4 0 0 }

   # Top left corner.
   TextureCoordinate2 { point [ 0 0.5, 0.5 0.5, 0.5 0, 0 0 ] }

   IndexedFaceSet {
      coordIndex [ 0, 1, 2, 3, -1 ]
      textureCoordIndex [ 0, 1, 2, 3, -1 ]
   }
}

\endcode

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    TextureCoordinate2 {
        point [  ]
    }
  \endcode

  \sa SoTextureCoordinateFunction, SoTextureCoordinateBinding
*/

// *************************************************************************

#include <Inventor/nodes/SoTextureCoordinate2.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/elements/SoTextureUnitElement.h>
#include <Inventor/elements/SoGLVBOElement.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/C/glue/gl.h>

#include "nodes/SoSubNodeP.h"
#include "misc/SoVBO.h"

// *************************************************************************

/*!
  \var SoMFVec2f SoTextureCoordinate2::point

  The set of 2D texture coordinates. Default value of field is an
  empty set.

  Texture coordinates are usually specified in normalized coordinates,
  ie in the range [0, 1]. (0, 0) is the lower left corner, while
  (1, 1) is the upper right corner of the texture image. Coordinates
  outside the [0, 1] range can be used to repeat the texture across a
  surface.

  \sa SoTexure2::wrapS, SoTexture2::wrapT
*/

// *************************************************************************

class SoTextureCoordinate2P {
 public:
  SoTextureCoordinate2P() : vbo(NULL) { }
  ~SoTextureCoordinate2P() { delete this->vbo; }
  SoVBO * vbo;
};

#define PRIVATE(obj) obj->pimpl

SO_NODE_SOURCE(SoTextureCoordinate2);

/*!
  Constructor.
*/
SoTextureCoordinate2::SoTextureCoordinate2(void)
{
  PRIVATE(this) = new SoTextureCoordinate2P;
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTextureCoordinate2);
  SO_NODE_ADD_FIELD(point, (NULL));
}

/*!
  Destructor.
*/
SoTextureCoordinate2::~SoTextureCoordinate2()
{
  delete PRIVATE(this);
}

// Documented in superclass.
void
SoTextureCoordinate2::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTextureCoordinate2, SO_FROM_INVENTOR_1|SoNode::VRML1);

  SO_ENABLE(SoGLRenderAction, SoGLTextureCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoTextureCoordinateElement);
  SO_ENABLE(SoGLRenderAction, SoGLMultiTextureCoordinateElement);
  SO_ENABLE(SoCallbackAction, SoMultiTextureCoordinateElement);

  SO_ENABLE(SoPickAction, SoTextureCoordinateElement);
  SO_ENABLE(SoPickAction, SoMultiTextureCoordinateElement);
}

// Documented in superclass.
void
SoTextureCoordinate2::doAction(SoAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);
  
  if (unit == 0) {
    SoTextureCoordinateElement::set2(action->getState(), this,
                                     this->point.getNum(),
                                     this->point.getValues(0));
  }
  else {
    SoMultiTextureCoordinateElement::set2(action->getState(), this, unit,
                                          this->point.getNum(),
                                          this->point.getValues(0));
  }
}

// Documented in superclass.
void
SoTextureCoordinate2::GLRender(SoGLRenderAction * action)
{
  SoState * state = action->getState();
  int unit = SoTextureUnitElement::get(state);

  if (unit == 0) {
    SoGLTextureCoordinateElement::setTexGen(action->getState(), this, NULL);
    SoTextureCoordinate2::doAction((SoAction *)action);
  }
  else {
    const cc_glglue * glue = cc_glglue_instance(SoGLCacheContextElement::get(state));
    int maxunits = cc_glglue_max_texture_units(glue);

    if (unit < maxunits) {
      SoGLMultiTextureCoordinateElement::setTexGen(action->getState(), this, unit, NULL);
      SoMultiTextureCoordinateElement::set2(action->getState(), this, unit,
                                            point.getNum(),
                                            point.getValues(0));
    }
  }

  SoBase::staticDataLock();
  const int num = this->point.getNum();
  SbBool setvbo = FALSE;
  if (SoGLVBOElement::shouldCreateVBO(state, num)) {
    setvbo = TRUE;
    SbBool dirty = FALSE;
    if (PRIVATE(this)->vbo == NULL) {
      PRIVATE(this)->vbo = new SoVBO(GL_ARRAY_BUFFER, GL_STATIC_DRAW); 
      dirty =  TRUE;
    }
    else if (PRIVATE(this)->vbo->getBufferDataId() != this->getNodeId()) {
      dirty = TRUE;
    }
    if (dirty) {
      PRIVATE(this)->vbo->setBufferData(this->point.getValues(0),
                                        num*sizeof(SbVec2f),
                                        this->getNodeId());
    }
  }
  else if (PRIVATE(this)->vbo && PRIVATE(this)->vbo->getBufferDataId()) {
    // clear buffers to deallocate VBO memory
    PRIVATE(this)->vbo->setBufferData(NULL, 0, 0);
  }
  SoBase::staticDataUnlock();
  SoGLVBOElement::setTexCoordVBO(state, 0, setvbo ? PRIVATE(this)->vbo : NULL);
}

// Documented in superclass.
void
SoTextureCoordinate2::callback(SoCallbackAction * action)
{
  SoTextureCoordinate2::doAction((SoAction *)action);
}

// Documented in superclass.
void
SoTextureCoordinate2::pick(SoPickAction * action)
{
  SoTextureCoordinate2::doAction((SoAction *)action);
}

#undef PRIVATE
