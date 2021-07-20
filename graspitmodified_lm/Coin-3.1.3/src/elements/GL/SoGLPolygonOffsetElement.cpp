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
  \class SoGLPolygonOffsetElement Inventor/elements/SoGLPolygonOffsetElement.h
  \brief The SoGLPolygonOffsetElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoGLPolygonOffsetElement.h>

#include <cassert>
#include <cstdlib>

#include <Inventor/C/tidbits.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/misc/SoGLDriverDatabase.h>

#include "misc/SoGL.h"

SO_ELEMENT_SOURCE(SoGLPolygonOffsetElement);

/*!
  This static method initializes static data for the
  SoGLPolygonOffsetElement class.
*/

void
SoGLPolygonOffsetElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoGLPolygonOffsetElement, inherited);
}

/*!
  The destructor.
*/

SoGLPolygonOffsetElement::~SoGLPolygonOffsetElement(void)
{
}

//! FIXME: write doc.

void
SoGLPolygonOffsetElement::init(SoState * stateptr)
{
  inherited::init(stateptr);
  this->state = stateptr;
}

//! FIXME: write doc.

void
SoGLPolygonOffsetElement::push(SoState * stateptr)
{
  SoGLPolygonOffsetElement * prev = (SoGLPolygonOffsetElement*)this->getNextInStack();

  this->style = prev->style;
  this->active = prev->active;
  this->offsetfactor = prev->offsetfactor;
  this->offsetunits = prev->offsetunits;
  this->state = stateptr;
  // capture previous element since we might or might not change the
  // GL state in set/pop
  prev->capture(stateptr);
}

//! FIXME: write doc.

void
SoGLPolygonOffsetElement::pop(SoState * stateptr, const SoElement * prevTopElement)
{
  const SoGLPolygonOffsetElement * prev =
    (const SoGLPolygonOffsetElement*)prevTopElement;

  if (this->style != prev->style ||
      this->active != prev->active ||
      this->offsetfactor != prev->offsetfactor ||
      this->offsetunits != prev->offsetunits) {
    this->updategl();
  }
}

//! FIXME: write doc.

void
SoGLPolygonOffsetElement::setElt(float factor, float units,
                                 Style styles, SbBool on)
{
  if (on != this->active ||
      styles != this->style ||
      factor != this->offsetfactor ||
      units != this->offsetunits) {
    this->active = on;
    this->style = styles;
    this->offsetfactor = factor;
    this->offsetunits = units;
    this->updategl();
  }
}


//! FIXME: write doc.

void
SoGLPolygonOffsetElement::updategl(void)
{
  const cc_glglue * w = sogl_glue_instance(this->state);
  if (!SoGLDriverDatabase::isSupported(w, SO_GL_POLYGON_OFFSET)) {
    static SbBool first = TRUE;
    if (first) {
      SoDebugError::postWarning("SoGLPolygonOffsetElement::updategl",
                                "OpenGL driver doesn't support z-buffer "
                                "offsetting");
      first = FALSE;
    }
    return;
  }


  if (this->active) {
    int styles = 0;
    styles |= (this->style & FILLED) ? cc_glglue_FILLED : 0;
    styles |= (this->style & LINES) ? cc_glglue_LINES : 0;
    styles |= (this->style & POINTS) ? cc_glglue_POINTS : 0;
    cc_glglue_glPolygonOffsetEnable(w, TRUE, styles);

    cc_glglue_glPolygonOffset(w, this->offsetfactor, this->offsetunits);
  }
  else { // ! active
    int all = cc_glglue_FILLED | cc_glglue_LINES | cc_glglue_POINTS;
    cc_glglue_glPolygonOffsetEnable(w, FALSE, all);
  }
}
