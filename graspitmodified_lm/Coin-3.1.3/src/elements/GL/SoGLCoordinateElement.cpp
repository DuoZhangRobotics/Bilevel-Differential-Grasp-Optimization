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
  \class SoGLCoordinateElement Inventor/elements/SoGLCoordinateElement.h
  \brief The SoGLCoordinateElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoGLCoordinateElement.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/system/gl.h>

#include <cassert>

SO_ELEMENT_SOURCE(SoGLCoordinateElement);

/*!
  This static method initializes static data for the SoGLCoordinateElement
  class.
*/

void
SoGLCoordinateElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoGLCoordinateElement, inherited);
}

/*!
  The destructor.
*/

SoGLCoordinateElement::~SoGLCoordinateElement(void)
{
}

/*!
  Send coordinates \a index to GL. Handles both 3D and 4D coordinates.
*/
void
SoGLCoordinateElement::send(const int index) const
{
  if (this->areCoords3D) glVertex3fv((const GLfloat*)(this->coords3D + index));
  else glVertex4fv((const GLfloat*)(this->coords4D + index));
}

//! FIXME: write doc.

//$ EXPORT INLINE
const SbVec3f *
SoGLCoordinateElement::getPtr3() const
{
  return this->coords3D;
}

//! FIXME: write doc.

//$ EXPORT INLINE
const SbVec4f *
SoGLCoordinateElement::getPtr4() const
{
  return this->coords4D;
}
