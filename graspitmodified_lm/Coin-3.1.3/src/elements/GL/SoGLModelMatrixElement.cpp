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
  \class SoGLModelMatrixElement Inventor/elements/SoGLModelMatrixElement.h
  \brief The SoGLModelMatrixElement class is yet to be documented.
  \ingroup elements

  FIXME: write doc.
*/

#include <Inventor/elements/SoGLModelMatrixElement.h>
#include <Inventor/elements/SoGLViewingMatrixElement.h>
#include <Inventor/SbRotation.h>
#include <Inventor/SbVec3f.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/system/gl.h>
#include <Inventor/C/tidbits.h>

#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

SO_ELEMENT_SOURCE(SoGLModelMatrixElement);

static int COIN_HANDLE_STACK_OVERFLOW = 0;

/*!
  This static method initializes static data for the
  SoGLModelMatrixElement class.
*/

void
SoGLModelMatrixElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoGLModelMatrixElement, inherited);

  const char * env = coin_getenv("COIN_HANDLE_STACK_OVERFLOW");
  if (env && atoi(env) > 0) COIN_HANDLE_STACK_OVERFLOW = 1;
  else COIN_HANDLE_STACK_OVERFLOW = 0;
}

/*!
  The destructor.
*/

SoGLModelMatrixElement::~SoGLModelMatrixElement(void)
{
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::init(SoState * stateptr)
{
  this->state = stateptr;
  this->viewEltNodeId = 0;
  this->stackoverflow = FALSE;
  inherited::init(stateptr);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::push(SoState * stateptr)
{
  SoGLModelMatrixElement * prev = (SoGLModelMatrixElement*)
    this->getNextInStack();

  // the stackoverflow test makes it possible to have scene graphs
  // with virtually unlimited depth and with transformations inside
  // each separator. If a GL_STACK_OVERFLOW error is encountered,
  // a glPopMatrix() will not be called in the pop() method, but
  // the GL matrix will be read from SoModelMatrixElement instead.
  //                                          pederb, 2000-12-20
  this->stackoverflow = prev->stackoverflow;
  this->state = prev->state;
  this->viewEltNodeId = prev->viewEltNodeId;

  if (COIN_HANDLE_STACK_OVERFLOW > 0) {
    if (!this->stackoverflow) {
      glPushMatrix();
      if (glGetError() == GL_STACK_OVERFLOW) {
        this->stackoverflow = TRUE;
      }
    }
  }
  else {
    glPushMatrix();
  }
  inherited::push(stateptr);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::pop(SoState * stateptr,
                            const SoElement * prevTopElement)
{
  inherited::pop(stateptr, prevTopElement);

  SoGLModelMatrixElement * prev = (SoGLModelMatrixElement*)
    prevTopElement;

  if (prev->stackoverflow) {
    SbMatrix mat = SoGLViewingMatrixElement::getResetMatrix(this->state);
    mat.multLeft(this->modelMatrix);
    glLoadMatrixf(mat[0]);
  }
  else {
    glPopMatrix();
  }
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::makeEltIdentity()
{
  SbMatrix mat = SoGLViewingMatrixElement::getResetMatrix(this->state);
  glLoadMatrixf(mat[0]);
  inherited::makeEltIdentity();
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::setElt(const SbMatrix & matrix)
{
  SbMatrix mat = SoGLViewingMatrixElement::getResetMatrix(this->state);
  mat.multLeft(matrix);
  glLoadMatrixf(mat[0]);
  inherited::setElt(matrix);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::multElt(const SbMatrix &matrix)
{
  glMultMatrixf(matrix[0]);
  inherited::multElt(matrix);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::translateEltBy(const SbVec3f &translation)
{
  glTranslatef(translation[0], translation[1], translation[2]);
  inherited::translateEltBy(translation);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::rotateEltBy(const SbRotation &rotation)
{
  SbVec3f axis;
  float angle;
  rotation.getValue(axis, angle);
  glRotatef(angle*180.0f/float(M_PI), axis[0], axis[1], axis[2]);
  inherited::rotateEltBy(rotation);
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::scaleEltBy(const SbVec3f &scaleFactor)
{
  glScalef(scaleFactor[0], scaleFactor[1], scaleFactor[2]);
  inherited::scaleEltBy(scaleFactor);
}

//! FIXME: write doc.

SbMatrix
SoGLModelMatrixElement::pushMatrixElt()
{
  this->viewEltNodeId = SoGLViewingMatrixElement::getNodeId(this->state);
  glPushMatrix();
  return inherited::pushMatrixElt();
}

//! FIXME: write doc.

void
SoGLModelMatrixElement::popMatrixElt(const SbMatrix &matrix)
{
  glPopMatrix();
  if (this->viewEltNodeId != SoGLViewingMatrixElement::getNodeId(this->state)) {
    this->setElt(matrix);
  }
  inherited::popMatrixElt(matrix);
}
