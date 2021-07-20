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
  \class SoTransform SoTransform.h Inventor/nodes/SoTransform.h
  \brief The SoTransform class is the "all-purpose" transformation node type.
  \ingroup nodes

  Like SoMatrixTransform, nodes of this type gives the application
  programmer maximum flexibility when specifying geometry
  transformations in a scene graph. If you want to set and keep the
  various components of the transformation matrix in separate
  entities, this node type is preferable, though.

  The order of operations is: first scaling is done, then rotation,
  then translation.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    Transform {
        translation 0 0 0
        rotation 0 0 1  0
        scaleFactor 1 1 1
        scaleOrientation 0 0 1  0
        center 0 0 0
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoTransform.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetMatrixAction.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoSFVec3f SoTransform::translation
  The translation part of the transformation.
*/
/*!
  \var SoSFRotation SoTransform::rotation
  The rotation part of the transformation.


  Note that there is one \e very common mistake that is easy to make
  when setting the value of a an SoSFRotation field, and that is to
  inadvertently use the wrong SbRotation constructor. This example
  should clarify the problem:

  \code
  mytransformnode->rotation.setValue(0, 0, 1, 1.5707963f);
  \endcode

  The programmer clearly tries to set a PI/2 rotation around the Z
  axis, but this will fail, as the SbRotation constructor invoked
  above is the one that takes as arguments the 4 floats of a \e
  quaternion. What the programmer almost certainly wanted to do was to
  use the SbRotation constructor that takes a rotation vector and a
  rotation angle, which is invoked like this:

  \code
  mytransformnode->rotation.setValue(SbVec3f(0, 0, 1), 1.5707963f);
  \endcode
*/
/*!
  \var SoSFVec3f SoTransform::scaleFactor
  The scaling part of the transformation.
*/
/*!
  \var SoSFRotation SoTransform::scaleOrientation
  The orientation the object is set to before scaling.
*/
/*!
  \var SoSFVec3f SoTransform::center
  The center point for the rotation.
*/

// *************************************************************************

SO_NODE_SOURCE(SoTransform);

/*!
  Constructor.
*/
SoTransform::SoTransform(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTransform);

  SO_NODE_ADD_FIELD(translation, (0.0f, 0.0f, 0.0f));
  SO_NODE_ADD_FIELD(rotation, (SbRotation(SbVec3f(0.0f, 0.0f, 1.0f), 0.0f)));
  SO_NODE_ADD_FIELD(scaleFactor, (1.0f, 1.0f, 1.0f));
  SO_NODE_ADD_FIELD(scaleOrientation, (SbRotation(SbVec3f(0.0f, 0.0f, 1.0f), 0.0f)));
  SO_NODE_ADD_FIELD(center, (0.0f, 0.0f, 0.0f));
}

/*!
  Destructor.
*/
SoTransform::~SoTransform()
{
}

// Doc from superclass.
void
SoTransform::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTransform, SO_FROM_INVENTOR_1|SoNode::VRML1);
} 

/*!
  Sets the transformation to translate to \a frompoint, with a rotation
  so that the (0,0,-1) vector is rotated into the vector from \a frompoint
  to \a topoint.
*/
void
SoTransform::pointAt(const SbVec3f & frompoint, const SbVec3f & topoint)
{
  this->scaleFactor = SbVec3f(1.0f, 1.0f, 1.0f);
  this->center = SbVec3f(0.0f, 0.0f, 0.0f);
  this->scaleOrientation = SbRotation(SbVec3f(0.0f, 0.0f, 1.0f), 0.0f);
  
  this->translation = frompoint;  
  SbVec3f dir = topoint - frompoint;
  if (dir.normalize() != 0.0f) {
    SbRotation rot(SbVec3f(0.0f, 0.0f, -1.0f), dir);
    this->rotation = rot;
  }
#if COIN_DEBUG
  else {
    SoDebugError::postWarning("SoTransform::pointAt",
                              "frompt == topoint");

  }
#endif // COIN_DEBUG
}

/*!
  Calculates the matrices to/from scale space.
*/
void
SoTransform::getScaleSpaceMatrix(SbMatrix & mat, SbMatrix & inv) const
{  
  SbMatrix tmp;
  mat.setTranslate(-center.getValue());
  tmp.setRotate(scaleOrientation.getValue().inverse());
  mat.multRight(tmp);
  tmp.setScale(scaleFactor.getValue());
  mat.multRight(tmp);
  inv = mat.inverse();
}

/*!
  Calculates the matrices to/from rotation space.
*/
void
SoTransform::getRotationSpaceMatrix(SbMatrix & mat, SbMatrix & inv) const
{
  SbMatrix tmp;
  mat.setTranslate(-this->center.getValue());
  tmp.setRotate(this->scaleOrientation.getValue().inverse());
  mat.multRight(tmp);
  tmp.setScale(this->scaleFactor.getValue());
  mat.multRight(tmp);
  tmp.setRotate(this->scaleOrientation.getValue());
  mat.multRight(tmp);
  tmp.setRotate(this->rotation.getValue());
  mat.multRight(tmp);
  inv = mat.inverse();
}

/*!
  Calculates the matrices to/from translation space.
*/
void
SoTransform::getTranslationSpaceMatrix(SbMatrix & mat, SbMatrix & inv) const
{
  SbMatrix tmp;
  mat.setTranslate(-this->center.getValue());
  tmp.setRotate(this->scaleOrientation.getValue().inverse());
  mat.multRight(tmp);
  tmp.setScale(this->scaleFactor.getValue());
  mat.multRight(tmp);
  tmp.setRotate(this->scaleOrientation.getValue());
  mat.multRight(tmp);
  tmp.setRotate(this->rotation.getValue());
  mat.multRight(tmp);
  tmp.setTranslate(this->translation.getValue());
  mat.multRight(tmp);
  inv = mat.inverse();
}

/*!
  Premultiplies this transformation by \a mat.
*/
void
SoTransform::multLeft(const SbMatrix & mat)
{
  SbMatrix matrix;
  matrix.setTransform(this->translation.getValue(),
                      this->rotation.getValue(),
                      this->scaleFactor.getValue(),
                      this->scaleOrientation.getValue(),
                      this->center.getValue());

  matrix.multLeft(mat);
  this->setMatrix(matrix);
}

/*!
  Postmultiplies this transformation by \a mat.
*/
void
SoTransform::multRight(const SbMatrix & mat)
{
  SbMatrix matrix;
  matrix.setTransform(this->translation.getValue(),
                      this->rotation.getValue(),
                      this->scaleFactor.getValue(),
                      this->scaleOrientation.getValue(),
                      this->center.getValue());
  matrix.multRight(mat);
  this->setMatrix(matrix);
}

/*!
  Premultiplies this transformation by the transformation in \a nodeonright.
*/
void
SoTransform::combineLeft(SoTransformation * nodeonright)
{
  SoGetMatrixAction ma(SbViewportRegion(100,100));
  ma.apply(nodeonright);
  this->multLeft(ma.getMatrix());
}

/*!
  Postmultiplies this transformation by the transformation in \a nodeonleft.
*/
void
SoTransform::combineRight(SoTransformation * nodeonleft)
{
  SoGetMatrixAction ma(SbViewportRegion(100,100));
  ma.apply(nodeonleft);
  this->multRight(ma.getMatrix());
}

/*!
  Sets the fields to create a transformation equal to \a mat.
*/
void
SoTransform::setMatrix(const SbMatrix & mat)
{  
  SbVec3f t,s,c = this->center.getValue();
  SbRotation r, so;
  mat.getTransform(t,r,s,so,c);
  
  this->translation = t;
  this->rotation = r;
  this->scaleFactor = s;
  this->scaleOrientation = so;
}

/*!  
  Sets the \e center field to \a newcenter. This might affect one
  or more of the other fields.  
*/
void
SoTransform::recenter(const SbVec3f & newcenter)
{
  SbMatrix matrix;
  matrix.setTransform(this->translation.getValue(),
                      this->rotation.getValue(),
                      this->scaleFactor.getValue(),
                      this->scaleOrientation.getValue(),
                      this->center.getValue());
  SbVec3f t,s;
  SbRotation r, so;
  matrix.getTransform(t, r, s, so, newcenter);
  this->translation = t;
  this->rotation = r;
  this->scaleFactor = s;
  this->scaleOrientation = so;
  this->center = newcenter;
}


// Doc from superclass.
void
SoTransform::doAction(SoAction * action)
{
  SbMatrix matrix;
  matrix.setTransform(this->translation.getValue(),
                      this->rotation.getValue(),
                      this->scaleFactor.getValue(),
                      this->scaleOrientation.getValue(),
                      this->center.getValue());
  
  SoModelMatrixElement::mult(action->getState(), this, matrix);
}

// Doc from superclass.
void
SoTransform::GLRender(SoGLRenderAction * action)
{
  SoTransform::doAction((SoAction *)action);
}

// Doc from superclass.
void
SoTransform::getBoundingBox(SoGetBoundingBoxAction * action)
{
  SoTransform::doAction((SoAction *)action);
}

// Doc from superclass.
void
SoTransform::getMatrix(SoGetMatrixAction * action)
{
  SbMatrix m;
  m.setTransform(this->translation.getValue(),
                 this->rotation.getValue(),
                 this->scaleFactor.getValue(),
                 this->scaleOrientation.getValue(),
                 this->center.getValue());
  action->getMatrix().multLeft(m);
  
  SbMatrix mi = m.inverse();
  action->getInverse().multRight(mi);
}

// Doc from superclass.
void
SoTransform::callback(SoCallbackAction * action)
{
  SoTransform::doAction((SoAction *)action);
}

// Doc from superclass.
void
SoTransform::pick(SoPickAction * action)
{
  SoTransform::doAction((SoAction *)action);
}

// Doc from superclass. Overrides the traversal method in this class for
// the SoGetPrimitiveCountAction because the number of primitives can
// be different depending on scene location (and thereby distance to
// camera) if there are e.g. SoLOD nodes in the scene.
void
SoTransform::getPrimitiveCount(SoGetPrimitiveCountAction * action)
{
  SoTransform::doAction((SoAction *)action);
}
