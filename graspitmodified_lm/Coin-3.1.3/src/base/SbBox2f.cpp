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
  \class SbBox2f SbBox.h Inventor/SbBox.h
  \brief The SbBox2f class is a 2 dimensional box with floating
  point corner coordinates.
  \ingroup base

  This box class is used by many other classes in Coin for data
  exchange and storage. It provides two box corners with floating
  point coordinates, which is among other things useful for
  representing screen or canvas dimensions in normalized coordinates.

  \sa SbBox2s, SbBox2d, SbBox3s, SbBox3f, SbBox3d, SbXfBox3f.
*/

// *************************************************************************

#include <Inventor/SbBox2f.h>

#include <limits>

#include <Inventor/SbBox2d.h>
#include <Inventor/SbBox2s.h>
#include <Inventor/SbBox2i32.h>
#include <Inventor/errors/SoDebugError.h>

// *************************************************************************

/*!
  \fn SbBox2f::SbBox2f(void)

  The default constructor makes an empty box.
*/

/*!
  \fn SbBox2f::SbBox2f(float xmin, float ymin, float xmax, float ymax)

  Constructs a box with the given corners.

  \a xmin should be less than \a xmax and \a ymin should be less than
  \a ymax if you want to make a valid box.
*/

/*!
  \fn SbBox2f::SbBox2f(const SbVec2f & min, const SbVec2f & max)

  Constructs a box with the given lower left and upper right corners.

  The coordinates of \a min should be less than the coordinates of
  \a max if you want to make a valid box.
*/

/*!
  \fn SbBox2f & SbBox2f::setBounds(float xmin, float ymin, float xmax, float ymax)
  Reset the boundaries of the box.

  \a xmin should be less than \a xmax and \a ymin should be less than
  \a ymax if you want to make a valid box.

  Returns reference to self.

  \sa getBounds().
*/

/*!
  \fn SbBox2f & SbBox2f::setBounds(const SbVec2f & min, const SbVec2f & max)

  Reset the boundaries of the box with the given corners.

  The coordinates of \a min should be less than the coordinates of
  \a max if you want to make a valid box.

  Returns reference to self.

  \sa getBounds().
*/

/*!
  Reset the boundaries of the box to the boundaries of the given \a box.

  Returns reference to self.

  \sa getBounds()
*/
SbBox2f &
SbBox2f::setBounds(const SbBox2d & box)
{
  if (box.isEmpty()) {
    makeEmpty();
  } else {
    minpt.setValue(box.getMin());
    maxpt.setValue(box.getMax());
  }
  return *this;
}

/*!
  Reset the boundaries of the box to the boundaries of the given \a box.

  Returns reference to self.

  \sa getBounds()
*/
SbBox2f &
SbBox2f::setBounds(const SbBox2s & box)
{
  if (box.isEmpty()) {
    makeEmpty();
  } else {
    minpt.setValue(box.getMin());
    maxpt.setValue(box.getMax());
  }
  return *this;
}

/*!
  Reset the boundaries of the box to the boundaries of the given \a box.

  Returns reference to self.

  \sa getBounds()
*/
SbBox2f &
SbBox2f::setBounds(const SbBox2i32 & box)
{
  if (box.isEmpty()) {
    makeEmpty();
  } else {
    minpt.setValue(box.getMin());
    maxpt.setValue(box.getMax());
  }
  return *this;
}

/*!
  Marks this as an empty box.

  \sa isEmpty().
*/
void
SbBox2f::makeEmpty(void)
{
  minpt.setValue(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
  maxpt.setValue(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
}

/*!
  \fn SbBool SbBox2f::isEmpty(void) const

  Check if this has been marked as an empty box.

  \sa makeEmpty().
*/

/*!
  \fn SbBool SbBox2f::hasArea(void) const

  Check if the box has "positive" area, i.e. the lower left corner is
  actually lower and more to the left than the other corner point.
*/

/*!
  \fn const SbVec2f & SbBox2f::getMin(void) const

  Returns the lower left corner of the box.

  \sa getOrigin(), getMax().
*/

/*!
  \fn SbVec2f & SbBox2f::getMin(void)

  Returns a modifiable reference ot the lower left corner of the box.

  \sa getOrigin(), getMax().
*/

/*!
  \fn const SbVec2f & SbBox2f::getMax(void) const

  Returns the upper right corner of the box.

  \sa getMin().
*/

/*!
  \fn SbVec2f SbBox2f::getCenter(void) const

  Returns the center point of the box.
*/

/*!
  Extend the boundaries of the box by the given point, i.e. make the
  box fit around the \a point if it isn't already situated within it.
*/
void
SbBox2f::extendBy(const SbVec2f & point)
{
  // The explicit cast to float is done to humour the HPUX aCC
  // compiler, which will otherwise say ``Template deduction failed to
  // find a match for the call to 'SbMin'''. mortene.
  this->minpt.setValue(SbMin(static_cast<float>(point[0]), static_cast<float>(this->minpt[0])),
                       SbMin(static_cast<float>(point[1]), static_cast<float>(this->minpt[1])));
  this->maxpt.setValue(SbMax(static_cast<float>(point[0]), static_cast<float>(this->maxpt[0])),
                       SbMax(static_cast<float>(point[1]), static_cast<float>(this->maxpt[1])));
}

/*!
  Extend the boundaries of the box by the given \a box parameter. This
  is equal to calling the above method twice with the corner points.
*/
void
SbBox2f::extendBy(const SbBox2f & box)
{
  if (box.isEmpty()) { return; }

  this->extendBy(box.getMin());
  this->extendBy(box.getMax());
}

/*!
  Check if \a point lies within the boundaries of this box.
 */
SbBool
SbBox2f::intersect(const SbVec2f & point) const
{
  if ((point[0] >= this->minpt[0]) && (point[0] <= this->maxpt[0]) &&
     (point[1] >= this->minpt[1]) && (point[1] <= this->maxpt[1])) return TRUE;
  return FALSE;
}

/*!
  Check if \a box lies wholly or partly within the boundaries
  of this box.
 */
SbBool
SbBox2f::intersect(const SbBox2f & box) const
{
  if ((box.getMax()[0] < this->getMin()[0]) ||
     (box.getMax()[1] < this->getMin()[1]) ||
     (box.getMin()[0] > this->getMax()[0]) ||
     (box.getMin()[1] > this->getMax()[1])) return FALSE;
  return TRUE;
}

/*!
  Return the point on the box closest to the given point \a p.
 */
SbVec2f
SbBox2f::getClosestPoint(const SbVec2f & p) const
{
  SbVec2f closest = p;

  SbVec2f center = this->getCenter();
  float devx = closest[0] - center[0];
  float devy = closest[1] - center[1];
  float halfwidth = (maxpt[0] - minpt[0]) / 2.0f;
  float halfheight = (maxpt[1] - minpt[1]) / 2.0f;

  // Move point to be on the nearest line of the box.
  if (fabs(devx) > fabs(devy))
    closest[0] = center[0] + halfwidth * ((devx < 0.0f) ? -1.0f : 1.0f);
  else
    closest[1] = center[1] + halfheight * ((devy < 0.0f) ? -1.0f : 1.0f);

  // Clamp to be inside box.
  closest[0] = SbMin(SbMax(closest[0], this->minpt[0]), this->maxpt[0]);
  closest[1] = SbMin(SbMax(closest[1], this->minpt[1]), this->maxpt[1]);

  return closest;
}

/*!
  \fn void SbBox2f::getBounds(float & xmin, float & ymin, float & xmax, float & ymax) const

  Returns the box boundaries.

  \sa setBounds(), getMin(), getMax().
*/

/*!
  \fn void SbBox2f::getBounds(SbVec2f & min, SbVec2f & max) const

  Returns the box corner points.

  \sa setBounds(), getMin(), getMax().
*/

/*!
  \fn void SbBox2f::getOrigin(float & originX, float & originY) const

  Returns the coordinates of the box origin (i.e. the lower left corner).

  \sa getMin().
*/

/*!
  \fn void SbBox2f::getSize(float & sizeX, float & sizeY) const

  Returns width and height of box.
*/

/*!
  \fn float SbBox2f::getAspectRatio(void) const

  Returns aspect ratio of box, which is defined as box width divided on
  box height.
*/

/*!
  \fn int operator == (const SbBox2f & b1, const SbBox2f & b2)
  \relates SbBox2f

  Check \a b1 and \a b2 for equality.
*/

/*!
  \fn int operator != (const SbBox2f & b1, const SbBox2f & b2)
  \relates SbBox2f

  Check \a b1 and \a b2 for inequality.
*/
