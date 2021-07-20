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
  \class SbBox2s SbBox.h Inventor/SbBox.h
  \brief The SbBox2s class is a 2 dimensional box with short
  integer coordinates.
  \ingroup base

  This box class is used by other classes in Coin for data
  exchange. It provides storage for two box corners with short integer
  coordinates, which is among other things useful for representing
  screen or canvas areas in absolute window coordinates.

  \sa SbBox2f, SbBox2d, SbBox3s, SbBox3f, SbBox3d, SbXfBox3f.
*/

// *************************************************************************

#include <Inventor/SbBox2s.h>

#include <limits>
#include <cassert>

#include <Inventor/SbBox2i32.h>
#include <Inventor/SbBox2f.h>
#include <Inventor/SbBox2d.h>
#include <Inventor/errors/SoDebugError.h>

// *************************************************************************

/*!
  \fn SbBox2s::SbBox2s(void)

  The default constructor makes an empty box.
*/

/*!
  \fn SbBox2s::SbBox2s(short xmin, short ymin, short xmax, short ymax)

  Constructs a box with the given corner coordinates.

  \a xmin should be less than \a xmax and \a ymin should be less than
  \a ymax if you want to make a valid box.
*/

/*!
  \fn SbBox2s::SbBox2s(const SbVec2s & boxmin, const SbVec2s & boxmax)

  Constructs a box with the given corners.

  The coordinates of \a min should be less than the coordinates of
  \a max if you want to make a valid box.
*/

/*!
  \fn SbBox2s & SbBox2s::setBounds(short xmin, short ymin, short xmax, short ymax)

  Reset the boundaries of the box.

  \a xmin should be less than \a xmax and \a ymin should be less than
  \a ymax if you want to make a valid box.

  Returns reference to self.

  \sa getBounds().
*/

/*!
  \fn SbBox2s & SbBox2s::setBounds(const SbVec2s & boxmin, const SbVec2s & boxmax)

  Reset the boundaries of the box with the given corners.

  The coordinates of \a min should be less than the coordinates of
  \a max if you want to make a valid box.

  Returns reference to self.

  \sa getBounds().
*/

/*!
  Reset the boundaries with the boundaries of the given \a box.

  Returns reference to self.

  \sa setBounds()
*/

SbBox2s &
SbBox2s::setBounds(const SbBox2i32 & box)
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
  Reset the boundaries with the boundaries of the given \a box.

  Returns reference to self.

  \sa setBounds()
*/

SbBox2s &
SbBox2s::setBounds(const SbBox2f & box)
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
  Reset the boundaries with the boundaries of the given \a box.

  Returns reference to self.

  \sa setBounds()
*/

SbBox2s &
SbBox2s::setBounds(const SbBox2d & box)
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
SbBox2s::makeEmpty(void)
{
  minpt.setValue(std::numeric_limits<short>::max(), std::numeric_limits<short>::max());
  maxpt.setValue(-std::numeric_limits<short>::max(), -std::numeric_limits<short>::max());
}

/*!
  \fn const SbVec2s & SbBox2s::getMin(void) const

  Returns the minimum point. This should usually be the lower left corner
  point of the box.

  \sa getOrigin(), getMax().
*/

/*!
  \fn const SbVec2s & SbBox2s::getMax(void) const

  Returns the maximum point. This should usually be the upper right corner
  point of the box.

  \sa getMin().
*/

/*!
  Extend the boundaries of the box by the given point, i.e. make the
  point fit inside the box if it isn't already within it.
*/
void
SbBox2s::extendBy(const SbVec2s & point)
{
  // The explicit casts are done to humour the HPUX aCC compiler,
  // which will otherwise say ``Template deduction failed to find a
  // match for the call to 'SbMin'''. mortene.
  this->minpt.setValue(SbMin(static_cast<short>(point[0]), static_cast<short>(this->minpt[0])),
                       SbMin(static_cast<short>(point[1]), static_cast<short>(this->minpt[1])));
  this->maxpt.setValue(SbMax(static_cast<short>(point[0]), static_cast<short>(this->maxpt[0])),
                       SbMax(static_cast<short>(point[1]), static_cast<short>(this->maxpt[1])));
}

/*!
  Extend the boundaries of the box by the given \a box parameter. This
  is equal to calling extendBy() twice with the corner points.
*/
void
SbBox2s::extendBy(const SbBox2s & box)
{
  if (box.isEmpty()) { return; }

  this->extendBy(box.getMin());
  this->extendBy(box.getMax());
}

/*!
  Check if the given point lies within the boundaries of this box.
*/
SbBool
SbBox2s::intersect(const SbVec2s & point) const
{
  if((point[0] >= this->minpt[0]) && (point[0] <= this->maxpt[0]) &&
     (point[1] >= this->minpt[1]) && (point[1] <= this->maxpt[1])) return TRUE;
  return FALSE;
}

/*!
  Check if \a box lies wholly or partly within the boundaries
  of this box.
*/
SbBool
SbBox2s::intersect(const SbBox2s & box) const
{
  if((box.getMax()[0] < this->getMin()[0]) ||
     (box.getMax()[1] < this->getMin()[1]) ||
     (box.getMin()[0] > this->getMax()[0]) ||
     (box.getMin()[1] > this->getMax()[1])) return FALSE;
  return TRUE;
}

/*!
  \fn void SbBox2s::getBounds(short & xmin, short & ymin, short & xmax, short & ymax) const

  Returns the box boundary coordinates.

  \sa setBounds(), getMin(), getMax().
*/

/*!
  \fn void SbBox2s::getBounds(SbVec2s & boxmin, SbVec2s & boxmax) const

  Returns the box corner points.

  \sa setBounds(), getMin(), getMax().
*/

/*!
  \fn void SbBox2s::getOrigin(short & originX, short & originY) const

  Returns the coordinates of the box origin (i.e. the lower left corner).

  \sa getMin().
*/

/*!
  \fn void SbBox2s::getSize(short & sizeX, short & sizeY) const

  Returns width and height of box.
*/

/*!
  \fn float SbBox2s::getAspectRatio(void) const

  Returns aspect ratio of box, which is defined as box width divided on
  box height.
*/

/*!
  \fn int operator == (const SbBox2s & b1, const SbBox2s & b2)
  \relates SbBox2s

  Check \a b1 and \a b2 for equality.
*/

/*!
  \fn int operator != (const SbBox2s & b1, const SbBox2s & b2)
  \relates SbBox2s

  Check \a b1 and \a b2 for inequality.
*/

/*!
  \fn SbBool SbBox2s::hasArea(void) const
*/
