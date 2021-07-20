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
  \class SoDetail SoDetail.h Inventor/details/SoDetail.h
  \brief The SoDetail class is the superclass for all classes storing
  detailed information about particular shapes.
  \ingroup details

  Detail information about shapes is used in relation to picking
  actions in Coin.  They typically contain the relevant information
  about what particular part of the shape a pick ray intersected with.

*/

// *************************************************************************

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <Inventor/SbName.h>
#include <Inventor/details/SoConeDetail.h>
#include <Inventor/details/SoCubeDetail.h>
#include <Inventor/details/SoCylinderDetail.h>
#include <Inventor/details/SoFaceDetail.h>
#include <Inventor/details/SoLineDetail.h>
#include <Inventor/details/SoPointDetail.h>
#include <Inventor/details/SoTextDetail.h>
#ifdef HAVE_NODEKITS
#include <Inventor/details/SoNodeKitDetail.h>
#endif // HAVE_NODEKITS

// *************************************************************************

/*!
  \fn SoDetail * SoDetail::copy(void) const
  Return a deep copy of ourself.

  \DANGEROUS_ALLOC_RETURN
*/

// Note: the following documentation for getTypeId() will also be
// visible for subclasses, so keep it general.  If you write any
// additional documentation for this method, check the other top-level
// classes with getTypeId() documentation to see if it is applicable
// to update those also (it probably is).
/*!
  \fn SoType SoDetail::getTypeId(void) const

  Returns the type identification of a detail derived from a class
  inheriting SoDetail.  This is used for run-time type checking and
  "downward" casting.

  Usage example:

  \code
  void fuhbear(SoDetail * detail)
  {
    if (detail->getTypeId() == SoFaceDetail::getClassTypeId()) {
      // safe downward cast, know the type
      SoFaceDetail * facedetail = (SoFaceDetail *)detail;
      /// [then something] ///
    }
    return; // ignore if not a SoFaceDetail
  }
  \endcode


  For application programmers wanting to extend the library with new
  detail classes: this method needs to be overridden in \e all
  subclasses. This is typically done as part of setting up the full
  type system for extension classes, which is usually accomplished by
  using the pre-defined macros available through
  Inventor/nodes/SoSubDetail.h: SO_DETAIL_SOURCE and
  SO_DETAIL_INIT_CLASS.
*/

// *************************************************************************

SoType SoDetail::classTypeId STATIC_SOTYPE_INIT;

// *************************************************************************

/*!
  Default constructor.
*/
SoDetail::SoDetail(void)
{
}

/*!
  Destructor.
*/
SoDetail::~SoDetail()
{
}

/*!
  Returns \c TRUE if \a type is derived from (or \e is) this class.
*/
SbBool
SoDetail::isOfType(const SoType type) const
{
  return this->getTypeId().isDerivedFrom(type);
}

/*!
  Returns the type for this class.
*/
SoType
SoDetail::getClassTypeId(void)
{
  return SoDetail::classTypeId;
}

// Note: the following documentation for initClass() will also be
// visible for subclasses, so keep it general.
/*!
  Initialize relevant common data for all instances, like the type
  system.
 */
void
SoDetail::initClass(void)
{
  SoDetail::classTypeId =
    SoType::createType(SoType::badType(), SbName("SoDetail"));
  SoDetail::initClasses();
}

/*!
  Call the initClass() methods of all built-in detail classes.

  (The initClass() method of user extension detail classes -- if any
  -- must be called explicitly by the application programmer in the
  application initialization code.)
 */
void
SoDetail::initClasses(void)
{
  SoConeDetail::initClass();
  SoCubeDetail::initClass();
  SoCylinderDetail::initClass();
  SoFaceDetail::initClass();
  SoLineDetail::initClass();
  SoPointDetail::initClass();
  SoTextDetail::initClass();
#ifdef HAVE_NODEKITS
  SoNodeKitDetail::initClass();
#endif // HAVE_NODEKITS
}
