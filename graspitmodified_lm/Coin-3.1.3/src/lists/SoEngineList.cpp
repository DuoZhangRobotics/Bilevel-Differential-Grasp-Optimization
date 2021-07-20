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
  \class SoEngineList SoEngineList.h Inventor/lists/SoEngineList.h
  \brief The SoEngineList class is a container for SoEngine objects.
  \ingroup engines

  As this class inherits SoBaseList, referencing and dereferencing
  will default be done on the objects at append(), remove(), insert()
  etc.
*/

#include <Inventor/lists/SoEngineList.h>


/*!
  Default constructor.
*/
SoEngineList::SoEngineList(void)
{
}

/*!
  Constructor with a hint about the number of elements the list will
  hold.

  \sa SoBaseList::SoBaseList(const int)
*/
SoEngineList::SoEngineList(const int size)
  : SoBaseList(size)
{
}

/*!
  Copy constructor.

  \sa SoBaseList::SoBaseList(const SoBaseList &)
*/
SoEngineList::SoEngineList(const SoEngineList & el)
  : SoBaseList(el)
{
}

/*!
  Destructor.

  \sa SoBaseList::~SoBaseList()
*/
SoEngineList::~SoEngineList()
{
}

/*!
  Append \a ptr to the list.

  \sa SoBaseList::append()
*/
void
SoEngineList::append(SoEngine * const ptr)
{
  SoBaseList::append((SoBase *)ptr);
}

/*!
  Return engine pointer at index \a i.

  \sa SoBaseList::operator[]()
*/
SoEngine *
SoEngineList::operator[](const int i) const
{
  return (SoEngine *)SoBaseList::operator[](i);
}

/*!
  Copy contents of list \a el to this list.

  \sa SoBaseList::operator=()
*/
SoEngineList &
SoEngineList::operator=(const SoEngineList & el)
{
  this->copy(el);
  return *this;
}
