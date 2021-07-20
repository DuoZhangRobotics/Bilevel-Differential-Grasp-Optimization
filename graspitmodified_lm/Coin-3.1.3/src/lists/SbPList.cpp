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
  \class SbPList SbPList.h Inventor/lists/SbPList.h
  \brief The SbPList class is a container class for void pointers.
  \ingroup base

*/


#include <Inventor/lists/SbPList.h>

/*!
  \fn SbPList::SbPList(const int sizehint)

  This constructor initializes the internal allocated size for the
  list to \a sizehint. Note that the list will still initially contain
  zero items.

*/

/*!
  \fn void SbPList::append(void * item)

  Append \a item to the end of the list.

  Automatically allocates more items internally if needed.
*/

/*!
  \fn void * SbPList::get(const int index) const

  Returns element at \a index. Does \e not expand array bounds if \a
  index is outside the list.
*/

/*!
  \fn void SbPList::set(const int index, void * item)

  Index operator to set element at \a index. Does \e not expand array
  bounds if \a index is outside the list.
*/

/*!
  \fn void SbPList::removeFast(const int index)

  Remove the item at \a index, moving the last item into its place and
  truncating the list.
*/

/*!
  \fn int SbPList::getLength(void) const

  Returns number of items in the list.
*/

/*!
 \fn void SbPList::truncate(const int length, const int fit)

 Shorten the list to contain \a length elements, removing items from
 \e index \a length and onwards.

 If \a fit is non-zero, will also shrink the internal size of the
 allocated array. Note that this is much less efficient than not
 re-fitting the array size.
*/

/*!
  \fn void ** SbPList::getArrayPtr(const int start = 0) const

  Returns pointer to a non-modifiable array of the lists elements.
  \a start specifies an index into the array.

  The caller is \e not responsible for freeing up the array, as it is
  just a pointer into the internal array used by the list.
*/

/*!
  \fn void *& SbPList::operator[](const int index) const

  Returns element at \a index.

  Will automatically expand the size of the internal array if \a index
  is outside the current bounds of the list. The values of any
  additional pointers are then set to \c NULL.
*/

/*!
  \fn SbBool SbPList::operator!=(const SbPList & l) const

  Inequality operator. Returns \c TRUE if this list and \a l are not
  equal.
*/

/*!
  \fn void SbPList::expand(const int size)

  Expand the list to contain \a size items. The new items added at the
  end have undefined value.
*/

/*!
  \fn int SbPList::getArraySize(void) const

  Return number of items there's allocated space for in the array.

  \sa getLength()
*/

/*!
  Default constructor.
*/
SbPList::SbPList(const int sizehint)
  : itembuffersize(DEFAULTSIZE), numitems(0), itembuffer(builtinbuffer)
{
  if (sizehint > DEFAULTSIZE) this->grow(sizehint);
}

/*!
  Copy constructor.
*/
SbPList::SbPList(const SbPList & l)
  : itembuffersize(DEFAULTSIZE), numitems(0), itembuffer(builtinbuffer)
{
  this->copy(l);
}

/*!
  Destructor.
*/
SbPList::~SbPList()
{
  if (this->itembuffer != builtinbuffer) delete[] this->itembuffer;
}

/*!
  Make this list a copy of \a l.
*/
void
SbPList::copy(const SbPList & l)
{
  if (this == &l) return;
  const int n = l.numitems;
  this->expand(n);
  for (int i = 0; i < n; i++) this->itembuffer[i] = l.itembuffer[i];
}

/*!
  Assignment operator
*/
SbPList &
SbPList::operator=(const SbPList & l)
{
  this->copy(l);
  return *this;
}

/*!
  Fit the allocated array exactly around the length of the list,
  descarding memory spent on unused pre-allocated array cells.

  You should normally not need or want to call this method, and it is
  only available for the sake of having the option to optimize memory
  usage for the unlikely event that you should throw around huge
  SbList objects within your application.
*/
void
SbPList::fit(void)
{
  const int items = this->numitems;

  if (items < this->itembuffersize) {
    void ** newitembuffer = this->builtinbuffer;
    if (items > DEFAULTSIZE) newitembuffer = new void*[items];

    if (newitembuffer != this->itembuffer) {
      for (int i = 0; i < items; i++) newitembuffer[i] = this->itembuffer[i];
    }

    if (this->itembuffer != this->builtinbuffer) delete[] this->itembuffer;
    this->itembuffer = newitembuffer;
    this->itembuffersize = items > DEFAULTSIZE ? items : DEFAULTSIZE;
  }
}

/*!
  Return index of first occurrence of \a item in the list, or -1 if \a
  item is not present.
*/
int
SbPList::find(void * item) const
{
  for (int i = 0; i < this->numitems; i++)
    if (this->itembuffer[i] == item) return i;
  return -1;
}

/*!
  Insert \a item at index \a insertbefore.

  \a insertbefore should not be larger than the current number of
  items in the list.
*/
void
SbPList::insert(void * item, const int insertbefore) {
#ifdef COIN_EXTRA_DEBUG
  assert(insertbefore >= 0 && insertbefore <= this->numitems);
#endif // COIN_EXTRA_DEBUG
  if (this->numitems == this->itembuffersize) this->grow();

  for (int i = this->numitems; i > insertbefore; i--)
    this->itembuffer[i] = this->itembuffer[i-1];
  this->itembuffer[insertbefore] = item;
  this->numitems++;
}

/*!
  Removes an \a item from the list. If there are several items with
  the same value, removes the \a item with the lowest index.
*/
void
SbPList::removeItem(void * item)
{
  int idx = this->find(item);
#ifdef COIN_EXTRA_DEBUG
  assert(idx != -1);
#endif // COIN_EXTRA_DEBUG
  if (idx >= 0) {
    this->remove(idx);
  }
}

/*!
  Remove the item at \a index, moving all subsequent items downwards
  one place in the list.
*/
void
SbPList::remove(const int index)
{
#ifdef COIN_EXTRA_DEBUG
  assert(index >= 0 && index < this->numitems);
#endif // COIN_EXTRA_DEBUG
  this->numitems--;
  for (int i = index; i < this->numitems; i++)
    this->itembuffer[i] = this->itembuffer[i + 1];
}

/*!
  Equality operator. Returns \c TRUE if this list and \a l are
  identical, containing the exact same ordered set of elements.
*/
int
SbPList::operator==(const SbPList & l) const
{
  if (this == &l) return TRUE;
  if (this->numitems != l.numitems) return FALSE;
  for (int i = 0; i < this->numitems; i++)
    if (this->itembuffer[i] != l.itembuffer[i]) return FALSE;
  return TRUE;
}

// Expand list to the given size, filling in with NULL pointers.
void
SbPList::expandlist(const int size) const
{
  const int oldsize = this->getLength();
  SbPList * thisp = (SbPList *)this;
  thisp->expand(size);
  for (int i = oldsize; i < size; i++) (*thisp)[i] = NULL;
}

// grow allocated array, not number of items
void
SbPList::grow(const int size)
{
  // Default behavior is to double array size.
  if (size == -1) this->itembuffersize <<= 1;
  else if (size <= this->itembuffersize) return;
  else { this->itembuffersize = size; }

  void ** newbuffer = new void*[this->itembuffersize];
  const int n = this->numitems;
  for (int i = 0; i < n; i++) newbuffer[i] = this->itembuffer[i];
  if (this->itembuffer != this->builtinbuffer) delete[] this->itembuffer;
  this->itembuffer = newbuffer;
}
