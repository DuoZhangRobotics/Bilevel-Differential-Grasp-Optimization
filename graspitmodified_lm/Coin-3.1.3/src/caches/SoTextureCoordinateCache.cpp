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
  \class SoTextureCoordinateCache SoTextureCoordinateCache.h Inventor/caches/SoTextureCoordinateCache.h
  The SoTextureCoordinateClass is used to cache generated texture coordinates.
*/

#include <Inventor/caches/SoTextureCoordinateCache.h>
#include <Inventor/lists/SbList.h>
#include <Inventor/SbBox3f.h>
#include <Inventor/SbVec2f.h>

#ifndef DOXYGEN_SKIP_THIS

class SoTextureCoordinateCacheP {
public:
  SbList <SbVec2f> texCoords;
};

#endif // DOXYGEN_SKIP_THIS

#define PRIVATE(obj) ((obj)->pimpl)

/*!
  Constructor.
*/
SoTextureCoordinateCache::SoTextureCoordinateCache(SoState * const state)
  : SoCache(state)
{
  PRIVATE(this) = new SoTextureCoordinateCacheP;
}

/*!
  Destructor.
*/
SoTextureCoordinateCache::~SoTextureCoordinateCache()
{
  delete PRIVATE(this);
}

/*!
  Generates texture coordinates based on the bounding box of the
  geometry. This is usually called default texture coordinates
  in OIV.
*/
void
SoTextureCoordinateCache::generate(const SbBox3f & bbox,
                                   const SbVec3f * vertices,
                                   const int numvertices)
{
  // FIXME: Support 3D texture coordinates. This functionality
  // is more or less a duplicate of that in SoTextureCoordinateBundle
  // (kintel 20020203)
  float sizes[3];
  float minvalues[3];
  int   offsets[3];

  bbox.getSize(sizes[0], sizes[1], sizes[2]);
  minvalues[0] = bbox.getMin()[0];
  minvalues[1] = bbox.getMin()[1];
  minvalues[2] = bbox.getMin()[2];
  offsets[0] = 0;
  offsets[1] = 1;
  offsets[2] = 2;

  // bubblesort :-)
  if (sizes[0] < sizes[1]) {
    SbSwap(sizes[0], sizes[1]);
    SbSwap(minvalues[0], minvalues[1]);
    SbSwap(offsets[0], offsets[1]);
  }
  if (sizes[1] < sizes[2]) {
    SbSwap(sizes[1], sizes[2]);
    SbSwap(minvalues[1], minvalues[2]);
    SbSwap(offsets[1], offsets[2]);
  }
  if (sizes[0] < sizes[1]) {
    SbSwap(sizes[0], sizes[1]);
    SbSwap(minvalues[0], minvalues[1]);
    SbSwap(offsets[0], offsets[1]);
  }

  float s, t;
  for (int i = 0; i < numvertices; i++) {
    s = vertices[i][offsets[0]];
    t = vertices[i][offsets[1]];
    s -= minvalues[0];
    t -= minvalues[1];
    s /= sizes[0];
    t /= sizes[1];
    // expand list array as needed
    if (i >= PRIVATE(this)->texCoords.getLength()) PRIVATE(this)->texCoords.append(SbVec2f());
    PRIVATE(this)->texCoords[i].setValue(s, t);
  }

  // fit list array in case we used to have more items than now
  PRIVATE(this)->texCoords.truncate(numvertices);
}

/*!
  Returns the generated texture coordinates.
*/
const SbVec2f *
SoTextureCoordinateCache::get(void) const
{
  return PRIVATE(this)->texCoords.getArrayPtr();
}

/*!
  Returns the number of generated texture coordinates.
*/
int
SoTextureCoordinateCache::getNum(void) const
{
  return PRIVATE(this)->texCoords.getLength();
}

#undef PRIVATE
