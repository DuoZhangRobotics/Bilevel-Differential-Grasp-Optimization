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
  \class SoIndexedShape SoIndexedShape.h Inventor/nodes/SoIndexedShape.h
  \brief The SoIndexedShape class is the superclass for all indexed vertex shapes.
  \ingroup nodes

  This is an abstract class which contains storage for four fields for
  indices to coordinates, normals, materials and texture coordinates
  for it's subclasses.
*/

#include <Inventor/nodes/SoIndexedShape.h>

#include <Inventor/actions/SoAction.h>
#include <Inventor/elements/SoCoordinateElement.h>
#include <Inventor/elements/SoShapeHintsElement.h>
#include <Inventor/elements/SoNormalBindingElement.h>
#include <Inventor/elements/SoTextureCoordinateBindingElement.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/caches/SoNormalCache.h>
#include <Inventor/nodes/SoVertexProperty.h>

#include "nodes/SoSubNodeP.h"
#include "coindefs.h" // COIN_OBSOLETED()

/*!
  \var SoMFInt32 SoIndexedShape::coordIndex
  Coordinate indices.
*/
/*!
  \var SoMFInt32 SoIndexedShape::materialIndex
  Material indices.
*/
/*!
  \var SoMFInt32 SoIndexedShape::normalIndex
  Normal indices.
*/
/*!
  \var SoMFInt32 SoIndexedShape::textureCoordIndex
  Texture coordinate indices.
*/


SO_NODE_ABSTRACT_SOURCE(SoIndexedShape);

/*!
  Constructor.
*/
SoIndexedShape::SoIndexedShape()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoIndexedShape);

  SO_NODE_ADD_FIELD(coordIndex,(0));
  SO_NODE_ADD_FIELD(materialIndex,(-1));
  SO_NODE_ADD_FIELD(normalIndex,(-1));
  SO_NODE_ADD_FIELD(textureCoordIndex,(-1));
}

/*!
  Destructor.
*/
SoIndexedShape::~SoIndexedShape()
{
}

// Documented in superclass.
void
SoIndexedShape::initClass(void)
{
  SO_NODE_INTERNAL_INIT_ABSTRACT_CLASS(SoIndexedShape, SO_FROM_INVENTOR_1);
}

// Collects common error msg code for computeBBox() for both
// Coordinate3 and Coordinate4 loops.
static void
error_idx_out_of_bounds(const SoIndexedShape * node, int idxidx, int upper)
{
  SbName n = node->getName();
  SbString ns(" ");
  ns += n;

  SbName t = node->getTypeId().getName();

  SbString bounds;
  if (upper > 0) {
    // (Note: the if-expression above should have been "upper >= 0",
    // if it weren't for the default SoCoordinateElement containing
    // one default coordinate.)
    bounds.sprintf("should be within [0, %d]", upper);
  }
  else {
    bounds = "but no coordinates are available from the state!";
  }

  SoDebugError::post("SoIndexedShape::computeBBox",
                     "coordinate index nr %d for %snode%s of type %s is "
                     "out of bounds: is %d, %s",
                     idxidx,
                     (n == "") ? "<unnamed> " : "",
                     (n != "") ? ns.getString() : "",
                     t.getString(),
                     node->coordIndex[idxidx],
                     bounds.getString());
}

// Documented in superclass. Overridden to calculate bounding box of
// all indexed coordinates, using the coordIndex field.
void
SoIndexedShape::computeBBox(SoAction * action, SbBox3f & box, SbVec3f & center)
{
  assert(box.isEmpty());
  SoState * state = action->getState();

  const SoCoordinateElement * coordelem = NULL;
  SoNode *vpnode = this->vertexProperty.getValue();
  SoVertexProperty *vp = 
    (vpnode && vpnode->isOfType(SoVertexProperty::getClassTypeId())) ?
    (SoVertexProperty *)vpnode : NULL;
  SbBool vpvtx = vp && (vp->vertex.getNum() > 0);
  if (!vpvtx) {
    coordelem = SoCoordinateElement::getInstance(state);
  }

  const int numcoords = vpvtx ? vp->vertex.getNum() : coordelem->getNum();
  int numacc = 0; // to calculate weighted center point
  center.setValue(0.0f, 0.0f, 0.0f);

  if (vpvtx || coordelem->is3D()) {
    const SbVec3f * coords = vpvtx ?
      vp->vertex.getValues(0) :
      coordelem->getArrayPtr3();

    const int32_t * ptr = this->coordIndex.getValues(0);
    const int32_t * endptr = ptr + this->coordIndex.getNum();
    while (ptr < endptr) {
      const int idx = *ptr++;
      if (idx < numcoords) {
        if (idx >= 0) {
          box.extendBy(coords[idx]);
          center += coords[idx];
        }
      }
#if COIN_DEBUG
      else {
        const ptrdiff_t faultyidxpos = (ptr - this->coordIndex.getValues(0)) - 1;
        error_idx_out_of_bounds(this, (int)faultyidxpos, numcoords - 1);
        if (numcoords <= 1) break; // give only one error msg on missing coords
        // (the default state is that there's a default
        // SoCoordinateElement element with a single default
        // coordinate point setup)
      }
#endif // COIN_DEBUG
    }
  }
  else {
    SbVec3f tmp;
    const SbVec4f * coords = coordelem->getArrayPtr4();
    const int32_t * ptr = this->coordIndex.getValues(0);
    const int32_t * endptr = ptr + this->coordIndex.getNum();
    while (ptr < endptr) {
      int idx = *ptr++;
      if (idx < numcoords) {
        if (idx >= 0) {
          SbVec4f h = coords[idx];
          h.getReal(tmp);
          box.extendBy(tmp);
          center += tmp;
          numacc++;
        }
      }
#if COIN_DEBUG
      else {
        const ptrdiff_t faultyidxpos = (ptr - this->coordIndex.getValues(0)) - 1;
        error_idx_out_of_bounds(this, (int)faultyidxpos, numcoords - 1);
        if (numcoords <= 1) break; // give only one error msg on missing coords
        // (the default state is that there's a default
        // SoCoordinateElement element with a single default
        // coordinate point setup)
      }
#endif // COIN_DEBUG
    }
  }
  if (numacc) center /= (float) numacc;
}

/*!
  Returns whether texture coordinates should be indexed or not.

  \sa SoTextureCoordinateBinding
*/
SbBool
SoIndexedShape::areTexCoordsIndexed(SoAction * action)
{
  return SoTextureCoordinateBindingElement::get(action->getState()) ==
    SoTextureCoordinateBindingElement::PER_VERTEX_INDEXED;
}

/*!  
  Starting at index position \a startCoord, returns the number of
  indices until either the end of index array or a separator index (-1).  
*/
int
SoIndexedShape::getNumVerts(const int startCoord)
{
  const int32_t * ptr = this->coordIndex.getValues(0);
  const int32_t * endptr = ptr + this->coordIndex.getNum();
  ptr += startCoord;
  int cnt = 0;
  while (ptr < endptr && *ptr >= 0) cnt++, ptr++;
  return cnt;
}

/*!
  Not implemented. Probably only used for internal purposes in SGI's
  original Open Inventor, which means it should have been private.

  Let us know if you need this method for any code you are porting and
  we'll look into implement it properly.
*/
void
SoIndexedShape::setupIndices(const int /* numParts */,
                             const int /* numFaces */,
                             const SbBool /* needNormals */,
                             const SbBool /* needTexCoords */)
{
  COIN_OBSOLETED();
}

/*!
  Not implemented. Probably only used for internal purposes in SGI's
  original Open Inventor, which means it should have been private.

  Let us know if you need this method for any code you are porting and
  we'll look into implement it properly.
*/
const int32_t *
SoIndexedShape::getNormalIndices()
{
  COIN_OBSOLETED();
  return NULL;
}

/*!
  Not implemented. Probably only used for internal purposes in SGI's
  original Open Inventor, which means it should have been private.

  Let us know if you need this method for any code you are porting and
  we'll look into implement it properly.
*/
const int32_t *
SoIndexedShape::getColorIndices()
{
  COIN_OBSOLETED();
  return NULL;
}

/*!
  Not implemented. Probably only used for internal purposes in SGI's
  original Open Inventor, which means it should have been private.

  Let us know if you need this method for any code you are porting and
  we'll look into implement it properly.
*/
const int32_t *
SoIndexedShape::getTexCoordIndices()
{
  COIN_OBSOLETED();
  return NULL;
}

/*!
  Convenience method that will fetch data needed for rendering or
  generating primitives. Takes care of normal cache.

  This method was not part of the original SGI Open Inventor API, and
  is an extension specific for Coin.
*/
SbBool
SoIndexedShape::getVertexData(SoState * state,
                              const SoCoordinateElement *& coords,
                              const SbVec3f *& normals,
                              const int32_t *& cindices,
                              const int32_t *& nindices,
                              const int32_t *& tindices,
                              const int32_t *& mindices,
                              int & numcindices,
                              const SbBool needNormals,
                              SbBool & normalCacheUsed)
{
  SoVertexShape::getVertexData(state, coords, normals, needNormals);
  
  cindices = this->coordIndex.getValues(0);
  numcindices = this->coordIndex.getNum();

  mindices = this->materialIndex.getValues(0);
  if (this->materialIndex.getNum() <= 0 || mindices[0] < 0) mindices = NULL;

  tindices = this->textureCoordIndex.getValues(0);
  if (this->textureCoordIndex.getNum() <= 0 || tindices[0] < 0) tindices = NULL;

  normalCacheUsed = FALSE;
  nindices = NULL;
  if (needNormals) {
    nindices = this->normalIndex.getValues(0);
    if (this->normalIndex.getNum() <= 0 || nindices[0] < 0) nindices = NULL;

    if (normals == NULL) {
      SoNormalCache * nc = this->generateAndReadLockNormalCache(state);
      normals = nc->getNormals();
      nindices = nc->getIndices();
      normalCacheUsed = TRUE;
     
      // if no normals were generated, unlock normal cache before
      // returning
      if (normals == NULL) {
        this->readUnlockNormalCache();
        normalCacheUsed = FALSE;
      }
    }
  }
  return TRUE;
}
