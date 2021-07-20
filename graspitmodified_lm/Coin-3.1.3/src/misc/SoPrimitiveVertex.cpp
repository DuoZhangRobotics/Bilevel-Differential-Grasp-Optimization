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
  \class SoPrimitiveVertex SoPrimitiveVertex.h Inventor/SoPrimitiveVertex.h
  \brief The SoPrimitiveVertex class represents a single vertex of a generated primitive.
  \ingroup general

  Instances of SoPrimitiveVertex are constructed when generating
  primitive data, primarily during an SoCallbackAction traversal.
  Depending on the context the vertex could represent a single 3D
  point, one of the two vertices in a line or one of the three
  vertices in a triangle.
*/

#include <Inventor/SoPrimitiveVertex.h>
#include <stdlib.h>

/*!
  Default constructor, sets up a "void" instance.
*/
SoPrimitiveVertex::SoPrimitiveVertex(void)
  : point(0.0f, 0.0f, 0.0f),
    normal(0.0f, 0.0f, 1.0f),
    textureCoords(0.0f, 0.0f, 0.0f, 1.0f),
    materialIndex(0),
    detail(NULL),
    packedColor(0)
{
}

/*!
  Copy operator.

  When \a pv is copied into this instance, a \e shallow copy is
  made. Ie, only the reference to the detail instance is copied (if
  any), not the detail itself.
*/

SoPrimitiveVertex &
SoPrimitiveVertex::operator = (const SoPrimitiveVertex & pv)
{
  this->point = pv.point;
  this->normal = pv.normal;
  this->textureCoords = pv.textureCoords;
  this->materialIndex = pv.materialIndex;
  this->detail = pv.detail;
  this->packedColor = pv.packedColor;
  return *this;
}

/*!
  Copy constructor. Does a shallow copy.

  \sa SoPrimitiveVertex::operator=()
*/

SoPrimitiveVertex::SoPrimitiveVertex(const SoPrimitiveVertex & pv)
{
  *this = pv;
}

/*!
  Destructor. The detail instance is owned by client code and will not
  be destructed here.
*/

SoPrimitiveVertex::~SoPrimitiveVertex()
{
}

/*!
  \fn const SbVec3f & SoPrimitiveVertex::getPoint(void) const

  Returns vertex coordinates, positioned in object space.
*/

/*!
  \fn const SbVec3f & SoPrimitiveVertex::getNormal(void) const

  Returns normal vector, oriented in object space.
*/

/*!
  \fn const SbVec4f & SoPrimitiveVertex::getTextureCoords(void) const

  Returns texture coordinates for vertex, specified in object space.
*/

/*!
  \fn int SoPrimitiveVertex::getMaterialIndex(void) const

  Returns index of the vertex into the currently active material, if
  any.
*/

/*!
  \fn uint32_t SoPrimitiveVertex::getPackedColor(void) const

  Returns the rgba packed color for the given vertex.

  \since Coin 3.0
*/

/*!
  \fn const SoDetail * SoPrimitiveVertex::getDetail(void) const

  Returns pointer to detail instance with more information about the
  vertex. A detail instance may not be available, and if so \c NULL is
  returned.
*/

/*!
  \fn void SoPrimitiveVertex::setPoint(const SbVec3f & pt)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setPoint(float x, float y, float z)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setNormal(const SbVec3f & n)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setNormal(float nx, float ny, float nz)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setTextureCoords(const SbVec2f & texcoords)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.

  Covenience function. Will fill in 0 and 1 in the last two texture
  coords in the internal SbVec4f texture coordinate instance.
*/

/*!
  \fn void SoPrimitiveVertex::setTextureCoords(float tx, float ty)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.

  Covenience function. Will fill in 0 and 1 in the last two texture
  coords in the internal SbVec4f texture coordinate instance.
*/

/*!
  \fn void SoPrimitiveVertex::setTextureCoords(const SbVec3f & texcoords)

  Covenience function. Will fill in 1 in the last coord.

  \COIN_FUNCTION_EXTENSION

  \since Coin 2.0
*/

/*!
  \fn void SoPrimitiveVertex::setTextureCoords(float tx, float ty, float tz)

  Covenience function. Will fill in 1 in the last coord.

  \COIN_FUNCTION_EXTENSION
  \since Coin 2.5
*/

/*!
  \fn void SoPrimitiveVertex::setTextureCoords(const SbVec4f & texcoords)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setMaterialIndex(int index)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setPackedColor(uint32_t rgba)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.
*/

/*!
  \fn void SoPrimitiveVertex::setDetail(SoDetail * detail)

  Used internally from library client code setting up an
  SoPrimitiveVertex instance.

  Note that it's the client's responsibility to do the deallocation of
  the detail instance after the SoPrimitiveVertex instance has gone
  out of scope.
*/
