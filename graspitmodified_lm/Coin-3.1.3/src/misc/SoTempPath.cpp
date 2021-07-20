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
  \class SoTempPath Inventor/misc/SoTempPath.h
  \brief The SoTempPath class is used to store temporary paths.
  \ingroup general

  The path simply turns off auditing in the constructor, and leaves
  the user with the responsibility of keeping the path valid.
*/

#include <Inventor/misc/SoTempPath.h>

/*!
  Constructor.
*/
SoTempPath::SoTempPath(const int approxlength)
  : SoFullPath(approxlength)
{
  this->auditPath(FALSE);
  this->nodes.addReferences(FALSE);
}

/*!
  Append a node (specified by \a node and parent child \a index) to the path.
  This method is only available in SoTempPath, since it will not
  consider auditing or hidden children.
*/
void
SoTempPath::simpleAppend(SoNode * const node, const int index)
{
  // this will make SoPath rescan the path for hidden children the
  // next time getLength() is called.
  this->firsthiddendirty = TRUE;

  // just append node and index
  this->nodes.append(node);
  this->indices.append(index);
}

/*!  
  Replace the tail of this path. The node is specified by \a node
  and parent child \a index. This method is only available in
  SoTempPath,, since it will not consider auditing or hidden children.  
*/
void 
SoTempPath::replaceTail(SoNode * const node, const int index)
{
  // this will make SoPath rescan the path for hidden children the
  // next time getLength() is called.
  this->firsthiddendirty = TRUE;

  // just replace the last node and index
  const int i = this->nodes.getLength() - 1;
  this->nodes.set(i, (SoBase*) node);
  this->indices[i] = index;
}
