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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_MANIPULATORS

/*!
  \class SoTabBoxManip SoTabBoxManip.h Inventor/manips/SoTabBoxManip.h
  \brief The SoTabBoxManip class wraps an SoTabBoxDragger.
  \ingroup manips

  <center>
  <img src="http://doc.coin3d.org/images/Coin/draggers/tabbox.png">
  </center>

  The SoTabBoxManip provides a convenient mechanism for the
  application programmer for setting up an SoTabBoxDragger in the
  scene connected to the relevant fields of an SoTransform node.

  The interaction from the end-user with the manipulator will then
  automatically influence the transformation matrix for the geometry
  following it in the scenegraph.
*/

#include <Inventor/manips/SoTabBoxManip.h>

#include <Inventor/nodes/SoSurroundScale.h>
#include <Inventor/draggers/SoTabBoxDragger.h>

#include "nodes/SoSubNodeP.h"

class SoTabBoxManipP {
public:
};

SO_NODE_SOURCE(SoTabBoxManip);

// Doc in superclass.
void
SoTabBoxManip::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTabBoxManip, SO_FROM_INVENTOR_1);
}

/*!
  Constructor sets us up with an SoTabBoxDragger for manipulating a
  transformation.
*/
SoTabBoxManip::SoTabBoxManip(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTabBoxManip);

  SoTabBoxDragger *dragger = new SoTabBoxDragger;
  this->setDragger(dragger);

  SoSurroundScale * ss =
    (SoSurroundScale *)dragger->getPart("surroundScale", TRUE);
  ss->numNodesUpToContainer = 4;
  ss->numNodesUpToReset = 3;
}


/*!
  Protected destructor. (SoHandleBoxManip is automatically destructed
  when it's reference count goes to 0.)
 */
SoTabBoxManip::~SoTabBoxManip()
{
}

#endif // HAVE_MANIPULATORS
