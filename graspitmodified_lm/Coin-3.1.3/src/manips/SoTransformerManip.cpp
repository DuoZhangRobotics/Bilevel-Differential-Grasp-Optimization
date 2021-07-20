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
  \class SoTransformerManip SoTransformerManip.h Inventor/manips/SoTransformerManip.h
  \brief The SoTransformerManip wraps an SoTransformerDragger for convenience.
  \ingroup manips

  <center>
  <img src="http://doc.coin3d.org/images/Coin/draggers/transformer.png">
  </center>

  The manipulator class takes care of wrapping up the
  SoTransformerDragger in a simple and convenient API for the
  application programmer, making it automatically surround the
  geometry it influences and taking care of the book-keeping routines
  for it's interaction with the relevant fields of an SoTransformation
  node.
*/

#include <Inventor/manips/SoTransformerManip.h>

#include <Inventor/nodes/SoSurroundScale.h>
#include <Inventor/draggers/SoTransformerDragger.h>

#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

#include "nodes/SoSubNodeP.h"

class SoTransformerManipP {
public:
};

SO_NODE_SOURCE(SoTransformerManip);


// doc in super
void
SoTransformerManip::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoTransformerManip, SO_FROM_INVENTOR_1);
}

/*!
  Default constructor. Allocates an SoTransformerDragger and an
  SoSurroundScale node to surround the geometry within our part of the
  scenegraph.
*/
SoTransformerManip::SoTransformerManip(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoTransformerManip);

  SoTransformerDragger *dragger = new SoTransformerDragger;
  this->setDragger(dragger);

  SoSurroundScale *ss = (SoSurroundScale*) dragger->getPart("surroundScale", TRUE);
  ss->numNodesUpToContainer = 4;
  ss->numNodesUpToReset = 3;
}

/*!
  Destructor.
*/
SoTransformerManip::~SoTransformerManip()
{
}

/*!
  Convenience function to use the
  SoTransformerDragger::isLocateHighlighting() method of the embedded
  dragger. See documentation of that method.
*/
SbBool
SoTransformerManip::isLocateHighlighting(void)
{
  SoDragger *dragger = this->getDragger();
  if (dragger && dragger->isOfType(SoTransformerDragger::getClassTypeId())) {
    return ((SoTransformerDragger*)dragger)->isLocateHighlighting();
  }
#if COIN_DEBUG
  SoDebugError::postWarning("SoTransformerManip::isLocateHighlighting",
                            "Not a valid dragger in manipulator");
#endif // debug
  return FALSE;
}

/*!
  Convenience function to use the
  SoTransformerDragger::setLocateHighlighting() method of the embedded
  dragger. See documentation of that method.
*/
void
SoTransformerManip::setLocateHighlighting(SbBool onoff)
{
  SoDragger *dragger = this->getDragger();
  if (dragger && dragger->isOfType(SoTransformerDragger::getClassTypeId())) {
    ((SoTransformerDragger*)dragger)->setLocateHighlighting(onoff);
  }
  else {
#if COIN_DEBUG
    SoDebugError::postWarning("SoTransformerManip::setLocateHighlighting",
                              "Not a valid dragger in manipulator");
#endif // debug
  }
}

/*!
  Convenience function to use the SoTransformerDragger::unsquishKnobs()
  method of the embedded dragger. See documentation of that method.
*/
void
SoTransformerManip::unsquishKnobs(void)
{
  SoDragger *dragger = this->getDragger();
  if (dragger && dragger->isOfType(SoTransformerDragger::getClassTypeId())) {
    ((SoTransformerDragger*)dragger)->unsquishKnobs();
  }
  else {
#if COIN_DEBUG
    SoDebugError::postWarning("SoTransformerManip::setLocateHighlighting",
                              "Not a valid dragger in manipulator");
#endif // debug
  }
}

#endif // HAVE_MANIPULATORS
