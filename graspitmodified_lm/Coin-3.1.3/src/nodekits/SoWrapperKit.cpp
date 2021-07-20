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

#ifdef HAVE_NODEKITS

/*!
  \class SoWrapperKit SoWrapperKit.h Inventor/nodekits/SoWrapperKit.h
  \brief The SoWrapperKit class is a simple kit for wrapping a transform and a sub-graph.
  \ingroup nodekits

  
  \NODEKIT_PRE_DIAGRAM

  \verbatim
  CLASS SoWrapperKit
  -->"this"
        "callbackList"
        "topSeparator"
           "pickStyle"
           "appearance"
           "units"
           "transform"
           "texture2Transform"
           "childList"
  -->      "localTransform"
  -->      "contents"
  \endverbatim

  \NODEKIT_POST_DIAGRAM


  \NODEKIT_PRE_TABLE

  \verbatim
  CLASS SoWrapperKit
  PVT   "this",  SoWrapperKit  --- 
        "callbackList",  SoNodeKitListPart [ SoCallback, SoEventCallback ] 
  PVT   "topSeparator",  SoSeparator  --- 
        "pickStyle",  SoPickStyle  --- 
        "appearance",  SoAppearanceKit  --- 
        "units",  SoUnits  --- 
        "transform",  SoTransform  --- 
        "texture2Transform",  SoTexture2Transform  --- 
        "childList",  SoNodeKitListPart [ SoShapeKit, SoSeparatorKit ] 
        "localTransform",  SoTransform  --- 
        "contents",  SoSeparator  --- 
  \endverbatim

  \NODEKIT_POST_TABLE
*/

#include <Inventor/nodekits/SoWrapperKit.h>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTransform.h>

#include "nodekits/SoSubKitP.h"

SO_KIT_SOURCE(SoWrapperKit);


/*!
  Constructor.
*/
SoWrapperKit::SoWrapperKit(void)
{
  SO_KIT_INTERNAL_CONSTRUCTOR(SoWrapperKit);

  // Note: we must use "" instead of , , to humour MS VisualC++ 6.

  SO_KIT_ADD_CATALOG_ENTRY(localTransform, SoTransform, TRUE, topSeparator, contents, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(contents, SoSeparator, TRUE, topSeparator, "", TRUE);

  SO_KIT_INIT_INSTANCE();
}

/*!
  Destructor.
*/
SoWrapperKit::~SoWrapperKit()
{
}

// Documented in superclass.
void
SoWrapperKit::initClass(void)
{
  SO_KIT_INTERNAL_INIT_CLASS(SoWrapperKit, SO_FROM_INVENTOR_1);
}

#endif // HAVE_NODEKITS
