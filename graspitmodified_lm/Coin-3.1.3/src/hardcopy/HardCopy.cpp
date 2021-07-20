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

// *************************************************************************

/*!
  \page hardcopy_overview A HardCopy Overview

  The main API for HardCopy support in Coin is the abstract class
  SoVectorizeAction. SoVectorizeAction will extract geometry from an
  Inventor scene graph, and project the geometry onto a specified
  page.  Since postscript and other vector based file formats do not
  support z-buffer or depth clipping, all geometry is rendered using a
  simple painter's algorithm (geometry is sorted based on distance to
  camera).

  SoVectorizePSAction inherits SoVectorizeAction, and will output a
  Postscript file.

  Texture-mapped polygons are not supported, since this is not
  supported by the vector file formats, at least it's not supported in
  Postscript. Gouraud shading is not supported in the Postscript
  language (at least not for V2.0), but an approximation is
  implemeting using an algorithm that divides the triangle into
  several small (flat-shaded) triangles. The gouraud shading quality
  (the number of sub-triangles) is controlled by an epsilon value. The
  gouraud shading function is written by Frederic Delhoume
  (delhoume (at) ilog.fr), and is free (public domain) software.

  Typical use of SoVectorizePSAction is shown in the following piece
  of code:

  \code

  SoVectorizePSAction * ps = new SoVectorizePSAction;
  SoVectorOutput * out = ps->getOutput();

  if (!out->openFile("output.ps")) {
    return -1; // unable to open output file
  }

  // to enable gouraud shading. 0.1 is a nice epsilon value
  // ps->setGouraudThreshold(0.1f);

  // clear to white background. Not really necessary if you
  // want a white background
  ps->setBackgroundColor(TRUE, SbColor(1.0f, 1.0f, 1.0f));

  // select LANDSCAPE or PORTRAIT orientation
  ps->setOrientation(SoVectorizeAction::LANDSCAPE);

  // start creating a new page (A4 page, with 10mm border).
  ps->beginPage(SbVec2f(10.0f, 10.0f), SbVec2f(190.0f, 277.0f));

  // There are also enums for A0-A10. Example:
  // ps->beginStandardPage(SoVectorizeAction::A4, 10.0f);

  // calibrate so that text, lines, points and images will have the
  // same size in the postscript file as on the monitor.
  ps->calibrate(viewer->getViewportRegion());

  // apply action on the viewer scenegraph. Remember to use
  // SoSceneManager's scene graph so that the camera is included.
  ps->apply(viewer->getSceneManager()->getSceneGraph());

  // this will create the postscript file
  ps->endPage();

  // close file
  out->closeFile();

  delete ps;

  \endcode

  It is also possible to have several viewports and/or layers on a
  page. This is useful if your application has several layers of
  geometry, for instance some annotations in 2D on top of a 3D scene
  graph. To create several layers, the beginViewport() and
  endViewport() functions can be used.

  \ingroup hardcopy
  \since Coin 2.1
  \since TGS provides HardCopy support as a separate extension for TGS Inventor.
*/

// *************************************************************************

/*!
  \class SoHardCopy SoHardCopy.h HardCopy/SoHardCopy.h
  \brief The SoHardCopy class is a static class for initializing the hardcopy support.
  \ingroup hardcopy
*/

// *************************************************************************

#include <Inventor/annex/HardCopy/SoHardCopy.h>

#include <Inventor/annex/HardCopy/SoVectorizePSAction.h>

#include "tidbitsp.h"

// *************************************************************************

static SbBool hardcopy_isinitialized = FALSE;

static void hardcopy_cleanup(void)
{
  hardcopy_isinitialized = FALSE;
}

/*!
  Initialization of the hardcopy support happens automatically from
  SoDB::init(), so the application programmer should normally not have
  to worry about it.
*/
void
SoHardCopy::init(void)
{
  if (hardcopy_isinitialized) { return; }

  SoVectorizeAction::initClass();
  SoVectorizePSAction::initClass();

  hardcopy_isinitialized = TRUE;
  coin_atexit((coin_atexit_f*)hardcopy_cleanup, CC_ATEXIT_NORMAL);
}

/*!
  Returns name of the hardcopy extension.
*/
const char *
SoHardCopy::getProductName(void)
{
  return "Kongsberg Oil & Gas Technologies hard-copy support for Coin.";
}

/*!
  Version number will always match Coin version number.
*/
const char *
SoHardCopy::getVersion(void)
{
  return COIN_VERSION;
}

// *************************************************************************
