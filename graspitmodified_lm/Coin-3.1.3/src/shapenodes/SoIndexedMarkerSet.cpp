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
  \class SoIndexedMarkerSet SoIndexedMarkerSet.h Inventor/nodes/SoIndexedMarkerSet.h
  \brief The SoIndexedMarkerSet class is used to display a set of bitmap markers at 3D positions.
  \ingroup nodes

  This node either uses the coordinates currently on the state
  (typically set up by a leading SoCoordinate3 node in the scenegraph)
  or from a SoVertexProperty node attached to this node to render a
  set of 3D points.

  To add new markers, use the static functions in SoMarkerSet.

  Here's a simple usage example of SoIndexedMarkerSet in a scenegraph:

  \verbatim
  #Inventor V2.1 ascii

  Separator {
     Material {
        diffuseColor [
         1 0 0, 0 1 0, 0 0 1, 1 1 0, 1 0 1, 1 1 1, 1 0.8 0.6, 0.6 0.8 1
        ]
     }
     MaterialBinding { value PER_VERTEX_INDEXED }

     Coordinate3 {
        point [
         -1 1 0, -1 -1 0, 1 -1 0, 1 1 0, 0 2 -1, -2 0 -1, 0 -2 -1, 2 0 -1
        ]
     }

     IndexedMarkerSet { 
        coordIndex [0, 1, 2, 3, 4, 5, 6, 7]
        markerIndex [0, 1, 0, 1, 0, 1, 0, 1]
     }
  }

  \endverbatim

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
  IndexedMarkerSet { 
    vertexProperty 	NULL 
    coordIndex 	0 
    materialIndex 	-1 
    normalIndex 	-1 
    textureCoordIndex 	-1
    markerIndex   -1
  }
  \endcode

  \since TGS Inventor 6.0, Coin 3.1
*/

#include <Inventor/nodes/SoIndexedMarkerSet.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/misc/SoState.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/system/gl.h>
#include <Inventor/nodes/SoVertexProperty.h>
#include <Inventor/nodes/SoMarkerSet.h>
#include <Inventor/elements/SoGLCoordinateElement.h>
#include <Inventor/elements/SoMaterialBindingElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoGLTexture3EnabledElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoViewVolumeElement.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/elements/SoProjectionMatrixElement.h>
#include <Inventor/elements/SoGLLazyElement.h>
#include <Inventor/elements/SoGLVBOElement.h>
#include <Inventor/elements/SoCullElement.h>
#include <Inventor/elements/SoGLCacheContextElement.h>

#include <Inventor/bundles/SoMaterialBundle.h>
#include <Inventor/SbViewVolume.h>

#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

#include "nodes/SoSubNodeP.h"
#include "misc/SoGL.h"

SO_NODE_SOURCE(SoIndexedMarkerSet);

/*!
  Constructor.
*/
SoIndexedMarkerSet::SoIndexedMarkerSet() : SoIndexedPointSet()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoIndexedMarkerSet);
  SO_NODE_ADD_FIELD(markerIndex, (0));
}

/*!
  Destructor.
*/
SoIndexedMarkerSet::~SoIndexedMarkerSet()
{
}

// doc from parent
void
SoIndexedMarkerSet::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoIndexedMarkerSet, SO_FROM_INVENTOR_6_0);
}

// doc from parent
void
SoIndexedMarkerSet::GLRender(SoGLRenderAction * action)
{
  int32_t numpts = this->coordIndex.getNum();
  if (numpts == 0) return;

  SoState * state = action->getState();
  state->push();

  if (this->vertexProperty.getValue()) {
    this->vertexProperty.getValue()->GLRender(action);
  }

  if (!this->shouldGLRender(action)){
    state->pop();
    return;
  }

  SoGLCacheContextElement::shouldAutoCache(state, SoGLCacheContextElement::DONT_AUTO_CACHE);

  // We just disable lighting and texturing for markers, since we
  // can't see any reason this should ever be enabled.  send an angry
  // email to <pederb@coin3d.org> if you disagree.

  SoLazyElement::setLightModel(state, SoLazyElement::BASE_COLOR);
  SoGLTextureEnabledElement::set(state, this, FALSE);
  SoGLTexture3EnabledElement::set(state, this, FALSE);

  SoMaterialBundle mb(action);

  const SbVec3f * normals;
  int numindices;
  const int32_t * cindices;
  const int32_t * nindices;
  const int32_t * tindices;
  const int32_t * mindices;
  const SoCoordinateElement * coords;
  Binding mbind;

  SbBool normalCacheUsed;

  this->getVertexData(state, coords, normals, cindices,
                      nindices, tindices, mindices, numindices,
                      false, normalCacheUsed);

  if (numindices == 0){//nothing to render
    state->pop();
    return;
  }

  mbind = this->findMaterialBinding(state);
  //if we don't have explicit material indices, use coord indices:
  if (mbind == PER_VERTEX_INDEXED && mindices == NULL) mindices = cindices;

  const SoGLCoordinateElement * glcoords = dynamic_cast<const SoGLCoordinateElement *>(coords);
  assert(glcoords && "could not cast to SoGLCoordinateElement");

  mb.sendFirst(); // always do this, even if mbind != OVERALL

  const SbMatrix & mat = SoModelMatrixElement::get(state);
  const SbViewVolume & vv = SoViewVolumeElement::get(state);
  const SbViewportRegion & vp = SoViewportRegionElement::get(state);
  const SbMatrix & projmatrix = (mat * SoViewingMatrixElement::get(state) *
                                 SoProjectionMatrixElement::get(state));
  SbVec2s vpsize = vp.getViewportSizePixels();

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, vpsize[0], 0, vpsize[1], -1.0f, 1.0f);

  for (int i = 0; i < numindices; i++) {
    int32_t idx = cindices[i];

    int midx = i;
#if COIN_DEBUG
      if (midx < 0 || midx > this->markerIndex.getNum() - 1) {
        static SbBool firsterror = TRUE;
        if (firsterror) {
          SoDebugError::postWarning("SoMarkerSet::GLRender",
                                    "no markerIndex for coordinate %d",
                                    i);
          firsterror = FALSE;
        }
        // Don't render, jump back to top of for-loop and continue with
        // next index.
        midx = this->markerIndex.getNum() - 1;
      }
#endif // COIN_DEBUG

    int marker = this->markerIndex[midx];
    if (marker == SoMarkerSet::NONE) { continue; }//no marker to render

    SbVec2s size;
    const unsigned char * bytes;
    SbBool isLSBFirst;

    if (marker >= SoMarkerSet::getNumDefinedMarkers()) continue;

    SbBool validMarker = SoMarkerSet::getMarker(marker, size, bytes, isLSBFirst);
    if (!validMarker) continue;

    if (mbind == PER_VERTEX_INDEXED) mb.send(mindices[i], TRUE);
    else if (mbind == PER_VERTEX) mb.send(i, TRUE);

    SbVec3f point = glcoords->get3(idx);

    // OpenGL's glBitmap() will not be clipped against anything but
    // the near and far planes. We want markers to also be clipped
    // against other clipping planes, to behave like the SoPointSet
    // superclass.
    const SbBox3f bbox(point, point);
    // FIXME: if there are *heaps* of markers, this next line will
    // probably become a bottleneck. Should really partition marker
    // positions in a oct-tree data structure and cull several at
    // the same time.  20031219 mortene.
    if (SoCullElement::cullTest(state, bbox, TRUE)) { continue; }

    projmatrix.multVecMatrix(point, point);
    point[0] = (point[0] + 1.0f) * 0.5f * vpsize[0];
    point[1] = (point[1] + 1.0f) * 0.5f * vpsize[1];      

    // To have the exact center point of the marker drawn at the
    // projected 3D position.  (FIXME: I haven't actually checked that
    // this is what TGS' implementation of the SoMarkerSet node does
    // when rendering, but it seems likely. 20010823 mortene.)

    point[0] = point[0] - (size[0] - 1) / 2;
    point[1] = point[1] - (size[1] - 1) / 2;

    //FIXME: this will probably fail if someone has overwritten one of the 
    //built-in markers. Currently there is no way of fetching a marker's 
    //alignment from outside the SoMarkerSet class though. 20090424 wiesener
    int align = (marker >= SoMarkerSet::NUM_MARKERS) ? 1 : 4;
    glPixelStorei(GL_UNPACK_ALIGNMENT, align);
    glRasterPos3f(point[0], point[1], -point[2]);
    glBitmap(size[0], size[1], 0, 0, 0, 0, bytes);
  }

  // FIXME: this looks wrong, shouldn't we rather reset the alignment
  // value to what it was previously?  20010824 mortene.
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  state->pop();

  // send approx number of points for autocache handling. Divide
  // by three so that three points is the same as one triangle.
  sogl_autocache_update(state, numindices/3, FALSE);
}

// Documented in superclass.
void
SoIndexedMarkerSet::notify(SoNotList * list)
{
  inherited::notify(list);
}


