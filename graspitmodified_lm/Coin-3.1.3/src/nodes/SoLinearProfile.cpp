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
  \class SoLinearProfile SoLinearProfile.h Inventor/nodes/SoLinearProfile.h
  \brief The SoLinearProfile class is a node for specifying linear profile curves.
  \ingroup nodes

  Use nodes of this type if you want to set up profiles that are
  simply straight lines connected by control points.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    LinearProfile {
        index 0
        linkage START_FIRST
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoLinearProfile.h>

#include <stdlib.h>

#include <Inventor/elements/SoProfileCoordinateElement.h>
#include <Inventor/lists/SbList.h>
#include <Inventor/threads/SbStorage.h>
#include <Inventor/errors/SoDebugError.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H
#include "tidbitsp.h"
#include "nodes/SoSubNodeP.h"

// *************************************************************************

typedef struct {
  SbList <float> * coordlist;
} so_linearprofile_data;

static void
so_linearprofile_construct_data(void * closure)
{
  so_linearprofile_data * data = (so_linearprofile_data*) closure;
  data->coordlist = NULL;
}

static void
so_linearprofile_destruct_data(void * closure)
{
  so_linearprofile_data * data = (so_linearprofile_data*) closure;
  delete data->coordlist;
}

static SbStorage * so_linearprofile_storage;

static void
so_linearprofile_cleanup(void)
{
  delete so_linearprofile_storage;
}

static SbList <float> *
so_linearprofile_get_coordlist(void)
{
  so_linearprofile_data * data = NULL;
  data = (so_linearprofile_data*) so_linearprofile_storage->get();
  if (data->coordlist == NULL) {
    data->coordlist = new SbList<float>;
  }
  return data->coordlist;
}

// *************************************************************************

SO_NODE_SOURCE(SoLinearProfile);

// *************************************************************************

/*!
  Constructor.
*/
SoLinearProfile::SoLinearProfile(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoLinearProfile);
}

/*!
  Destructor.
*/
SoLinearProfile::~SoLinearProfile()
{
}

// Doc from superclass.
void
SoLinearProfile::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoLinearProfile, SO_FROM_INVENTOR_1);

  so_linearprofile_storage = new SbStorage(sizeof(so_linearprofile_data),
                                           so_linearprofile_construct_data,
                                           so_linearprofile_destruct_data);
  coin_atexit((coin_atexit_f*) so_linearprofile_cleanup, CC_ATEXIT_NORMAL);
}

// Doc from superclass.
void
SoLinearProfile::getTrimCurve(SoState * state, int32_t & numpoints,
                              float *& points, int & floatspervec,
                              int32_t & numknots, float *& knotvector)
{
  SbList <float> * coordListLinearProfile = 
    so_linearprofile_get_coordlist();

  const SoProfileCoordinateElement * elem = (const SoProfileCoordinateElement*)
    SoProfileCoordinateElement::getInstance(state);

  coordListLinearProfile->truncate(0);

  // Get the number of SoProfileCoordinate2/3 points
  int32_t numcoords = elem->getNum();
  // Get the number of profile coordinate indices
  int n = this->index.getNum();

  if (numcoords) {
    // Both 2D or 3D profile coordinates might have been specified, so
    // get the appropriate coordinates and save the number of floats
    // per vector for later usage.
    if (elem->is2D()) {
      points = (float*) elem->getArrayPtr2();
      floatspervec = 2;
    }
    else {
      points = (float*) elem->getArrayPtr3();
      floatspervec = 3;
    }

    assert(points);
  }

  // Append the coordinates to a list over the profile coordinates.
  for (int i = 0; i < n; i++) {
    int idx = this->index[i];
    
    // If valid profile coordinates have been specified
    if (idx >= 0 && idx < numcoords) {
      for (int j = 0; j < floatspervec; j++) {
        coordListLinearProfile->append(points[idx * floatspervec + j]);
      }
    }
    // If invalid profile coordinates have been specified
    else {
      // Add dummy coordinate for robustness
      for (int j = 0; j < floatspervec; j++) {
        coordListLinearProfile->append(0.0f);
      }
      
      // Print errormessage
      static uint32_t current_errors = 0;
      if (current_errors < 1) {
        SoDebugError::postWarning("SoLinearProfile::getVertices", "Illegal profile "
                                  "coordinate index specified: %d. Should be within "
                                  "[0, %d]", idx, numcoords - 1);
      }
      current_errors++;
    }
  }

  // Set  return variables
  points = (float*) coordListLinearProfile->getArrayPtr();
  numpoints = n;
  numknots = 0;
}

// Doc from superclass.
void
SoLinearProfile::getVertices(SoState * state, int32_t & numvertices,
                             SbVec2f *& vertices)
{
  SbList <float> * coordListLinearProfile = 
    so_linearprofile_get_coordlist();

  const SoProfileCoordinateElement * elem = (const SoProfileCoordinateElement*)
    SoProfileCoordinateElement::getInstance(state);

  coordListLinearProfile->truncate(0);

  // Get the number of SoProfileCoordinate2/3 points
  int32_t numcoords = elem->getNum();
  // Get the number of profile coordinate indices
  int n = this->index.getNum();

  if (numcoords) {
    float * points;
    int floatspervec;
    // Both 2D or 3D profile coordinates might have been specified, so
    // get the appropriate coordinates and save the number of floats
    // per vector for later usage.
    if (elem->is2D()) {
      points = (float *) elem->getArrayPtr2();
      floatspervec = 2;
    }
    else {
      points = (float *) elem->getArrayPtr3();
      floatspervec = 3;
    }

    assert(points);

    // Append the coordinates to a list over the profile coordinates.
    // When 3D profile coordinates have been specified, only the
    // 2D-part of the coordinates will be used.
    for (int i = 0; i < n; i++) {
      int idx = this->index[i];

      // If valid profile coordinates have been specified
      if (idx >= 0 && idx < numcoords) {
        coordListLinearProfile->append(points[(idx * floatspervec)]);
        coordListLinearProfile->append(points[(idx * floatspervec) + 1]);
      }
      // If invalid profile coordinates have been specified
      else {
      // Append dummy coordinate for robustness
        coordListLinearProfile->append(0.0f);
        coordListLinearProfile->append(0.0f);

        // Print errormessage
        static uint32_t current_errors = 0;
        if (current_errors < 1) {
          SoDebugError::postWarning("SoLinearProfile::getVertices", "Illegal profile "
                                    "coordinate index specified: %d. Should be within "
                                    "[0, %d]", idx, numcoords - 1);
        }
        current_errors++;
      }
    }

    vertices = (SbVec2f *) coordListLinearProfile->getArrayPtr();
    numvertices = n;
  }
  else {
    vertices = NULL;
    numvertices = 0;
  }
}
