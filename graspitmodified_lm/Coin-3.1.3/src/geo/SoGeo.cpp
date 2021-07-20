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

#include "SoGeo.h"
#include <Inventor/nodes/SoGeoOrigin.h>
#include <Inventor/nodes/SoGeoLocation.h>
#include <Inventor/nodes/SoGeoSeparator.h>
#include <Inventor/nodes/SoGeoCoordinate.h>
#include <Inventor/elements/SoGeoElement.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SbVec3d.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbDPMatrix.h>
#include <Inventor/errors/SoDebugError.h>

#include "SbGeoProjection.h"
#include "SbUTMProjection.h"
#include "SbGeoAngle.h"
#include "SbGeoEllipsoid.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void
SoGeo::init(void)
{
  SoGeoElement::initClass();
  SoGeoOrigin::initClass();
  SoGeoLocation::initClass();
  SoGeoSeparator::initClass();
  SoGeoCoordinate::initClass();
}

// return zone number stored in string, or 0 if parsing failed
static int find_utm_zone(const SbString & s)
{
  if (s.getLength() < 2) return 0;
  if (s[0] != 'Z' && s[1] != 'z') return 0;

  return (int) atoi(s.getString()+1);
}

static SbDPMatrix find_coordinate_system(const SbString * system,
                                         const int numsys,
                                         const SbVec3d & coords)
{
  SbVec3d p;
  if (system[0] == "GC") {
    p = coords;
  }
  else {
    double latitude, longitude, elev;

    if (system[0] == "UTM") {
      SbUTMProjection proj(find_utm_zone(system[1]), SbGeoEllipsoid("WGS84"));
      SbGeoAngle lat, lng;

      proj.unproject(coords[0], coords[1], &lat, &lng);

      latitude = lat.rad();
      longitude = lng.rad();
      elev = coords[2];
    }
    else if (system[0] == "GD") {
      latitude = coords[0] * M_PI / 180.0;
      longitude = coords[1] * M_PI / 180.0;
      elev = coords[2];
    }
    else {
      assert(0 && "not supported");
      latitude = longitude = elev = 0.0;
    }

    // formula based on http://en.wikipedia.org/wiki/Geodetic_system

    double a = 6378137.0; // earth semimajor axis in meters
    double f = 1.0/298.257223563; // reciprocal flattening
    double e2 = 2*f - f*f; // eccentricity squared

    double sinlat = sin(latitude);
    double coslat = cos(latitude);
    double chi = sqrt(1.0 - e2 * (sinlat*sinlat));


    p[0] = (a / chi + elev) * coslat * cos(longitude);
    p[1] = (a / chi + elev) * coslat * sin(longitude);
    p[2] = (a * (1.0-e2)/ chi + elev) * sinlat;


#if 0 // for debugging
    SbUTMProjection utm(17,  SbGeoEllipsoid("WGS84"));
    double east, north;
    SbGeoAngle lat(latitude);
    SbGeoAngle lng(longitude);

    utm.project(lat,lng, &east, &north);

    fprintf(stderr,"zone 17 coords: %g %g\n",
            east, north);
#endif // debugging

  }

  SbVec3d Z = p;
  (void) Z.normalize();

  // FIXME: handle the case when origin is at the north or south pole
  SbVec3d Y(0.0, 0.0, 1.0);
  SbVec3d X = Y.cross(Z);
  (void) X.normalize();
  Y = Z.cross(X);
  (void) Y.normalize();

  SbDPMatrix m = SbDPMatrix::identity();

  for (int i = 0; i < 3; i++) {
    m[0][i] = X[i];
    m[1][i] = Y[i];
    m[2][i] = Z[i];
    m[3][i] = p[i];
  }

#if 0 // for debugging
  SbGeoAngle lat(latitude);
  SbGeoAngle lng(longitude);

  fprintf(stderr,"coordinate system matrix: (%f %f) (%d %d %.2f, %d %d %.2f\n",
          latitude*180/M_PI, longitude*180.0/M_PI,
          lat.deg(), lat.minutes(), lat.seconds(),
          lng.deg(), lng.minutes(), lng.seconds());
  m.print(stderr);
  fprintf(stderr,"\n");
#endif // debugging

  return m;
}

//
// Currently not used. Kept here since it might be useful to find and
// UTM zone from lat/long
//
static SbUTMProjection find_utm_projection(const SbString * system,
                                           const int numsystem,
                                           const SbVec3d & coords,
                                           SbVec3d & projcoords)
{
  if (system[0] != "UTM") { // find an UTM zone to project into
    assert(system[0] == "GD");
    // project to an UTM zone

    double degree = coords[1];
    int zone = int ((degree + 180.0) / 6.0);

    SbGeoAngle lat(coords[0], 0.0, 0.0, 'N');
    SbGeoAngle lng(coords[1], 0.0, 0.0, 'N');

    SbUTMProjection proj(zone, SbGeoEllipsoid("WGS84"));

    double east, north;
    proj.project(lat, lng, &east, &north);

    projcoords[0] = east;
    projcoords[1] = north;
    projcoords[2] = coords[2];

    return proj;
  }

  // just return the original UTM Zone and coords
  projcoords = coords;
  return SbUTMProjection(find_utm_zone(system[1]), SbGeoEllipsoid("WGS84"));
}

SbMatrix
SoGeo::calculateTransform(const SbString * originsystem,
                          const int numoriginsys,
                          const SbVec3d & geocoords,
                          const SbString * localsystem,
                          const int numlocalsys,
                          const SbVec3d & localcoords)
{
  // start on 2; the first index is always the projection type, and if UTM the second should always be a zone
  for (int i = 2; i < numoriginsys; i++) {
    if (originsystem[i] == "FLAT") {
      SbMatrix m;
      m.makeIdentity();
      bool valid = (originsystem[0] == "UTM" && 
                    localsystem[0] == originsystem[0] && 
                    localsystem[1] == originsystem[1]);
      if (!valid){
        SoDebugError::post("SoGeo::calculateTransform", "FLAT projections only supported within the same UTM zone");
        return m;      
      }
      m.setTranslate(SbVec3f(localcoords - geocoords));
      return m;
    }
  }

  SbDPMatrix om = find_coordinate_system(originsystem, numoriginsys, geocoords);
  SbDPMatrix lm = find_coordinate_system(localsystem, numlocalsys, localcoords);
  SbDPMatrix r = lm * om.inverse();

  // transform to a single precision matrix.
  return SbMatrix(static_cast<float>(r[0][0]), static_cast<float>(r[0][1]),
                  static_cast<float>(r[0][2]), static_cast<float>(r[0][3]),
                  static_cast<float>(r[1][0]), static_cast<float>(r[1][1]),
                  static_cast<float>(r[1][2]), static_cast<float>(r[1][3]),
                  static_cast<float>(r[2][0]), static_cast<float>(r[2][1]),
                  static_cast<float>(r[2][2]), static_cast<float>(r[2][3]),
                  static_cast<float>(r[3][0]), static_cast<float>(r[3][1]),
                  static_cast<float>(r[3][2]), static_cast<float>(r[3][3]));
}
