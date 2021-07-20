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
  \class SoPolygonOffset SoPolygonOffset.h Inventor/nodes/SoPolygonOffset.h
  \brief The SoPolygonOffset class is a node type for "layering" rendering primitives.
  \ingroup nodes

  A common problem with realtime 3D rendering systems is that rendered
  primitives which are at approximately the same depth with regard to
  the camera viewpoint will appear to flicker. I.e.: from one angle
  one primitive will appear to be closer, while at another angle,
  another primitive will appear closer. When this happens, the
  rendered graphics at that part of the scene will of course look a
  lot less visually pleasing.

  One common situation where this problem often occurs is when you
  attempt to put a wireframe grid as an outline on top of filled
  polygons.

  The cause of the problem described above is that the Z-buffer of any
  render system has a limited resolution, often at 16, 24 or 32
  bits. Because of this, primitives which are close will sometimes get
  the \e same depth value in the Z-buffer, even though they are \a not
  actually at the same depth-coordinate.

  To rectify the flickering problem, this node can be inserted in the
  scene graph at the proper place(s) to explicitly define how
  polygons, lines and/or points should be offset with regard to other
  primitives.

  As for the details of how the SoPolygonOffset::factor and
  SoPolygonOffset::units should be set, we quote the OpenGL
  documentation:

  \verbatim

      The value of the offset is

          factor * DZ + r * units

      where DZ is a measurement of the change in depth relative to the
      screen area of the polygon, and r is the smallest value that is
      guaranteed to produce a resolvable offset for a given
      implementation. The offset is added before the depth test is
      performed and before the value is written into the depth buffer.

  \endverbatim

  One word of notice with regard to the above quote from the OpenGL
  documentation: it doesn't really make sense to set "factor" and
  "units" to values with different signs, i.e. "factor" to a negative
  value and "units" to a positive value, or vice versa.

  The pixels would then be "pushed back" in z-order by one part of the
  equation, but at the same time be "pushed forward" by the other part
  of the equation. This would most likely give very inconsistent
  results, but which may at first look ok.

  We mention this potential for making a mistake, as it seems to be a
  quite common error.


  Below is a simple, correct usage example:

  \verbatim
  #Inventor V2.1 ascii
  
  Separator {
     Coordinate3 { point [ -1 -1 0, 1 -1 0, 1 1 0, -1 1 0 ] }
     
     Separator {
        BaseColor { rgb 1 1 0 }
        # needs to draw polygons-as-line, and not "real" lines -- see
        # documentation below on why this is so:
        DrawStyle { style LINES }
        # draw two triangles, to get a line crossing the face of the
        # polygon:
        IndexedFaceSet { coordIndex [ 0, 1, 2, -1, 0, 2, 3 -1 ] }
     }
     
     PolygonOffset {
        styles FILLED
        factor 1.0
        units 1.0
     }
     
     BaseColor { rgb 0 0.5 0 }
     FaceSet { numVertices [ 4 ] }
  }
  \endverbatim

  Without the polygonoffset node in the above example, the lines may
  look irregularly stippled with some graphics card drivers, as parts
  of it will show through the faceset, others not. This happen on
  seemingly random parts, as the z-buffer floating point calculations
  will be fickle with regard to whether or not the polygon or the line
  will be closer to the camera.

  See the API documentation of the SoPolygonOffset::styles field below
  for a discussion of one important limitation of OpenGL's Z-buffer
  offset mechanism: it only works with polygons or polygons rendered
  in line or point mode, using the SoDrawStyle::style field.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    PolygonOffset {
        factor 1
        units 1
        styles FILLED
        on TRUE
    }
  \endcode

  \since TGS Inventor 2.5
  \since Coin 1.0
*/

// *************************************************************************

#include <Inventor/nodes/SoPolygonOffset.h>

#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLPolygonOffsetElement.h>
#include <Inventor/elements/SoOverrideElement.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \enum SoPolygonOffset::Style

  Enumeration of the rendering primitives which can be influenced by
  an SoPolygonOffset node.
*/

/*!
  \var SoSFFloat SoPolygonOffset::factor

  Offset multiplication factor. Scales the variable depth in the
  z-buffer of the rendered primitives.

  See SoPolygonOffset's main class documentation above for detailed
  information on how the factor value is used.

  Default value is 1.0.
*/
/*!
  \var SoSFFloat SoPolygonOffset::units

  Offset translation multiplication factor. Will be multiplied with
  the value which represents the smallest discrete step that can be
  distinguished with the underlying Z-buffer resolution.

  See SoPolygonOffset's main class documentation above for detailed
  information on how the units value is used.

  Note that positive values will push geometry "away" into the
  Z-buffer, while negative values will "move" geometry closer.

  Default value is 1.0.
*/
/*!
  \var SoSFBitMask SoPolygonOffset::styles

  The rendering primitive type to be influenced by this node. This is
  a bitmask variable, so you can select several primitive types (out
  of filled polygons, lines and points) be influenced by the offset at
  the same time.

  There is one very important OpenGL limitation to know about in this
  regard: Z-buffer offsetting can \e only be done for either polygons,
  or for \e polygons rendered \e as \e lines or \e as \e points.

  So attempts at using this node to offset e.g. SoLineSet /
  SoIndexedLineSet or SoPointSet primitives will \e not work.

  See the comments in the scene graph below for a detailed example on
  what SoPolygonOffset can and can not do:

  \code
  #Inventor V2.1 ascii
  
  Separator {
     # render polygon:
  
     Coordinate3 { point [ -1.1 -1.1 0, 1.1 -1.1 0, 1.1 1.1 0, -1.1 1.1 0 ] }
     BaseColor { rgb 0 0.5 0 }
     FaceSet { numVertices [ 4 ] }
  
     # offset polygon-as-lines to be in front of above polygon:
  
     PolygonOffset {
        styles LINES
        factor -2.0
        units 1.0
     }
     
     # render lines:
  
     Coordinate3 { point [ -1 -1 0, 1 -1 0, 1 1 0, -1 1 0 ] }
     BaseColor { rgb 1 1 0 }
  
     Switch {
        # change this to '0' to see how glPolygonOffset() does *not* work
        # with "true" lines
        whichChild 1
  
        DEF child0 Group {
           # can *not* be offset
           IndexedLineSet { coordIndex [ 0, 1, 2, 3, 0, 2, -1, 1, 3 -1 ]
           }
        }
  
        DEF child1 Group {
           # will be offset
           DrawStyle { style LINES }
           FaceSet { numVertices [ 4 ] }
        }
     }
  }
  \endcode

  Field default value is SoPolygonOffset::FILLED.
*/

/*!
  \var SoSFBool SoPolygonOffset::on

  Whether the offset is on or off. Default is for SoPolygonOffset::on
  to be \c TRUE.
*/


// *************************************************************************

SO_NODE_SOURCE(SoPolygonOffset);

/*!
  Constructor.
*/
SoPolygonOffset::SoPolygonOffset(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoPolygonOffset);

  SO_NODE_ADD_FIELD(factor, (1.0f));
  SO_NODE_ADD_FIELD(units, (1.0f));
  SO_NODE_ADD_FIELD(styles, (SoPolygonOffset::FILLED));
  SO_NODE_ADD_FIELD(on, (TRUE));

  SO_NODE_DEFINE_ENUM_VALUE(Style, FILLED);
  SO_NODE_DEFINE_ENUM_VALUE(Style, LINES);
  SO_NODE_DEFINE_ENUM_VALUE(Style, POINTS);
  SO_NODE_SET_SF_ENUM_TYPE(styles, Style);
}

/*!
  Destructor.
*/
SoPolygonOffset::~SoPolygonOffset()
{
}

// doc in super
void
SoPolygonOffset::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoPolygonOffset, SO_FROM_INVENTOR_2_5|SO_FROM_COIN_1_0);

  SO_ENABLE(SoCallbackAction, SoPolygonOffsetElement);
  SO_ENABLE(SoGLRenderAction, SoGLPolygonOffsetElement);
}


void
SoPolygonOffset::doAction(SoAction * action)
{
  SoState * state = action->getState();
  
  if (SoOverrideElement::getPolygonOffsetOverride(state)) return;
  
  float factorval, units_val;
  SoPolygonOffsetElement::Style styles_val;
  SbBool offset_on;
  
  factorval = this->factor.getValue();
  units_val = this->units.getValue();
  styles_val = (SoPolygonOffsetElement::Style)this->styles.getValue();
  offset_on = this->on.getValue();
  
  SoPolygonOffsetElement::set(action->getState(),
                              this,
                              factorval,
                              units_val,
                              styles_val,
                              offset_on);
  
  if (this->isOverride()) {
    SoOverrideElement::setPolygonOffsetOverride(state, this, TRUE);
  }
}

void
SoPolygonOffset::callback(SoCallbackAction * action)
{
  SoPolygonOffset::doAction((SoAction *)action);
}

void
SoPolygonOffset::GLRender(SoGLRenderAction * action)
{
  SoPolygonOffset::doAction((SoAction *)action);
}
