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
  \class SoWWWAnchor SoWWWAnchor.h Inventor/nodes/SoWWWAnchor.h
  \brief The SoWWWAnchor class adds URL callbacks to the highlighted geometry.
  \ingroup nodes

  In addition to highlighting geometry under the cursor, the application
  programmer can set callbacks. It is possible to set one callback for
  picking, the fetch callback, and one callback for highlighting.

  \verbatim
  #Inventor V2.1 ascii
  
  WWWAnchor {
     name "http://www.coin3d.org/Coin/egg.iv"
     description "Easter Egg"
  
     Separator {
        Transform { scaleFactor 0.8 1.2 0.8 }
        Sphere { }
     }
  }
  \endverbatim

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    WWWAnchor {
        renderCaching AUTO
        boundingBoxCaching AUTO
        renderCulling AUTO
        pickCulling AUTO
        color 0.3 0.3 0.3
        style EMISSIVE
        mode AUTO
        name "<Undefined URL>"
        description ""
        map NONE
    }
  \endcode

  \since Inventor 2.1
*/


#include <Inventor/nodes/SoWWWAnchor.h>

#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/SoPickedPoint.h>

#include "nodes/SoSubNodeP.h"
#include "coindefs.h"
#include "tidbitsp.h"


/*!
  \enum SoWWWAnchor::Mapping
  Enum that says how a picked node's position should be mapped to the URL.
*/
/*!
  \var SoWWWAnchor::Mapping SoWWWAnchor::NONE
  The position of the picked node is not mapped to the URL.
*/
/*!
  \var SoWWWAnchor::Mapping SoWWWAnchor::POINT
  The position of the picked node is mapped to the URL as object space
  coordinates, adding a parameter string to the end of the URL. To
  assure that the URL works with all browsers, the coordinates are
  divided by commas sent as the hex representation.

  If a model by the name of sim.wrl resided at www.coin3d.org and the
  picked point had the coordinates [1.5, 10, 6.77], the resulting URL
  would be "http://www.coin3d.org/sim.wrl?1.5%2c10%2c6.77".
*/


/*!
  \var SoSFString SoWWWAnchor::name
  The name of the URL which the anchor points to.
*/
/*!
  \var SoSFString SoWWWAnchor::description
  The description of the URL.
*/
/*!
  \var SoSFEnum SoWWWAnchor::map
  Enum describing how a node's position should be mapped to the URL.
*/

// *************************************************************************

#ifndef DOXYGEN_SKIP_THIS

class SoWWWAnchorP {
 public:
  SoWWWAnchorP(SoWWWAnchor * ownerptr) {
    this->owner = ownerptr;
    this->fullname = "";
  }
  SoWWWAnchor * owner;
  SbString fullname;

  static SoWWWAnchorCB * fetchfunc;
  static void * fetchdata;
  static SoWWWAnchorCB * highlightfunc;
  static void * highlightdata;

  static void atexit_cleanup() {
    SoWWWAnchorP::fetchfunc = NULL;
    SoWWWAnchorP::fetchdata = NULL;
    SoWWWAnchorP::highlightfunc = NULL;
    SoWWWAnchorP::highlightdata = NULL;
  }
};

// static members
SoWWWAnchorCB * SoWWWAnchorP::fetchfunc = NULL;
void * SoWWWAnchorP::fetchdata = NULL;
SoWWWAnchorCB * SoWWWAnchorP::highlightfunc = NULL;
void * SoWWWAnchorP::highlightdata = NULL;

#endif // DOXYGEN_SKIP_THIS

SO_NODE_SOURCE(SoWWWAnchor);

#define PRIVATE(p) (p->pimpl)

/*!
  Constructor.
*/
SoWWWAnchor::SoWWWAnchor()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoWWWAnchor);

  PRIVATE(this) = new SoWWWAnchorP(this);

  SO_NODE_ADD_FIELD(name, ("<Undefined URL>"));
  SO_NODE_ADD_FIELD(description, (""));
  SO_NODE_ADD_FIELD(map, (NONE));

  SO_NODE_DEFINE_ENUM_VALUE(Map, NONE);
  SO_NODE_DEFINE_ENUM_VALUE(Map, POINT);
  SO_NODE_SET_SF_ENUM_TYPE(map, Map);
}

/*!
  Destructor.
*/
SoWWWAnchor::~SoWWWAnchor()
{
  delete PRIVATE(this);
}

// doc in super
void
SoWWWAnchor::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoWWWAnchor, SO_FROM_INVENTOR_2_1|SoNode::VRML1);
  coin_atexit((coin_atexit_f*)SoWWWAnchorP::atexit_cleanup, CC_ATEXIT_NORMAL);
}


/*!
  Sets the full URL to \a url. If this is set, this URL will be used in
  callbacks instead of the URL set in SoWWWAnchor::name.

  \sa SoWWWAnchor::getFullURLName()
 */
void
SoWWWAnchor::setFullURLName(const SbString & url)
{
  PRIVATE(this)->fullname = url;
}

/*!
  Returns the full URL if it's set by
  SoWWWAnchor::setFullURLName(). Otherwise the contents of
  SoWWWAnchor::name is returned.

  \sa SoWWWAnchor::setFullURLName()
 */
const SbString &
SoWWWAnchor::getFullURLName(void)
{
  if (PRIVATE(this)->fullname.getLength() > 0) {
    return PRIVATE(this)->fullname;
  }

  return this->name.getValue();
}

// documented in superclass
void
SoWWWAnchor::handleEvent(SoHandleEventAction * action)
{
  const SoEvent * event = action->getEvent();
  if (event->isOfType(SoMouseButtonEvent::getClassTypeId()) &&
      SoWWWAnchorP::fetchfunc) {
    const SoMouseButtonEvent * mbevent = (SoMouseButtonEvent*)event;
    if (SoMouseButtonEvent::isButtonPressEvent(mbevent,
                                               SoMouseButtonEvent::BUTTON1)) {
      SbString s = this->getFullURLName();
      if (this->map.getValue() == POINT) {
        const SoPickedPoint * pp = action->getPickedPoint();
        const SbVec3f point = pp->getObjectPoint(NULL);
        SbString temp;
        temp.sprintf("?%g%%2c%g%%2c%g", point[0], point[1], point[2]);
        s.operator+=(temp);
      }

      SoWWWAnchorP::fetchfunc(s, SoWWWAnchorP::fetchdata, this);
    }
  }
  inherited::handleEvent(action);
}

/*!
  Sets the callback function \a f that is called when a SoWWWAnchor node is
  clicked on. This callback can among other things be used to provide a
  browser with the URL of this node.

  The callback will be called with the URL, \a userData and a pointer to
  this node as arguments.
 */
void
SoWWWAnchor::setFetchURLCallBack(SoWWWAnchorCB * f, void * userData)
{
  SoWWWAnchorP::fetchfunc = f;
  SoWWWAnchorP::fetchdata = userData;
}

/*!
  Sets the callback function \a f that is called when a SoWWWAnchor node
  is highlighted. This callback can among other things be used to provide
  the user with a visual clue on which URL the node points to, for example
  by showing the URL as a string.

  The callback will be called with the URL, \a userData and a pointer to
  this node as arguments.
 */
void
SoWWWAnchor::setHighlightURLCallBack(SoWWWAnchorCB * f, void * userData)
{
  SoWWWAnchorP::highlightfunc = f;
  SoWWWAnchorP::highlightdata = userData;
}

/*!
  Calls the highlight callback set up with
  SoWWWAnchor::setHighlightURLCallBack().
*/
void
SoWWWAnchor::redrawHighlighted(SoAction * act, SbBool isNowHighlighting)
{
  inherited::redrawHighlighted(act, isNowHighlighting);

  if (SoWWWAnchorP::highlightfunc) {
    SbString s = this->getFullURLName();
    SoWWWAnchorP::highlightfunc(s, SoWWWAnchorP::highlightdata, this);
  }
}

#undef PRIVATE
