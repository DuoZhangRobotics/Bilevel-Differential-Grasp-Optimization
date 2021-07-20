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

#include "scxml/ScXMLInitial.h"

#include <algorithm>
#include <cassert>
#include <cstring>

#include <boost/scoped_ptr.hpp>

#include <Inventor/scxml/ScXML.h>

#include "scxml/ScXMLTransition.h"
#include "scxml/ScXMLCommonP.h"

// *************************************************************************

/*!
  \class ScXMLInitial ScXMLInitial.h Inventor/scxml/ScXMLInitial.h
  \brief Implementation of the &lt;initial&gt; SCXML element.

  \since Coin 3.0
  \ingroup scxml
*/

class ScXMLInitial::PImpl {
public:
  boost::scoped_ptr<ScXMLTransition> transitionptr;

};

SCXML_OBJECT_SOURCE(ScXMLInitial);

#define PRIVATE(obj) ((obj)->pimpl)

void
ScXMLInitial::initClass(void)
{
  SCXML_OBJECT_INIT_CLASS(ScXMLInitial, ScXMLObject, SCXML_DEFAULT_NS, "initial");
}

ScXMLInitial::ScXMLInitial(void)
  : id(NULL)
{
}

ScXMLInitial::~ScXMLInitial(void)
{
  this->setIdAttribute(NULL);
}

void
ScXMLInitial::setIdAttribute(const char * idstr)
{
  if (this->id && this->id != this->getXMLAttribute("id")) {
    delete [] this->id;
  }
  this->id = NULL;
  if (idstr) {
    this->id = new char [ strlen(idstr) + 1 ];
    strcpy(this->id, idstr);
  }
}

// const char * ScXMLInitial::getIdAttrbute(void) const

SbBool
ScXMLInitial::handleXMLAttributes(void)
{
  if (!inherited::handleXMLAttributes()) return FALSE;

  this->id = const_cast<char *>(this->getXMLAttribute("id"));

  if (!this->id) { return FALSE; }

  return TRUE;
}

// *************************************************************************

SCXML_SINGLE_OBJECT_API_IMPL(ScXMLInitial, ScXMLTransition, PRIVATE(this)->transitionptr, Transition);

#undef PRIVATE
