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

#include "scxml/ScXMLHistory.h"

#include <assert.h>
#include <string.h>

#include <boost/scoped_ptr.hpp>

#include <Inventor/scxml/ScXML.h>

#include "scxml/ScXMLTransition.h"
#include "scxml/ScXMLCommonP.h"

// *************************************************************************

/*!
  \class ScXMLHistory ScXMLHistory.h Inventor/scxml/ScXMLHistory.h
  \brief Implementation of the &lt;history&gt; SCXML element.

  The &lt;history&gt; element is not supported, other than as a dummy
  state.
  
  \since Coin 3.0
  \ingroup scxml
*/

class ScXMLHistory::PImpl {
public:
  boost::scoped_ptr<ScXMLTransition> transitionptr;


};

#define PRIVATE(obj) ((obj)->pimpl)

SCXML_OBJECT_SOURCE(ScXMLHistory);

void
ScXMLHistory::initClass(void)
{
  SCXML_OBJECT_INIT_CLASS(ScXMLHistory, ScXMLObject, SCXML_DEFAULT_NS, "history");
}

ScXMLHistory::ScXMLHistory(void)
  : id(NULL), type(NULL)
{
}

ScXMLHistory::~ScXMLHistory(void)
{
  this->setIdAttribute(NULL);
  this->setTypeAttribute(NULL);
}

void
ScXMLHistory::setIdAttribute(const char * idstr)
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

void
ScXMLHistory::setTypeAttribute(const char * typestr)
{
  if (this->type && this->type != this->getXMLAttribute("type")) {
    delete [] this->type;
  }
  this->type = NULL;
  if (typestr) {
    this->type = new char [ strlen(typestr) + 1 ];
    strcpy(this->type, typestr);
  }
}

SbBool
ScXMLHistory::handleXMLAttributes(void)
{
  if (!inherited::handleXMLAttributes()) return FALSE;

  this->id = const_cast<char *>(this->getXMLAttribute("id"));
  this->type = const_cast<char *>(this->getXMLAttribute("type"));

  if (!this->id) { return FALSE; }

  if (this->type) {
    if ((strcmp(this->type, "deep") != 0) &&
        (strcmp(this->type, "shallow") != 0)) {
      return FALSE;
    }
  }

  return TRUE;
}

// *************************************************************************

SCXML_SINGLE_OBJECT_API_IMPL(ScXMLHistory, ScXMLTransition, PRIVATE(this)->transitionptr, Transition);

#undef PRIVATE
