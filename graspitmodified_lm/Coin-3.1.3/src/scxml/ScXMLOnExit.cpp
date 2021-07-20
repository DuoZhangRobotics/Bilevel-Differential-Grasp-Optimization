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

#include "scxml/ScXMLOnExit.h"

#include <assert.h>
#include <algorithm>
#include <vector>

#include <Inventor/scxml/ScXML.h>
#include <Inventor/scxml/ScXMLInvoke.h>

#include "scxml/ScXMLCommonP.h"

// *************************************************************************

/*!
  \class ScXMLOnExit ScXMLOnExit.h Inventor/scxml/ScXMLOnExit.h
  \brief Implementation of the &lt;onexit&gt; SCXML element.

  \since Coin 3.0
  \ingroup scxml
*/

class ScXMLOnExit::PImpl {
public:
  std::vector<ScXMLInvoke *> invokelist;

  ~PImpl(void)
  {
    std::vector<ScXMLInvoke *>::iterator it = this->invokelist.begin();
    while (it != this->invokelist.end()) {
      delete *it;
      ++it;
    }
    this->invokelist.clear();
  }

};

#define PRIVATE(obj) ((obj)->pimpl)

SCXML_OBJECT_SOURCE(ScXMLOnExit);

void
ScXMLOnExit::initClass(void)
{
  SCXML_OBJECT_INIT_CLASS(ScXMLOnExit, ScXMLObject, SCXML_DEFAULT_NS, "onexit");
}

ScXMLOnExit::ScXMLOnExit(void)
{
}

ScXMLOnExit::~ScXMLOnExit(void)
{
}

// *************************************************************************
// executable content

SCXML_LIST_OBJECT_API_IMPL(ScXMLOnExit, ScXMLInvoke, PRIVATE(this)->invokelist, Invoke, Invokes);

void
ScXMLOnExit::invoke(ScXMLStateMachine * statemachine)
{
  std::vector<ScXMLInvoke *>::const_iterator it =
    PRIVATE(this)->invokelist.begin();
  while (it != PRIVATE(this)->invokelist.end()) {
    (*it)->invoke(statemachine);
    ++it;
  }
}

#undef PRIVATE
