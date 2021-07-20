#ifndef COIN_SCXMLONENTRY_H
#define COIN_SCXMLONENTRY_H

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

#include <Inventor/scxml/ScXMLObject.h>
#include <Inventor/tools/SbPimplPtr.h>

class ScXMLInvoke;
class ScXMLEvent;
class ScXMLStateMachine;

class ScXMLOnEntry : public ScXMLObject {
  typedef ScXMLObject inherited;
  SCXML_OBJECT_HEADER(ScXMLOnEntry);

public:
  static void initClass(void);

  ScXMLOnEntry(void);
  virtual ~ScXMLOnEntry(void);

  // executable content
  virtual int getNumInvokes(void) const;
  virtual ScXMLInvoke * getInvoke(int idx) const;
  virtual void addInvoke(ScXMLInvoke * invoke);
  virtual void removeInvoke(ScXMLInvoke * invoke);
  virtual void clearAllInvokes(void);

  // invoke
  virtual void invoke(ScXMLStateMachine * statemachine);
  
private:
  ScXMLOnEntry(const ScXMLOnEntry & rhs); // N/A
  ScXMLOnEntry & operator = (const ScXMLOnEntry & rhs); // N/A

  class PImpl;
  SbPimplPtr<PImpl> pimpl;

}; // ScXMLOnEntry

#endif // !COIN_SCXMLONENTRY_H
