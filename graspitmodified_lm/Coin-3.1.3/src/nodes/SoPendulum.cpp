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
  \class SoPendulum SoPendulum.h Inventor/nodes/SoPendulum.h
  \brief The SoPendulum class is used to create oscillating rotations.
  \ingroup nodes

  A smooth transition between rotation0 and rotation1 is created using
  a cosine function. In the beginning of the cycle, rotation0 is
  used. Halfway through the cycle, the resulting rotation equals
  rotation1, and at the end of the cycle, we're at rotation0 again.

  <b>FILE FORMAT/DEFAULTS:</b>
  \code
    Pendulum {
        rotation 0 0 1  0
        rotation0 0 0 1  0
        rotation1 0 0 1  0
        speed 1
        on TRUE
    }
  \endcode
*/

// *************************************************************************

#include <Inventor/nodes/SoPendulum.h>

#include <Inventor/SbVec3f.h>
#include <Inventor/SoOutput.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/engines/SoCalculator.h>
#include <Inventor/engines/SoElapsedTime.h>
#include <Inventor/engines/SoInterpolateRotation.h>

#include "nodes/SoSubNodeP.h"

// *************************************************************************

/*!
  \var SoSFRotation SoPendulum::rotation0
  The first rotation limit of the interpolation.
*/
/*!
  \var SoSFRotation SoPendulum::rotation1
  The other rotation limit of the interpolation.
*/
/*!
  \var SoSFFloat SoPendulum::speed
  Speed in cycles per second. Defaults to 1.
*/
/*!
  \var SoSFBool SoPendulum::on
  Toggles animation on or off. Defaults to being \c on.
*/

// *************************************************************************

SO_NODE_SOURCE(SoPendulum);

/*!
  Constructor.
*/
SoPendulum::SoPendulum(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoPendulum);

  SO_NODE_ADD_FIELD(rotation0, (SbRotation(SbVec3f(0.0f, 0.0f, 1.0f), 0.0f)));
  SO_NODE_ADD_FIELD(rotation1, (SbRotation(SbVec3f(0.0f, 0.0f, 1.0f), 0.0f)));
  SO_NODE_ADD_FIELD(speed, (1.0f));
  SO_NODE_ADD_FIELD(on, (TRUE));

  this->interpolator = new SoInterpolateRotation;
  this->interpolator->ref();
  this->calculator = new SoCalculator;
  this->calculator->ref();
  this->timer = new SoElapsedTime;
  this->timer->ref();

  this->calculator->expression = "oa = (1.0 - cos(a*b*2*M_PI)) * 0.5";
  this->calculator->a.connectFrom(&this->timer->timeOut);
  this->timer->on.connectFrom(&this->on);
  this->calculator->b.connectFrom(&this->speed);
  this->interpolator->input0.connectFrom(&this->rotation0);
  this->interpolator->input1.connectFrom(&this->rotation1);
  this->interpolator->alpha.connectFrom(&this->calculator->oa);
  this->rotation.connectFrom(&this->interpolator->output, TRUE);
}

/*!
  Destructor.
*/
SoPendulum::~SoPendulum()
{
  this->interpolator->unref();
  this->calculator->unref();
  this->timer->unref();
}

// Doc from superclass.
void
SoPendulum::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoPendulum, SO_FROM_INVENTOR_1);
}

// Documented in superclass.
void
SoPendulum::write(SoWriteAction * action)
{
  // Overridden to not write out internal engine connections.

  SoOutput * out = action->getOutput();

  // Decouple connections to/from internal engines to avoid them being
  // written. (Only done at first pass.)
  if (out->getStage() == SoOutput::COUNT_REFS)
    this->deconnectInternalEngine();

  inherited::write(action);

  // Reenable all connections to/from internal engine. (Only done at
  // last pass.)
  if (out->getStage() == SoOutput::WRITE)
    this->reconnectInternalEngine();
}

// FIXME: I _think_ we made a mistake when overriding SoNode::copy()
// and making it virtual. See FIXME-comment above
// SoBlinker::copy(). 20011220 mortene.

// Overridden to decouple and reconnect engine around copy operation.
SoNode *
SoPendulum::copy(SbBool copyconnections) const
{
  // Decouple connections to/from internal engines to avoid them being
  // copied.
  ((SoPendulum *)this)->deconnectInternalEngine();

  SoPendulum * cp = (SoPendulum *)inherited::copy(copyconnections);

  // Reenable all connections to/from internal engines.
  ((SoPendulum *)this)->reconnectInternalEngine();

  return cp;
}

// Remove connections to and from internal engines.
void
SoPendulum::deconnectInternalEngine(void)
{
  // Do this first, to avoid field being set due to subsequent engine
  // input value change.
  this->rotation.disconnect(&this->interpolator->output);

  this->timer->on.disconnect(&this->on);
  this->timer->on = FALSE;
  this->calculator->b.disconnect(&this->speed);
  this->interpolator->input0.disconnect(&this->rotation0);
  this->interpolator->input1.disconnect(&this->rotation1);
}


// Reenable all connections to/from internal engines.
void
SoPendulum::reconnectInternalEngine(void)
{
  this->timer->on.connectFrom(&this->on);
  this->calculator->b.connectFrom(&this->speed);
  this->interpolator->input0.connectFrom(&this->rotation0);
  this->interpolator->input1.connectFrom(&this->rotation1);
  this->interpolator->alpha.connectFrom(&this->calculator->oa);

  this->rotation.connectFrom(&this->interpolator->output, TRUE);
}
