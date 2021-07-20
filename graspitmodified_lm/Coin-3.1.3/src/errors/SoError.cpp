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
  \class SoError SoError.h Inventor/errors/SoError.h
  \brief The SoError class is the base class for all the error handling classes.
  \ingroup errors

  The default error handler just prints messages on the standard error
  output channel, but this can be overridden by client applications.

  Being able to override the default handler is useful when you want
  to collect error messages upon e.g. model import for later
  presentation of the messages to the user in any custom manner (like
  for instance in a GUI messagebox).

  The SoError class is not designed to be particularly useful for
  "direct use". Within the Coin library it is only used through its
  subclasses.

*/

// *************************************************************************

#include <Inventor/errors/SoErrors.h>

#include "coindefs.h"

#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

#include <Inventor/C/tidbits.h>
#include <Inventor/SbName.h>
#include <Inventor/SoPath.h>
#include <Inventor/engines/SoEngine.h>
#include <Inventor/nodes/SoNode.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

// *************************************************************************

SoType SoError::classTypeId STATIC_SOTYPE_INIT;
SoErrorCB * SoError::callback = NULL; // make use of default cc_error handler
void * SoError::callbackData = NULL;

// *************************************************************************

// "Converter" constructor.
SoError::SoError(const cc_error * error)
{
  cc_error_init(&this->err);
  cc_error_copy(error, &this->err);
}

void
SoError::callbackForwarder(const cc_error * error, void * COIN_UNUSED_ARG(data))
{
  SoError wrappederr(error);

  assert(SoError::callback != NULL);
  (*SoError::callback)(&wrappederr, SoError::callbackData);
}


/*!
  \fn SoError::~SoError()

  The default destructor does nothing.
*/


/*!
  This method takes care of initializing all static data for the class.
*/
void
SoError::initClass(void)
{
  SoError::callback = defaultHandlerCB;
  SoError::callbackData = NULL;
  SoError::classTypeId =
    SoType::createType(SoType::badType(), SbName("Error"));
}

/*!
  This static method initializes all the SoError classes.
*/
void
SoError::initClasses(void)
{
  SoError::initClass();
  SoDebugError::initClass();
  SoMemoryError::initClass();
  SoReadError::initClass();
}

/*!
  This static method returns the SoType for this class.

  \sa getTypeId()
*/
SoType
SoError::getClassTypeId(void)
{
  return SoError::classTypeId;
}

// FIXME: grab better version of getTypeId() doc from SoBase, SoAction
// and / or SoDetail. 20010913 mortene.
/*!
  This method returns the SoType of a particular object instance.

  \sa getClassTypeId()
*/
SoType
SoError::getTypeId(void) const
{
  return SoError::classTypeId;
}

/*!
  This method returns \c TRUE if the error instance is of - or derived
  from - \a type, and \c FALSE otherwise.
*/
SbBool
SoError::isOfType(const SoType type) const
{
  const SoType myType = this->getTypeId();
  return (myType == type) ? TRUE : myType.isDerivedFrom(type);
}

/*!
  This method sets the error handler callback for messages posted via
  this class.

  Note that this will not override the error/debug message handler for
  subclasses, these will have to be overrided by calling the subclass'
  setHandlerCallback() method.

  \sa defaultHandlerCB()
*/
void
SoError::setHandlerCallback(SoErrorCB * const function, void * const data)
{
  if (SoError::callback == SoError::defaultHandlerCB ||
      SoError::callback == NULL) {
    // The user is overriding the default handler, so set up a
    // "converter" callback function that makes an SoError out of an
    // cc_error and forwards control to the callback function given as
    // an argument to setHandlerCallback().
    cc_error_set_handler_callback(
       reinterpret_cast<cc_error_cb *>(SoError::callbackForwarder),
       NULL);
  }

  SoError::callback = function;
  SoError::callbackData = data;
}

/*!
  Returns the error handler callback for messages posted via this
  class.
*/
SoErrorCB *
SoError::getHandlerCallback(void)
{
  return SoError::callback;
}

/*!
  This method returns the pointer used for passing data back to the
  callback handler method.
*/
void *
SoError::getHandlerData(void)
{
  return SoError::callbackData;
}

/*!
  This method returns an SbString containing error info from the given
  error instance.
*/
const SbString &
SoError::getDebugString(void) const
{
  // Cast away constness and fetch value from underlying cc_error instance.
  const_cast<SbString &>(this->debugstring) = cc_string_get_text(cc_error_get_debug_string(&this->err));
  return this->debugstring;
}


/*!
  This method posts an error message.  The \a format string and the
  trailing aguments should follow the printf() standard.
*/
void
SoError::post(const char * const format, ...)
{
  va_list args;
  va_start(args, format);
  cc_error_post_arglist(format, args);
  va_end(args);
}

/*!
  Constructs a string identifying the \a node with name (if available)
  and memory pointer.
*/
SbString
SoError::getString(const SoNode * const node)
{
  SbString str;
  SoError::generateBaseString(str, node, "node");
  return str;
}

/*!
  Constructs a string identifying the \a path with name (if available)
  and memory pointer.
*/
SbString
SoError::getString(const SoPath * const path)
{
  SbString str;
  SoError::generateBaseString(str, path, "path");
  return str;
}

/*!
  Constructs a string identifying the \a engine with name (if available)
  and memory pointer.
*/
SbString
SoError::getString(const SoEngine * const engine)
{
  SbString str;
  SoError::generateBaseString(str, engine, "engine");
  return str;
}

/*!
  Contains the default code for handling error strings.

  Default treatment of an error message is to print it out on the
  standard error file handle.
*/
void
SoError::defaultHandlerCB(const SoError * error, void * data)
{
  cc_error_default_handler_cb(&error->err, data);
}

/*!
  This is just a convenience wrapper around the getHandlerCallback() and
  getHandlerData() methods.
*/
SoErrorCB *
SoError::getHandler(void * & data) const
{
  data = SoError::callbackData;
  return SoError::callback;
}

/*!
  Replace the latest stored debug string with \a str.
*/
void
SoError::setDebugString(const char * const str)
{
  cc_error_set_debug_string(&this->err, str);
}

/*!
  Add \a str at the end of the currently stored debug string.
*/
void
SoError::appendToDebugString(const char * const str)
{
  cc_error_append_to_debug_string(&this->err, str);
}

/*!
  This method calls the appropriate handler for an error instance.  All
  error handling goes through this method, and is therefore a good
  candidate for a debugger breakpoint.
*/
void
SoError::handleError(void)
{
  // If this assert() statement were present earlier, it could have
  // saved me half the day of debugging today.. 20021015 mortene.
  assert((SoError::classTypeId != SoType::badType()) &&
         "SoError attempted used before class was initialized");

  cc_error_handle(&this->err);
}

/*!
  This method is used by the getString methods.  It just generates a
  '<what> named "<name>" at address <address>' string.
*/
void
SoError::generateBaseString(SbString & str,
                            const SoBase * const base,
                            const char * const what)
{
  str.sprintf("%s named \"%s\" at address %p",
              what, base->getName().getString(), base);
}
