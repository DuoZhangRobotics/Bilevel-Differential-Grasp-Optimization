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
  \class SoDebugError SoDebugError.h Inventor/errors/SoDebugError.h
  \brief The SoDebugError class is the internal debugging message passing mechanism.
  \ingroup errors

  This class basically serves two purposes:

  1) It is the message interface through which error and warning
  conditions are passed to the programmer which is using the Coin
  library API for building applications.  These messages are generated
  when API methods are used in an incorrect manner, or if actions on
  the library has caused it to enter an inconsistent state. Coin
  programmers should then quickly be able to trace down errors in
  their application code.

  For this service to be available, the Coin library must be compiled
  with debugging information. For release builds, most of this
  "programmer's safety net" is removed to gain maximum speed and to
  use minimum space (build the smallest possible library).

  2) Coin application programmers can call the SoDebugError methods
  within their own code as a convenient way of debugging application
  code, replacing the usual cascades of fprintf() calls.


  Coin supports an environment variable to set conditional
  breakpoints.  The COIN_DEBUG_BREAK environment variable can be set
  to any number of functions in the form of a list separated by commas
  or spaces.  The functionnames must be given as
  "classname::functioname" (ie without return type, parenthesis or
  argument types or names).  If a debug message is posted from one of
  those functions, your program will be stopped (using assert(0)).

  This can be useful if you want to get core-dumps or enter a debugger
  whenever for instance a warning or an error is posted from some
  function (e.g. "SbVec3f::normalize"). This feature is only enabled
  in the debug version of Coin.
*/

// *************************************************************************

#include <Inventor/errors/SoDebugError.h>

#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <Inventor/C/errors/debugerror.h>
#include <Inventor/C/tidbits.h>
#include <Inventor/SbName.h>
#include <Inventor/SoType.h>
#include <Inventor/lists/SbList.h>

#include "coindefs.h"
#include "tidbitsp.h"

#ifndef COIN_WORKAROUND_NO_USING_STD_FUNCS
using std::strcmp;
using std::strlen;
using std::strcpy;
using std::strchr;
using std::memcpy;
#endif // !COIN_WORKAROUND_NO_USING_STD_FUNCS

#ifdef HAVE_UNISTD_H
  #include <sys/types.h>
  #include <unistd.h>
#endif

#if defined(_MSC_VER)
  #if _MSC_VER >= COIN_MSVC_8_0_VERSION
    #include <intrin.h>
  #endif
  #if _MSC_VER >= COIN_MSVC_7_0_VERSION
    #define COIN_DEBUGGER_BREAK(x) __debugbreak()
  #else
    #define COIN_DEBUGGER_BREAK(x) __asm { int 3 }
  #endif
#else
  #if defined (MINGW32) || defined (CYGWIN)
    #include <windows.h>
    #define COIN_DEBUGGER_BREAK(x) ::DebugBreak()
  #else
    #if defined (_POSIX_VERSION)
      #include <signal.h>
      #define COIN_DEBUGGER_BREAK(x) ::raise(SIGINT)
    #else
      #define COIN_DEBUGGER_BREAK(x) assert(0 && x)
    #endif
  #endif
#endif

// *************************************************************************

SoType SoDebugError::classTypeId STATIC_SOTYPE_INIT;
SoErrorCB * SoDebugError::callback = SoError::defaultHandlerCB;
void * SoDebugError::callbackData = NULL;

// *************************************************************************

#if COIN_DEBUG

// variables for run-time breakpoints
static int num_breakpoints = 0;
static char ** breakpoints = NULL;

extern "C" {

static void
debug_break_cleanup(void)
{
  for (int i = 0; i < num_breakpoints; i++) {
    delete[] breakpoints[i];
  }
  delete[] breakpoints;
  breakpoints = NULL;
  num_breakpoints = 0;
}

} // extern "C"

#endif // COIN_DEBUG

/*!
  \enum SoDebugError::Severity

  Specifies the available severity levels of the debug messages.
*/


// Documented for parent class.
void
SoDebugError::initClass(void)
{
  SoDebugError::callback = SoError::defaultHandlerCB;
  SoDebugError::callbackData = NULL;
  SoDebugError::classTypeId =
    SoType::createType(SoError::getClassTypeId(), "DebugError");

  // FIXME: the following "tokenizer" code could perhaps be factored
  // out and moved to e.g. SbStringList? 20030820 mortene.

#if COIN_DEBUG
  const char * env = coin_getenv("COIN_DEBUG_BREAK");
  if (env) {
    num_breakpoints = 1;
    const char * ptr = env;
    while (*ptr) {
      if (*ptr == ' ' || *ptr == ',') num_breakpoints++;
      ptr++;
    }
    breakpoints = new char*[num_breakpoints];
    coin_atexit(debug_break_cleanup, CC_ATEXIT_MSG_SUBSYSTEM);
    const size_t envstrlen = strlen(env);
    char * cpy = new char[envstrlen + 1];
    (void)strcpy(cpy, env);
    ptr = cpy;
    const char * end = strchr(ptr, ' ');
    const char * tst = strchr(ptr, ',');
    if (end == NULL || (tst && tst < end)) end = tst;
    if (end == NULL) end = strchr(ptr, '\0');
    int i = 0;
    while (end && i < num_breakpoints) {
      const ptrdiff_t len = end - ptr;
      breakpoints[i] = new char[len + 1];
      (void)memcpy(breakpoints[i], ptr, len);
      breakpoints[i][len] = '\0';
      i++;
      ptr = end+1;
      if (ptr < (cpy + envstrlen)) {
        end = strchr(ptr, ' ');
        tst = strchr(ptr, ',');
        if (end == NULL || (tst && tst < end)) end = tst;
        if (end == NULL) end = strchr(ptr, 0);
      }
    }
    num_breakpoints = i; // just in case parsing failed for some reason
    delete[] cpy;
  }
#endif // COIN_DEBUG
}

void
SoDebugError::callbackForwarder(const struct cc_debugerror * error,
                            void * COIN_UNUSED_ARG(data)
                            )
{
  SoDebugError wrappederr;

  switch (cc_debugerror_get_severity(error)) {
  case CC_DEBUGERROR_ERROR:
    wrappederr.severity = SoDebugError::ERROR;
    break;
  case CC_DEBUGERROR_WARNING:
    wrappederr.severity = SoDebugError::WARNING;
    break;
  case CC_DEBUGERROR_INFO:
    wrappederr.severity = SoDebugError::INFO;
    break;
  default:
    assert(FALSE);
    break;
  }

  const cc_string * dbgstr =
    cc_error_get_debug_string(reinterpret_cast<const cc_error *>(error));
  const char * dbgstrc = cc_string_get_text(dbgstr);
  wrappederr.setDebugString(dbgstrc);

  assert(SoDebugError::callback != NULL);
  SoDebugError::callback(&wrappederr, SoDebugError::callbackData);
}

// Documented for parent class.
void
SoDebugError::setHandlerCallback(SoErrorCB * const function, void * const data)
{
  if (SoDebugError::callback == SoError::defaultHandlerCB) {
    // The user is overriding the default handler, so set up a
    // "converter" callback function that makes an SoDebugError out of
    // an cc_debugerror and forwards control to the callback function
    // given as an argument to setHandlerCallback().
    cc_debugerror_set_handler_callback(
       reinterpret_cast<cc_debugerror_cb *>(SoDebugError::callbackForwarder), NULL);
  }

  SoDebugError::callback = function;
  SoDebugError::callbackData = data;
}

// Documented for parent class.
SoErrorCB *
SoDebugError::getHandlerCallback(void)
{
  return SoDebugError::callback;
}

// Documented for parent class.
void *
SoDebugError::getHandlerData(void)
{
  return SoDebugError::callbackData;
}

// Documented for parent class.
SoType
SoDebugError::getClassTypeId(void)
{
  return SoDebugError::classTypeId;
}

// Documented for parent class.
SoType
SoDebugError::getTypeId(void) const
{
  return SoDebugError::classTypeId;
}

/*!
  Returns severity level of current message.

  You can use this to filter out debug messages and only let warnings
  and/or errors pass through to the end-user if you have set your own
  handler callback.

  \sa setHandlerCallback()
*/
SoDebugError::Severity
SoDebugError::getSeverity(void) const
{
  return this->severity;
}

#if COIN_DEBUG

static inline void
check_breakpoints(const char * source)
{
  for (int i = 0; i < num_breakpoints; i++) {
    if (strcmp(breakpoints[i], source) == 0) {
      COIN_DEBUGGER_BREAK("Coin debug break");
    }
  }
}

#else // COIN_DEBUG

static inline void
check_breakpoints(const char *)
{
  // should be empty
}

#endif // ! COIN_DEBUG


void
SoDebugError::commonPostHandling(Severity severity, const char * type,
                                 const char * source, const SbString & s)
{
  SoDebugError error;
  error.severity = severity;
  error.setDebugString("Coin ");
  error.appendToDebugString(type);
  error.appendToDebugString(" in ");
  error.appendToDebugString(source);
  error.appendToDebugString("(): ");
  error.appendToDebugString(s.getString());

  if (SoDebugError::callback != SoError::defaultHandlerCB) {
    SoDebugError::callback(&error, SoDebugError::callbackData);
  }
  else {
    error.handleError();
  }

  check_breakpoints(source);
}

#define SODEBUGERROR_POST(SEVERITY, TYPE) \
  va_list args; \
  va_start(args, format); \
  SbString s; \
  s.vsprintf(format, args); \
  va_end(args); \
 \
  SoDebugError::commonPostHandling(SEVERITY, TYPE, source, s)


/*!
  This method posts a message with severity level "ERROR". \a source should
  be the name of the calling function.
*/
void
SoDebugError::post(const char * const source, const char * const format, ...)
{
  SODEBUGERROR_POST(ERROR, "error");
}

/*!
  This method posts a message with severity level "WARNING". \a source should
  be the name of the calling function.
*/
void
SoDebugError::postWarning(const char * const source, const char * const format, ...)
{
  SODEBUGERROR_POST(WARNING, "warning");
}

/*!
  This method posts a message with severity level "INFO". \a source should
  be the name of the calling function.
*/
void
SoDebugError::postInfo(const char * const source, const char * const format, ...)
{
  SODEBUGERROR_POST(INFO, "info");
}

// Documented for parent class.
SoErrorCB *
SoDebugError::getHandler(void * & data) const
{
  // FIXME: why the heck isn't this a static method? Looks to me like
  // an API bug. 20020526 mortene.

  data = SoDebugError::callbackData;
  return SoDebugError::callback;
}

#undef SODEBUGERROR_POST
#undef COIN_DEBUGGER_BREAK
