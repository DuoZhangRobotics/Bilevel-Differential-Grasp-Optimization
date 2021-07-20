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

// The purpose of the class in this file is to try to aid in detecting
// when multiple Coin DLLs are loaded into the same run-time process
// image under Windows.
//
// This is useful because linking with multiple instances of Coin, for
// instance both a release build (e.g. indirectly through SoWin or
// SoQt) and a debug build (e.g. from the application executable
// linking) at the same time, will cause hard to find errors.
//
// A typical symptom of this problem would be an assert-error from the
// constructor of an SoBase-derived class because of missing library
// initialization for one of the Coin-libraries in memory.

// *************************************************************************

#include "CoinStaticObjectInDLL.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#ifdef HAVE_WINDOWS_H
#include <windows.h>
#endif // HAVE_WINDOWS_H

#include <Inventor/SbTime.h>

// *************************************************************************

void * CoinStaticObjectInDLL::mutexhandle = NULL;
CoinStaticObjectInDLL * CoinStaticObjectInDLL::singleton = NULL;

// *************************************************************************

#ifndef HAVE_WIN32_API

CoinStaticObjectInDLL::CoinStaticObjectInDLL(void) { assert(FALSE); }
CoinStaticObjectInDLL::~CoinStaticObjectInDLL() { assert(FALSE); }
void CoinStaticObjectInDLL::init(void) { /* will be called, so don't assert */ }
// dummy, avoid linker errors
SbString CoinStaticObjectInDLL::mutexName(void) { return SbString(""); }
SbBool CoinStaticObjectInDLL::activateMutex(void) {  assert(FALSE); return TRUE; }
void CoinStaticObjectInDLL::deactivateMutex(void) { assert(FALSE); }

#else // HAVE_WIN32_API

// FIXME: this check should be possible to turn off with an envvar, as
// it is sometimes not an indication of a problem to have multiple
// Coin instances in the same process -- for instance if Coin is part
// of a browser plug-in (or several plug-ins).  20041021 mortene.
//
// UPDATE 20041108 mortene: an idea for improvement, from kyrah: if
// each built Coin library had a unique key / ID, we could incorporate
// that into the check. That would be an improvement, as most errors
// of "multiple-instance-inclusion" is very likely that a debug build
// and a release build of Coin is loaded at the same time, typically
// because e.g. SoWin or SoQt was linked against coin2d.dll, while the
// app links against coin2.dll.

// *************************************************************************

// Only one of these objects should be allocated for the process. If
// two or more of them are it means that multiple instances of the
// Coin library is loaded for the same process image -- and we'll
// throw up the error message box.

#if 1 // enable / disable the check -- provides a quick workaround in
      // case someone runs into trouble with it.
static CoinStaticObjectInDLL dllobject;
#endif // enable / disable

// *************************************************************************

// Must use the C library getenv() calls below, not Coin's
// coin_getenv(), since this code is executed before SoDB::init().

static SbBool
debug(void)
{
  static int dbg = -1;
  if (dbg == -1) {
    const char * env = getenv("COIN_DEBUG_STATICOBJECT");
    dbg = (env && (atoi(env) > 0)) ? 1 : 0;
  }
  return dbg;
}

static SbBool
runtime_disabled(void)
{
  static int val = -1;
  if (val == -1) {
    const char * env = getenv("COIN_INSTANCEMUTEX_DISABLE");
    val = (env && (atoi(env) > 0)) ? 1 : 0;
  }
  return val;
}

// *************************************************************************

CoinStaticObjectInDLL::CoinStaticObjectInDLL(void)
{
  if (runtime_disabled()) { return; }

  // Can't use SoDebugError, as nothing has been initialized yet.
  if (debug()) {
    printf("%p %f CoinStaticObjectInDLL constructor\n",
           this, SbTime::getTimeOfDay().getValue());
  }

  assert(CoinStaticObjectInDLL::singleton == NULL);
  CoinStaticObjectInDLL::singleton = this;

  if (!CoinStaticObjectInDLL::activateMutex()) {
    MessageBox(NULL,
               "Detected two instances of the Coin library in the same\n"
               "process image!!\n\n"

               "Application can not continue without errors, and\n"
               "will exit when you quit this dialog box.\n\n"

               "This is an indication of a serious problem with the\n"
               "settings in your project configuration.\n\n"

               "One likely cause of this error is that a different\n"
               "configuration of the Coin library was linked with\n"
               "the GUI-binding library (e.g. SoWin or SoQt) than\n"
               "for the application's linker settings. Try for instance\n"
               "to check if there is a mismatch between debug and release\n"
               "libraries.\n\n"

               "The depends.exe program that comes with VisualC++ is a\n"
               "good tool for tracking things like this down. Make sure\n"
               "you inspect the complete path of each loaded dll.\n\n"

               "If you are completely lost as how to find and fix\n"
               "this problem on your own, try the\n"
               "<coin-discuss@coin3d.org> mailing list (or the support\n"
               "address <coin-support@coin3d.org> if you hold a Coin\n"
               "Professional Edition License).\n",


               "Fatal error!", MB_OK | MB_ICONERROR | MB_TASKMODAL);
    exit(1);
  }
}

CoinStaticObjectInDLL::~CoinStaticObjectInDLL()
{
  if (runtime_disabled()) { return; }

  if (debug()) {
    printf("%p %f CoinStaticObjectInDLL destructor\n",
           this, SbTime::getTimeOfDay().getValue());
  }

  // Wrap in if-check, since the mutex would usually be deallocated
  // from CoinStaticObjectInDLL::init(). This doesn't always happen,
  // though, as a Coin DLL can be loaded and then unloaded without
  // SoDB::init() being called (as from
  // src/glue/dl.c:cc_dl_coin_handle()).

  if (CoinStaticObjectInDLL::mutexhandle) {
    CoinStaticObjectInDLL::deactivateMutex();
  }

  assert(CoinStaticObjectInDLL::singleton);
  CoinStaticObjectInDLL::singleton = NULL;
}

// Called from SoDB::init().
void
CoinStaticObjectInDLL::init(void)
{
  if (runtime_disabled()) { return; }

  if (debug()) {
    printf("%p %f CoinStaticObjectInDLL::init()\n",
           CoinStaticObjectInDLL::singleton,
           SbTime::getTimeOfDay().getValue());
  }

  // Program control has been handed over to actual program code at
  // this point, so we take away the mutex again -- the check should
  // have hit by now if two Coin DLLs were loaded.

  // (Doing this removes the potential for a rather obscure problem
  // we've seen: we sometimes "manually" load a Coin DLL as an
  // indirect means to get at the values of OpenGL extension function
  // symbols (in src/glue/dl.c:cc_dl_coin_handle()). That would
  // sometimes trigger the instantiation of two CoinStaticObjectInDLL
  // objects, which made the above check hit.)

  if (CoinStaticObjectInDLL::singleton) { // in case the static object
                                          // has been disabled in the
                                          // code
    CoinStaticObjectInDLL::deactivateMutex();
  }
}

// Returns TRUE if mutex was created ok, or FALSE if the mutex was
// already created.
SbBool
CoinStaticObjectInDLL::activateMutex(void)
{
  if (debug()) {
    printf("%p %f CoinStaticObjectInDLL::activateMutex(), mutexname=='%s'\n",
           CoinStaticObjectInDLL::singleton,
           SbTime::getTimeOfDay().getValue(),
           CoinStaticObjectInDLL::mutexName().getString());
  }

  assert(CoinStaticObjectInDLL::mutexhandle == NULL);

  SetLastError(0); // so we don't react to an old error for the check below

  CoinStaticObjectInDLL::mutexhandle = (HANDLE)
    CreateMutex(NULL, TRUE, CoinStaticObjectInDLL::mutexName().getString());
  // (The mutex is automatically destructed by the operating system
  // when the process exits.)

  return (GetLastError() == ERROR_ALREADY_EXISTS) ? FALSE : TRUE;
}

void
CoinStaticObjectInDLL::deactivateMutex(void)
{
  if (debug()) {
    printf("%p %f CoinStaticObjectInDLL::deactivateMutex(), handle==%p\n",
           CoinStaticObjectInDLL::singleton,
           SbTime::getTimeOfDay().getValue(),
           CoinStaticObjectInDLL::mutexhandle);
  }

  // it's only necessary to close the mutex handle the first time 
  // CoinStaticObjectInDLL::init() is called. In subsequent calls the handle
  // will be NULL.
  if (CoinStaticObjectInDLL::mutexhandle != NULL) {
    const BOOL ok = CloseHandle((HANDLE)CoinStaticObjectInDLL::mutexhandle);
    if (!ok) { // just in case
      if (debug()) {
        printf("%p %f CoinStaticObjectInDLL::deactivateMutex(), "
               "ERROR: CloseHandle(%p) failed!\n",
               CoinStaticObjectInDLL::singleton,
               SbTime::getTimeOfDay().getValue(),
               CoinStaticObjectInDLL::mutexhandle);
      }

      MessageBox(NULL,
                 "CloseHandle() in CoinStaticObjectInDLL::deactivateMutex()\n"
                 "failed! Please report to <coin-support@coin3d.org>.\n",
                 "Warning!", MB_OK | MB_ICONERROR | MB_TASKMODAL);
    }
    CoinStaticObjectInDLL::mutexhandle = (HANDLE)NULL;
  }
}

SbString
CoinStaticObjectInDLL::mutexName(void)
{
  SbString s;
  s.sprintf("COIN_LIBRARY_PROCESS_%d", GetCurrentProcessId());
  return s;
}

#endif // HAVE_WIN32_API

// *************************************************************************
