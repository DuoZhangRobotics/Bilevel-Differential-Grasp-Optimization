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

#include "misc/AudioTools.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <Inventor/SbString.h>
#include <Inventor/C/tidbits.h>

#include "tidbitsp.h"
#include "glue/openal_wrapper.h"

// *************************************************************************

// There are some differences in error defines between different
// versions of OpenAL, and we try to fix them up below.
//
// (It even looks like some of these defines have been set up
// exclusively for MSWindows usage, while others are for other
// platforms. Really great interface design there from the OpenAL
// developers...)
#ifndef AL_ILLEGAL_ENUM
#define AL_ILLEGAL_ENUM AL_INVALID_ENUM
#endif // !AL_ILLEGAL_ENUM
#ifndef AL_ILLEGAL_COMMAND
#define AL_ILLEGAL_COMMAND AL_INVALID_OPERATION
#endif // !AL_ILLEGAL_COMMAND

// *************************************************************************

const char *
coin_get_openal_error(int errcode)
{
  switch (errcode) {
  case AL_INVALID_NAME:
    return "AL_INVALID_NAME - Illegal name passed as an argument to an AL call";

  case AL_INVALID_ENUM:
    return "AL_INVALID_ENUM - Illegal enum passed as an argument to an AL call";

  case AL_INVALID_VALUE:
    return "AL_INVALID_VALUE - Illegal value passed as an argument to an AL call";

  case AL_INVALID_OPERATION:
    return "AL_INVALID_OPERATION - A function was called at an inappropriate time or in an inappropriate way, causing an illegal state. This can be an incompatible ALenum, object ID, and/or function";

  case AL_OUT_OF_MEMORY:
    return "AL_OUT_OF_MEMORY - A function could not be completed, because there is not enough memory available.";

  default:
    return "UNDEFINED ERROR";
  }
}

/* Return value of COIN_DEBUG_AUDIO environment variable. */
int
coin_debug_audio(void)
{
  static int d = -1;
  if (d == -1) {
    const char * val = coin_getenv("COIN_DEBUG_AUDIO");
    d = val ? atoi(val) : 0;
  }
  return (d > 0) ? 1 : 0;
}

// *************************************************************************

// The following data and functions are used to avoid
// SoAudioRenderAction being applied to the scene graph unless there
// has been at least one sound-emitting node instantiated.

static SbBool should_traverse = FALSE;

static void audio_cleanup(void) 
{
  should_traverse = FALSE;
}

// Called from the constructor of SoVRMLAudioClip and SoVRMLSound.
// After this has been called once, SoSceneManager instances will
// start applying an SoAudioRenderAction to its scene graph before
// rendering, at each frame.
//
// This is far from perfect design, but it is better than how it used
// to be (SoAudioRenderAction _always_ applied, even when no sound
// nodes had even been instantiated by the application -- which I
// presume is the common case, actually).
//
// 20050627 mortene.
void
coin_sound_enable_traverse(void)
{
  should_traverse = TRUE;
  coin_atexit((coin_atexit_f*) audio_cleanup, CC_ATEXIT_NORMAL);
}

// Called from SoSceneManager::render() to decide whether or not to
// invoke its SoAudioRenderAction on its scene graph.
SbBool
coin_sound_should_traverse(void)
{
  return should_traverse;
}

// *************************************************************************
