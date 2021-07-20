#ifndef COIN_SIMAGEWRAPPER_H
#define COIN_SIMAGEWRAPPER_H

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

#ifndef COIN_INTERNAL
#error this is a private header file
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#endif /* HAVE_WINDOWS_H */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* Typedefinitions of function signatures for simage calls we
     use. We need these for casting from the void-pointer return of
     dlsym().

     Note specifically for MSWindows that we do _not_ use the APIENTRY
     keyword in the function typedefs, as that would set them up with
     the __stdcall calling convention -- and simage functions are
     built with the __cdecl calling convention.
  */

  typedef void (*simage_version_t)(int *, int *, int *);
  typedef int (*simage_check_supported_t)(const char * filename);
  typedef unsigned char * (*simage_read_image_t)(const char * filename,
                                                          int * w, int * h,
                                                          int * numcomponents);
  typedef int (*simage_check_save_supported_t)(const char * filenameextension);
  typedef int (*simage_save_image_t)(const char * filename,
                                              const unsigned char * bytes,
                                              int w, int h, int numcomponents,
                                              const char * filenameextension);
  typedef const char * (*simage_get_last_error_t)(void);
  typedef unsigned char * (*simage_resize_t)(unsigned char * imagedata,
                                                      int width, int height,
                                                      int numcomponents,
                                                      int newwidth, int newheight);
  typedef unsigned char * (*simage_resize3d_t)(unsigned char * imagedata,
                                                        int width, int height,
                                                        int numcomponents,
                                                        int layers,
                                                        int newwidth, 
                                                        int newheight,
                                                        int newlayers);
  typedef void (*simage_free_image_t)(unsigned char * imagedata);
  typedef int (*simage_next_power_of_two_t)(int val);

  typedef int (*simage_get_num_savers_t)(void);
  typedef void * (*simage_get_saver_handle_t)(int idx);
  typedef const char * (*simage_get_saver_extensions_t)(void * handle);
  typedef const char * (*simage_get_saver_fullname_t)(void * handle);
  typedef const char * (*simage_get_saver_description_t)(void * handle);

  /* This define is set up in the simage_wrapper.cpp file, according to
     whether or not we link static at compile-time or dynamic at
     run-time to the simage library. */
#if !defined(SIMAGEWRAPPER_ASSUME_SIMAGE)
  /* This wrapping of the enum and typedefs is necessary to avoid
     multiple definitions (they are copy'n'pasted from the simage.h
     header file). */
  enum {
    S_INTEGER_PARAM_TYPE,
    S_BOOL_PARAM_TYPE = S_INTEGER_PARAM_TYPE,
    S_FLOAT_PARAM_TYPE,
    S_DOUBLE_PARAM_TYPE,
    S_STRING_PARAM_TYPE,
    S_POINTER_PARAM_TYPE,
    S_FUNCTION_PARAM_TYPE
  };

  typedef struct simage_parameters_s s_params;
#endif

  /* Streams implementation was added for simage v1.4. */
#if !defined(SIMAGEWRAPPER_ASSUME_SIMAGE) || !defined(SIMAGE_VERSION_1_4)
  typedef struct simage_stream_s s_stream;
#endif

  typedef s_params * (*s_params_create_t)(void);
  typedef void (*s_params_destroy_t)(s_params * params);
  typedef void (*s_params_set_t)(s_params * params, ...);
  typedef int (*s_params_get_t)(s_params * params, ...);

  typedef s_stream * (*s_stream_open_t)(const char * filename, 
                                          s_params * params /* | NULL */);
  typedef void * (*s_stream_get_buffer_t)(s_stream * stream, 
                                           void * prealloc /* | NULL */,
                                           int *size /* | NULL */,
                                           s_params * params /* | NULL */);
  typedef void (*s_stream_close_t)(s_stream * stream);
  typedef void (*s_stream_destroy_t)(s_stream * stream);
  typedef s_params * (*s_stream_params_t)(s_stream * stream);

  typedef struct {
    /* Is the simage library at all available? */
    int available;

    /* simage versioning. */
    struct {
      int major, minor, micro;
    } version;
    int (*versionMatchesAtLeast)(int major,
                                 int minor,
                                 int micro);

    simage_version_t simage_version;
    simage_check_supported_t simage_check_supported;
    simage_read_image_t simage_read_image;
    simage_check_save_supported_t simage_check_save_supported;
    simage_save_image_t simage_save_image;
    simage_get_last_error_t simage_get_last_error;
    simage_resize_t simage_resize;
    simage_free_image_t simage_free_image;
    simage_next_power_of_two_t simage_next_power_of_two;
    simage_get_num_savers_t simage_get_num_savers;
    simage_get_saver_handle_t simage_get_saver_handle;
    simage_get_saver_extensions_t simage_get_saver_extensions;
    simage_get_saver_fullname_t simage_get_saver_fullname;
    simage_get_saver_description_t simage_get_saver_description;

    s_params_create_t s_params_create;
    s_params_destroy_t s_params_destroy;
    s_params_set_t s_params_set;
    s_params_get_t s_params_get;

    simage_resize3d_t simage_resize3d;

    s_stream_open_t s_stream_open;
    s_stream_get_buffer_t s_stream_get_buffer;
    s_stream_close_t s_stream_close;
    s_stream_destroy_t s_stream_destroy;
    s_stream_params_t s_stream_params;

  } simage_wrapper_t;

  const simage_wrapper_t * simage_wrapper(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* COIN_SIMAGEWRAPPER_H */
