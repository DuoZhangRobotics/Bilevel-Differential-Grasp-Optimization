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

#include "SoGLSLShaderParameter.h"
#include "SoGLSLShaderObject.h"

#include <Inventor/errors/SoDebugError.h>

// *************************************************************************

SoGLSLShaderParameter::SoGLSLShaderParameter(void)
{
  this->location  = -1;
  this->cacheType = GL_FLOAT;
  this->cacheName = "";
  this->cacheSize =  0;
  this->isActive = TRUE;
  this->programid = 0;
}

SoGLSLShaderParameter::~SoGLSLShaderParameter()
{
}

SoShader::Type
SoGLSLShaderParameter::shaderType(void) const
{ 
  return SoShader::GLSL_SHADER;
}

void
SoGLSLShaderParameter::set1f(const SoGLShaderObject * shader,
                             const float value, const char *name, const int)
{
  if (this->isValid(shader, name, GL_FLOAT))
    shader->GLContext()->glUniform1fARB(this->location, value);
}

void
SoGLSLShaderParameter::set2f(const SoGLShaderObject * shader,
                             const float * value, const char *name, const int)
{
  if (this->isValid(shader, name, GL_FLOAT_VEC2_ARB))
    shader->GLContext()->glUniform2fARB(this->location, value[0], value[1]);
}

void
SoGLSLShaderParameter::set3f(const SoGLShaderObject * shader,
                             const float * v, const char *name, const int)
{
  if (this->isValid(shader, name, GL_FLOAT_VEC3_ARB))
    shader->GLContext()->glUniform3fARB(this->location, v[0], v[1], v[2]);
}

void
SoGLSLShaderParameter::set4f(const SoGLShaderObject * shader,
                             const float * v, const char *name, const int)
{
  if (this->isValid(shader, name, GL_FLOAT_VEC4_ARB))
    shader->GLContext()->glUniform4fARB(this->location, v[0], v[1], v[2], v[3]);
}


void
SoGLSLShaderParameter::set1fv(const SoGLShaderObject * shader, const int num,
                              const float *value, const char * name, const int)
{
  int cnt = num;
  if (this->isValid(shader, name, GL_FLOAT, &cnt))
    shader->GLContext()->glUniform1fvARB(this->location, cnt, value);
}

void
SoGLSLShaderParameter::set2fv(const SoGLShaderObject * shader, const int num,
                              const float* value, const char* name, const int)
{
  int cnt = num;
  if (this->isValid(shader, name, GL_FLOAT_VEC2_ARB, &cnt))
    shader->GLContext()->glUniform2fvARB(this->location, cnt, value);
}

void
SoGLSLShaderParameter::set3fv(const SoGLShaderObject * shader, const int num,
                              const float* value, const char * name, const int)
{
  int cnt = num;
  if (this->isValid(shader, name, GL_FLOAT_VEC3_ARB, &cnt))
    shader->GLContext()->glUniform3fvARB(this->location, cnt, value);
}

void
SoGLSLShaderParameter::set4fv(const SoGLShaderObject * shader, const int num,
                              const float* value, const char * name, const int)
{
  int cnt = num;
  if (this->isValid(shader, name, GL_FLOAT_VEC4_ARB, &cnt))
    shader->GLContext()->glUniform4fvARB(this->location, cnt, value);
}

void
SoGLSLShaderParameter::setMatrix(const SoGLShaderObject *shader,
                                 const float * value, const char * name,
                                 const int)
{
  if (this->isValid(shader, name, GL_FLOAT_MAT4_ARB))
    shader->GLContext()->glUniformMatrix4fvARB(this->location,1,FALSE,value);
}

  
void
SoGLSLShaderParameter::setMatrixArray(const SoGLShaderObject *shader,
                                      const int num, const float *value,
                                      const char *name, const int)
{
  int cnt = num;
  if (this->isValid(shader, name, GL_FLOAT_MAT4_ARB, &cnt))
    shader->GLContext()->glUniformMatrix4fvARB(this->location,cnt,FALSE,value);
}


void
SoGLSLShaderParameter::set1i(const SoGLShaderObject * shader,
                             const int32_t value, const char * name, const int)
{
  if (this->isValid(shader, name, GL_INT))
    shader->GLContext()->glUniform1iARB(this->location, value);
}

void
SoGLSLShaderParameter::set2i(const SoGLShaderObject * shader,
                             const int32_t * value, const char * name,
                             const int)
{
  if (this->isValid(shader, name, GL_INT_VEC2_ARB))
    shader->GLContext()->glUniform2iARB(this->location, value[0], value[1]);
}

void
SoGLSLShaderParameter::set3i(const SoGLShaderObject * shader,
                             const int32_t * v, const char * name,
                             const int)
{
  if (this->isValid(shader, name, GL_INT_VEC3_ARB))
    shader->GLContext()->glUniform3iARB(this->location, v[0], v[1], v[2]);
}

void
SoGLSLShaderParameter::set4i(const SoGLShaderObject * shader,
                             const int32_t * v, const char * name,
                             const int)
{
  if (this->isValid(shader, name, GL_INT_VEC4_ARB))
    shader->GLContext()->glUniform4iARB(this->location, v[0], v[1], v[2], v[3]);
}

void
SoGLSLShaderParameter::set1iv(const SoGLShaderObject * shader,
                              const int num,
                              const int32_t * value, const char * name,
                              const int)
{
  if (this->isValid(shader, name, GL_INT))
    shader->GLContext()->glUniform1ivARB(this->location, num, (const GLint*) value);
}

void
SoGLSLShaderParameter::set2iv(const SoGLShaderObject * shader,
                              const int num,
                              const int32_t * value, const char * name,
                              const int)
{
  if (this->isValid(shader, name, GL_INT_VEC2_ARB))
    shader->GLContext()->glUniform2ivARB(this->location, num, (const GLint*)value);
}

void
SoGLSLShaderParameter::set3iv(const SoGLShaderObject * shader,
                              const int num,
                              const int32_t * v, const char * name,
                              const int)
{
  if (this->isValid(shader, name, GL_INT_VEC3_ARB))
    shader->GLContext()->glUniform3ivARB(this->location, num, (const GLint*)v);
}

void
SoGLSLShaderParameter::set4iv(const SoGLShaderObject * shader,
                              const int num,
                              const int32_t * v, const char * name,
                              const int)
{
  if (this->isValid(shader, name, GL_INT_VEC4_ARB))
    shader->GLContext()->glUniform4ivARB(this->location, num, (const GLint*)v);
}

#include <stdio.h>
SbBool
SoGLSLShaderParameter::isValid(const SoGLShaderObject * shader,
                               const char * name, GLenum type,
                               int * num)
{
  assert(shader);
  assert(shader->shaderType() == SoShader::GLSL_SHADER);
  
  COIN_GLhandle pHandle = ((SoGLSLShaderObject*)shader)->programHandle;
  int32_t pId = ((SoGLSLShaderObject*)shader)->programid;
  
  // return TRUE if uniform isn't active. We warned the user about
  // this when we found it to be inactive.
  if ((pId == this->programid) && (this->location > -1) && !this->isActive) return TRUE;
  
  if ((pId == this->programid) && (this->location > -1) && 
      (this->cacheName == name) && (this->cacheType == type)) {
    if (num) { // assume: ARRAY
      if (this->cacheSize < *num) {
        // FIXME: better error handling - 20050128 martin
        SoDebugError::postWarning("SoGLSLShaderParameter::isValid",
                                  "parameter %s[%d] < input[%d]!",
                                  this->cacheName.getString(),
                                  this->cacheSize, *num);
        *num = this->cacheSize;
      }
      return (*num > 0);
    }
    return TRUE;
  }
  
  const cc_glglue * g = shader->GLContext();
  
  this->cacheSize = 0;  
  this->location = g->glGetUniformLocationARB(pHandle,
                                              (const COIN_GLchar *)name);
  this->programid = pId;
  
  if (this->location == -1)  {
#if COIN_DEBUG
    SoDebugError::postWarning("SoGLSLShaderParameter::isValid",
                              "parameter '%s' not found in program.",
                              name);
#endif // COIN_DEBUG
    return FALSE;
  }
  GLint activeUniforms = 0;
  g->glGetObjectParameterivARB(pHandle, GL_OBJECT_ACTIVE_UNIFORMS_ARB, &activeUniforms);

  GLint i;
  GLint tmpSize = 0;
  GLenum tmpType;
  GLsizei length;
  COIN_GLchar myName[256];
  
  this->cacheName = name;
  this->isActive = FALSE; // set uniform to inactive while searching
  
  // this will only happen once after the variable has been added so
  // it's not a performance issue that we have to search for it here.
  for (i = 0; i < activeUniforms; i++) {
    g->glGetActiveUniformARB(pHandle, i, 128, &length, &tmpSize, 
                             &tmpType, myName);
    if (this->cacheName == myName) {
      this->cacheSize = tmpSize;
      this->cacheType = tmpType;
      this->isActive = TRUE;
      break;
    }
  }
  if (!this->isActive) {
    // not critical, but warn user so they can remove the unused parameter
#if COIN_DEBUG
    SoDebugError::postWarning("SoGLSLShaderParameter::isValid",
                              "parameter '%s' not active.",
                              this->cacheName.getString());
#endif // COIN_DEBUG
    // return here since cacheSize and cacheType will not be properly initialized
    return TRUE;
  }

  if (type == GL_INT) {
    switch (this->cacheType) {
    case GL_INT:
    case GL_SAMPLER_1D_ARB:
    case GL_SAMPLER_2D_ARB:
    case GL_SAMPLER_3D_ARB:
    case GL_SAMPLER_CUBE_ARB:
    case GL_SAMPLER_1D_SHADOW_ARB:
    case GL_SAMPLER_2D_SHADOW_ARB:
    case GL_SAMPLER_2D_RECT_ARB:
    case GL_SAMPLER_2D_RECT_SHADOW_ARB: 
      break;
    default: 
      return FALSE;
    }
  }
  else if (this->cacheType != type)
    return FALSE;

  if (num) { // assume: ARRAY
    if (this->cacheSize < *num) {
      // FIXME: better error handling - 20050128 martin
      SoDebugError::postWarning("SoGLSLShaderParameter::isValid",
                                "parameter %s[%d] < input[%d]!",
                                this->cacheName.getString(),
                                this->cacheSize, *num);
      *num = this->cacheSize;
    }
    return (*num > 0);
  }
  return TRUE;
}
