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
  \class SoShaderParameter SoShaderParameter.h Inventor/nodes/SoShaderParameter.h

  The SoShaderParameter class is the base class for all shader parameter classes.

  In addition to the \a name and \a identifier field, all subclasses have a
  \a value field which is used for specifying the parameter value.
*/

/*!
  SoSFString SoShaderParameter::name

  The shader parameter name. Used for Cg and GLSL programs.
*/

/*!
  SoSFInt32 SoShaderParameter::identifier

  The shader parameter identifier. Used for ARB shader programs.
*/

#include <Inventor/nodes/SoShaderParameter.h>

#include <assert.h>

#include "nodes/SoSubNodeP.h"
#include "misc/SbHash.h"
#include "glue/cg.h"
#include "shaders/SoGLShaderObject.h"
#include "shaders/SoGLShaderParameter.h"
#include "shaders/SoGLCgShaderParameter.h"

/* **************************************************************************
 * *** SoShaderParameter ***
 * **************************************************************************/

SO_NODE_ABSTRACT_SOURCE(SoShaderParameter);

// doc from parent
void
SoShaderParameter::initClass(void)
{
  SO_NODE_INTERNAL_INIT_ABSTRACT_CLASS(SoShaderParameter,
                                       SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

/*!
  Constructor.
*/
SoShaderParameter::SoShaderParameter(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter);

  SO_NODE_ADD_FIELD(name, (""));
  SO_NODE_ADD_FIELD(identifier, (0));
}

/*!
  Destructor.
*/
SoShaderParameter::~SoShaderParameter()
{
}

#define PRIVATE(obj) obj->pimpl

/* **************************************************************************
 * *** SoUniformShaderParameter ***
 * **************************************************************************/

class SoUniformShaderParameterP {
public:
  SoUniformShaderParameterP() { }
  ~SoUniformShaderParameterP() {
    SbList <uint32_t> keylist;
    this->glparams.makeKeyList(keylist);
    for (int i = 0; i < keylist.getLength(); i++) {
      SoGLShaderParameter * param;
      (void) this->glparams.get(keylist[i], param);
      deleteGLParameter(param);
    }
  }
  static void deleteGLParameter(SoGLShaderParameter * param) {
    // FIXME: schedule for delete, pederb 2005-11-30
  }
  // FIXME: add a cache context destruction callback, pederb 2005-11-30
  SbHash <SoGLShaderParameter *, uint32_t> glparams;
};

SO_NODE_ABSTRACT_SOURCE(SoUniformShaderParameter);

void
SoUniformShaderParameter::initClass(void)
{
  SO_NODE_INTERNAL_INIT_ABSTRACT_CLASS(SoUniformShaderParameter,
                                       SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoUniformShaderParameter::SoUniformShaderParameter(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoUniformShaderParameter);

  PRIVATE(this) = new SoUniformShaderParameterP;
}

SoUniformShaderParameter::~SoUniformShaderParameter()
{
  delete PRIVATE(this);
}

void
SoUniformShaderParameter::ensureParameter(SoGLShaderObject * shader)
{
  assert(shader);
  const uint32_t context = shader->getCacheContext();
  SoGLShaderParameter * param;
  if (!PRIVATE(this)->glparams.get(context, param)) {
    param = shader->getNewParameter();
    (void) PRIVATE(this)->glparams.put(context, param);
  }
  if (param->shaderType() != shader->shaderType()) {
    SoUniformShaderParameterP::deleteGLParameter(param);
    param = shader->getNewParameter();
    (void) PRIVATE(this)->glparams.put(context, param);
  }
}

SoGLShaderParameter *
SoUniformShaderParameter::getGLShaderParameter(const uint32_t cachecontext)
{
  SoGLShaderParameter * glparam;
  if (PRIVATE(this)->glparams.get(cachecontext, glparam)) return glparam;
  return NULL;
}


/* **************************************************************************
 * *** SoShaderParameter1f ***
 * **************************************************************************/


SO_NODE_SOURCE(SoShaderParameter1f);

void SoShaderParameter1f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter1f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameter1f::SoShaderParameter1f(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter1f);
  SO_NODE_ADD_FIELD(value, (0));
}

SoShaderParameter1f::~SoShaderParameter1f()
{
}

void
SoShaderParameter1f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set1f(shader,
                                                               this->value.getValue(),
                                                               this->name.getValue().getString(),
                                                               this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter2f ***
 * **************************************************************************/


SO_NODE_SOURCE(SoShaderParameter2f);

void SoShaderParameter2f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter2f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameter2f::SoShaderParameter2f(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter2f);
  SO_NODE_ADD_FIELD(value, (0,0));
}

SoShaderParameter2f::~SoShaderParameter2f()
{
}

void SoShaderParameter2f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set2f(shader,
                         this->value.getValue().getValue(),
                         this->name.getValue().getString(),
                         this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter3f ***
 * **************************************************************************/


SO_NODE_SOURCE(SoShaderParameter3f);

void SoShaderParameter3f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter3f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameter3f::SoShaderParameter3f(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter3f);
  SO_NODE_ADD_FIELD(value, (0,0,0));
}

SoShaderParameter3f::~SoShaderParameter3f()
{
}

void SoShaderParameter3f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set3f(shader,
                         this->value.getValue().getValue(),
                         this->name.getValue().getString(),
                         this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter4f ***
 * **************************************************************************/


SO_NODE_SOURCE(SoShaderParameter4f);

void SoShaderParameter4f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter4f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameter4f::SoShaderParameter4f(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter4f);
  SO_NODE_ADD_FIELD(value, (0,0,0,0));
}

SoShaderParameter4f::~SoShaderParameter4f()
{
}

void SoShaderParameter4f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set4f(shader,
                         this->value.getValue().getValue(),
                         this->name.getValue().getString(),
                         this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderStateMatrixParameter ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderStateMatrixParameter);

void SoShaderStateMatrixParameter::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderStateMatrixParameter,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderStateMatrixParameter::SoShaderStateMatrixParameter(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderStateMatrixParameter);

  SO_NODE_DEFINE_ENUM_VALUE(MatrixType, MODELVIEW);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixType, PROJECTION);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixType, TEXTURE);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixType, MODELVIEW_PROJECTION);

  SO_NODE_ADD_FIELD(matrixType, (MODELVIEW));
  SO_NODE_SET_SF_ENUM_TYPE(matrixType, MatrixType);


  SO_NODE_DEFINE_ENUM_VALUE(MatrixTransform, IDENTITY);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixTransform, TRANSPOSE);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixTransform, INVERSE);
  SO_NODE_DEFINE_ENUM_VALUE(MatrixTransform, INVERSE_TRANSPOSE);

  SO_NODE_ADD_FIELD(matrixTransform, (IDENTITY));
  SO_NODE_SET_SF_ENUM_TYPE(matrixTransform, MatrixTransform);
}

SoShaderStateMatrixParameter::~SoShaderStateMatrixParameter()
{

}

// State matrices only work with CG!!!
// FIXME: check up why. 20050125 mortene.
// COMMENT: Because they are only defined in CG (and not in ARB or GLSL)
//          a state matrix uniform delivers the current GL_MODELVIEW,
//          GL_PROJECTION,... matrices, which can be also accessed via
//          glstate.matrix.modelview, glstate.matrix.projection,...
//          since CG 1.2 (or earlier)
//                                           -- 20050126 martin.
void
SoShaderStateMatrixParameter::updateParameter(SoGLShaderObject *shader)
{
  if (shader->shaderType() != SoShader::CG_SHADER) return;
  if (this->name.isDefault()) return;

  this->ensureParameter(shader);

  CGGLenum type;
  switch (this->matrixType.getValue()) {
  case MODELVIEW: type = CG_GL_MODELVIEW_MATRIX; break;
  case PROJECTION: type = CG_GL_PROJECTION_MATRIX; break;
  case TEXTURE: type = CG_GL_TEXTURE_MATRIX; break;
  case MODELVIEW_PROJECTION: type = CG_GL_MODELVIEW_PROJECTION_MATRIX; break;
  default: assert(FALSE); break;
  }

  CGGLenum tform;
  switch (this->matrixTransform.getValue()) {
  case IDENTITY: tform = CG_GL_MATRIX_IDENTITY; break;
  case TRANSPOSE: tform = CG_GL_MATRIX_TRANSPOSE; break;
  case INVERSE: tform = CG_GL_MATRIX_INVERSE; break;
  case INVERSE_TRANSPOSE: tform = CG_GL_MATRIX_INVERSE_TRANSPOSE; break;
  default: assert(FALSE); break;
  }

  SoGLCgShaderParameter * param = (SoGLCgShaderParameter *)
    this->getGLShaderParameter(shader->getCacheContext());
  param->setState(shader, type, tform, this->name.getValue().getString());
}

/* **************************************************************************
 * ***                      SoShaderParameterArray1f                      ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray1f);

void SoShaderParameterArray1f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray1f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterArray1f::SoShaderParameterArray1f(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterArray1f);
  SO_NODE_ADD_FIELD(value, (0));
}

SoShaderParameterArray1f::~SoShaderParameterArray1f()
{
}

void SoShaderParameterArray1f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set1fv(shader,
                                                                this->value.getNum(),
                                                                this->value.getValues(0),
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
}

/* **************************************************************************
 * ***                      SoShaderParameterArray2f                      ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray2f);

void SoShaderParameterArray2f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray2f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterArray2f::SoShaderParameterArray2f(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterArray2f);
  SO_NODE_ADD_FIELD(value, (SbVec2f(0,0)));
}

SoShaderParameterArray2f::~SoShaderParameterArray2f()
{
}

void SoShaderParameterArray2f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);

  int     num    = this->value.getNum();
  float * buffer = NULL;

  if (num > 0) {
    buffer = new float[2*num];
    for (int i=0; i<num; i++) {
      buffer[2*i+0] = this->value[i][0];
      buffer[2*i+1] = this->value[i][1];
    }
  }

  this->getGLShaderParameter(shader->getCacheContext())
        ->set2fv(shader, num, buffer,
                 this->name.getValue().getString(),
                 this->identifier.getValue());
  if (buffer) delete[] buffer;
}

/* **************************************************************************
 * ***                      SoShaderParameterArray3f                      ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray3f);

void SoShaderParameterArray3f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray3f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterArray3f::SoShaderParameterArray3f(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterArray3f);
  SO_NODE_ADD_FIELD(value, (SbVec3f(0,0,0)));
}

SoShaderParameterArray3f::~SoShaderParameterArray3f()
{
}

void SoShaderParameterArray3f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);

  int     num    = this->value.getNum();
  float * buffer = NULL;

  if (num > 0) {
    buffer = new float[3*num];
    for (int i=0; i<num; i++) {
      buffer[3*i+0] = this->value[i][0];
      buffer[3*i+1] = this->value[i][1];
      buffer[3*i+2] = this->value[i][2];
    }
  }

  this->getGLShaderParameter(shader->getCacheContext())->set3fv(shader, num, buffer,
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
  if (buffer) delete[] buffer;
}

/* **************************************************************************
 * ***                      SoShaderParameterArray4f                      ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray4f);

void SoShaderParameterArray4f::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray4f,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterArray4f::SoShaderParameterArray4f(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterArray4f);
  SO_NODE_ADD_FIELD(value, (SbVec4f(0,0,0,0)));
}

SoShaderParameterArray4f::~SoShaderParameterArray4f()
{
}

void SoShaderParameterArray4f::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);

  int     num    = this->value.getNum();
  float * buffer = NULL;

  if (num > 0) {
    buffer = new float[4*num];
    for (int i=0; i<num; i++) {
      buffer[4*i+0] = this->value[i][0];
      buffer[4*i+1] = this->value[i][1];
      buffer[4*i+2] = this->value[i][2];
      buffer[4*i+3] = this->value[i][2];
    }
  }

  this->getGLShaderParameter(shader->getCacheContext())
        ->set4fv(shader, num, buffer,
                 this->name.getValue().getString(),
                 this->identifier.getValue());
  if (buffer) delete[] buffer;
}

/* **************************************************************************
 * ***                       SoShaderParameterMatrix                      ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterMatrix);

void SoShaderParameterMatrix::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterMatrix,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterMatrix::SoShaderParameterMatrix(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterMatrix);
  SO_NODE_ADD_FIELD(value, (SbMatrix(1,0,0,0,
                                     0,1,0,0,
                                     0,0,1,0,
                                     0,0,0,1)));
}

SoShaderParameterMatrix::~SoShaderParameterMatrix()
{
}

void SoShaderParameterMatrix::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);

  this->getGLShaderParameter(shader->getCacheContext())->setMatrix(shader,  this->value.getValue()[0],
                                                                   this->name.getValue().getString(),
                                                                   this->identifier.getValue());
}

/* **************************************************************************
 * ***                       SoShaderParameterMatrixArray                 ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterMatrixArray);

void SoShaderParameterMatrixArray::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterMatrixArray,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameterMatrixArray::SoShaderParameterMatrixArray(void)
{
  SO_NODE_CONSTRUCTOR(SoShaderParameterMatrixArray);
  SO_NODE_ADD_FIELD(value, (SbMatrix(1,0,0,0,
                                     0,1,0,0,
                                     0,0,1,0,
                                     0,0,0,1)));
}

SoShaderParameterMatrixArray::~SoShaderParameterMatrixArray()
{
}

void SoShaderParameterMatrixArray::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);

  int     num    = this->value.getNum();
  float * buffer = NULL;

  if (num > 0) {
    buffer = new float[16*num];
    float * ptr = buffer;
    for (int i=0; i<num; i++) {
      const float * matrix = this->value[i].getValue()[0];
      for (int j=0; j<16; j++) *(ptr++) = *(matrix++);
    }
  }

  this->getGLShaderParameter(shader->getCacheContext())
        ->setMatrixArray(shader, num, buffer,
                         this->name.getValue().getString(),
                         this->identifier.getValue());

  if (buffer) delete[] buffer;
}


/* **************************************************************************
 * *** SoShaderParameter1i ***
 * **************************************************************************/


SO_NODE_SOURCE(SoShaderParameter1i);

void SoShaderParameter1i::initClass(void)
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter1i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

SoShaderParameter1i::SoShaderParameter1i(void)
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter1i);
  SO_NODE_ADD_FIELD(value, (0));
}

SoShaderParameter1i::~SoShaderParameter1i()
{
}

void
SoShaderParameter1i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set1i(shader,
                                                               this->value.getValue(),
                                                               this->name.getValue().getString(),
                                                               this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter2i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameter2i);

SoShaderParameter2i::SoShaderParameter2i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter2i);
  SO_NODE_ADD_FIELD(value, (0,0));
}

SoShaderParameter2i::~SoShaderParameter2i()
{
}

void
SoShaderParameter2i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter2i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameter2i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set2i(shader,
                                                               this->value.getValue().getValue(),
                                                               this->name.getValue().getString(),
                                                               this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter3i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameter3i);

SoShaderParameter3i::SoShaderParameter3i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter3i);
  SO_NODE_ADD_FIELD(value, (0,0,0));
}

SoShaderParameter3i::~SoShaderParameter3i()
{
}

void
SoShaderParameter3i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter3i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameter3i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set3i(shader,
                                                               this->value.getValue().getValue(),
                                                               this->name.getValue().getString(),
                                                               this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameter4i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameter4i);

SoShaderParameter4i::SoShaderParameter4i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameter4i);
  SO_NODE_ADD_FIELD(value, (0,0,0,0));
}

SoShaderParameter4i::~SoShaderParameter4i()
{
}

void
SoShaderParameter4i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameter4i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameter4i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set4i(shader,
                                                               this->value.getValue().getValue(),
                                                               this->name.getValue().getString(),
                                                               this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameterArray1i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray1i);

SoShaderParameterArray1i::SoShaderParameterArray1i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameterArray1i);
  SO_NODE_ADD_FIELD(value, (0));
}

SoShaderParameterArray1i::~SoShaderParameterArray1i()
{
}

void
SoShaderParameterArray1i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray1i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameterArray1i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set1iv(shader,
                                                                this->value.getNum(),
                                                                (const int32_t*) this->value.getValues(0),
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameterArray2i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray2i);

SoShaderParameterArray2i::SoShaderParameterArray2i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameterArray2i);
  SO_NODE_ADD_FIELD(value, (0,0));
}

SoShaderParameterArray2i::~SoShaderParameterArray2i()
{
}

void
SoShaderParameterArray2i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray2i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameterArray2i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set2iv(shader,
                                                                this->value.getNum(),
                                                                (const int32_t*) this->value.getValues(0),
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameterArray3i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray3i);

SoShaderParameterArray3i::SoShaderParameterArray3i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameterArray3i);
  SO_NODE_ADD_FIELD(value, (0,0,0));
}

SoShaderParameterArray3i::~SoShaderParameterArray3i()
{
}

void
SoShaderParameterArray3i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray3i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameterArray3i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set3iv(shader,
                                                                this->value.getNum(),
                                                                (const int32_t*) this->value.getValues(0),
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
}

/* **************************************************************************
 * *** SoShaderParameterArray4i ***
 * **************************************************************************/

SO_NODE_SOURCE(SoShaderParameterArray4i);

SoShaderParameterArray4i::SoShaderParameterArray4i()
{
  SO_NODE_INTERNAL_CONSTRUCTOR(SoShaderParameterArray4i);
  SO_NODE_ADD_FIELD(value, (0,0,0,0));
}

SoShaderParameterArray4i::~SoShaderParameterArray4i()
{
}

void
SoShaderParameterArray4i::initClass()
{
  SO_NODE_INTERNAL_INIT_CLASS(SoShaderParameterArray4i,
                              SO_FROM_COIN_2_5|SO_FROM_INVENTOR_5_0);
}

void
SoShaderParameterArray4i::updateParameter(SoGLShaderObject *shader)
{
  this->ensureParameter(shader);
  this->getGLShaderParameter(shader->getCacheContext())->set4iv(shader,
                                                                this->value.getNum(),
                                                                (const int32_t*) this->value.getValues(0),
                                                                this->name.getValue().getString(),
                                                                this->identifier.getValue());
}

#undef PRIVATE

// *************************************************************************

#ifdef COIN_TEST_SUITE

BOOST_AUTO_TEST_CASE(initialized)
{
  {
    SoShaderParameter1f * parameter1f = new SoShaderParameter1f;
    assert(parameter1f);
    parameter1f->ref();
    BOOST_CHECK_MESSAGE(parameter1f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter1f->unref();
  }
  {
    SoShaderParameter1i * parameter1i = new SoShaderParameter1i;
    assert(parameter1i);
    parameter1i->ref();
    BOOST_CHECK_MESSAGE(parameter1i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter1i->unref();
  }
  {
    SoShaderParameter2f * parameter2f = new SoShaderParameter2f;
    assert(parameter2f);
    parameter2f->ref();
    BOOST_CHECK_MESSAGE(parameter2f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter2f->unref();
  }
  {
    SoShaderParameter2i * parameter2i = new SoShaderParameter2i;
    assert(parameter2i);
    parameter2i->ref();
    BOOST_CHECK_MESSAGE(parameter2i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter2i->unref();
  }
  {
    SoShaderParameter3f * parameter3f = new SoShaderParameter3f;
    assert(parameter3f);
    parameter3f->ref();
    BOOST_CHECK_MESSAGE(parameter3f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter3f->unref();
  }
  {
    SoShaderParameter3i * parameter3i = new SoShaderParameter3i;
    assert(parameter3i);
    parameter3i->ref();
    BOOST_CHECK_MESSAGE(parameter3i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter3i->unref();
  }
  {
    SoShaderParameter4f * parameter4f = new SoShaderParameter4f;
    assert(parameter4f);
    parameter4f->ref();
    BOOST_CHECK_MESSAGE(parameter4f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter4f->unref();
  }
  {
    SoShaderParameter4i * parameter4i = new SoShaderParameter4i;
    assert(parameter4i);
    parameter4i->ref();
    BOOST_CHECK_MESSAGE(parameter4i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parameter4i->unref();
  }

  {
    SoShaderParameterArray1f * parametera1f = new SoShaderParameterArray1f;
    assert(parametera1f);
    parametera1f->ref();
    BOOST_CHECK_MESSAGE(parametera1f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera1f->unref();
  }
  {
    SoShaderParameterArray1i * parametera1i = new SoShaderParameterArray1i;
    assert(parametera1i);
    parametera1i->ref();
    BOOST_CHECK_MESSAGE(parametera1i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera1i->unref();
  }
  {
    SoShaderParameterArray2f * parametera2f = new SoShaderParameterArray2f;
    assert(parametera2f);
    parametera2f->ref();
    BOOST_CHECK_MESSAGE(parametera2f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera2f->unref();
  }
  {
    SoShaderParameterArray2i * parametera2i = new SoShaderParameterArray2i;
    assert(parametera2i);
    parametera2i->ref();
    BOOST_CHECK_MESSAGE(parametera2i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera2i->unref();
  }
  {
    SoShaderParameterArray3f * parametera3f = new SoShaderParameterArray3f;
    assert(parametera3f);
    parametera3f->ref();
    BOOST_CHECK_MESSAGE(parametera3f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera3f->unref();
  }
  {
    SoShaderParameterArray3i * parametera3i = new SoShaderParameterArray3i;
    assert(parametera3i);
    parametera3i->ref();
    BOOST_CHECK_MESSAGE(parametera3i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera3i->unref();
  }
  {
    SoShaderParameterArray4f * parametera4f = new SoShaderParameterArray4f;
    assert(parametera4f);
    parametera4f->ref();
    BOOST_CHECK_MESSAGE(parametera4f->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera4f->unref();
  }
  {
    SoShaderParameterArray4i * parametera4i = new SoShaderParameterArray4i;
    assert(parametera4i);
    parametera4i->ref();
    BOOST_CHECK_MESSAGE(parametera4i->getTypeId() != SoType::badType(),
                        "missing class initialization");
    parametera4i->unref();
  }

  {
    SoShaderParameterMatrix * matrix = new SoShaderParameterMatrix;
    assert(matrix);
    matrix->ref();
    BOOST_CHECK_MESSAGE(matrix->getTypeId() != SoType::badType(),
                        "missing class initialization");
    matrix->unref();
  }
  {
    SoShaderParameterMatrixArray * matrixarray = new SoShaderParameterMatrixArray;
    assert(matrixarray);
    matrixarray->ref();
    BOOST_CHECK_MESSAGE(matrixarray->getTypeId() != SoType::badType(),
                        "missing class initialization");
    matrixarray->unref();
  }

  {
    SoShaderStateMatrixParameter * statematrix = new SoShaderStateMatrixParameter;
    assert(statematrix);
    statematrix->ref();
    BOOST_CHECK_MESSAGE(statematrix->getTypeId() != SoType::badType(),
                        "missing class initialization");
    statematrix->unref();
  }
}

#endif // COIN_TEST_SUITE
