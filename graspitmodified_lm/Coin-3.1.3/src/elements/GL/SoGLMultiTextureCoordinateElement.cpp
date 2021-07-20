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
  \class SoGLMultiTextureCoordinateElement Inventor/elements/SoGLMultiTextureCoordinateElement.h
  \brief The SoGLMultiTextureCoordinateElement class stores the current gltexture coordinates for several units.
  \ingroup elements
*/

#include <Inventor/elements/SoGLMultiTextureCoordinateElement.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/misc/SoState.h>

#include <cassert>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/system/gl.h>
#include <Inventor/C/glue/gl.h>

#define MAX_UNITS 16 // FIXME: make dynamic? pederb, 2003-11-03

class SoGLMultiTextureCoordinateElementP {
public:
  SoGLMultiTextureCoordinateElement::GLUnitData unitdata[MAX_UNITS];
  int contextid;

  // switch/case table for faster rendering.
  enum SendLookup {
    NONE,
    FUNCTION,
    TEXCOORD2,
    TEXCOORD3,
    TEXCOORD4
  };
  SendLookup sendlookup[MAX_UNITS];
};

#define PRIVATE(obj) obj->pimpl

SO_ELEMENT_CUSTOM_CONSTRUCTOR_SOURCE(SoGLMultiTextureCoordinateElement);

/*!
  This static method initializes static data for the
  SoGLMultiTextureCoordinateElement class.
*/

void
SoGLMultiTextureCoordinateElement::initClass(void)
{
  SO_ELEMENT_INIT_CLASS(SoGLMultiTextureCoordinateElement, inherited);
}

/*!
  The constructor.
*/
SoGLMultiTextureCoordinateElement::SoGLMultiTextureCoordinateElement(void)
{
  PRIVATE(this) = new SoGLMultiTextureCoordinateElementP;

  this->setTypeId(SoGLMultiTextureCoordinateElement::classTypeId);
  this->setStackIndex(SoGLMultiTextureCoordinateElement::classStackIndex);
}

/*!
  The destructor.
*/

SoGLMultiTextureCoordinateElement::~SoGLMultiTextureCoordinateElement()
{
  delete PRIVATE(this);
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::init(SoState * state)
{
  SoAction * action = state->getAction();
  assert(action->isOfType(SoGLRenderAction::getClassTypeId()));
  // fetch cache context id from action since SoGLCacheContextElement
  // might not be initialized yet.
  SoGLRenderAction * glaction = (SoGLRenderAction*) action;
  PRIVATE(this)->contextid = glaction->getCacheContext();

  inherited::init(state);
  for (int i = 0; i < MAX_UNITS; i++) {
    GLUnitData & ud = PRIVATE(this)->unitdata[i];
    ud.texgenCB = NULL;
    ud.texgenData = NULL;
  }
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::push(SoState * state)
{
  inherited::push(state);
  SoGLMultiTextureCoordinateElement * prev = (SoGLMultiTextureCoordinateElement*)this->getNextInStack();

  PRIVATE(this)->contextid = PRIVATE(prev)->contextid;
  for (int i = 0; i < MAX_UNITS; i++) {
    PRIVATE(this)->unitdata[i] = PRIVATE(prev)->unitdata[i];
  }
  // capture previous element since we might or might not change the
  // GL state in set/pop
  prev->capture(state);
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::pop(SoState * state,
                                       const SoElement * prevTopElement)
{
  inherited::pop(state, prevTopElement);
  SoGLMultiTextureCoordinateElement * prev = (SoGLMultiTextureCoordinateElement*) prevTopElement;

  const cc_glglue * glue = cc_glglue_instance(PRIVATE(this)->contextid);

  for (int i = 0; i < MAX_UNITS; i++) {
    const GLUnitData & thisud = PRIVATE(this)->unitdata[i];
    const GLUnitData & prevud = PRIVATE(prev)->unitdata[i];

    SbBool enablegen = FALSE;
    SbBool disablegen = FALSE;
    SbBool docallback = FALSE;

    if (thisud.texgenCB && !prevud.texgenCB) {enablegen = TRUE; docallback = TRUE;}
    else if (!thisud.texgenCB && prevud.texgenCB) disablegen = TRUE;
    else if (thisud.texgenCB != prevud.texgenCB) docallback = TRUE;

    if (enablegen || disablegen || docallback) {
      // must change texture unit while updating OpenGL
      cc_glglue_glActiveTexture(glue, (GLenum) (int(GL_TEXTURE0) + i));
    }
    if (enablegen) {
      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
      glEnable(GL_TEXTURE_GEN_R);
      glEnable(GL_TEXTURE_GEN_Q);
    }
    if (disablegen) {
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_GEN_R);
      glDisable(GL_TEXTURE_GEN_Q);
    }
    if (docallback) {
      this->doCallback(i);
    }
    // restore default unit
    if (enablegen || disablegen || docallback) {
      cc_glglue_glActiveTexture(glue, (GLenum) GL_TEXTURE0);
    }
  }
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::setTexGen(SoState * const state,
                                             SoNode * const node,
                                             const int unit,
                                             SoTexCoordTexgenCB * const texgenFunc,
                                             void * const texgenData,
                                             SoTextureCoordinateFunctionCB * const func,
                                             void * const funcData)
{
  SoMultiTextureCoordinateElement::setFunction(state, node, unit, func, funcData);

  SoGLMultiTextureCoordinateElement *element = (SoGLMultiTextureCoordinateElement *)
    SoElement::getElement(state, classStackIndex);
  if (element) {
    element->setElt(unit, texgenFunc, texgenData);
  }
}

//!  FIXME: write doc.

SoTextureCoordinateElement::CoordType
SoGLMultiTextureCoordinateElement::getType(const int unit) const
{
  assert(unit >= 0 && unit < MAX_UNITS);
  if (PRIVATE(this)->unitdata[unit].texgenCB) return SoTextureCoordinateElement::NONE;
  return inherited::getType(unit);
}

//!  FIXME: write doc.

const SoGLMultiTextureCoordinateElement *
SoGLMultiTextureCoordinateElement::getInstance(SoState * const state)
{
  return (SoGLMultiTextureCoordinateElement*)
    SoElement::getConstElement(state, classStackIndex);
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::send(const int unit, const int index) const
{
  const UnitData & ud = this->getUnitData(unit);
  GLenum glunit = (GLenum) (int(GL_TEXTURE0) + unit);
  const cc_glglue * glue = cc_glglue_instance(PRIVATE(this)->contextid);

  switch (PRIVATE(this)->sendlookup[unit]) {
  case SoGLMultiTextureCoordinateElementP::NONE:
    break;
  case SoGLMultiTextureCoordinateElementP::FUNCTION:
    assert(0 && "should not happen");
    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD2:
    assert(index < ud.numCoords);
    cc_glglue_glMultiTexCoord2fv(glue, glunit, ud.coords2[index].getValue());
    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD3:
    cc_glglue_glMultiTexCoord3fv(glue, glunit, ud.coords3[index].getValue());
    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD4:
    cc_glglue_glMultiTexCoord4fv(glue, glunit, ud.coords4[index].getValue());
    break;
  default:
    assert(0 && "should not happen");
    break;
  }
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::send(const int unit,
                                        const int index,
                                        const SbVec3f &c,
                                        const SbVec3f &n) const
{
  const UnitData & ud = this->getUnitData(unit);
  GLenum glunit = (GLenum) (int(GL_TEXTURE0) + unit);
  const cc_glglue * glue = cc_glglue_instance(PRIVATE(this)->contextid);

  switch (PRIVATE(this)->sendlookup[unit]) {
  case SoGLMultiTextureCoordinateElementP::NONE:
    break;
  case SoGLMultiTextureCoordinateElementP::FUNCTION:
    assert(ud.funcCB);
    cc_glglue_glMultiTexCoord4fv(glue, glunit,
                                 ud.funcCB(ud.funcCBData, c, n).getValue());

    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD2:
    cc_glglue_glMultiTexCoord2fv(glue, glunit, ud.coords2[index].getValue());
    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD3:
    cc_glglue_glMultiTexCoord3fv(glue, glunit, ud.coords3[index].getValue());
    break;
  case SoGLMultiTextureCoordinateElementP::TEXCOORD4:
    cc_glglue_glMultiTexCoord4fv(glue, glunit, ud.coords4[index].getValue());
    break;
  default:
    assert(0 && "should not happen");
    break;
  }
}

//!  FIXME: write doc.

void
SoGLMultiTextureCoordinateElement::setElt(const int unit,
                                          SoTexCoordTexgenCB * func,
                                          void *data)
{
  assert(unit >= 0 && unit < MAX_UNITS);
  GLUnitData & ud = PRIVATE(this)->unitdata[unit];

  SbBool enablegen = FALSE;
  SbBool disablegen = FALSE;
  SbBool docallback = FALSE;

  if (func && !ud.texgenCB) {enablegen = TRUE; docallback = TRUE;}
  else if (!func && ud.texgenCB) disablegen = TRUE;
  else if (func && func != ud.texgenCB) docallback = TRUE;

  if (func) {
    // update SoMultiTextureCoordinateElement type
    this->getUnitData(unit).whatKind = SoTextureCoordinateElement::FUNCTION;
  }
  ud.texgenCB = func;
  ud.texgenData = data;

  const cc_glglue * glue = cc_glglue_instance(PRIVATE(this)->contextid);

  if (enablegen || disablegen || docallback) {
    cc_glglue_glActiveTexture(glue, (GLenum) (int(GL_TEXTURE0) + unit));
  }

  if (enablegen) {
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glEnable(GL_TEXTURE_GEN_R);
    glEnable(GL_TEXTURE_GEN_Q);
  }
  if (disablegen) {
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
    glDisable(GL_TEXTURE_GEN_R);
    glDisable(GL_TEXTURE_GEN_Q);
  }
  if (docallback) this->doCallback(unit);

  if (enablegen || disablegen || docallback) {
    cc_glglue_glActiveTexture(glue, (GLenum) GL_TEXTURE0);
  }
}

void
SoGLMultiTextureCoordinateElement::doCallback(const int unit) const
{
  if (PRIVATE(this)->unitdata[unit].texgenCB) {
    PRIVATE(this)->unitdata[unit].texgenCB(PRIVATE(this)->unitdata[unit].texgenData);
  }
}

/*!
  Internal method that is called from SoGLTextureCoordinateElement to
  set up optimized rendering.
*/
void
SoGLMultiTextureCoordinateElement::initRender(const SbBool * enabled, const int maxenabled) const
{
  // need writeable instance
  SoGLMultiTextureCoordinateElement * elem = ( SoGLMultiTextureCoordinateElement *) this;
  for (int i = 1; i <= maxenabled; i++) {
    PRIVATE(elem)->sendlookup[i] = SoGLMultiTextureCoordinateElementP::NONE;
    // init the sendloopup variable
    if (enabled[i]) {
      const UnitData & ud = this->getUnitData(i);
      switch (ud.whatKind) {
      case SoTextureCoordinateElement::DEFAULT:
        assert(0 && "should not happen");
        break;
      case SoTextureCoordinateElement::FUNCTION:
        if (ud.funcCB) {
          PRIVATE(elem)->sendlookup[i] = SoGLMultiTextureCoordinateElementP::FUNCTION;
        }
        break;
      case SoTextureCoordinateElement::NONE:
        break;
      case SoTextureCoordinateElement::EXPLICIT:
        {
          switch (ud.coordsDimension) {
          case 2:
            PRIVATE(elem)->sendlookup[i] = SoGLMultiTextureCoordinateElementP::TEXCOORD2;
            break;
          case 3:
            PRIVATE(elem)->sendlookup[i] = SoGLMultiTextureCoordinateElementP::TEXCOORD3;
            break;
          case 4:
            PRIVATE(elem)->sendlookup[i] = SoGLMultiTextureCoordinateElementP::TEXCOORD4;
            break;
          default:
            assert(0 && "should not happen");
            break;
          }
        }
        break;
      default:
        assert(0 && "should not happen");
        break;
      }
    }
  }
}



#undef PRIVATE
#undef MAX_UNITS
