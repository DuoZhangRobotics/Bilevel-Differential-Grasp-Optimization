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
  \class SoSField SoSField.h Inventor/fields/SoSField.h
  \brief The SoSField class is the base class for fields which contains
  only a single value.
  \ingroup fields

  All field types which should always contain only a single member
  value inherits this class. SoSField is an abstract class.

  You use methods setValue() and getValue() to store or fetch the
  value of single-value fields. Example:

  \code
    SoSpotLight * lightnode = new SoSpotLight;
    lightnode->on.setValue(TRUE); // The "on" field of SoSpotLight is
                                  // a single value field of type SoSFBool.
    ...
    ...
    // Change lightswitch.
    if (lightnode->on.getValue() == FALSE)
      lightnode->on = TRUE; // We can use operator = instead of setValue().
    else
      lightnode->on = FALSE;
    ...
  \endcode

  When nodes, engines or other types of field containers are written
  to file, their single-value fields are written to file in this
  format:

  \code
    containerclass {
      fieldname value
      fieldname value
      ...
    }
  \endcode

  ..like this, for instance, a SpotLight node from a scene
  graph which will be default \e off when read back from file:

  \code
    SpotLight {
      on FALSE
    }
  \endcode

  \sa SoMField
*/


// FIXME: document the SO_SFIELD_HEADER, SO_SFIELD_SOURCE and
// SO_SFIELD_INIT_CLASS macros. Also write a note (with link to an
// online bookstore?) about the Inventor Toolmaker, if the application
// programmer should need to extend Coin with new single-value
// fields. 20010913 mortene.

// *************************************************************************

#include <Inventor/fields/SoSField.h>

#include <cassert>

#include <Inventor/fields/SoSubField.h>
#include <Inventor/errors/SoDebugError.h>

#include "tidbitsp.h"

// *************************************************************************

SoType SoSField::classTypeId STATIC_SOTYPE_INIT;

// *************************************************************************

// Overridden from parent class.
SoType
SoSField::getClassTypeId(void)
{
  return SoSField::classTypeId;
}

/*!
  The SoSField constructor is protected, as this is an abstract
  class.
*/
SoSField::SoSField(void)
{
}

/*!
  The SoSField destructor is empty, and is only defined so we could
  make it virtual.
*/
SoSField::~SoSField()
{
}

// Overridden from parent class.
void
SoSField::initClass(void)
{
  PRIVATE_FIELD_INIT_CLASS(SoSField, "SField", inherited, NULL);
}

void
SoSField::atexit_cleanup(void)
{
  SoType::removeType(SoSField::classTypeId.getName());
  SoSField::classTypeId STATIC_SOTYPE_INIT;
}
