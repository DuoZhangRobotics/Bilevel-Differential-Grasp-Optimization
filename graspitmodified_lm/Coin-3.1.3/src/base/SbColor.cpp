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
  \class SbColor SbColor.h Inventor/SbColor.h
  \brief The SbColor class contains the red, green and blue components
  which make up a color value.
  \ingroup base

  This class is used within other classes in Coin.  It inherits the
  SbVec3f class, interpreting the 3 component vector as a vector in
  the RGB cube where the red, green and blue components corresponds to
  x, y and z respectively.

  SbColor also adds a few extra methods for convenient handling of
  setting and getting color values as 32 bit packed values or as HSV
  values.

  \sa SbVec3f, SbColor4f */

#include <Inventor/SbColor.h>
#if COIN_DEBUG
#include <Inventor/errors/SoDebugError.h>
#endif // COIN_DEBUG

/*!
  Default constructor. The color value will be uninitialized.
 */
SbColor::SbColor(void)
  : SbVec3f()
{
}

/*!
  Construct and initialize an SbColor with the red, green and blue
  values given by the \c v vector.
 */
SbColor::SbColor(const SbVec3f& v)
  : SbVec3f(v)
{
}

/*!
  Construct and initialize an SbColor with the red, green and blue
  taken from given \c rgb array.
 */
SbColor::SbColor(const float* const rgb)
  : SbVec3f(rgb)
{
}

/*!
  Construct and initialize an SbColor with the given red, green and blue
  values.
 */
SbColor::SbColor(const float r, const float g, const float b)
  : SbVec3f(r, g, b)
{
}

/*!
  Set the color value as a 32 bit combined red/green/blue/alpha value.
  Each component is 8 bit wide (i.e. from 0x00 to 0xff), and the red
  value should be stored leftmost, like this: 0xRRGGBBAA.

  The transparency value is not stored internally in SbColor, just
  converted to a transparency value in [0, 1] and returned in the
  \c transparency field. A value of 1.0 means completely transparent
  and a value of 0.0 is completely opaque.

  \sa getPackedValue().
 */
SbColor&
SbColor::setPackedValue(const uint32_t rgba, float& transparency)
{
  // Should work regardless of endianness on run-time architecture.
  this->setValue((rgba >> 24)/255.0f,
                 ((rgba >> 16)&0xff)/255.0f,
                 ((rgba >> 8)&0xff)/255.0f);
  transparency = 1.0f - (rgba&0xff)/255.0f;
  return *this;
}

/*!
  Return color as a 32 bit packed integer in the form 0xRRGGBBAA. The
  transparency part of the return value is taken from the supplied
  \c transparency argument.

  \sa setPackedValue().
 */
uint32_t
SbColor::getPackedValue(const float transparency) const
{
  // Should work regardless of endianness on run-time architecture.
  return ((static_cast<uint32_t>(red()*255.0f + 0.5f) << 24) |
          (static_cast<uint32_t>(green()*255.0f + 0.5f) << 16) |
          (static_cast<uint32_t>(blue()*255.0f + 0.5f) << 8) |
          static_cast<uint32_t>((1.0f - transparency)*255.0f + 0.5f));
}

/*!
  Set the color as a \c hue, \c saturation, \c value triplet.
  The hue component should be normalized to within [0, 1] before you
  call this method, where 0 is equal to 0� and 1 is equal to 360�.

  \sa getHSVValue().
 */
SbColor&
SbColor::setHSVValue(float hue, float saturation,
                     float value)
{
#if COIN_DEBUG
  if (!(hue>=0.0f && hue<=1.0f) ||
      !(saturation>=0.0f && saturation<=1.0f) ||
      !(value>=0.0f && value<=1.0f)) {
    SoDebugError::postWarning("SbColor::setHSVValue",
                              "One or more values out of range, clamping.");
    hue = SbClamp(hue, 0.f, 1.f);
    saturation = SbClamp(saturation, 0.f, 1.f);
    value = SbClamp(value, 0.f, 1.f);
  }
#endif // COIN_DEBUG

  // HSV to RGB conversion routine based on the one presented in
  // "Computer Graphics: Principles and Practice", 2nd ed., by Foley et
  // al.
  //
  // Valid input range for all values are [0, 1].
  //                                                       19980809 mortene.

  float h = hue;
  float s = saturation;
  float v = value;

  if(h == 1.0f) h = 0.0f;
  h *= 6.0f;
  int i = static_cast<int>(floor(h));
  float fraction = h - static_cast<float>(i);
  float p = v * (1.0f - s);
  float q = v * (1.0f - (s * fraction));
  float t = v * (1.0f - (s * (1.0f - fraction)));

  switch(i) {
  case 0: this->setValue(v, t, p); break;
  case 1: this->setValue(q, v, p); break;
  case 2: this->setValue(p, v, t); break;
  case 3: this->setValue(p, q, v); break;
  case 4: this->setValue(t, p, v); break;
  case 5: this->setValue(v, p, q); break;
  }

  return *this;
}

/*!
  Return the color as a \c hue, \c saturation, \c value triplet.

  \sa setHSVValue().
 */
void
SbColor::getHSVValue(float &h, float &s, float &v) const
{
  // RGB to HSV conversion routine based on the one presented in
  // "Computer Graphics: Principles and Practice", 2nd ed., by Foley et
  // al.
  //
  // Hue is normalized to [0, 1]. Undefined hue and value is set to 0.0
  // on valid <r, g, b> triplets. Results are undefined when any of the
  // red, green or blue values are invalid (i.e. outside [0, 1]).
  //
  //                                                     19980809 mortene.
  h = s = 0.0f;
  v = SbMax(SbMax(red(), green()), blue());
  float min = SbMin(SbMin(red(), green()), blue());

  float delta = v - min;
  if(v > 0.0f) s = delta/v;
  if(delta > 0.0f) {
    if(red() == v) h = (green() - blue())/delta;
    else if(green() == v) h = 2 + (blue() - red())/delta;
    else if(blue() == v) h = 4 + (red() - green())/delta;

    h *= 60.0f;
    if(h < 0.0f) h += 360.0f;
    h /= 360.0f;
  }
}

/*!
  Set the color as a \c hue, \c saturation, \c value triplet.
  The hue component should be normalized to within [0, 1] before you
  call this method, where 0 is equal to 0� and 1 is equal to 360�.

  \sa getHSVValue().
 */
SbColor &
SbColor::setHSVValue(const float hsv[3])
{
  return this->setHSVValue(hsv[0], hsv[1], hsv[2]);
}

/*!
  Return the color as a \c hue, \c saturation, \c value triplet.

  \sa setHSVValue().
 */
void
SbColor::getHSVValue(float hsv[3]) const
{
  this->getHSVValue(hsv[0], hsv[1], hsv[2]);
}
