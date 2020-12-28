#ifndef HASH_H
#define HASH_H

#include <CommonFile/MathBasic.h>

PRJ_BEGIN

struct Hash
{
  size_t operator()(const Vec2i& key) const;
  size_t operator()(const Vec3i& key) const;
  size_t operator()(const Vec4i& key) const;
  size_t operator()(const std::pair<Vec3i,Vec3i>& key) const;
};

PRJ_END

#endif
