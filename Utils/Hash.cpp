#include "Hash.h"

USE_PRJ_NAMESPACE

template <class T>
FORCE_INLINE void hash_combine(std::size_t& seed,const T& v)
{
  std::hash<T> hasher;
  seed^=hasher(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}
size_t Hash::operator()(const Vec2i& key) const
{
  size_t seed=0;
  hash_combine<sizeType>(seed,key[0]);
  hash_combine<sizeType>(seed,key[1]);
  return seed;
}
size_t Hash::operator()(const Vec3i& key) const
{
  size_t seed=0;
  hash_combine<sizeType>(seed,key[0]);
  hash_combine<sizeType>(seed,key[1]);
  hash_combine<sizeType>(seed,key[2]);
  return seed;
}
size_t Hash::operator()(const Vec4i& key) const
{
  size_t seed=0;
  hash_combine<sizeType>(seed,key[0]);
  hash_combine<sizeType>(seed,key[1]);
  hash_combine<sizeType>(seed,key[2]);
  hash_combine<sizeType>(seed,key[3]);
  return seed;
}
size_t Hash::operator()(const std::pair<Vec3i,Vec3i>& key) const
{
  size_t seed=0;
  hash_combine<sizeType>(seed,key.first[0]);
  hash_combine<sizeType>(seed,key.first[1]);
  hash_combine<sizeType>(seed,key.first[2]);
  hash_combine<sizeType>(seed,key.second[0]);
  hash_combine<sizeType>(seed,key.second[1]);
  hash_combine<sizeType>(seed,key.second[2]);
  return seed;
}
