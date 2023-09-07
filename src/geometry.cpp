#include "geometry.h"

// 从浮点类型的三维向量构造整型类型的三维向量
template <> template <> vec<3, int>  ::vec(const vec<3, float>& v) : x(int(v.x + .5f)), y(int(v.y + .5f)), z(int(v.z + .5f)) {}

// 从整型类型的三维向量构造浮点类型的三维向量
template <> template <> vec<3, float>::vec(const vec<3, int>& v) : x(v.x), y(v.y), z(v.z) {}

// 从浮点类型的二维向量构造整型类型的二维向量
template <> template <> vec<2, int>  ::vec(const vec<2, float>& v) : x(int(v.x + .5f)), y(int(v.y + .5f)) {}

// 从整型类型的二维向量构造浮点类型的二维向量
template <> template <> vec<2, float>::vec(const vec<2, int>& v) : x(v.x), y(v.y) {}