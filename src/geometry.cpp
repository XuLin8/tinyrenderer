#include "geometry.h"

// �Ӹ������͵���ά���������������͵���ά����
template <> template <> vec<3, int>  ::vec(const vec<3, float>& v) : x(int(v.x + .5f)), y(int(v.y + .5f)), z(int(v.z + .5f)) {}

// ���������͵���ά�������측�����͵���ά����
template <> template <> vec<3, float>::vec(const vec<3, int>& v) : x(v.x), y(v.y), z(v.z) {}

// �Ӹ������͵Ķ�ά���������������͵Ķ�ά����
template <> template <> vec<2, int>  ::vec(const vec<2, float>& v) : x(int(v.x + .5f)), y(int(v.y + .5f)) {}

// ���������͵Ķ�ά�������측�����͵Ķ�ά����
template <> template <> vec<2, float>::vec(const vec<2, int>& v) : x(v.x), y(v.y) {}