#include <math.h>
#include "vec3.h"


struct vec3 _vec3(scalar x, scalar y, scalar z)
{
	struct vec3 v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

struct vec3 vec3_zero() 
{
	struct vec3 v;
	v.x = 0;
	v.y = 0;
	v.z = 0;
	return v;
}

struct vec3 vec3_add(struct vec3 v, struct vec3 u) 
{
	struct vec3 t;
	t.x = v.x + u.x;
	t.y = v.y + u.y;
	t.z = v.z + u.z;
	return t;
}

struct vec3 vec3_sub(struct vec3 v, struct vec3 u) 
{
	struct vec3 t;
	t.x = v.x - u.x;
	t.y = v.y - u.y;
	t.z = v.z - u.z;
	return t;
}

struct vec3 vec3_mul(struct vec3 v, struct vec3 u) 
{
	struct vec3 t;
	t.x = v.x * u.x;
	t.y = v.y * u.y;
	t.z = v.z * u.z;
	return t;
}

struct vec3 vec3_muls(scalar s, struct vec3 v) 
{
	struct vec3 t;
	t.x = s * v.x;
	t.y = s * v.y;
	t.z = s * v.z;
	return t;
}

double vec3_norm(struct vec3 v)
{
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

double vec3_dot(struct vec3 v, struct vec3 u)
{
    return v.x*u.x + v.y*u.y + v.z*u.z;
}
