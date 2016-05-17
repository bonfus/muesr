#include <math.h>
#include "mat3.h"


struct mat3 mat3_zero() 
{
	struct mat3 m;
	m.a = vec3_zero();
	m.b = vec3_zero();
	m.c = vec3_zero();

	return m;
}

struct mat3 mat3_identity() 
{
	struct mat3 m;
	m.a = vec3_zero();
	m.b = vec3_zero();
	m.c = vec3_zero();

	m.a.x = 1;
	m.b.y = 1;
	m.c.z = 1;

	return m;
}

struct mat3 mat3_diag(scalar a, scalar b, scalar c) 
{
	struct mat3 m;
	m.a = vec3_zero();
	m.b = vec3_zero();
	m.c = vec3_zero();

	m.a.x = a;
	m.b.y = b;
	m.c.z = c;

	return m;
}

struct mat3 mat3_add(struct mat3 m, struct mat3 n) 
{
	struct mat3 o;
	o.a.x = m.a.x + n.a.x;	
	o.a.y = m.a.y + n.a.y;
	o.a.z = m.a.z + n.a.z;
	o.b.x = m.b.x + n.b.x;	
	o.b.y = m.b.y + n.b.y;
	o.b.z = m.b.z + n.b.z;
	o.c.x = m.c.x + n.c.x;	
	o.c.y = m.c.y + n.c.y;
	o.c.z = m.c.z + n.c.z;

	return o;
}

struct mat3 mat3_mul(struct mat3 m, struct mat3 n) 
{
	struct mat3 o;
	o.a.x = m.a.x * n.a.x + m.a.y * n.b.x + m.a.z * n.c.x;	
	o.a.y = m.a.x * n.a.y + m.a.y * n.b.y + m.a.z * n.c.y;
	o.a.z = m.a.x * n.a.z + m.a.y * n.b.z + m.a.z * n.c.z;
	o.b.x = m.b.x * n.a.x + m.b.y * n.b.x + m.b.z * n.c.x;	
	o.b.y = m.b.x * n.a.y + m.b.y * n.b.y + m.b.z * n.c.y;
	o.b.z = m.b.x * n.a.z + m.b.y * n.b.z + m.b.z * n.c.z;
	o.c.x = m.c.x * n.a.x + m.c.y * n.b.x + m.c.z * n.c.x;	
	o.c.y = m.c.x * n.a.y + m.c.y * n.b.y + m.c.z * n.c.y;
	o.c.z = m.c.x * n.a.z + m.c.y * n.b.z + m.c.z * n.c.z;
	return o;
}

struct vec3 mat3_mulv(struct mat3 m, struct vec3 v)
{
	struct vec3 u;
	u.x = m.a.x * v.x +	m.a.y * v.y + m.a.z * v.z;
	u.y = m.b.x * v.x +	m.b.y * v.y + m.b.z * v.z;
	u.z = m.c.x * v.x +	m.c.y * v.y + m.c.z * v.z;
	return u;	
}

struct vec3 mat3_vmul(struct vec3 v,struct mat3 m)
{
	struct vec3 u;
	u.x = m.a.x * v.x +	m.b.x * v.y + m.c.x * v.z;
	u.y = m.a.y * v.x +	m.b.y * v.y + m.c.y * v.z;
	u.z = m.a.z * v.x +	m.b.z * v.y + m.c.z * v.z;
	return u;	
}

// Get the rotation matrix representation a rotation 'r' around axis 'v'
// make sure 'v' is a unit vector
struct mat3 mat3_aangle(struct vec3 v, scalar r)
{
	struct mat3 m;
	scalar c, s, C;

	c = cos(r);
	s = -sin(r);
	C = 1 - c;
	m = mat3_identity();

	m.a.x = v.x * v.x * C + c;
	m.a.y = v.y * v.x * C + v.z * s;
	m.a.z = v.z * v.x * C - v.y * s;
	m.b.x = v.x * v.y * C - v.z * s;
	m.b.y = v.y * v.y * C + c;
	m.b.z = v.z * v.y * C + v.x * s;
	m.c.x = v.x * v.z * C + v.y * s;
	m.c.y = v.y * v.z * C - v.x * s;
	m.c.z = v.z * v.z * C + c;
	return m;
}

struct mat3 mat3_inv(struct mat3 i)
{
	struct mat3 minv;
	// computes the inverse of a matrix m
	double det = i.a.x * (i.b.y * i.c.z - i.c.y * i.b.z) -
				i.a.y * (i.b.x * i.c.z - i.b.z * i.c.x) +
				i.a.z * (i.b.x * i.c.y - i.b.y * i.c.x);
	
	double invdet = 1.0 / det;
	
	minv.a.x = (i.b.y * i.c.z - i.c.y * i.b.z) * invdet;
	minv.a.y = (i.a.z * i.c.y - i.a.y * i.c.z) * invdet;
	minv.a.z = (i.a.y * i.b.z - i.a.z * i.b.y) * invdet;
	minv.b.x = (i.b.z * i.c.x - i.b.x * i.c.z) * invdet;
	minv.b.y = (i.a.x * i.c.z - i.a.z * i.c.x) * invdet;
	minv.b.z = (i.b.x * i.a.z - i.a.x * i.b.z) * invdet;
	minv.c.x = (i.b.x * i.c.y - i.c.x * i.b.y) * invdet;
	minv.c.y = (i.c.x * i.a.y - i.a.x * i.c.y) * invdet;
	minv.c.z = (i.a.x * i.b.y - i.b.x * i.a.y) * invdet;
	
	return minv;
}
