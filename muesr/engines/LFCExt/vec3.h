#ifndef VEC3_H
#define VEC3_H

typedef double scalar;

struct vec3 {
	scalar x, y, z;
};

struct vec3 _vec3(scalar, scalar, scalar);
struct vec3 vec3_zero(void);
struct vec3 vec3_add(struct vec3, struct vec3);
struct vec3 vec3_sub(struct vec3, struct vec3);
struct vec3 vec3_mul(struct vec3, struct vec3);
struct vec3 vec3_muls(scalar, struct vec3);
double vec3_norm(struct vec3 v);
double vec3_dot(struct vec3 v, struct vec3 u);
#endif
