#ifndef MAT3_H
#define MAT3_H

#include "vec3.h"

struct mat3 {
	struct vec3 a, b, c;
};

struct mat3 mat3_zero(void);
struct mat3 mat3_identity(void);
struct mat3 mat3_diag(scalar a, scalar b, scalar c) ;
struct mat3 mat3_mul(struct mat3, struct mat3); 
struct mat3 mat3_add(struct mat3, struct mat3); 
struct vec3 mat3_mulv(struct mat3, struct vec3);
struct vec3 mat3_vmul(struct vec3, struct mat3);

struct mat3 mat3_aangle(struct vec3, scalar);
struct mat3 mat3_inv(struct mat3 i);

#endif
