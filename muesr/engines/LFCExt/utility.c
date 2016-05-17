#include <math.h>
#include "utility.h"

void print_vec(const char * msg, struct vec3 i)
{
    printf(msg);
    printf(" %e %e %e\n", i.x, i.y, i.z);
}

void print_mat(const char * msg, struct mat3 i)
{
    printf(msg);
    printf(" %e %e %e\n", i.a.x, i.a.y, i.a.z);
    printf(" %e %e %e\n", i.b.x, i.b.y, i.b.z);
    printf(" %e %e %e\n", i.c.x, i.c.y, i.c.z);
}
