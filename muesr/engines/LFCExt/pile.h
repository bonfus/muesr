#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>
#include "vec3.h"

typedef struct {
	unsigned int nElements;
	double * ranks; // ,ust be positive
	struct vec3 * elements;
} pile;

void pile_init(pile * p, unsigned int nElements);

void pile_add_element(pile * p, double rank, struct vec3 v);

void pile_move_elements_from_position(pile * p, unsigned int pos);

void pile_free(pile * p);
