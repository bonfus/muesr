#include <math.h>
#include "pile.h"
#include "utility.h"

int main(){
	pile p;
	
	unsigned int s = 10;
	
	pile_init(&p, s);
	
	struct vec3 v = vec3_zero() ;
	
	v.x = 1.;
	
	pile_add_element(&p,8.291562e+00, v);
	
	struct vec3 v2 = vec3_zero() ;

	v2.x = 2.;
	
	pile_add_element(&p, 4.330127e-04, v2);
	
	v.x = 3.;
	
	pile_add_element(&p, 1.0, v);
	
	v.x = 4.;
	
	pile_add_element(&p, 0.1, v);
	
	v.x = 5.;
	
	pile_add_element(&p, 0.4, v);
	
	for (unsigned int i = 0; i<s; i++)
	{
		printf("Element %d has rank %f and is %e %e %e\n", i, p.ranks[i], p.elements[i].x,p.elements[i].y,p.elements[i].z);
	}
	
	pile_free(&p);
     
}
