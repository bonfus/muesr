#include <math.h>
#include "mat3.h"
#include "utility.h"

int main(){
   struct mat3 d;
   struct mat3 dinv;
   d.a.x = 1.0; d.a.y = 0.0; d.a.z = 0.0;
   d.b.x = 0.0; d.b.y = 1.0; d.b.z = 0.0;
   d.c.x = 0.0; d.c.y = 0.0; d.c.z = 1.0;
   
   print_mat("Initial matrix\n", d);
   
   dinv = mat3_inv(d);
   
   print_mat("Inverted matrix\n", dinv);
   
   // second test


   d.a.x =0.26436857; d.a.y = 0.85196823; d.a.z = 0.07217137;
   d.b.x =0.55844147; d.b.y = 0.99388193; d.b.z = 0.07394236;
   d.c.x =0.20257205; d.c.y = 0.63541224; d.c.z = 0.36302912;

	print_mat("Initial matrix\n", d);
	// determinat -0.065914490444935595
	// expected output
	//[[-4.76108061,  3.99655223,  0.1324941 ],
    // [ 2.84841554, -1.23422918, -0.31488448],
    // [-2.32889377, -0.06981654,  3.23181197]])
    
   dinv = mat3_inv(d);
   
   print_mat("Inverted matrix\n", dinv);
    
     
}
