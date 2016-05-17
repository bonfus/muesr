#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "mat3.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

//Dipolar Tensor
//
// This function constructs the dipolar tensor with the positions given in in_positions
void DT(const double *in_positions, 
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, unsigned int size,
          double *out_field) 
{

    unsigned int scx, scy, scz = 10; //supercell sizes
    unsigned int i,j,k; // counters for supercells
    
    struct vec3 atmpos;
    struct vec3 muonpos;
    struct vec3 r;
//    struct vec3 m;   //magnetic moment of atom
//    struct vec3 u;   // unit vector
        
    struct mat3 sc_lat;
    
    double n;
//    double c,s; //cosine and sine of K.R
    double onebrcube; // 1/r^3
    double onebrfive; // 1/r^5
    
//    struct vec3 sk  ;
//    struct vec3 isk ;
//    double  phi ;
    
    struct mat3 A, D;
#ifdef _OPENMP
    double Bxx=0.0; double Bxy=0.0; double Bxz=0.0;
    double Byx=0.0; double Byy=0.0; double Byz=0.0;
    double Bzx=0.0; double Bzy=0.0; double Bzz=0.0;
#endif    
    
    unsigned int atom;     // counter for atoms
    
    
    
    
    
    
    // define dupercell size
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

#ifdef _DEBUG    
    printf("I use: %i %i %i\n",scx, scy, scz);    
    printf("Size is: %i\n",size);
#endif 
    
    sc_lat.a.x = in_cell[0];
    sc_lat.a.y = in_cell[1];
    sc_lat.a.z = in_cell[2];
    sc_lat.b.x = in_cell[3];
    sc_lat.b.y = in_cell[4];
    sc_lat.b.z = in_cell[5];
    sc_lat.c.x = in_cell[6];
    sc_lat.c.y = in_cell[7];
    sc_lat.c.z = in_cell[8];
 
#ifdef _DEBUG      
    for (i=0;i<3;i++)
        printf("Cell is: %i %e %e %e\n",i,in_cell[i*3],in_cell[i*3+1],in_cell[i*3+2]);
        
    //printf("a %e %e %e\n", sc_lat.a.x, sc_lat.a.y, sc_lat.a.z);
#endif     

    sc_lat = mat3_mul(
                        mat3_diag((float) scx, (float) scy, (float) scz),
                        sc_lat);
    
    
    // muon position in reduced coordinates
    muonpos.x =  (in_muonpos[0] + (scx/2) ) / (float) scx;
    muonpos.y =  (in_muonpos[1] + (scy/2) ) / (float) scy;
    muonpos.z =  (in_muonpos[2] + (scz/2) ) / (float) scz;

#ifdef _DEBUG
    printf("Muon pos (frac): %e %e %e\n",muonpos.x,muonpos.y,muonpos.z);
#endif

    muonpos = mat3_vmul(muonpos,sc_lat);

#ifdef _DEBUG
    printf("Muon pos (cart): %e %e %e\n",muonpos.x,muonpos.y,muonpos.z);
#endif


#ifdef _DEBUG
    for (atom = 0; atom < size; ++atom)
    {
                    
        // atom position in reduced coordinates
        atmpos.x =  in_positions[3*atom] ;
        atmpos.y =  in_positions[3*atom+1] ;
        atmpos.z =  in_positions[3*atom+2] ;
        
        printf("Atom pos (crys): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
        
        // go to cartesian coordinates (in Angstrom!)
        atmpos = mat3_vmul(atmpos,sc_lat);  
        
        printf("Atom pos (cart): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
    }
        
#endif


A = mat3_zero();
D = mat3_zero();

#pragma omp parallel
{    
#pragma omp for collapse(3) private(i,j,k,atom,r,n,atmpos,D,onebrcube,onebrfive) reduction(+:Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz)
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                // loop over atoms
                for (atom = 0; atom < size; ++atom)
                {
                    
                    // atom position in reduced coordinates
                    atmpos.x = ( in_positions[3*atom] + (float) i) / (float) scx;
                    atmpos.y = ( in_positions[3*atom+1] + (float) j) / (float) scy;
                    atmpos.z = ( in_positions[3*atom+2] + (float) k) / (float) scz;
                    

                    
                    // go to cartesian coordinates (in Angstrom!)
                    atmpos = mat3_vmul(atmpos,sc_lat);
                    
                    //printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z);
                    // difference between atom pos and muon pos (cart coordinates)
                    
                    r = vec3_sub(atmpos,muonpos);
                    
                    n = vec3_norm(r);
                    if (n < radius)
                    {


                        // vector
                        onebrcube = 1.0/pow(n,3);
                        onebrfive = 1.0/pow(n,5);
                        
                        
                        // See uSR bible (Yaouanc Dalmas De Reotier, page 81)
                        // alpha = x
                        D.a.x = -onebrcube+3.0*r.x*r.x*onebrfive;
                        D.a.y = 3.0*r.x*r.y*onebrfive;
                        D.a.z = 3.0*r.x*r.z*onebrfive;
                        
                        // alpha = y
                        D.b.x = D.a.y;
                        D.b.y = -onebrcube+3.0*r.y*r.y*onebrfive;
                        D.b.z = 3.0*r.y*r.z*onebrfive;
                        
                        // alpha = z
                        D.c.x = D.a.z;
                        D.c.y = D.b.z;
                        D.c.z = -onebrcube+3.0*r.z*r.z*onebrfive;
                        
#ifdef _OPENMP                        
                        Bxx += D.a.x; Bxy += D.a.y; Bxz += D.a.z;
                        Byx += D.b.x; Byy += D.b.y; Byz += D.b.z;
                        Bzx += D.c.x; Bzy += D.c.y; Bzz += D.c.z;
                        
#else
                        A = mat3_add( A,D );
#endif                                    
                        
                    }                    

                }
                
            }
        }
    }
}
#ifdef _OPENMP
    A.a.x = Bxx; A.a.y = Bxy; A.a.z = Bxz;
    A.b.x = Byx; A.b.y = Byy; A.b.z = Byz;
    A.c.x = Bzx; A.c.y = Bzy; A.c.z = Bzz;
#endif     
    // B = vec3_muls(0.9274009, B); // we should multiply for a volume
    out_field[0] = A.a.x; out_field[1] = A.a.y; out_field[2] = A.a.z;
    out_field[3] = A.b.x; out_field[4] = A.b.y; out_field[5] = A.b.z;
    out_field[6] = A.c.x; out_field[7] = A.c.y; out_field[8] = A.c.z;

}


