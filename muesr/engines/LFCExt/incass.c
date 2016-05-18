#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat3.h"
#include "pile.h"
#include "config.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
#include <omp.h>
#endif 


//Arbitrary size sum for incommensurate magnetic orders
void INCASS(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius, 
          unsigned int size, unsigned int in_nangles,
          double *out_field_cont, double *out_field_dip, double *out_field_lor) 
{

    unsigned int scx, scy, scz = 10; //supercell sizes
    unsigned int i,j,k; // counters for supercells
    
    struct vec3 atmpos;
    struct vec3 muonpos;
    
    struct vec3 r;
    struct vec3 u;   // unit vector
    
    struct mat3 sc_lat;
    struct mat3 inv_sc_lat;
    
    double n;
    double c,s; //cosine and sine of K.R
    double onebrcube; // 1/r^3
    
    double  phi ;
    
    // tmp value for speed and clearness
    struct vec3 crysvec;
    
    double stagmom[size];                    // this is m_0
    struct vec3 refatmpos[size];             // reference atom used to produce C and S    
    struct vec3 Ahelix[size], Bhelix[size]; // two unit vectors describing the helix in the m_0 (cos(phi).a +/- sin(phi).b)
    struct vec3 SDip[size], CDip[size]; // sums of contribution providing cosine and sine prefactors
    struct vec3 SLor[size], CLor[size]; // sums of contribution providing cosine and sine prefactors
    
    pile CCont, SCont;
    
    struct vec3 K;
    

    unsigned int a, angn;     // counter for atoms

    // initialize variables

    
    pile_init(&CCont, nnn_for_cont);
    pile_init(&SCont, nnn_for_cont);

    for (a = 0; a < size; ++a)
    {
        Ahelix[a] = vec3_zero();
        Bhelix[a] = vec3_zero();
        CDip[a] = vec3_zero();
        SDip[a] = vec3_zero();
        CLor[a] = vec3_zero();
        SLor[a] = vec3_zero();
    }
    
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

    inv_sc_lat = mat3_inv(sc_lat);

#ifdef _DEBUG      
    for (i=0;i<3;i++)
        printf("Cell is: %i %e %e %e\n",i,in_cell[i*3],in_cell[i*3+1],in_cell[i*3+2]);
      
    printf("Inverse cell is %e %e %e\n", inv_sc_lat.a.x, inv_sc_lat.a.y, inv_sc_lat.a.z);
    printf("Inverse cell is %e %e %e\n", inv_sc_lat.b.x, inv_sc_lat.b.y, inv_sc_lat.b.z);
    printf("Inverse cell is %e %e %e\n", inv_sc_lat.c.x, inv_sc_lat.c.y, inv_sc_lat.c.z);

    for (a = 0; a < size; ++a)
    {
                    
        // atom position in reduced coordinates
        atmpos.x =  in_positions[3*a] ;
        atmpos.y =  in_positions[3*a+1] ;
        atmpos.z =  in_positions[3*a+2] ;
        
        printf("Atom pos (crys): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
        
        // go to cartesian coordinates (in Angstrom!)
        atmpos = mat3_vmul(atmpos,sc_lat);  
        
        printf("Atom pos (cart): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
#ifndef _EXTENSION
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+1],in_fc[6*a+2],in_fc[6*a+3],in_fc[6*a+4],in_fc[6*a+5]);
#else
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+2],in_fc[6*a+4],in_fc[6*a+1],in_fc[6*a+3],in_fc[6*a+5]);
#endif
        printf("phi: %e\n",in_phi[a]);
    }
        

#endif     
    
    K.x = in_K[0];
    K.y = in_K[1];
    K.z = in_K[2];

#ifdef _DEBUG     
    printf("K is: %e %e %e \n",K.x,K.y,K.z);
    printf("Radius is: %e\n",radius);
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



    for (a = 0; a < size; ++a)
    {
        // reference atom in reduced coordinates
        //   the first atom is chosen as reference
        refatmpos[a].x =  (in_positions[3*a+0] + (scx/2) ) / (float) scx;
        refatmpos[a].y =  (in_positions[3*a+1] + (scy/2) ) / (float) scy;
        refatmpos[a].z =  (in_positions[3*a+2] + (scz/2) ) / (float) scz;


#ifdef _DEBUG
        printf("Reference atom pos (frac): %e %e %e\n",refatmpos[a].x,refatmpos[a].y,refatmpos[a].z);
#endif

        refatmpos[a] = mat3_vmul(refatmpos[a],sc_lat);

#ifdef _DEBUG
        printf("Reference atom pos (cart): %e %e %e\n",refatmpos[a].x,refatmpos[a].y,refatmpos[a].z);
#endif

        
        
        // now take care of magntism
        struct vec3 tmp;
#ifndef _EXTENSION        
        tmp.x = in_fc[6*a]; 
        tmp.y = in_fc[6*a+1]; 
        tmp.z = in_fc[6*a+2];
#else
        tmp.x = in_fc[6*a]; 
        tmp.y = in_fc[6*a+2]; 
        tmp.z = in_fc[6*a+4];
#endif
        stagmom[a] = vec3_norm(tmp);
        Ahelix[a] = vec3_muls(1.0/stagmom[a],tmp);

        
        // now B
#ifndef _EXTENSION
        printf("ERROR!!! If you see this the extension compilation went wrong!\n");
        tmp.x = in_fc[6*a+3]; 
        tmp.y = in_fc[6*a+4]; 
        tmp.z = in_fc[6*a+5];
#else
        tmp.x = in_fc[6*a+1]; 
        tmp.y = in_fc[6*a+3]; 
        tmp.z = in_fc[6*a+5];
#endif      
        // check if they are the same
        if (fabs(stagmom[a] - vec3_norm(tmp))>EPS)
        {
            printf("ERROR!!! Staggered moment is different in real and imag parts of atom %u\n Use another routine!\n",a);
        }
        Bhelix[a] =  vec3_muls(1.0/vec3_norm(tmp),tmp);
        
        if (fabs(vec3_dot(Ahelix[a],Bhelix[a])) > EPS)
        {
            printf("ERROR!!! Real and imaginary part of atom %u are not orthogonal by %e!\n",a,vec3_dot(Ahelix[a],Bhelix[a]));
        }

#ifdef _DEBUG
        printf("Unit vector a (cart): %e %e %e\n",Ahelix[a].x,Ahelix[a].y,Ahelix[a].z);
        printf("Unit vector b (cart): %e %e %e\n",Bhelix[a].x,Bhelix[a].y,Bhelix[a].z);
        printf("Stag mom: %e\n", stagmom[a]);
#endif        
        
        if (fabs(in_phi[a]) > EPS)
        {
            printf("ERROR!!! Phi not supported!!!!\n");
        }
        
    }

// parallel execution starts here
// the shared variables are listed just to remember about data races!
// other variable shaed by default: refatmpos,atmpos,phi,Ahelix,Bhelix
#pragma omp parallel shared(SDip,CDip,SLor,CLor,SCont,CCont) 
{
#pragma omp for collapse(3) private(i,j,k,a,r,n,c,s,u,crysvec,onebrcube,atmpos)
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                // loop over atoms
                for (a = 0; a < size; ++a)
                {
                    
                    // atom position in reduced coordinates
                    atmpos.x = ( in_positions[3*a] + (float) i) / (float) scx;
                    atmpos.y = ( in_positions[3*a+1] + (float) j) / (float) scy;
                    atmpos.z = ( in_positions[3*a+2] + (float) k) / (float) scz;
                    
                    
                    
                    // go to cartesian coordinates (in Angstrom!)
                    atmpos = mat3_vmul(atmpos,sc_lat);
                    
                    //printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z);
                    // difference between atom pos and muon pos (cart coordinates)
                    
                    r = vec3_sub(atmpos,muonpos);
                    
                    n = vec3_norm(r);
                    if (n < radius)
                    {

                        phi = in_phi[a];
                        
                        // unit vector
                        u = vec3_muls(1.0/n,r);
                        onebrcube = 1.0/pow(n,3);
                        
                        // go back to fractional (crystal) definition!
                        // !!!!! CHECK THIS DEFINITION !!!!
                        crysvec = mat3_vmul(vec3_sub(vec3_sub(atmpos,muonpos),refatmpos[a]),inv_sc_lat);
                        //
                        c = cos ( 2.0*M_PI * (vec3_dot(K,crysvec) + phi ) );
                        s = sin ( 2.0*M_PI * (vec3_dot(K,crysvec) + phi ) );
#ifdef _DEBUG
                        printf("crysvec : %e %e %e\n", crysvec.x, crysvec.y,crysvec.z);
                        struct vec3 tmp = vec3_sub(atmpos,muonpos);
                        printf("vec3_sub(atmpos,muonpos) : %e %e %e\n", tmp.x, tmp.y,tmp.z);
                        tmp = vec3_sub(vec3_sub(atmpos,muonpos),refatmpos[a]);
                        printf("vec3_sub(vec3_sub(atmpos,muonpos),refatmpos[a]) : %e %e %e\n",  tmp.x, tmp.y,tmp.z);
                        
                        printf("vec3_dot(K,vec3_sub(vec3_sub(atmpos,muonpos),refatmpos[a])) : %e \n",  vec3_dot(K,vec3_sub(vec3_sub(atmpos,muonpos),refatmpos[a])));
                        
                        
                        printf("cos is this time : %e \n", c);
                        printf("sin is this time : %e \n", s);

                        
                        printf("u is : %e %e %e\n", u.x, u.y, u.z);
                        tmp = vec3_muls(onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Ahelix[a],u),u), Ahelix[a]));
                        printf("A part %d is : %e %e %e\n", a, tmp.x, tmp.y, tmp.z);
                        tmp = vec3_muls(onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Bhelix[a],u),u), Bhelix[a]));
                        printf("B part %d is : %e %e %e\n", a, tmp.x, tmp.y, tmp.z);
                        
                        tmp = vec3_add (
                                            vec3_muls( c * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Ahelix[a],u),u), Ahelix[a])),
                                            vec3_muls( s * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Bhelix[a],u),u), Bhelix[a]))
                                       );
                                       
                        printf("CDip %d to be added : %e %e %e\n", a, tmp.x, tmp.y, tmp.z);
                        tmp =  vec3_sub(
                                            vec3_muls( s * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Ahelix[a],u),u), Ahelix[a])),
                                            vec3_muls( c * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Bhelix[a],u),u), Bhelix[a]))
                                        );
                        printf("SDip %d to be added : %e %e %e\n", a, tmp.x, tmp.y, tmp.z);
#endif
                        // sum all data
                        #pragma omp critical(dipolar)
                        {
                            // Dipolar
                            CDip[a] = vec3_add(
                                            CDip[a],
                                            vec3_add(
                                                vec3_muls( c * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Ahelix[a],u),u), Ahelix[a])),
                                                vec3_muls( s * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Bhelix[a],u),u), Bhelix[a]))
                                            )
                                        );
                            
                            SDip[a] = vec3_add(
                                            SDip[a],
                                            vec3_sub(
                                                vec3_muls( s * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Ahelix[a],u),u), Ahelix[a])),
                                                vec3_muls( c * onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(Bhelix[a],u),u), Bhelix[a]))
                                            )
                                        );
                        }
                        #pragma omp critical(lorentz)
                        {                                        
                            // Lorentz
                            CLor[a] = vec3_add(
                                            CLor[a],
                                            vec3_add(
                                                vec3_muls( c , Ahelix[a]),
                                                vec3_muls( s , Bhelix[a])
                                            )
                                        );
                            SLor[a] = vec3_add(
                                            SLor[a],
                                            vec3_sub(
                                                vec3_muls( s ,Ahelix[a]),
                                                vec3_muls( c ,Bhelix[a])
                                            )
                                        );
                        }
						// Contact
						if (n < cont_radius) {
                            #pragma omp critical(contact)
                            {
                                pile_add_element(&CCont, pow(n,CONT_SCALING_POWER), 
							  								vec3_add(
							  									vec3_muls( stagmom[a] * c , Ahelix[a]),
							  									vec3_muls( stagmom[a] * s , Bhelix[a])
							  								)							
                                                );
                                pile_add_element(&SCont, pow(n,CONT_SCALING_POWER), 
							  								vec3_sub(
							  									vec3_muls( stagmom[a]* s , Ahelix[a]),
							  									vec3_muls( stagmom[a]* c , Bhelix[a])
							  								)
                                                );
						    }
                        }
#ifdef _DEBUG                      
                        printf("CDip %d is now : %e %e %e\n", a, CDip[a].x, CDip[a].y, CDip[a].z);
                        printf("SDip %d is now : %e %e %e\n", a, SDip[a].x, SDip[a].y, SDip[a].z);
#endif
                    }                    
                }
            }
        }
    }
}
    
    double angle=0;
    struct vec3 BDip;
    struct vec3 BLor;
    struct vec3 BCont;

    // for contact field evaluation
    struct vec3 CBCont = vec3_zero();
    struct vec3 SBCont = vec3_zero();
    int NofM = 0; // Number of moments considered
    double SumOfWeights = 0;    
    
#pragma omp parallel sections private(angn,angle,BDip,BLor,BCont,i) firstprivate(CBCont,SBCont,NofM,SumOfWeights)
{
    // first portion, dipolar fields and Lorentz
    #pragma omp section
    {
        for (angn = 0; angn < in_nangles; ++angn)
        {
            angle = 2*M_PI*((float) angn / (float) in_nangles);
            
            //  === Dipolar Field ===
            BDip = vec3_zero();
            // loop over atoms
            for (a = 0; a < size; ++a)
            {
                BDip =  vec3_add(BDip,
                                vec3_muls(
                                    stagmom[a],
                                    vec3_sub(
                                        vec3_muls(cos(angle) , CDip[a]),
                                        vec3_muls(sin(angle) , SDip[a])
                                    )
                                )
                            );
            }
            
            BDip = vec3_muls(0.9274009, BDip); // to tesla units
            out_field_dip[3*angn+0] = BDip.x;
            out_field_dip[3*angn+1] = BDip.y;
            out_field_dip[3*angn+2] = BDip.z;        
            
            
            //  === Lorentz Field ===
            BLor = vec3_zero();
            // loop over atoms
            for (a = 0; a < size; ++a)
            {
                BLor =  vec3_add(BLor,
                                vec3_muls(
                                    stagmom[a],
                                    vec3_sub(
                                        vec3_muls(cos(angle) , CLor[a]),
                                        vec3_muls(sin(angle) , SLor[a])
                                    )
                                )
                            );
            }
            
            BLor = vec3_muls(0.33333333333*11.654064, vec3_muls(3./(4.*M_PI*pow(radius,3)),BLor));
            out_field_lor[3*angn+0] = BLor.x;
            out_field_lor[3*angn+1] = BLor.y;
            out_field_lor[3*angn+2] = BLor.z;        
            
        }
    }


    // second portion, contact fields
    #pragma omp section
    {
        //  === Contact Field ===
        BCont = vec3_zero();
        
        for (i=0; i < nnn_for_cont; i++) {
            if ((CCont.ranks[i] >= 0.0) && (fabs(CCont.ranks[i] - SCont.ranks[i])<EPS)) {
                CBCont = vec3_add(CBCont, vec3_muls(1./CCont.ranks[i],
                                                    CCont.elements[i])
                                    );
                SBCont = vec3_add(SBCont, vec3_muls(1./SCont.ranks[i],
                                                    SCont.elements[i])
                                    );
                SumOfWeights += 1./CCont.ranks[i];
                NofM++;
            } else {
                printf("Something VERY odd ! ranks 1: %e ranks 2: %e \n", CCont.ranks[i] , SCont.ranks[i] );
            }
        }
            
        // (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
        //   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
        //   ≈ 7.769376 T⋅Å^3
        
        for (angn = 0; angn < in_nangles; ++angn) {
            
            angle = 2*M_PI*((float) angn / (float) in_nangles);
            
            if (NofM >0) {
                BCont = vec3_muls((1./SumOfWeights) * 7.769376 , 
                                    vec3_sub(
                                        vec3_muls(cos(angle) , CBCont),
                                        vec3_muls(sin(angle) , SBCont)
                                    )
                                );
            } // otherwise is zero anyway!
            
            out_field_cont[3*angn+0] = BCont.x;
            out_field_cont[3*angn+1] = BCont.y;
            out_field_cont[3*angn+2] = BCont.z;        
        }
    }
} // end of omp parallel sections

    // free stuff used for contact field
    pile_free(&CCont);
    pile_free(&SCont);

}


