
#include <assert.h>

#include "params.h"
#include "steps.h"
#include "uint.h"
#include "fp.h"

#include "mont.h"
#include <stdio.h>
#include <string.h>




void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A24)
{
   fp tmp0, tmp1, tmp2;        //requires precomputation of A24=(A+2C:4C)

    fp_add3(&tmp0, &P->x, &P->z);
    fp_sub3(&tmp1, &P->x, &P->z);
    fp_sq2(&R->x, &tmp0);
    fp_sub3(&tmp2, &Q->x, &Q->z);
    fp_add3(&S->x, &Q->x, &Q->z);
    fp_mul2(&tmp0, &tmp2);
    fp_sq2(&R->z, &tmp1);
    fp_mul2(&tmp1, &S->x);
    fp_sub3(&tmp2, &R->x, &R->z);
    fp_mul2(&R->z, &A24->z);
    fp_mul2(&R->x, &R->z);
    fp_mul3(&S->x, &A24->x, &tmp2);
    fp_sub3(&S->z, &tmp0, &tmp1);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &tmp0, &tmp1);
    fp_mul2(&R->z, &tmp2);
    fp_sq1(&S->z);
    fp_sq1(&S->x);
    fp_mul2(&S->z, &PQ->x);
    fp_mul2(&S->x, &PQ->z);

}

// A->x = (C-D)^2
// A->z 4CD
void xDBLADD_Huff(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{

    fp t0, t1, t2, t3, t4, t5, t6;
    
    fp_add3(&t0, &P->x, &P->z);                         // t0 = Px+Pz
    fp_sub3(&t1, &P->x, &P->z);                         // t1 = Px-Pz
    fp_add3(&t2, &Q->x, &Q->z);                         // t2 = Qx+Qz
    fp_sub3(&t3, &Q->x, &Q->z);                         // t3 = Qx-Qz
    fp_sq2(&t4, &t0);                                   // t4 = (Px+Pz)^2
    fp_sq2(&t5, &t1);                                   // t5 = (Px-Pz)^2
    fp_mul3(&t2, &t1, &t2);                        // t2 = (Px-Pz)(Qx+Qz) = PxQx+PxQz-PzQx-PzQz
    fp_mul3(&t3, &t3, &t0);                        // t3 = (Px+Pz)(Qx-Qz) = PxQx-PxQz+PzQx-PzQz
    fp_add3(&t0, &t2, &t3);                             // t0 = 2(PxQx-PzQz) => Sx
    fp_sub3(&t1, &t2, &t3);                             // t1 = 2(PxQz-PzQx) => Sz
    fp_mul3(&t6, &t4, &A->z);                        // t6 = 4CD(Px+Pz)^2
    fp_mul3(&R->x, &t6, &t5);                      // Rx = 4CD(Px+Pz)^2(Px-Pz)^2
    fp_sub3(&t5, &t4, &t5);                             // t5 = 4PxPz
    fp_mul3(&R->z, &A->x, &t5);                           // X2 = (C-D)^24PxPz
    fp_add3(&R->z, &t6, &R->z);                           // X2 =  4CD(Px+Pz)^2+(C-D)^24PxPz
    fp_mul3(&R->z, &R->z, &t5);      
    fp_sq2(&t0, &t0);                               // t0 = (PxQx-PzQz)^2 => Sx
    fp_sq2(&t1, &t1);
    fp_mul3(&S->z, &PQ->x, &t1);
    fp_mul3(&S->x, &PQ->z, &t0);     
}



void xDBL(proj *Q, proj const *A, proj const *P)
{

    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);

}


// A->x = (C-D)^2
// A->z = 4CD
void xDBL_Huff(proj *Q, proj const *A, proj const *P)
{
    fp t0, t1, t2;
    
    fp_add3(&t0, &P->x, &P->z);                         // t0 = X1+Z1
    fp_sq2(&t0, &t0);                                   // t0 = (X1+Z1)^2
    fp_sub3(&t1, &P->x, &P->z);                         // t1 = X1-Z1
    fp_sq2(&t1, &t1);                                   // t1 = (X1-Z1)^2
    fp_mul3(&t2, &A->z, &t0);                           // t2 = 4CD(X1+Z1)^2
    fp_mul3(&Q->x, &t2, &t1);                           // X2 = 4CD(X1+Z1)^2(X1-Z1)^2

    fp_sub3(&t0, &t0, &t1);                             // t0 = 4X1Z1
    fp_mul3(&Q->z, &t0, &A->x);                         // Z2 = 4X1Z1(C-D)^2
    fp_add3(&Q->z, &Q->z, &t2);                         // Z2 = 4X1Z1(C-D)^2^2+4CD(X1+Z1)^2
    fp_mul3(&Q->z, &Q->z, &t0);

}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
}


void xADD_Huff(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;

    fp_add3(&a, &P->x, &P->z);                          // a = X+Z
    fp_sub3(&b, &P->x, &P->z);                          // b = X-Z
    fp_add3(&c, &Q->x, &Q->z);                          // c = x+z
    fp_sub3(&d, &Q->x, &Q->z);                          // d = x-z
    fp_mul2(&a, &d);                                    // a = (X+Z)(x-z)=Xx-Xz+xZ-Zz
    fp_mul2(&b, &c);                                    // b = (X-Z)(x+z)=Xx+Xz-Zx-Zz
    fp_add3(&c, &a, &b);                                // c = Xx-Zz //x
    fp_sub3(&d, &a, &b);                                // d = Xz-xZ
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
}


/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    proj A24;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    fp_add3(&A24.x, &A->z, &A->z);    //precomputation of A24=(A+2C:4C)
    fp_add3(&A24.z, &A24.x, &A24.x);
    fp_add2(&A24.x, &A->x);

    unsigned long i = 64 * LIMBS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, &A24);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}


void xMUL_Huff(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    unsigned long i = 64 * LIMBS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD_Huff(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}


void exp_by_squaring_(fp* x, fp* y, uint64_t exp)
{
	fp result1, result2;
	fp_set(&result1, 1);
	fp_set(&result2, 1);

    while (exp)
    {
        if (exp & 1){
          fp_mul2(&result1, x);
          fp_mul2(&result2, y);
	}
	
        fp_sq1(x);
	fp_sq1(y);
        exp >>= 1;
    }

    fp_cswap(&result1, x, 1);
    fp_cswap(&result2, y, 1);

}



void biquad_both_opt(fp *out,fp *outinv, fp *coeff, fp *coeffQ, const proj *A)
{
  fp t0, t1, t2;
  // coeffQ[0]=Qx+Qz
  // coeffQ[1]=Qx-Qz
  // coeffQ[2]= 4QxQz
 // coeffQ[3]=2(Qx^2+Qz^2)

  fp_mul3(&t0, &coeff[3], &coeffQ[1]); // t0 = (Px+Pz)*(Qx-Qz)
  fp_mul3(&t1, &coeff[4], &coeffQ[0]); // t1 = (Px-Pz)*(Qx+Qz)

  fp_add3(&out[0], &t0, &t1); // t2 = 2(PxQx-PzQz)
  fp_sq1(&out[0]); //4(PxQx-PzQz)^2
  fp_mul2(&out[0],&A->z); //4C(PxQx-PzQz)^2
  

  fp_sub3(&out[2], &t1, &t0); // 2(PxQz-PzQx)
  fp_sq1(&out[2]); // 4(PxQz-PzQx)^2
  fp_mul2(&out[2],&A->z);

  fp_mul3(&t0, &coeff[1], &coeffQ[3]); //t0= 4CPxPz*2(Qx^2+Qz^2)
  fp_mul3(&t1, &coeff[0], &coeffQ[2]); // t1 = 4QxQz*2C(Px^2+Pz^2)
  fp_mul3(&t2, &coeff[2], &coeffQ[2]); // t2 = 4APxPz* 4QxQz

  fp_add3(&t0,&t0,&t1);
  fp_add3(&t0,&t0,&t2);
  fp_sub3(&out[1], &fp_0, &t0);

 
  outinv[1] = out[1];
  outinv[2] = out[0];
  outinv[0] = out[2];
}

// For isogeny evaluation
void biquad_both_Huff_Ed(fp *out,fp *outinv, fp *coeff, fp *coeffQ, const proj *A)
{

  // coeffQ[0]=Qx+Qz
  // coeffQ[1]=Qx-Qz
  // coeffQ[2]= 4QxQz
  // coeffQ[3]=2(Qx^2+Qz^2)


  fp t0, t1, t2;
  
  fp_mul3(&t0, &coeff[3], &coeffQ[1]); // t0 = (Px+Pz)*(Qx-Qz)
  fp_mul3(&t1, &coeff[4], &coeffQ[0]); // t1 = (Px-Pz)*(Qx+Qz)

  fp_add3(&out[0], &t0, &t1); // t2 = 2(PxQx-PzQz)
  fp_sq1(&out[0]); //4(PxQx-PzQz)^2
  fp_mul2(&out[0],&A->z); //4C(PxQx-PzQz)^2
  

  fp_sub3(&out[2], &t1, &t0); // 2(PxQz-PzQx)
  fp_sq1(&out[2]); // 4(PxQz-PzQx)^2
  fp_mul2(&out[2],&A->z);

  fp_mul3(&t0, &coeff[1], &coeffQ[3]); //t0= 4CPxPz*2(Qx^2+Qz^2)
  fp_mul3(&t1, &coeff[0], &coeffQ[2]); // t1 = 4QxQz*2C(Px^2+Pz^2)
  fp_mul3(&t2, &coeff[2], &coeffQ[2]); // t2 = 4APxPz* 4QxQz

  fp_add3(&t0,&t0,&t1);
  fp_add3(&t0,&t0,&t2);
  fp_sub3(&out[1], &fp_0, &t0);

 
  outinv[1] = out[1];
  outinv[2] = out[0];
  outinv[0] = out[2];


}

void biquad_both_Huff_opt(fp *out,fp *outinv, fp *out2, fp *outinv2, const proj *P, const proj *A, fp *coeff)
{
  // coeff[0]=c(Qx+Qz)
  // coeff[1]=c(Qx-Qz)
  // coeff[2] =c(Qx2+Qz2)
  // coeff[3] =c(Qx2-Qz2)
  // coeff[4] =4CQxQz
  // coeff[5] =2A+4C
   // coeff[6] =4CQx2Qz2
 // coeff[7]=2c(Qx^2+Qz^2)
// coeff[8]=2c(Qx2^2+Qz2^2)



  fp Pplus; fp_add3(&Pplus, &P->x, &P->z);
  fp Pminus; fp_sub3(&Pminus, &P->x, &P->z);

  fp t0, t1, t2, t3, t4, t5;
  // for evaluation
  fp_mul3(&t0, &Pplus, &coeff[1]); // t0 = (Px+Pz0 )c(Qx-Qz)=c(PxQx-PxQz+PzQx-PzQz)
  fp_mul3(&t1, &Pminus, &coeff[0]); // t1 = (Px-Pz)c(Qx+Qz)=c(PxQx+PxQz-PzQx-PzQz)

  fp_add3(&out[0], &t0, &t1); //out0 = 2C(PxQx-PzQz)
  fp_sq1(&out[0]); //out0 = 4C^2(PxQx-PzQz)^2

  fp_sub3(&out[2], &t1, &t0); // 2c(PxQz-PzQx)
  fp_sq1(&out[2]); // 4c^2(PxQz-PzQx)^2


  // for coefficient
  fp_mul3(&t0, &Pplus, &coeff[3]); // t0 = (Px+Pz)(Qx-Qz)=PxQx-PxQz+PzQx-PzQz
  fp_mul3(&t1, &Pminus, &coeff[2]); // t1 = (Px-Pz)(Qx+Qz)=PxQx+PxQz-PzQx-PzQz

  fp_add3(&out2[0], &t0, &t1); //out0 = 2(PxQx-PzQz)
  fp_sq1(&out2[0]);

  fp_sub3(&out2[2], &t1, &t0); // 2(PxQz-PzQx)
  fp_sq1(&out2[2]); // 4(PxQz-PzQx)^2



  // computing out[1], out2[1]
    fp_sq1(&Pplus);  //(Px+Pz)^2
    fp_sq1(&Pminus);  //(Px-Pz)^2
    fp_add3(&t0,&Pplus,&Pminus); // t0 = 2(Px^2+Pz^2)
    fp_sub3(&t5,&Pplus,&Pminus); // t5 = 4PxPz
     fp_mul3(&t1, &t5, &A->z); // t1 = 4cPxPz
    fp_mul3(&t0, &t0, &A->z); // t0 = 2c(Px^2+Pz^2)


  fp_mul3(&t2, &t1, &coeff[7]); //t2= 4cPxPz*2C(Qx^2+Qz^2)
  fp_mul3(&t3, &t0, &coeff[4]); // t3 = 4cQxQz*2c(Px^2+Pz^2)

  fp_mul3(&t4, &t5, &coeff[4]); // t4 = 4PxPz* 4cQxQz
  fp_mul3(&t4, &t4, &coeff[5]); // t4 = 4PxPz* 4cQxQz *(c+2d)
  fp_add3(&t4, &t4, &t4);

  fp_add3(&t2,&t2,&t3);
  fp_add3(&t2,&t2,&t4);
  fp_sub3(&out[1],&fp_0, &t2);

  //compute out2

  fp_mul3(&t2, &t1, &coeff[8]); //t2= 4PcxPz*2C(Qx^2+Qz^2)
  fp_mul3(&t3, &t0, &coeff[6]); // t3 = 4cQxQz*2(Px^2+Pz^2)
  fp_mul3(&t4, &t5, &coeff[6]); // t4 = 4PxPz* 4cQxQz
  fp_mul3(&t4, &t4, &coeff[5]); // t4 = 4(2C+4D)PxPz* 4cQxQz
  fp_add3(&t4, &t4, &t4);

  fp_add3(&t2,&t2,&t3);
  fp_add3(&t2,&t2,&t4);
  fp_sub3(&out2[1], &fp_0, &t2);


  outinv[1] = out[1];
  outinv[2] = out[0];
  outinv[0] = out[2];


  
  outinv2[1] = out2[1];
  outinv2[2] = out2[0];
  outinv2[0] = out2[2];

}

void biquad_pm1(fp *outplus,fp *outminus,const proj *P,const proj *A)
{
  fp Pplus; fp_add3(&Pplus,&P->x,&P->z); //plus=x+z
  fp Pminus; fp_sub3(&Pminus,&P->x,&P->z); //pminus=x-z
  fp_sq1(&Pplus); //pplus=(x+z)^2
  fp_sq1(&Pminus); //pminus=(x-z)^2

  fp_mul3(&outplus[0],&Pminus,&A->z); // outp =C(x-z)^2
  outplus[2] = outplus[0];
  fp_mul3(&outminus[0],&Pplus,&A->z); // outm=C(x+z)^2
  outminus[2] = outminus[0];

  fp u;
  fp_sub3(&u,&Pminus,&Pplus); // u=-4xz
  fp_mul2(&u,&A->x); // u=-4xzA

  fp t;
  fp_add3(&t,&outminus[0],&outminus[0]); // t=2C(x^2+2xz+z^2)
  fp_sub3(&outplus[1],&u,&t); //outplus1=-4Axz*2C(x^2+2xz+z^2)

  fp_add3(&t,&outplus[0],&outplus[0]);
  fp_sub3(&outminus[1],&t,&u);
}
void biquad_pm1_opt(fp *coeff, fp *outplus,fp *outminus,const proj *P,const proj *A)
{
  fp Pplus; fp_add3(&coeff[3],&P->x,&P->z); //plus=x+z
  fp Pminus; fp_sub3(&coeff[4],&P->x,&P->z); //pminus=x-z
  fp_sq2(&Pplus, &coeff[3]); //pplus=(x+z)^2
  fp_sq2(&Pminus, &coeff[4]); //pminus=(x-z)^2

  fp_mul3(&outplus[0],&Pminus,&A->z); // outp =C(x-z)^2
  outplus[2] = outplus[0];
  fp_mul3(&outminus[0],&Pplus,&A->z); // outm=C(x+z)^2
  outminus[2] = outminus[0];

  fp_add3(&coeff[0],&outplus[0],&outminus[0] ); // coeff[0]=2C(X^2+Z^2
  fp_sub3(&coeff[1], &outminus[0],&outplus[0] ); //4CXZ

  fp u;
  fp_sub3(&coeff[2], &Pplus, &Pminus); // coeff[2]=4XZ
  fp_mul2(&coeff[2],&A->x); // coeff[2=4xzA
  fp_sub3(&u, &fp_0,  &coeff[2]);

 // fp_sub3(&u,&Pminus,&Pplus); // u=-4xz
//  fp_mul2(&u,&A->x); // u=-4xzA

  fp t;
  fp_add3(&t,&outminus[0],&outminus[0]); // t=2C(x^2+2xz+z^2)
  fp_sub3(&outplus[1],&u,&t); //outplus1=-4Axz-2C(x^2+2xz+z^2)

  fp_add3(&t,&outplus[0],&outplus[0]);
  fp_sub3(&outminus[1],&t,&u);
}


void biquad_pm1_opt_Huff(fp *coeff, fp *outplus,fp *outminus,const proj *P,const proj *A)
{
  fp Pplus; fp_add3(&coeff[3],&P->x,&P->z); //plus=x+z
  fp Pminus; fp_sub3(&coeff[4],&P->x,&P->z); //pminus=x-z
  fp_sq2(&Pplus, &coeff[3]); //pplus=(x+z)^2
  fp_sq2(&Pminus, &coeff[4]); //pminus=(x-z)^2

  fp_mul3(&outplus[0],&Pminus,&A->z); // outp =C(x-z)^2
  outplus[2] = outplus[0];
  fp_mul3(&outminus[0],&Pplus,&A->z); // outm=C(x+z)^2
  outminus[2] = outminus[0];

  fp_add3(&coeff[0],&outplus[0],&outminus[0] ); // coeff[0]=2C(X^2+Z^2)
  fp_sub3(&coeff[1], &outminus[0],&outplus[0] ); //4CXZ

  fp u, t;


  fp_sub3(&coeff[2], &Pplus, &Pminus); // coeff[2]=4XZ
  fp_mul2(&coeff[2],&A->x); // coeff[2=4xzA
  fp_sub3(&u, &fp_0,  &coeff[2]);


  fp_add3(&t,&outminus[0],&outminus[0]); // t=2C(x^2+2xz+z^2)
  fp_sub3(&outplus[1],&u,&t); //outplus1=-4Axz-2C(x^2+2xz+z^2)

  fp_add3(&t,&outplus[0],&outplus[0]);
  fp_sub3(&outminus[1],&t,&u);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
void xISOG(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1, tmp2, Psum, Pdif;
    fp T[4] = {K->z, K->x, K->x, K->z};
    proj Q;

    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);

    fp_sub3(&tmp1, &K->x, &K->z);
    fp_add3(&tmp0, &K->x, &K->z);
    
    fp_mul3(&tmp1, &tmp1, &Psum);
    fp_mul3(&tmp0, &tmp0, &Pdif);
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);


    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);


      	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
    	  fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);

    	fp_mul2(&tmp1, &Psum);
    	fp_mul2(&tmp0, &Pdif);
    	fp_add3(&tmp2, &tmp0, &tmp1);
	    fp_mul2(&Q.x, &tmp2);
    	fp_sub3(&tmp2, &tmp0, &tmp1);
	    fp_mul2(&Q.z, &tmp2);

    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);
}

// Isogeny eval : Montgomery curve
// coefficient : Edwards
void xISOG_hy(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);


    fp tmp0, tmp1, tmp2, Psum, Pdif;
    proj Q, Aed, prod;

    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients // C'=2C
    fp_add3(&Aed.x, &A->x, &Aed.z); // D'=A+2C
    fp_sub3(&Aed.z, &A->x, &Aed.z); //C' = A+2C-2C A
   
    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);

    fp_sub3(&prod.x, &K->x, &K->z);
    fp_add3(&prod.z, &K->x, &K->z);
    
    fp_mul3(&tmp1, &prod.x, &Psum);
    fp_mul3(&tmp0, &prod.z, &Pdif);
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
           xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
    	fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
	fp_mul2(&prod.x, &tmp1);
        fp_mul2(&prod.z, &tmp0);
    	fp_mul2(&tmp1, &Psum);
    	fp_mul2(&tmp0, &Pdif);
    	fp_add3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.x, &tmp2);
    	fp_sub3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.z, &tmp2);

    }


    // point evaluation
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);

}

// Huff isogney : standard Velu
void xISOG_Huff(proj *A, proj *Aw, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1, tmp2, tmp3, tmp4, t0, t1;
    fp Psum, Pdif;
    fp T[2]={fp_1, fp_1};
    proj Q;


    fp_add3(&Psum, &P->x, &P->z);       // Psum = (X+Z)
    fp_sub3(&Pdif, &P->x, &P->z);       // Pdif = (X-Z)
    fp_add3(&tmp2, &K->x, &K->z);       // tmp2 = Xi+Zi
    fp_sub3(&tmp3, &K->x, &K->z);       // tmp3 = Xi-Zi
    fp_mul3(&tmp0, &tmp2, &Pdif);       // tmp0 = (Xi+Zi)(X-Z)=XXi-ZXi+XZi-ZZi
    fp_mul3(&tmp1, &tmp3, &Psum);       // tmp1 = (Xi-Zi)(X+Z)=XXi+ZXi-XZi-ZZi
    fp_add3(&Q.x, &tmp0, &tmp1);        // Qx = XXi-ZZi
    fp_sub3(&Q.z, &tmp1, &tmp0);

    fp_add3(&t0, &A->x, &A->z);        //t0 = (C+D)
    fp_sub3(&t1, &A->x, &A->z);        //t1 = (C-D)


    fp_mul3(&tmp0, &t0, &tmp2);        // tmp0= (C+D)(X+Z)=CX+CZ+DX+DZ
    fp_mul3(&tmp1, &t1, &tmp3);        // tmp1= (C-D)(X-Z)=CX-CZ-DX+DZ
    fp_add3(&T[0], &tmp0, &tmp1);      // tmp3 = 2(CX+DZ)
    //fp_mul3(&T[1], &T[1], &tmp3);

    fp_sub3(&T[1], &tmp1, &tmp0);
    //fp_mul3(&T[0], &T[0], &tmp3);



    proj M[3] = {*K};
    xDBL_Huff(&M[1], Aw, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD_Huff(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

        fp_add3(&tmp2, &M[i % 3].x, &M[i % 3].z);  // tmp2 = x+z
        fp_sub3(&tmp3, &M[i % 3].x, &M[i % 3].z);  // tmp3 = x-z

        fp_mul3(&tmp0, &t0, &tmp2);        // tmp0= (C+D)(X+Z)=CX+CZ+DX+DZ
        fp_mul3(&tmp1, &t1, &tmp3);        // tmp0= (C-D)(X-Z)=CX-CZ-DX+DZ
        fp_add3(&tmp4, &tmp0, &tmp1);     // tmp3 = CX+DZ   
        fp_mul3(&T[0], &T[0], &tmp4);

        fp_sub3(&tmp4, &tmp1, &tmp0);       // tmp3 = CZ+DX
        fp_mul3(&T[1], &T[1], &tmp4);


     //   fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);       // tmp0 = Xi+Zi
    //    fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);       // tmp1 = Xi-Zi
        fp_mul3(&tmp0, &tmp2, &Pdif);       // tmp0 = (Xi+Zi)(X-Z)=XXi-ZXi+XZi-ZZi
        fp_mul3(&tmp1, &tmp3, &Psum);       // tmp1 = (Xi-Zi)(X+Z)=XXi+ZXi-XZi-ZZi
        fp_add3(&tmp3, &tmp0, &tmp1);        // Qx = XXi-ZZi
        fp_mul2(&Q.x, &tmp3);

        fp_sub3(&tmp3, &tmp1, &tmp0);
        fp_mul2(&Q.z,  &tmp3);
    }
    fp_sq1(&T[0]);
    fp_sq1(&T[1]);
    fp_mul3(&A->x, &A->x, &T[0]);
    fp_mul3(&A->z, &A->z, &T[1]);

    
    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

                        // transform 
    fp_add3(&tmp0, &A->x, &A->z);       // t0 =
    fp_sub3(&tmp1, &A->x, &A->z);
    fp_sq2(&Aw->x, &tmp1);             // Aw.x=(C-D)^2
    fp_sq1(&tmp0);
    fp_sub3(&Aw->z, &tmp0, &Aw->x);   // Aw.z=4AB


    
}

// Huff isogeny, standard Velu
// Isogeny evaluation : Huff curve
// Coefficients : Edwards curve
void xISOG_Huff_Edwards(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1, tmp2;
    fp Psum, Pdif;
  
    proj Q, prod;


    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);

    fp_sub3(&prod.x, &K->x, &K->z); //prodx=(Xi-Zi)
    fp_add3(&prod.z, &K->x, &K->z); //prodz=(Xi+Zi)
    
    fp_mul3(&tmp1, &prod.x, &Psum); // tmp1=(Xi-Zi)(X+Z) (tmp1)
    fp_mul3(&tmp0, &prod.z, &Pdif); // tmp0=(Xi+Zi)(X-Z) (tmp0)
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);

    proj Aed; //Huff coefficient to Edwards
    
    fp_add3(&Aed.z, &A->x, &A->x);
    fp_add3(&Aed.z, &Aed.z, &Aed.z); // Aedz = 4Ax

    fp_add3(&Aed.x, &A->z, &A->z);
    fp_add3(&Aed.x, &Aed.x, &Aed.x);

    fp_add3(&Aed.x, &Aed.x, &Aed.z);


    proj M[3] = {*K};
    xDBL_Huff(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD_Huff(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
    	fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
	fp_mul2(&prod.x, &tmp1);
        fp_mul2(&prod.z, &tmp0);
    	fp_mul2(&tmp1, &Psum);
    	fp_mul2(&tmp0, &Pdif);
    	fp_add3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.x, &tmp2);
    	fp_sub3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.z, &tmp2);
    }


    
    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    // Edwards to Huff
    fp_add3(&A->x, &Aed.z, &Aed.z); //
    fp_add3(&A->x, &A->x, &A->x); // 

    fp_add3(&A->z, &Aed.x, &Aed.x); //
    fp_add3(&A->z, &A->z, &A->z); //
    fp_sub3(&A->z, &A->z, &A->x);


    
}

// Montgomery isogeny - square-root Velu
void xISOG_sqrt(proj *A, proj *P, proj const *K, long long k)
{
   assert (k >= 3);
    assert (k % 2 == 1);


    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);
 

    fp tmp0, tmp1, tmp2;

    proj Aed;
    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp_add3(&Aed.x, &A->x, &Aed.z);
    fp_sub3(&Aed.z, &A->x, &Aed.z);
   
    fp Psum, Pdif;
    fp coeffQ[4];
    
    fp_add3(&Psum, &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&Pdif, &P->x, &P->z);   // Pdif= x-z
    fp_add3(&coeffQ[0], &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&coeffQ[1], &P->x, &P->z);   // Pdif= x-z
    fp_sq2(&tmp0, &coeffQ[0]); //(x+z)^2
    fp_sq2(&tmp1, &coeffQ[1]); //(x-z)^2
    fp_add3(&coeffQ[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp_sub3(&coeffQ[2], &tmp0, &tmp1); // 



    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL(&M[2], A, K); Minit[2] = 1;

    
      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    

    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        fp_sub3(&TI[2*i],&fp_0, &M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp T1[3*gs];
      fp Tminus1[3*gs];
      fp coeff[5];

      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
        biquad_pm1_opt(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],A);
        biquad_both_opt(TP+3*j,TPinv+3*j, coeff, coeffQ, A);

      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2_selfreciprocal(T1,gs); // t1=1=c
      poly_multiprod2_selfreciprocal(Tminus1,gs); //tminus1=-1=1/c

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
        fp_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
        fp_mul2(&Abatch.x,&tmp1);
        fp_mul2(&Abatch.z,&tmp0);
        fp_mul2(&tmp1, &Psum);
        fp_mul2(&tmp0, &Pdif);
        fp_add3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.x, &tmp2);
        fp_sub3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.z, &tmp2);
      }
    

    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &Abatch.x);
    fp_mul2(&Aed.x, &Abatch.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);
}


// Huff isogeny - square-root Velu
void xISOG_Huff_sqrt_opt(proj *A, proj *Aw, proj *P, proj const *K, long long k)
{
       assert (k >= 3);
    assert (k % 2 == 1);


    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);


    fp tmp0, tmp1, tmp2, a0;
    proj Am;
    fp coeff[9];

  fp_add3(&a0, &Aw->x, &Aw->x); 
    fp_add3(&coeff[5], &a0, &Aw->z); //2C+D


   fp_sub3(&Am.x, &fp_0, &A->x);
    memcpy(&Am.z, &A->z, sizeof(fp));


fp Psum, Pdif;
    fp Csum, Cdif;
          fp t0, t1, t2, t3;
    fp_add3(&Psum, &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&Pdif, &P->x, &P->z);   // Pdif= x-z
    fp_add3(&Csum, &A->x, &A->z);   // Csum C+D
    fp_sub3(&Cdif, &A->x, &A->z);   // Cdif= C-D

    fp_mul3(&coeff[0], &Psum, &Aw->z);
    fp_mul3(&coeff[1], &Pdif, &Aw->z);
    fp_add3(&t0, &Am.x, &Am.z);
    fp_sub3(&t1, &Am.x, &Am.z);


    fp_mul3(&coeff[2], &t0, &Aw->z);
    fp_mul3(&coeff[3], &t1, &Aw->z);
    fp_sq2(&t2, &Psum); //(Qx+Qz)^2
    fp_sq2(&t3, &Pdif); //(Qx-Qz)^2
    fp_sub3(&coeff[4], &t2, &t3);   // 4QxQz
    fp_mul3(&coeff[4], &coeff[4], &Aw->z);
    fp_add3(&coeff[7], &t2, &t3); //C^2
    fp_mul3(&coeff[7], &coeff[7], &Aw->z);

    fp_sq2(&t0, &t0); //(Qx+Qz)^2
    fp_sq2(&t1, &t1); //(Qx-Qz)^2
    fp_sub3(&coeff[6], &t0, &t1);   // 4QxQz
    fp_mul3(&coeff[6], &coeff[6], &Aw->z);
    fp_add3(&coeff[8], &t0, &t1);
    fp_mul3(&coeff[8], &coeff[8], &Aw->z);
    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL_Huff(&M[2], Aw, K); Minit[2] = 1;

    
      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD_Huff(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL_Huff(&M[s],Aw,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL_Huff(&M[s],Aw,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD_Huff(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD_Huff(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    

    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

   
      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        fp_sub3(&TI[2*i],&fp_0,&M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp TC[3*gs]; //-c
      fp TCinv[3*gs];

      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);

        biquad_both_Huff_opt(TP+3*j,TPinv+3*j, TC+3*j,TCinv+3*j, &M[2*bs*(2*j+1)], Aw, coeff);


      
      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2(TC,gs);
      poly_multiprod2(TCinv,gs);
    //  poly_multiprod2_selfreciprocal(TC,gs); // t1=1=c
   //   poly_multiprod2_selfreciprocal(TCinv,gs); //tminus1=-1=1/c

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,TC,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,TCinv,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z); //tmp1 = xi-zi
        fp_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z); // tmp 0 = xi+zi
        fp_mul3(&t0, &tmp0, &Csum);   // t0 =(C+D)(xi+zi)
        fp_mul3(&t1, &tmp1, &Cdif);
        fp_add3(&t2, &t0, &t1); //Ax
        fp_mul2(&Abatch.z, &t2);
        fp_sub3(&t2, &t0, &t1);
        fp_mul2(&Abatch.x, &t2);
        //fp_mul2(&Abatch.x,&tmp1); //*(x-z)
        //fp_mul2(&Abatch.z,&tmp0); // (x+z)

        fp_mul2(&tmp1, &Psum); // tmp1 (xi-zi)(x+z) xxi+xiz-zix-zzi
        fp_mul2(&tmp0, &Pdif); // tmp0 (xi+zi)(x-z) xxi-zxi
        fp_add3(&tmp2, &tmp0, &tmp1); //tmp2 = xxi-zzi
        fp_mul2(&Qbatch.x, &tmp2);
        fp_sub3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.z, &tmp2);
      }
    

    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);


    //compute prod.x^8, prod.z^8
    fp_sq1(&Abatch.x); //tc
    fp_sq1(&Abatch.z);

    //compute image curve parameters
    fp_mul2(&A->z, &Abatch.x);
    fp_mul2(&A->x, &Abatch.z);

                            // transform 
    fp_add3(&tmp0, &A->x, &A->z);       // t0 =(C+D)
    fp_sub3(&tmp1, &A->x, &A->z);        // t1 = (C-D)
    fp_sq2(&Aw->x, &tmp1);             // Aw.x=(C-D)^2
    fp_sq1(&tmp0);                      // (C+D)^2
    fp_sub3(&Aw->z, &tmp0, &Aw->x);   // Aw.z=4AB

}

// Huff isogeny - square-root Velu
// Isogeny evaluation - Huff curve
// Coefficients - Edwards curve
void xISOG_sqrt_Huff_Edwards(proj *A, proj *P, proj const *K, long long k)
{
   assert (k >= 3);
    assert (k % 2 == 1);


    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);
 

    fp tmp0, tmp1, tmp2;

    proj Aed, Am; //Huff coefficient to Edwards
    
    fp_add3(&Aed.z, &A->x, &A->x); //2(C-D)^2
    fp_add3(&Am.x, &Aed.z, &A->z); // Anx= 2(C-D)^2+4CD
    fp_add3(&Am.x, &Am.x, &Am.x); // Anx= 2(C-D)^2+4CD
    fp_add3(&Aed.z, &Aed.z, &Aed.z); // Aedz = 4(C-D)^2

  

    fp_add3(&Aed.x, &A->z, &A->z);
    fp_add3(&Aed.x, &Aed.x, &Aed.x);

    fp_add3(&Aed.x, &Aed.x, &Aed.z);
    memcpy(&Am.z,&A->z, sizeof(fp));
   
    fp Psum, Pdif;
    fp coeffQ[4];
    
    fp_add3(&Psum, &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&Pdif, &P->x, &P->z);   // Pdif= x-z
    fp_add3(&coeffQ[0], &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&coeffQ[1], &P->x, &P->z);   // Pdif= x-z
    fp_sq2(&tmp0, &coeffQ[0]); //(x+z)^2
    fp_sq2(&tmp1, &coeffQ[1]); //(x-z)^2
    fp_add3(&coeffQ[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp_sub3(&coeffQ[2], &tmp0, &tmp1); // 



    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL_Huff(&M[2], A, K); Minit[2] = 1;

    
      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD_Huff(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL_Huff(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD_Huff(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL_Huff(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD_Huff(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD_Huff(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    

    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        fp_sub3(&TI[2*i],&fp_0, &M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp T1[3*gs];
      fp Tminus1[3*gs];
      fp coeff[5];

      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
        biquad_pm1_opt_Huff(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],&Am);

        biquad_both_Huff_Ed(TP+3*j,TPinv+3*j, coeff, coeffQ, A);

      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2_selfreciprocal(T1,gs); // t1=1=c
      poly_multiprod2_selfreciprocal(Tminus1,gs); //tminus1=-1=1/c

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
        fp_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
        fp_mul2(&Abatch.x,&tmp1);
        fp_mul2(&Abatch.z,&tmp0);
        fp_mul2(&tmp1, &Psum);
        fp_mul2(&tmp0, &Pdif);
        fp_add3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.x, &tmp2);
        fp_sub3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.z, &tmp2);
      }
    

    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &Abatch.x);
    fp_mul2(&Aed.x, &Abatch.z);

    // Edwards to Huff
    fp_add3(&A->x, &Aed.z, &Aed.z); //
    fp_add3(&A->x, &A->x, &A->x); // 

    fp_add3(&A->z, &Aed.x, &Aed.x); //
    fp_add3(&A->z, &A->z, &A->z); //
    fp_sub3(&A->z, &A->z, &A->x);
}


