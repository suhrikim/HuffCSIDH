#ifndef MONT_H
#define MONT_H

#include "params.h"
#include "poly.h"

void xDBL(proj *Q, proj const *A, proj const *P);
void xDBL_Huff(proj *Q, proj const *A, proj const *P);

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xADD_Huff(proj *S, proj const *P, proj const *Q, proj const *PQ);

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xDBLADD_Huff(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);

void xMUL(proj *Q, proj const *A, proj const *P, uint const *k);
void xMUL_Huff(proj *Q, proj const *A, proj const *P, uint const *k);

void xISOG(proj *A, proj *P, proj const *K, uint64_t k);
void xISOG_hy(proj *A, proj *P, proj const *K, uint64_t k);
void xISOG_Huff(proj *A, proj *Aw, proj *P, proj const *K, uint64_t k);
void xISOG_sqrt(proj *A, proj *P, proj const *K, long long k);

void xISOG_Huff_sqrt_opt(proj *A, proj *Aw, proj *P, proj const *K, long long k);



#endif
