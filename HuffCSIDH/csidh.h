#ifndef CSIDH_H
#define CSIDH_H

#include <stdbool.h>

#include "params.h"

typedef struct private_key {
    int8_t e[(NUM_PRIMES + 1) / 2]; /* packed int4_t */
} private_key;

typedef struct public_key {
    fp A; /* Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;

extern const public_key base;

void csidh_private(private_key *priv);
bool csidh(public_key *out, public_key const *in, private_key const *priv);
bool csidh_hy(public_key *out, public_key const *in, private_key const *priv);
bool csidh_Huff(public_key *out, public_key const *in, private_key const *priv);

bool csidh_hy_sqrt(public_key *out, public_key const *in, private_key const *priv);
bool csidh_Huff_sqrt_opt(public_key *out, public_key const *in, private_key const *priv);


bool csidh_Huff_Edwards(public_key *out, public_key const *in, private_key const *priv);
bool csidh_Huff_sqrt_Edwards(public_key *out, public_key const *in, private_key const *priv);


#endif
