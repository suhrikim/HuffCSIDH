
#include <string.h>
#include <assert.h>

#include "uint.h"
#include "fp.h"
#include "mont.h"
#include "csidh.h"
#include "rng.h"

const public_key base = {0}; /* A = 0 */

/* TODO allow different encodings depending on parameters */
/* TODO waste less randomness */
void csidh_private(private_key *priv)
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < NUM_PRIMES; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= NUM_PRIMES)
                    break;
            }
        }
    }
}

static bool validate_rec(proj *P, proj const *A, size_t lower, size_t upper, uint *order, bool *is_supersingular)
{
    assert(lower < upper);

    if (upper - lower == 1) {

        /* now P is [(p+1) / l_lower] times the original random point */
        /* we only gain information if this multiple is non-zero */

        if (memcmp(&P->z, &fp_0, sizeof(fp))) {

            uint tmp;
            uint_set(&tmp, primes[lower]);
            xMUL(P, A, P, &tmp);

            if (memcmp(&P->z, &fp_0, sizeof(fp))) {
                /* order does not divide p+1. */
                *is_supersingular = false;
                return true;
            }

            uint_mul3_64(order, order, primes[lower]);

            if (uint_sub3(&tmp, &four_sqrt_p, order)) { /* returns borrow */
                /* order > 4 sqrt(p), hence definitely supersingular */
                *is_supersingular = true;
                return true;
            }
        }

        /* inconclusive */
        return false;
    }

    size_t mid = lower + (upper - lower + 1) / 2;

    uint cl = uint_1, cu = uint_1;
    for (size_t i = lower; i < mid; ++i)
        uint_mul3_64(&cu, &cu, primes[i]);
    for (size_t i = mid; i < upper; ++i)
        uint_mul3_64(&cl, &cl, primes[i]);

    proj Q;

    xMUL(&Q, A, P, &cu);
    xMUL(P, A, P, &cl);

    /* start with the right half; bigger primes help more */
    return validate_rec(&Q, A, mid, upper, order, is_supersingular)
        || validate_rec(P, A, lower, mid, order, is_supersingular);
}

/* never accepts invalid keys. */
bool validate(public_key const *in)
{
    /* make sure the curve is nonsingular: A^2-4 != 0 */
    {
        uint dummy;
        if (!uint_sub3(&dummy, (uint *) &in->A, &p)) /* returns borrow */
            /* A >= p */
            return false;

        fp fp_pm2;
        fp_set(&fp_pm2, 2);
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = 2 */
            return false;

        fp_sub3(&fp_pm2, &fp_0, &fp_pm2);
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = -2 */
            return false;
    }

    const proj A = {in->A, fp_1};

    do {
        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        /* maximal 2-power in p+1 */
        xDBL(&P, &A, &P);
        xDBL(&P, &A, &P);

        bool is_supersingular;
        uint order = uint_1;

        if (validate_rec(&P, &A, 0, NUM_PRIMES, &order, &is_supersingular))
            return is_supersingular;

    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}

/* compute w(aw+b)(a+bw) */
static void Huff_rhs(fp *rhs, fp const *A,  fp const *x)
{
    fp t0, t1;
    fp_sq2(&t0, x);         // t0 = w^2
    fp_mul3(&t1, x, A);     // t1 = cw
    fp_add3(&t0, &t0, &t1); // t1 =w^2+cw
    fp_add3(&t1, &t1, &fp_1);   // t0 = cw+1
    fp_mul3(rhs, &t0, &t1);
 
}

static void Huff_rhs_ed(fp *rhs, fp const *A,  fp const *x)
{
    fp t0, t1;

    fp_add3(&t0, x, &fp_1); 
    fp_sq2(&t0, &t0);         // t0 = (w+1)^2
    fp_mul3(&t1, x, A);     // t1 = cw
    fp_add3(&t1, &t1, &t1);
    fp_add3(&t1, &t1, &t1);

    fp_add3(&t0, &t0, &t1); // t0=(w+1)^2+cw
    fp_mul3(rhs, &t0, x);
 
}


/* totally not constant-time. */
void action(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
}


void action_hy(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */


    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1}; // A=A+2c, 4C


    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);
        done[sign] = true;

        for (size_t i = NUM_PRIMES-1; i < NUM_PRIMES; --i) {  //changed loop direction

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = i - 1; j < NUM_PRIMES; --j)   //changed loop direction
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
               xMUL(&K, &A, &P, &cof);
              
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                   xISOG_hy(&A, &P, &K, primes[i]);
   

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
    
}


void action_Huff(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;

            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};
    proj Aw;
    fp t0, t1;

    bool done[2] = {false, false};
    
    // transform coeff
    fp_add3(&t0, &A.x, &A.z);       // t0 =
    fp_sub3(&t1, &A.x, &A.z);
    fp_sq2(&Aw.x, &t1);             // Aw.x=(C-D)^2
    fp_sq1(&t0);
    fp_sub3(&Aw.z, &t0, &Aw.x);   // Aw.z=4AB


    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        Huff_rhs(&rhs, &A.x,  &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_Huff(&P, &Aw, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_Huff(&K, &Aw, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {


                    xISOG_Huff(&A, &Aw, &P, &K, primes[i]);


                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;
        fp_add3(&t0, &A.x, &A.z);       // t0 =
        fp_sub3(&t1, &A.x, &A.z);
        fp_sq2(&Aw.x, &t1);             // Aw.x=(C-D)^2
        fp_sq1(&t0);
        fp_sub3(&Aw.z, &t0, &Aw.x);   // Aw.z=4AB


    } while (!(done[0] && done[1]));

    out->A = A.x;

}


void action_hy_sqrt(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */


    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);
        done[sign] = true;

        for (size_t i = NUM_PRIMES-1; i < NUM_PRIMES; --i) {  //changed loop direction

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = i - 1; j < NUM_PRIMES; --j)   //changed loop direction
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
               xMUL(&K, &A, &P, &cof);
              
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    if(primes[i]>100)
                    {
                        xISOG_sqrt(&A, &P, &K, primes[i]);
                    }                        
                    else
                    {
                        xISOG_hy(&A, &P, &K, primes[i]);
                    }
   

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
    
}





void action_Huff_sqrt_opt(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};
    proj Aw;
    fp t0, t1;

    bool done[2] = {false, false};
    
    // transform coeff
    fp_add3(&t0, &A.x, &A.z);       // t0 =
    fp_sub3(&t1, &A.x, &A.z);
    fp_sq2(&Aw.x, &t1);             // Aw.x=(C-D)^2
    fp_sq1(&t0);
    fp_sub3(&Aw.z, &t0, &Aw.x);   // Aw.z=4AB


    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        Huff_rhs(&rhs, &A.x,  &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_Huff(&P, &Aw, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_Huff(&K, &Aw, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    if(primes[i]>100)
                    {
                        xISOG_Huff_sqrt_opt(&A, &Aw, &P, &K, primes[i]);
                    }                        
                    else
                    {
                        xISOG_Huff(&A, &Aw, &P, &K, primes[i]);
                    }

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;
        fp_add3(&t0, &A.x, &A.z);       // t0 =
        fp_sub3(&t1, &A.x, &A.z);
        fp_sq2(&Aw.x, &t1);             // Aw.x=(C-D)^2
        fp_sq1(&t0);
        fp_sub3(&Aw.z, &t0, &Aw.x);   // Aw.z=4AB


    } while (!(done[0] && done[1]));

    out->A = A.x;

}





void action_Huff_Edwards(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;

            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};


    bool done[2] = {false, false};
    
    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        Huff_rhs_ed(&rhs, &A.x,  &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_Huff(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_Huff(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG_Huff_Edwards(&A, &P, &K, primes[i]);


                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1; // (c+c)


    } while (!(done[0] && done[1]));

    out->A = A.x;

}




void action_Huff_sqrt_Edwards(public_key *out, public_key const *in, private_key const *priv)
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;

            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};
 
    bool done[2] = {false, false};
    
    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        Huff_rhs_ed(&rhs, &A.x,  &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_Huff(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_Huff(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {


                    if(primes[i]>100)
                    {
                        xISOG_sqrt_Huff_Edwards(&A, &P, &K, primes[i]);
                    }                        
                    else
                    {
                        xISOG_Huff_Edwards(&A, &P, &K, primes[i]);
                    }


                    


                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1; // (c+c)


    } while (!(done[0] && done[1]));

    out->A = A.x;

}
bool csidh(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action(out, in, priv);
    return true;
}


bool csidh_hy(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_hy(out, in, priv);
    return true;
}




bool csidh_Huff(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_Huff(out, in, priv);
    return true;
}


bool csidh_hy_sqrt(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_hy_sqrt(out, in, priv);
    return true;
}




bool csidh_Huff_sqrt_opt(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_Huff_sqrt_opt(out, in, priv);
    return true;
}

bool csidh_Huff_sqrt_Edwards(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_Huff_sqrt_Edwards(out, in, priv);
    return true;
}



bool csidh_Huff_Edwards(public_key *out, public_key const *in, private_key const *priv)
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */

    action_Huff_Edwards(out, in, priv);
    return true;
}