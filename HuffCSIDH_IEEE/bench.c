
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <inttypes.h>

#include "rng.h"
#include "csidh.h"


#define OS_WIN       1
#define OS_LINUX     2

#if defined(__WINDOWS__)        // Microsoft Windows OS
    #define OS_TARGET OS_WIN
#else defined(__LINUX__)        // Linux OS
    #define OS_TARGET OS_LINUX 
#endif

int64_t cpucycles(void)
{ // Access system counter for benchmarking
#if (OS_TARGET == OS_WIN) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    return __rdtsc();
#elif (OS_TARGET == OS_WIN) && (TARGET == TARGET_ARM)
    return __rdpmccntr64();
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    unsigned int hi, lo;

    __asm__ __volatile__ ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    struct timespec time;

    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
#else
    return 0;            
#endif
}

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

#define loop 100000
/* defaults */

#ifndef BENCH_ITS
    #define BENCH_ITS 1000
#endif

#if !defined(BENCH_VAL) && !defined(BENCH_ACT)
    #define BENCH_VAL 1
    #define BENCH_ACT 1
#endif
#ifndef BENCH_VAL
    #define BENCH_VAL 0
#endif
#ifndef BENCH_ACT
    #define BENCH_ACT 0
#endif


const unsigned long its = BENCH_ITS;
const bool val = BENCH_VAL, act = BENCH_ACT;


const size_t stacksz = 0x8000;  /* 32k */

/* csidh.c */
bool validate(public_key const *in);
void action(public_key *out, public_key const *in, private_key const *priv);


void csidh_bench()
{

    unsigned long long  alice_keygen=0;
    int i;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;


    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 

    for(i=0; i<loop; i++)
    {
        csidh_private(&priv_alice);
        csidh_private(&priv_bob);

        cycles1 = cpucycles();
        csidh(&pub_alice, &base, &priv_alice);
        csidh(&pub_bob, &base, &priv_bob);
        csidh(&shared_alice, &pub_bob, &priv_alice);
        csidh(&shared_bob, &pub_alice, &priv_bob);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");
    printf("\n\n");

}

void csidh_hy_bench()
{

    unsigned long long  alice_keygen=0;
    int i;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;


    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH hybrid original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 

    for(i=0; i<loop; i++)
    {
        csidh_private(&priv_alice);
        csidh_private(&priv_bob);

        cycles1 = cpucycles();
        csidh_hy(&pub_alice, &base, &priv_alice);
        csidh_hy(&pub_bob, &base, &priv_bob);
        csidh_hy(&shared_alice, &pub_bob, &priv_alice);
        csidh_hy(&shared_bob, &pub_alice, &priv_bob);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");
    printf("\n\n");

}


void csidh_Huff_bench()
{

    unsigned long long  alice_keygen=0;
    int i;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    public_key base_Huff = {base_a}; 
    unsigned long long cycles1, cycles2;

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH HUFF  \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
  
    for(i=0; i<loop; i++)
    {
        csidh_private(&priv_alice);
        csidh_private(&priv_bob);

        cycles1 = cpucycles();
        csidh_Huff(&pub_alice, &base_Huff, &priv_alice);
        csidh_Huff(&pub_bob, &base_Huff, &priv_bob);
        csidh_Huff(&shared_alice, &pub_bob, &priv_alice);
        csidh_Huff(&shared_bob, &pub_alice, &priv_bob);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");
    printf("\n\n");

}



void csidh_hy_sqrt_bench()
{

    unsigned long long  alice_keygen=0;
    int i;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH HYBRID SQRT \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 

    for(i=0; i<loop; i++)
    {
        csidh_private(&priv_alice);
        csidh_private(&priv_bob);

        cycles1 = cpucycles();
        csidh_hy_sqrt(&pub_alice, &base, &priv_alice);
        csidh_hy_sqrt(&pub_bob, &base, &priv_bob);
        csidh_hy_sqrt(&shared_alice, &pub_bob, &priv_alice);
        csidh_hy_sqrt(&shared_bob, &pub_alice, &priv_bob);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}


void csidh_Huff_opt_sqrt_bench()
{

    unsigned long long  alice_keygen=0;
    int i;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    public_key base_Huff = {base_a}; 
    unsigned long long cycles1, cycles2;


    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH HUFF OPT_SQRT \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 

    for(i=0; i<loop; i++)
    {
        csidh_private(&priv_alice);
        csidh_private(&priv_bob);

        cycles1 = cpucycles();
        csidh_Huff_sqrt_opt(&pub_alice, &base_Huff, &priv_alice);
        csidh_Huff_sqrt_opt(&pub_bob, &base_Huff, &priv_bob);
        csidh_Huff_sqrt_opt(&shared_alice, &pub_bob, &priv_alice);
        csidh_Huff_sqrt_opt(&shared_bob, &pub_alice, &priv_bob);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");
    printf("\n\n");

}


int main()
{
    csidh_bench(); 
    csidh_hy_bench();
    csidh_Huff_bench();
   csidh_hy_sqrt_bench();
   csidh_Huff_opt_sqrt_bench();
}

