
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "fp.h"
#include "csidh.h"

void uint_print(uint const *x)
{
    for (size_t i = 8*LIMBS-1; i < 8*LIMBS; --i)
        printf("%02hhx", i[(unsigned char *) x->c]);
}

void fp_print(fp const *x)
{
    uint y;
    fp_dec(&y, x);
    uint_print(&y);
}

void test_csidh()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh(&pub_alice, &base, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh(&pub_bob, &base, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_Hybrid()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Hybrid (Mont. + Ed.) \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hy(&pub_alice, &base, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hy(&pub_bob, &base, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hy(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hy(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}

void test_Huff()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
     public_key base_Huff = {base_a}; /* A = 0 */
    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH using Huff curve\n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff(&pub_alice, &base_Huff, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff(&pub_bob, &base_Huff, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_Hybrid_sqrt()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Hybrid (Mont. + Ed.) using SQRT Velu\n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hy_sqrt(&pub_alice, &base, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hy_sqrt(&pub_bob, &base, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hy_sqrt(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hy_sqrt(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}



void test_Huff_sqrt_opt()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
public_key base_Huff = {base_a}; /* A = 0 */
    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Huff using SQRT Velu\n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_sqrt_opt(&pub_alice, &base_Huff, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_sqrt_opt(&pub_bob, &base_Huff, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_sqrt_opt(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_sqrt_opt(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_Huff_edwards()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
public_key base_Huff = {base_a}; /* A = 0 */
    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Huff using Edwards curves\n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    memcpy(&base_Huff, &fp_1, sizeof(fp));


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_Edwards(&pub_alice, &base_Huff, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_Edwards(&pub_bob, &base_Huff, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_Edwards(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_Edwards(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_Huff_sqrt_edwards()
{

    clock_t t0, t1;

    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
public_key base_Huff = {base_a}; /* A = 0 */
    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Huff using Edwards curves + Sqrt_Velu\n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    memcpy(&base_Huff, &fp_1, sizeof(fp));


    t0 = clock();
    csidh_private(&priv_alice);
    t1 = clock();

    printf("Alice's private key   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_alice); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_alice]);
    printf("\n\n");

    t0 = clock();
    csidh_private(&priv_bob);
    t1 = clock();

    printf("Bob's private key     (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    for (size_t i = 0; i < sizeof(priv_bob); ++i)
        printf("%02hhx", i[(uint8_t *) &priv_bob]);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_sqrt_Edwards(&pub_alice, &base_Huff, &priv_alice));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_sqrt_Edwards(&pub_bob, &base_Huff, &priv_bob));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_Huff_sqrt_Edwards(&shared_alice, &pub_bob, &priv_alice));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_Huff_sqrt_Edwards(&shared_bob, &pub_alice, &priv_bob));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}
int main()
{

    
   test_csidh();
   test_Hybrid();
   test_Huff();
   test_Hybrid_sqrt();

   test_Huff_sqrt_opt();

   test_Huff_edwards();
   test_Huff_sqrt_edwards();
}

