/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#define _XOPEN_SOURCE 500 // M_PI

#include "mtwister.h"
#include "params.h"
#include "wtime.h"

#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";


// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];


// useful for checking if a vector mask has its 8 elements equal to 0xFFFFFFFF

unsigned int check_t_correct(__m256 t)
{
    for (unsigned int i = 0; i < 8; ++i) {
        if (1.0f < t[i]) {
            return 0;
        }
    }
    return 1;
}



/***
 * Photon
 ***/

static void photon(MTRand r)
{
    
    const __m256 albedo = _mm256_set1_ps(MU_S * (1.0f / (MU_S + MU_A)));
    const __m256 shells_per_mfp = _mm256_set1_ps(1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S)));


    // Random arrays.
    float array_rnd[8];
    float array_rnd2[8];
    float index_array[8];

    /* launch */


    // helper vectors
    const __m256 zeros_vector = _mm256_set1_ps(0.0f);
    const __m256 ones_vector = _mm256_set1_ps(1.0f);
    const __m256 twos_vector = _mm256_set1_ps(2.0f);
    const __m256 tens_vector = _mm256_set1_ps(10.0f);
    const __m256 l_vector = _mm256_set1_ps(0.001f);
    const __m256 tl_vector = _mm256_set1_ps(0.1f);

    // seteando todos los vectores con el valor dado
    __m256 x = _mm256_set1_ps(0.0f);
    __m256 y = _mm256_set1_ps(0.0f);
    __m256 z = _mm256_set1_ps(0.0f);
    __m256 u = _mm256_set1_ps(0.0f);
    __m256 v = _mm256_set1_ps(0.0f);
    __m256 w = _mm256_set1_ps(1.0f);
    __m256 weight = _mm256_set1_ps(1.0f);

    unsigned int photon_count = 0;

    while (photon_count < PHOTONS) {

        for (unsigned int i = 0; i < 8; ++i) {
            array_rnd[i] = -logf((float)genRand(&r));
        }

        /* move */

        //float t = -logf((float)genRand(&r));
        __m256 t = _mm256_load_ps(array_rnd);


        // fmadd_ps(a, b, c) == (a * b) + c
        x = _mm256_fmadd_ps(t, u, x);
        y = _mm256_fmadd_ps(t, v, y);
        z = _mm256_fmadd_ps(t, w, z);

        // cuadrados de números
        __m256 x_squared = _mm256_mul_ps(x, x);
        __m256 y_squared = _mm256_mul_ps(y, y);
        __m256 z_squared = _mm256_mul_ps(z, z);

        __m256 sum_cord = _mm256_add_ps(_mm256_add_ps(x_squared, y_squared), z_squared);


        /* absorb */
        __m256 shell_vector = _mm256_mul_ps(_mm256_sqrt_ps(sum_cord), shells_per_mfp);
        __m256 max_shell_vector = _mm256_set1_ps(SHELLS - 1);

        shell_vector = _mm256_min_ps(shell_vector, max_shell_vector);


        __m256 helper_vector = _mm256_mul_ps(_mm256_sub_ps(ones_vector, albedo), weight);
        __m256 helper_vector_squared = _mm256_mul_ps(helper_vector, helper_vector); /* add up squares */

       _mm256_store_ps(index_array, shell_vector);

        for (unsigned int i = 0; i < 8; ++i) {
            heat[(unsigned int)index_array[i]] += (float)helper_vector[i];
            heat2[(unsigned int)index_array[i]] += (float)helper_vector_squared[i];
        }

        weight = _mm256_mul_ps(weight, albedo);


        __m256 vec_mask = _mm256_set1_ps(0.0f);

        for (unsigned int i = 0; i < 8; ++i) {
            array_rnd[i] = 2.0f * (float)genRand(&r) - 1.0f;
            array_rnd2[i] = 2.0f * (float)genRand(&r) - 1.0f;
        }
        __m256 xi1 = _mm256_load_ps(array_rnd);
        __m256 xi2 = _mm256_load_ps(array_rnd2);


        t = _mm256_add_ps(_mm256_mul_ps(xi1, xi1), _mm256_mul_ps(xi2, xi2));
        int mask_int;
        do {


            for (unsigned int i = 0; i < 8; ++i) {
                array_rnd[i] = 2.0f * (float)genRand(&r) - 1.0f;
                array_rnd2[i] = 2.0f * (float)genRand(&r) - 1.0f;
            }
            __m256 xi1 = _mm256_load_ps(array_rnd);
            __m256 xi2 = _mm256_load_ps(array_rnd2);

            // 1 ==_CMP_LT_OS == <
            vec_mask = _mm256_cmp_ps(t, ones_vector, 1);
            mask_int = _mm256_movemask_ps(vec_mask);

            t = _mm256_blendv_ps(_mm256_add_ps(_mm256_mul_ps(xi1, xi1), _mm256_mul_ps(xi2, xi2)), t, vec_mask);
        } while (mask_int != 0xFF);
        assert(check_t_correct(t));


        // _mm256_fmsub_ps(a, b, c) == (a * b) - c
        __m256 u = _mm256_fmsub_ps(twos_vector, t, ones_vector);

        __m256 root = _mm256_sqrt_ps(_mm256_mul_ps(_mm256_sub_ps(ones_vector, _mm256_mul_ps(u, u)), _mm256_div_ps(ones_vector, t)));

        __m256 v = _mm256_mul_ps(xi1, root);

        __m256 w = _mm256_mul_ps(xi2, root);


        // weight < 0.001f
        __m256 weight_mask = _mm256_cmp_ps(weight, l_vector, _CMP_LT_OS);

        weight = _mm256_blendv_ps(weight, _mm256_mul_ps(weight, tens_vector), weight_mask);

        for (unsigned int i = 0; i < 8; ++i) {
            array_rnd[i] = genRand(&r);
        }
        __m256 rand_vec = _mm256_load_ps(array_rnd);

        __m256 roulette_mask = _mm256_cmp_ps(rand_vec, tl_vector, _CMP_GT_OS);


        x = _mm256_blendv_ps(x, _mm256_blendv_ps(x, zeros_vector, weight_mask), roulette_mask);
        y = _mm256_blendv_ps(y, _mm256_blendv_ps(y, zeros_vector, weight_mask), roulette_mask);
        z = _mm256_blendv_ps(z, _mm256_blendv_ps(z, zeros_vector, weight_mask), roulette_mask);
        u = _mm256_blendv_ps(u, _mm256_blendv_ps(u, zeros_vector, weight_mask), roulette_mask);
        v = _mm256_blendv_ps(v, _mm256_blendv_ps(v, zeros_vector, weight_mask), roulette_mask);
        w = _mm256_blendv_ps(w, _mm256_blendv_ps(w, ones_vector, weight_mask), roulette_mask);
        weight = _mm256_blendv_ps(weight, _mm256_blendv_ps(weight, ones_vector, weight_mask), roulette_mask);

        for (unsigned int i = 0; i < 8; ++i) {
            if (roulette_mask[i] && weight_mask[i]) {
                photon_count++;
            }
        }
    }
}


/***
 * Main matter
 ***/

int main(void)
{
    // heading
    /*
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n#\n", PHOTONS);
    */
    // configure RNG
    //srand(SEED);
    MTRand r = seedRand(SEED);

    // start timer
    double start = wtime();
    // simulation
    photon(r);

    // stop timer
    double end = wtime();
    assert(start <= end);
    double elapsed = end - start;

    FILE* fptr = fopen("photons_results.txt", "a");
    if (fptr == NULL) {
        exit(EXIT_FAILURE);
    }
    fprintf(fptr, "%lf\n", 1e-3 * PHOTONS / elapsed);

    fclose(fptr);

    /*
    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);

    
    printf("%lf\n", 1e-3 * PHOTONS / elapsed);


    printf("# Radius\tHeat\n");
    printf("# [microns]\t[W/cm^3]\tError\n");
    float t = 4.0f * M_PI * powf(MICRONS_PER_SHELL, 3.0f) * PHOTONS / 1e12;
    for (unsigned int i = 0; i < SHELLS - 1; ++i) {
        printf("%6.0f\t%12.5f\t%12.5f\n", i * (float)MICRONS_PER_SHELL,
               heat[i] / t / (i * i + i + 1.0 / 3.0),
               sqrt(heat2[i] - heat[i] * heat[i] / PHOTONS) / t / (i * i + i + 1.0f / 3.0f));
    }
    printf("# extra\t%12.5f\n", heat[SHELLS - 1] / PHOTONS);
    */

    return 0;
}
