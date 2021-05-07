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


/***
 * Photon
 ***/

static void photon(MTRand r)
{
    /*
    const float albedo = MU_S * (1.0f / (MU_S + MU_A));
    const float shells_per_mfp = 1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S));
    */
    const __m256 albedo = _mm256_set1_ps(MU_S * (1.0f / (MU_S + MU_A)));
    const __m256 shells_per_mfp = _mm256_set1_ps(1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S)));

    float array_rnd[8];
    /* launch */

    /*
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float w = 1.0f;
    float weight = 1.0f;
    */

    // seteando todos los vectores con el valor dado
    __m256 x = _mm256_set1_ps(0.0f);
    __m256 y = _mm256_set1_ps(0.0f);
    __m256 z = _mm256_set1_ps(0.0f);
    __m256 u = _mm256_set1_ps(0.0f);
    __m256 v = _mm256_set1_ps(0.0f);
    __m256 w = _mm256_set1_ps(1.0f);
    __m256 weight = _mm256_set1_ps(1.0f);


    // helper vectors
    __m256 ones_vector = _mm256_set1_ps(1.0f);

    for (;;) {

        for (unsigned int i = 0; i < 8; ++i) {
            array_rnd[i] = (float)genRand(&r);
        }

        /* move */

        //float t = -logf((float)genRand(&r));
        __m256 t = _mm256_set_ps(-logf(array_rnd[0]),
                                 -logf(array_rnd[1]),
                                 -logf(array_rnd[2]),
                                 -logf(array_rnd[3]),
                                 -logf(array_rnd[4]),
                                 -logf(array_rnd[5]),
                                 -logf(array_rnd[6]),
                                 -logf(array_rnd[7]));

        /*
        x += t * u;
        y += t * v;
        z += t * w;
        */
        //fmadd_ps(a, b, c) == (a * b) + c
        x = _mm256_fmadd_ps(t, u, x);
        y = _mm256_fmadd_ps(t, v, y);
        z = _mm256_fmadd_ps(t, w, z);

        // cuadrados de números
        __m256 x_squared = _mm256_mul_ps(x, x);
        __m256 y_squared = _mm256_mul_ps(y, y);
        __m256 z_squared = _mm256_mul_ps(z, z);

        __m256 sum_cord = _mm256_add_ps(_mm256_add_ps(x_squared, y_squared), z_squared);
        //unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp; 
        /* absorb */
        __m256i shell_vector = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_sqrt_ps(sum_cord), shells_per_mfp));
        __m256i max_shell_vector = _mm256_set1_epi32(SHELLS - 1);
        
        /*
        if (shell > SHELLS - 1) {
            shell = SHELLS - 1;
        }
        */
        shell_vector = _mm256_min_epi32(shell_vector, max_shell_vector);


        /*
        heat[shell] += (1.0f - albedo) * weight;
        heat2[shell] += (1.0f - albedo) * (1.0f - albedo) * weight * weight;  
        weight *= albedo;
        */
        __m256 helper_vector = _mm256_mul_ps(_mm256_sub_ps(ones_vector, albedo), weight);
        __m256 helper_vector_squared = _mm256_mul_ps(helper_vector, helper_vector); /* add up squares */

        for (unsigned int i=0; i < 8; ++i){
            heat[shell_vector[i]] += helper_vector[i];
            heat2[shell_vector[i]] += helper_vector_squared[i];
        }

        weight = _mm256_mul_ps(weight, albedo);



        /* New direction, rejection method */
        //float xi1, xi2;

        /*
        do {

            xi1 = 2.0f * genRand(&r) - 1.0f;
            xi2 = 2.0f * genRand(&r) - 1.0f;
            t = xi1 * xi1 + xi2 * xi2;

        } while (1.0f < t);
        */

       // Aux. DELETE WHEN THE DOWHILE IS UPDATED.

       __m256 xi1 = _mm256_set1_ps(1.0f);
       __m256 xi2 = _mm256_set1_ps(1.0f);

        /*
        u = 2.0f * t - 1.0f;
        v = xi1 * sqrtf((1.0f - u * u) * (1.0f / t));
        w = xi2 * sqrtf((1.0f - u * u) * (1.0f / t));
        */

        // Aux vector.
        __m256 twos_vector = _mm256_set1_ps(2.0f);

        __m256 u = _mm256_sub_ps (_mm256_mul_ps(twos_vector,t),ones_vector);

        __m256 root = _mm256_sqrt_ps(_mm256_mul_ps(_mm256_sub_ps(ones_vector,_mm256_mul_ps(u,u)),_mm256_div_ps(ones_vector,t)));

        __m256 v = _mm256_mul_ps(xi1,root);
            
        __m256 w = _mm256_mul_ps(xi2,root);

        
        /*
        if (weight < 0.001f) { 
            if ((float)genRand(&r) > 0.1f) {
                break;
            }
            weight *= 10.0f;
        }
        */
       
       // Vector de 0.001s.
       __m256 little_vector = _mm256_set1_ps(0.001f);

       // Resultado de comparar si "weight" es menor a 0.001.
       __m256 ifResult1 = _mm256_cmp_ps(weight,little_vector,_CMP_LT_OS);

       // Random vector.
       float array_rnd2[8];

       for (unsigned int i = 0; i < 8; ++i) {
            array_rnd2[i] = (float)genRand(&r);
        }

        __m256 Rdm = _mm256_set_ps(-logf(array_rnd2[0]),
                                 -logf(array_rnd2[1]),
                                 -logf(array_rnd2[2]),
                                 -logf(array_rnd2[3]),
                                 -logf(array_rnd2[4]),
                                 -logf(array_rnd2[5]),
                                 -logf(array_rnd2[6]),
                                 -logf(array_rnd2[7]));

        // Vector de 0.1s.
        __m256 not_that_little_vector = _mm256_set1_ps(0.1f);

        // Resultado de comparar si "Rdm" es mayor a 0.1.
        __m256 ifResult2 = _mm256_cmp_ps(Rdm,not_that_little_vector,_CMP_GT_OS);
        
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
    for (unsigned int i = 0; i < PHOTONS; ++i) {
        photon(r);
    }
    // stop timer
    double end = wtime();
    assert(start <= end);
    double elapsed = end - start;

    /*
    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);
    */
    //
    printf("%lf\n", 1e-3 * PHOTONS / elapsed);

    /*
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
