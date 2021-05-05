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

        x += t * u;
        y += t * v;
        z += t * w;

        unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp; /* absorb */
        if (shell > SHELLS - 1) {
            shell = SHELLS - 1;
        }
        heat[shell] += (1.0f - albedo) * weight;
        heat2[shell] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */
        weight *= albedo;

        /* New direction, rejection method */
        float xi1, xi2;

        do {

            xi1 = 2.0f * genRand(&r) - 1.0f;
            xi2 = 2.0f * genRand(&r) - 1.0f;
            t = xi1 * xi1 + xi2 * xi2;

        } while (1.0f < t);

        u = 2.0f * t - 1.0f;
        v = xi1 * sqrtf((1.0f - u * u) * (1.0f / t));
        w = xi2 * sqrtf((1.0f - u * u) * (1.0f / t));

        if (weight < 0.001f) { /* roulette */
            if ((float)genRand(&r) > 0.1f) {
                break;
            }
            weight *= 10.0f;
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
