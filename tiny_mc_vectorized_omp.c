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
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";


// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];
unsigned int photon_count = 0;

int check_t_greater_1(float* t)
{
    for (unsigned int i = 0; i < 8; ++i) {
        if (1.0f < t[i]) {
            return 1;
        }
    }
    return 0;
}

int check_t_correct(float * t) {
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
    float heat_pr[SHELLS];
    float heat2_pr[SHELLS];
    float x[8] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f };
    float y[8] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f };
    float z[8] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f };
    float u[8] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f };
    float v[8] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f };
    float w[8] = {
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
    };
    float weight[8] = {
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
    };

    float t[8];
    float albedo[8];
    float shells_per_mfp[8];
    unsigned int shell[8];

    #pragma omp simd aligned(albedo, shells_per_mfp:32)
    for (unsigned int i = 0; i < 8; ++i) {
        albedo[i] = MU_S * (1.0f / (MU_S + MU_A));
        shells_per_mfp[i] = 1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S));
    }


    while (photon_count < PHOTONS) {

/* launch */
        #pragma omp simd aligned(t, x, y, z, shell, weight:32)
        for (unsigned int i = 0; i < 8; ++i) {
            t[i] = -logf((float)genRand(&r)); /*move*/
            x[i] += t[i] * u[i];
            y[i] += t[i] * v[i];
            z[i] += t[i] * w[i];
            shell[i] = (unsigned int)sqrtf(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) * shells_per_mfp[i]; /*absorb*/

            if (shell[i] > SHELLS - 1) {
                shell[i] = SHELLS - 1;
            }

            heat_pr[shell[i]] += (1.0f - albedo[i]) * weight[i];
            heat2_pr[shell[i]] += (1.0f - albedo[i]) * (1.0f - albedo[i]) * weight[i] * weight[i]; /* add up squares */
            weight[i] *= albedo[i];
        }


        /* New direction, rejection method */
        float xi1[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        float xi2[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        int t_mask[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

        do {
            #pragma omp simd aligned(xi1, xi2:32)
            for (unsigned int i = 0; i < 8; ++i) {
                xi1[i] = 2.0f * genRand(&r) - 1.0f;
                xi2[i] = 2.0f * genRand(&r) - 1.0f;

                if (t_mask[i] == 0) {
                    t[i] = xi1[i] * xi1[i] + xi2[i] * xi2[i];
                }
            }
            #pragma omp simd aligned(t:32)
            for (unsigned int i = 0; i < 8; ++i) {
                if (t[i] <= 1.0f) {
                    t_mask[i] = 1;
                }
            }


        } while (check_t_greater_1(t));


        #pragma omp simd aligned(u, v, w, t, xi1, xi2, weight, x, y, z:32)
        for (unsigned int i = 0; i < 8; ++i) {

            u[i] = 2.0f * t[i] - 1.0f;
            v[i] = xi1[i] * sqrtf((1.0f - u[i] * u[i]) * (1.0f / t[i]));
            w[i] = xi2[i] * sqrtf((1.0f - u[i] * u[i]) * (1.0f / t[i]));


            if (weight[i] < 0.001f) { /* roulette */
                weight[i] *= 10.0f;

                if ((float)genRand(&r) > 0.1f) {
                    // reset lane
                    x[i] = 0.0f;
                    y[i] = 0.0f;
                    z[i] = 0.0f;
                    u[i] = 0.0f;
                    v[i] = 0.0f;
                    w[i] = 1.0f;
                    weight[i] = 1.0f;
                    #pragma omp atomic
                    photon_count++;
                    #pragma omp flush
                }
            }
        }
    }
    #pragma omp simd
    for (unsigned int i = 0; i < SHELLS; ++i) {
        #pragma omp atomic
            heat[i] += heat_pr[i];
        #pragma omp atomic
            heat2[i] += heat2_pr[i];
        
        
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
#pragma omp parallel shared(photon_count, heat, heat2)
    {
        photon(r);
    }

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


    //printf("%lf\n", 1e-3 * PHOTONS / elapsed);
    /*
    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);

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
