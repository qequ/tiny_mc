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


/***
 * Photon
 ***/

static void photon(MTRand r)
{
    const float albedo = MU_S * (1.0f / (MU_S + MU_A));
    const float shells_per_mfp = 1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S));

    /* launch */
    float xi1, xi2;
    float t;
    unsigned int shell;
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float w = 1.0f;
    float weight = 1.0f;

    float heat_pr[SHELLS];
    float heat2_pr[SHELLS];

    #pragma omp parallel for schedule(dynamic) private(heat_pr, heat2_pr, x, y, z, u, v, w, t, xi1, xi2, shell, weight)
    for (unsigned int i = 0; i < PHOTONS; ++i) {

        for (;;) {
            t = -logf((float)genRand(&r)); /* move */
            x += t * u;
            y += t * v;
            z += t * w;

            shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp; /* absorb */
            if (shell > SHELLS - 1) {
                shell = SHELLS - 1;
            }
            heat_pr[shell] += (1.0f - albedo) * weight;
            heat2_pr[shell] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */

            weight *= albedo;

            /* New direction, rejection method */

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
                    // cargando a los shared heats a partir de los privados
                    #pragma omp critical
                    {
                        for (int n = 0; n < SHELLS; ++n) {
                            heat[n] += heat_pr[n];
                            heat2[n] += heat2_pr[n];
                        }
                    }
                    break;
                }
                weight *= 10.0f;
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

    printf("%lf\n", 1e-3 * PHOTONS / elapsed);
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
