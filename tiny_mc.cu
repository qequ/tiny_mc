/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

//#define _XOPEN_SOURCE 500 // M_PI

#include "params.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// headers useful for cuda

#include <ctime>
#include <cuda_runtime.h>
#include <curand_kernel.h>


// global state, heat and heat square in each shell


/***
 * Photon
 ***/
__global__ void init_curand(curandState* rng_states)
{
    int gtid = blockDim.x * blockIdx.x + threadIdx.x;
    // setup a seed for every thread
    curand_init((unsigned long long)clock(), gtid, 0, &rng_states[gtid]);
}

__global__ void test_curand(curandState* rng_states)
{
    int gtid = blockDim.x * blockIdx.x + threadIdx.x;
    float rnd = curand_uniform(&rng_states[gtid]);
    printf("gtid = %i - rnd = %f\n", gtid, rnd);
}


__global__ void print_kernel()
{
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

/*
static void photon()
{
    const float albedo = MU_S * (1.0f / (MU_S + MU_A));
    const float shells_per_mfp = 1e4 * (1.0f / MICRONS_PER_SHELL) * (1.0f / (MU_A + MU_S));

    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float w = 1.0f;
    float weight = 1.0f;

    for (;;) {
        float t = -logf((float)rand()); 
        x += t * u;
        y += t * v;
        z += t * w;

        unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp; 
        if (shell > SHELLS - 1) {
            shell = SHELLS - 1;
        }
        heat[shell] += (1.0f - albedo) * weight;
        heat2[shell] += (1.0f - albedo) * (1.0f - albedo) * weight * weight;
        weight *= albedo;

        float xi1, xi2;

        do {

            xi1 = 2.0f * rand() - 1.0f;
            xi2 = 2.0f * rand() - 1.0f;
            t = xi1 * xi1 + xi2 * xi2;

        } while (1.0f < t);

        u = 2.0f * t - 1.0f;
        v = xi1 * sqrtf((1.0f - u * u) * (1.0f / t));
        w = xi2 * sqrtf((1.0f - u * u) * (1.0f / t));

        if (weight < 0.001f) { 
            if ((float)rand() > 0.1f) {
                break;
            }
            weight *= 10.0f;
        }
    }
}
*/

/***
 * Main matter
 ***/

int main(void)
{

    double block_count = ceil(PHOTONS / BLOCK_SIZE);
    unsigned int total_num_threads = block_count * BLOCK_SIZE;


    curandState* rng_states;
    cudaMallocManaged(&rng_states, total_num_threads * sizeof(curandState));

    init_curand<<<block_count, BLOCK_SIZE>>>(rng_states);
    test_curand<<<block_count, BLOCK_SIZE>>>(rng_states);
    cudaDeviceSynchronize();


    return 0;
}
