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

#include "helper_cuda.h"
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


__global__ void photon(float* global_heat, float* global_heat2, curandState* rng_states)
{
    int gtid = blockDim.x * blockIdx.x + threadIdx.x;

    // a photon per thread and if PHOTONS is not a multiple of 32 cape it
    if (gtid <= PHOTONS) {
        curandState thread_rng_state = rng_states[gtid];

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
            float t = -logf((float)curand_uniform(&thread_rng_state));
            x += t * u;
            y += t * v;
            z += t * w;

            unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;
            if (shell > SHELLS - 1) {
                shell = SHELLS - 1;
            }
            // atomic add
            atomicAdd(&global_heat[shell], (1.0f - albedo) * weight);
            atomicAdd(&global_heat2[shell], (1.0f - albedo) * (1.0f - albedo) * weight * weight);
            weight *= albedo;

            float xi1, xi2;

            do {

                xi1 = 2.0f * curand_uniform(&thread_rng_state) - 1.0f;
                xi2 = 2.0f * curand_uniform(&thread_rng_state) - 1.0f;
                t = xi1 * xi1 + xi2 * xi2;

            } while (1.0f < t);

            u = 2.0f * t - 1.0f;
            v = xi1 * sqrtf((1.0f - u * u) * (1.0f / t));
            w = xi2 * sqrtf((1.0f - u * u) * (1.0f / t));

            if (weight < 0.001f) {
                if ((float)curand_uniform(&thread_rng_state) > 0.1f) {
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
    //get block_according to number of threads per block and PHOTONS
    double block_count = ceil(PHOTONS / BLOCK_SIZE);
    unsigned int total_num_threads = block_count * BLOCK_SIZE;

    // initialize heat and heat2 to be shared between cpu and gpu
    float * heat;
    float * heat2;
    checkCudaCall(cudaMallocManaged(&heat, SHELLS * sizeof(float)));
    checkCudaCall(cudaMallocManaged(&heat2, SHELLS * sizeof(float)));

    for (int i = 0; i < SHELLS; i++) {
        heat[i] = 0;
        heat2[i] = 0;
    }

    // init curand
    curandState* rng_states;
    checkCudaCall(cudaMallocManaged(&rng_states, total_num_threads * sizeof(curandState)));

    init_curand<<<block_count, BLOCK_SIZE>>>(rng_states);
    //test_curand<<<block_count, BLOCK_SIZE>>>(rng_states);
    cudaDeviceSynchronize();

    checkCudaCall(cudaGetLastError());

    // variables to measure time
    float time;
    cudaEvent_t start, stop;

    checkCudaCall( cudaEventCreate(&start) );
    checkCudaCall( cudaEventCreate(&stop) );
    checkCudaCall( cudaEventRecord(start, 0) );


    //photon
    photon<<<block_count, BLOCK_SIZE>>>(heat, heat2, rng_states);
    checkCudaCall(cudaGetLastError());

    checkCudaCall( cudaEventRecord(stop, 0) );
    checkCudaCall( cudaEventSynchronize(stop) );
    checkCudaCall( cudaEventElapsedTime(&time, start, stop) );

    // the measure of time in cuda is ms
    float elapsed = time / 1000;
    
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
    

    return 0;
}
