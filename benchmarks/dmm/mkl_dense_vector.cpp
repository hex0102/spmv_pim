/* C source code is found in dgemm_example.c */
#define min(x,y) (((x) < (y)) ? (x) : (y))
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include <sys/time.h>
#include <mkl_blas.h>

//Calculate time elapsed in secs
inline float calc_elapsed (struct timeval t1, struct timeval t0, int number_of_runs) {
    float time = ((float) (t1.tv_sec-t0.tv_sec) + (float) (t1.tv_usec-t0.tv_usec) / 1000000.0)/ number_of_runs;
    printf (" %.10f\n", time);
    return time;
}

int main(int argc, char *argv[])
{
    float *A, *x, *y; // x and y are vectors
    int m, k, n, i, j;
    float alpha, beta;
    float s_initial, s_elapsed;

    struct timeval t0_beg, t0_end;

    // Number of Warm up cycles and benchmarking cycles
    int WARM_COUNT = 10;
    int LOOP_COUNT = 1000;

    printf ("\n This example computes real vector y=alpha*A*x+beta*y using \n"
            " Intel(R) MKL function dgemm, where A is a matrix, x and y are vectors, and \n"
            " alpha and beta are float precision scalars\n\n");

    if (argc == 2) {
        m = atoi(argv[1]);
        k = m;
    }
    else {
        m = 32, k = 32;
    }
    n = 1;

    printf (" Initializing data for matrix multiplication y=A*x for matrix \n"
            " A(%ix%i) and vector x(%ix%i)\n\n", m, k, k, n);

    alpha = 0.9; beta = 0.9;

    printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");

    A = (float *)mkl_malloc( m*k*sizeof( float ), 64 );
    x = (float *)mkl_malloc( k*n*sizeof( float ), 64 );
    y = (float *)mkl_malloc( m*n*sizeof( float ), 64 );

    if (A == NULL || x == NULL || y == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A);
      mkl_free(x);
      mkl_free(y);
      return 1;
    }
    printf (" Intializing matrix data \n\n");
    for (i = 0; i < (m*k); i++) {
        A[i] = (float)1/(i+1);
    }
    for (i = 0; i < (k*n); i++) {
        x[i] = (float)1/(-i-1);
    }
    for (i = 0; i < (m*n); i++) {
        y[i] = 0.0;
    }

    printf (" Finding max number of threads Intel(R) MKL can use for parallel runs \n\n");
    
    int max_threads = mkl_get_max_threads();

    printf (" Running Intel(R) MKL for %i threads \n\n", max_threads);

    mkl_set_num_threads(max_threads);

    printf (" Warm up %d matrix vector product using Intel(R) MKL dgemm function via CBLAS interface \n\n", WARM_COUNT);

    for (int r = 0; r < WARM_COUNT; r++) {
        cblas_sgemv(CblasRowMajor, CblasNoTrans, m, k, alpha, A, k, x, n, beta, y, n);
    }

    printf (" Computing %d matrix vector product using Intel(R) MKL dgemm function via CBLAS interface \n\n", LOOP_COUNT);

    gettimeofday (&t0_beg, 0);
    for (int r = 0; r < LOOP_COUNT; r++) {
        cblas_sgemv(CblasRowMajor, CblasNoTrans, m, k, alpha, A, k, x, n, beta, y, n);
    }
    gettimeofday (&t0_end, 0);
    
    printf (" == Matrix multiplication using Intel(R) MKL dgemv completed (s) == \n");
    calc_elapsed (t0_end, t0_beg, LOOP_COUNT);

    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
    //            m, n, k, alpha, A, k, B, n, beta, C, n);
    printf ("\n Computations completed.\n\n");
    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A[j+i*k]);
      }
      printf ("\n");
    }
    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(k,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.0f", x[j+i*n]);
      }
      printf ("\n");
    }
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", y[j+i*n]);
      }
      printf ("\n");
    }

    printf (" == Matrix multiplication using Intel(R) MKL dgemv completed (s) == \n");
    float time = calc_elapsed (t0_end, t0_beg, LOOP_COUNT);

    #if (POWER == 1)
    
    //Get time for number of runs
    int number_of_runs = (int) 50/time; //Make sure power scripts is measuring for time <= 10s
    printf("Number of runs: %d\n", number_of_runs); 

    // system("/home/aporvaa/research/outerspace/gemma_mkl/power_script.sh &");
    system("/home/fengsy/Documents/dmm_bench/power_script.sh &");
    for(int i = 0; i < number_of_runs; i++) {
        cblas_sgemv(CblasRowMajor, CblasNoTrans, m, k, alpha, A, k, x, n, beta, y, n);
    }

    #endif

    // printf ("\n Deallocating memory \n\n");
    mkl_free(A);
    mkl_free(x);
    mkl_free(y);
    // printf (" Example completed. \n\n");

    return 0;
}
