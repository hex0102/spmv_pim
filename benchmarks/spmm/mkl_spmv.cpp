/*******************************************************************************
* Copyright 2013-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/******************************************************************************
*
* Consider the matrix A
*
*                 |  10     11      0     0     0   |
*                 |   0      0     12    13     0   |
*   A    =        |  15      0      0     0    14   |,
*                 |   0     16     17     0     0   |
*                 |   0      0      0    18    19   |
*
* and diagonal matrix B
*
*                 |   5      0      0     0     0   |
*                 |   0      6      0     0     0   |
*   B    =        |   0      0      7     0     0   |.
*                 |   0      0      0     8     0   |
*                 |   0      0      0     0     9   |
*
*  Both matrices A and B are stored in a zero-based compressed sparse row (CSR) storage
*  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
*  Intel Math Kernel Library Developer Reference) as follows:
*
*           values_A = ( 10  11  12  13  15  14  16  17  18  19 )
*          columns_A = (  0   1   2   3   0   4   1   2   3   4 )
*         rowIndex_A = (  0       2       4       6       8      10 )
*
*           values_B = ( 5  6  7  8  9  )
*          columns_B = ( 0  1  2  3  4  )
*         rowIndex_B = ( 0  1  2  3  4  5 )
*
*  The example computes two scalar products :
*
*         < (A*B)*x ,       y > = left,   using MKL_SPARSE_D_SPMM and CBLAS_DDOT.
*         <     B*x , (A^t)*y > = right,  using MKL_SPARSE_D_MV and CBLAS_DDOT.
*
*         These products should result in the same value. To obtain matrix C,
*         use MKL_SPARSE_D_EXPORT_CSR and print the result.
*
******************************************************************************/

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include "matrix_market.hpp"
#include <sys/time.h>
#include <omp.h>

using namespace std;

// #define M 5
// #define NNZ 10
#define ALIGN 128
// #define DEBUG
/* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
   use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = 1;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

void readMatA (float ** vals, MKL_INT ** cols, MKL_INT ** rowPtr, MKL_INT * m, MKL_INT * nnz){

    vector<map<unsigned int, float> > stl_A;
    
    //Read Matrix
    if (!read_matrix_market_file(stl_A, "./matA.txt"))
    {
      cout << "Error reading Matrix file" << endl;
      exit(-1);
    }

    // convert to plain CSR arrays:
    *m = stl_A.size();
    *nnz = 0;
    for (size_t row=0; row<stl_A.size(); ++row)
      *nnz += stl_A[row].size();

    
    /* Allocation of memory */
    *rowPtr = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (stl_A.size()+1), ALIGN);
    *cols = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * *nnz, ALIGN);
    *vals = (float *)mkl_malloc(sizeof(float) * *nnz, ALIGN);
    
     
    int current_index = 0;
    for (size_t row=0; row<stl_A.size(); ++row) {
      (*rowPtr)[row] = current_index;
      for (std::map<unsigned int, float>::const_iterator it = stl_A[row].begin(); it != stl_A[row].end(); ++it) {
        (*cols)[current_index] = it-> first;
        (*vals)[current_index] = it-> second;
        ++ current_index;
      }   
    }
    (*rowPtr)[stl_A.size()] = current_index;
    // set_1based_ind(*rowPtr, *cols, stl_A.size(), *nnz);

    // printf("VALS: ");
    // for (MKL_INT i = 0; i < *nnz; i++) {
    //   printf(" %f", *(*vals+i));
    // }
    // printf("\n");

    // printf("COLS: ");
    // for (int i = 0; i < *nnz; i++) {
    //   printf(" %lld", *(*cols+i));
    // }
    // printf("\n");

    // printf("ROWS: ");
    // for (int i = 0; i < stl_A.size()+1; i++) {
    //   printf(" %lld", *(*rowPtr+i));
    // }
    // printf("\n");
}

int main(int argc, char *argv[]) {

    int WARM_COUNT = 10;
    int LOOP_COUNT = 1000;

    if (argc == 2) {
        LOOP_COUNT = atoi(argv[1]);
    }

    /* Declaration of values */
    // double  *values_A = NULL, *values_B = NULL, *values_C = NULL;
    // MKL_INT *columns_A = NULL, *columns_B = NULL, *columns_C = NULL;
    // MKL_INT *rowIndex_A = NULL, *rowIndex_B = NULL, *pointerB_C = NULL, *pointerE_C = NULL;
    float  *values_A = NULL;
    MKL_INT *columns_A = NULL;
    MKL_INT *rowIndex_A = NULL;
    // MKL_INT *pointerB_A = NULL, *pointerE_A = NULL;


    // double  *rslt_mv = NULL, *rslt_mv_trans = NULL, *x = NULL, *y = NULL;
    float  *rslt_mv = NULL, *x = NULL;

    // double   left, right, residual;
    MKL_INT  rows, cols, i, j, ii, status;

    MKL_INT M, NNZ;

    sparse_index_base_t    indexing;
    struct matrix_descr    descr_type_gen;
    // sparse_matrix_t        csrA = NULL, csrB = NULL, csrC = NULL;
    sparse_matrix_t        csrA = NULL;

    float s_initial, s_elapsed;

    readMatA(&values_A, &columns_A, &rowIndex_A, &M, &NNZ);
    printf("N = %lld, NNZ = %lld\n",M, NNZ);

    /* Setting MAX number of threads */
    MKL_INT max_threads = mkl_get_max_threads();
    printf("Max threads = %lld\n", max_threads);
    mkl_set_num_threads(max_threads);

    /* Allocation of memory */
    // values_A = (double *)mkl_malloc(sizeof(double) * NNZ, ALIGN);
    // columns_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * NNZ, ALIGN);
    // rowIndex_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (M + 1), ALIGN);

    // values_B = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    // columns_B = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * M, ALIGN);
    // rowIndex_B = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (M + 1), ALIGN);

    x = (float *)mkl_malloc(sizeof(float) * M, ALIGN);
    // y = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    rslt_mv = (float *)mkl_malloc(sizeof(float) * M, ALIGN);
    // rslt_mv_trans = (double *)mkl_malloc(sizeof(double) * M, ALIGN);

    /* Set values of the variables*/
    descr_type_gen.type = SPARSE_MATRIX_TYPE_GENERAL;
    status = 0, ii = 0;
    //Matrix A 
    // for( i = 0; i < NNZ; i++ )
    //       values_A[i] = i + 10;
    // for( i = 0; i < NNZ; i++ )
    //       columns_A[i] = i % 5;
    // rowIndex_A[0] = 0;
    // for( i = 1; i < M + 1; i++ )
    //       rowIndex_A[i] = rowIndex_A[i - 1] + 2;

    //Matrix B
    // for( i = 0; i < M; i++ )
    //       values_B[i] = i + 5;
    // for( i = 0; i < M; i++ )
    //       columns_B[i] = i % 5;
    // for( i = 0; i < M + 1; i++ )
    //       rowIndex_B[i] = i;

    //Vectors x and y
    for( i = 0; i < M; i++ )
    {
          x[i] = (float) i; //y[i] = 1.0;
    }

    /* Printing usable data */
    #ifdef DEBUG
    printf( "\n\n_______________Example program for MKL_SPARSE_D_SPMM_________________\n\n" );
    printf( " COMPUTE  A * x = y, where matrices are stored in CSR format\n" );

    float *correct_answer = (float *)mkl_malloc(sizeof(float) * M, ALIGN);
    printf( "\n MATRIX A:\nrow# : (value, column) (value, column)\n" );
    for( i = 0; i < M; i++ )
    {
        printf("row#%lld:", i + 1); fflush(0);
        float temp = 0;
        for( j = rowIndex_A[i]; j < rowIndex_A[i+1]; j++ )
        {
            printf(" (%5.5f, %6lld)", values_A[ii], columns_A[ii] ); fflush(0);
            temp += values_A[ii] * columns_A[ii];
            ii++;
        }
        correct_answer[i] = temp;
        printf( "\n" );
    }

    printf( "\n VECTOR X:\n" );
    for( i = 0; i < M; i++ )
    {
      printf("  %.5f ", x[i]);
        
    }
    printf( "\n" );
    #endif

    // ii = 0;
    // printf( "\n MATRIX B:\nrow# : (value, column)\n" );
    // for( i = 0; i < M; i++ )
    // {
    //     printf("row#%d:", i + 1); fflush(0);
    //     for( j = rowIndex_B[i]; j < rowIndex_B[i+1]; j++ )
    //     {
    //         printf(" (%5.0f, %6d)", values_B[ii], columns_B[ii] ); fflush(0);
    //         ii++;
    //     }
    //     printf( "\n" );
    // }
    // printf( "\n Check the resultant matrix C, using two scalar products\n" );
    // printf( " (values of these scalar products must match).\n" );

/* Prepare arrays, which are related to matrices.
   Create handles for matrices A and B stored in CSR format */
    CALL_AND_CHECK_STATUS(mkl_sparse_s_create_csr( &csrA, SPARSE_INDEX_BASE_ZERO, M, M, rowIndex_A, rowIndex_A+1, columns_A, values_A ),
                          "Error after MKL_SPARSE_D_CREATE_CSR, csrA \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &csrB, SPARSE_INDEX_BASE_ZERO, M, M, rowIndex_B, rowIndex_B+1, columns_B, values_B ),
    //                       "Error after MKL_SPARSE_D_CREATE_CSR, csrB \n");

/* Compute C = A * B  */
    // CALL_AND_CHECK_STATUS(mkl_sparse_spmm( SPARSE_OPERATION_NON_TRANSPOSE, csrA, csrB, &csrC ),
    //                       "Error after MKL_SPARSE_SPMM \n");

/* Analytic Routines for MKL_SPARSE_D_MV.
   HINTS: provides estimate of number and type of upcoming matrix-vector operations
   OPTIMIZE: analyze sparse matrix; choose proper kernels and workload balancing strategy */
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrA, SPARSE_OPERATION_TRANSPOSE,     descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrA \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrB, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrB \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrC, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrC \n");

    CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrA, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
                          "Error after MKL_SPARSE_SET_MV_HINT, csrA \n");

    CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrA ),
                          "Error after MKL_SPARSE_OPTIMIZE, csrA \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrB ),
    //                       "Error after MKL_SPARSE_OPTIMIZE, csrB \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrC ),
    //                       "Error after MKL_SPARSE_OPTIMIZE, csrC \n");

/* Execution Routines */
/* Step 1:
          Need to compute the following variables:
                 rslt_mv = C * x
                    left = <rslt_mv, y>              */
    CALL_AND_CHECK_STATUS(mkl_sparse_s_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr_type_gen, x, 0.0, rslt_mv ),
                          "Error after MKL_SPARSE_D_MV, csrC*x  \n");

    /* WARM UP USING ITERATIONS */
    for (int k = 0; k < WARM_COUNT; k++) {
        mkl_sparse_s_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr_type_gen, x, 0.0, rslt_mv );
    }

    /* TESTING */
    s_initial = dsecnd();
    for (int k = 0; k < LOOP_COUNT; k++) {
        mkl_sparse_s_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr_type_gen, x, 0.0, rslt_mv );
    }
    s_elapsed = (dsecnd() - s_initial) / LOOP_COUNT;

    printf (" == Sparse matrix dense vector multiplication using Intel(R) MKL dgemm completed == \n" " == at %.10f seconds == \n", (s_elapsed));
    printf ("%.10f\n", (s_elapsed));
    // left = cblas_ddot( M, rslt_mv, 1, y, 1 );

/* Step 2:
          Need to compute the following variables:
           rslt_mv       =     B * x
           rslt_mv_trans = (A)^t * y
                   right = <rslt_mv, rslt_mv_trans>  */

    // CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrB, descr_type_gen, x, 0.0, rslt_mv ),
    //                       "Error after MKL_SPARSE_D_MV, csrB*x  \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_TRANSPOSE,     1.0, csrA, descr_type_gen, y, 0.0, rslt_mv_trans),
    //                       "Error after MKL_SPARSE_D_MV, csrA*y  \n");
    // right = cblas_ddot( M, rslt_mv, 1, rslt_mv_trans, 1);

/* Step 3:
          Compare values obtained for left and right  */
    // residual = fabs(left - right)/(fabs(left)+1);

    // printf( "\n The difference between < C*x , y > and < B*x , (A^t)*y > = %g,\n", residual );
    // printf( " which means that MKL_SPARSE_D_SPMM arrived correct at a solution.\n" );
/* Printing OUTPUT DATA */
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( csrA, &indexing, &rows, &cols, &pointerB_A, &pointerE_A, &columns_A, &values_A ),
    //                       "Error after MKL_SPARSE_D_EXPORT_CSR  \n");

    // printf( "\n RESULTANT MATRIX A:\nrow# : (value, column) (value, column)\n" );
    // ii = 0;
    // for( i = 0; i < M; i++ )
    // {
    //     printf("row#%d:", i + 1); fflush(0);
    //     for( j = pointerB_A[i]; j < pointerE_A[i]; j++ )
    //     {
    //         printf(" (%5.0f, %6d)", values_C[ii], columns_C[ii] ); fflush(0);
    //         ii++;
    //     }
    //     printf( "\n" );
    // }
    // printf( "_____________________________________________________________________  \n" );

    #ifdef DEBUG
    printf( "\n RESULTANT VECTOR Y:\n" );
    for( i = 0; i < M; i++ )
    {
      printf("  %.5f ", rslt_mv[i]);
        
    }
    printf( "\n" );
    printf( "_____________________________________________________________________  \n" );
    #endif

    #ifdef DEBUG

    MKL_INT c;
    c = 0;

    printf( "\n RESULT CHECKING:\n" );

    for( i = 0; i < M; i++ )
    {
        if (correct_answer[i] - rslt_mv[i] < 0.001) {
            c++;
        }
        else {
            printf( "The %lld element is WRONG where y = %.5f and result = %.5f\n" , i, correct_answer[i], rslt_mv[i] );
        }
    }

    if (c == M) {
        printf( "*******************  RESULT CORRECT  *******************\n" );
    }
    else {
        printf( "*******************  RESULT WRONG!!  *******************\n" );
    }
    #endif

    /* Deallocate memory */
    memory_free:
    //Release matrix handle. Not necessary to deallocate arrays for which we don't allocate memory: values_C, columns_C, pointerB_C, and pointerE_C.
    //These arrays will be deallocated together with csrC structure.
      // if( mkl_sparse_destroy( csrC ) != SPARSE_STATUS_SUCCESS)
      // { printf(" Error after MKL_SPARSE_DESTROY, csrC \n");fflush(0); status = 1; }

    //Deallocate arrays for which we allocate memory ourselves.
      // mkl_free(rslt_mv_trans); mkl_free(rslt_mv); mkl_free(x); mkl_free(y);
    mkl_free(rslt_mv); mkl_free(x);

    //Release matrix handle and deallocate arrays for which we allocate memory ourselves.
    if( mkl_sparse_destroy( csrA ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, csrA \n");fflush(0); status = 1; }
    mkl_free(values_A); mkl_free(columns_A); mkl_free(rowIndex_A);

    // if( mkl_sparse_destroy( csrB ) != SPARSE_STATUS_SUCCESS)
    // { printf(" Error after MKL_SPARSE_DESTROY, csrB \n");fflush(0); status = 1; }
    // mkl_free(values_B); mkl_free(columns_B); mkl_free(rowIndex_B);


    return status;
}