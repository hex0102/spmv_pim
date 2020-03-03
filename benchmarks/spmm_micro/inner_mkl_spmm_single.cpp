#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <assert.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include "mkl.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include "matrix_market.hpp"
//#include "common.hpp"

#include <sys/stat.h>
#define ALIGN 128
#define NUMBER_RUNS 1

using namespace std;

//Calculate time elapsed in secs
inline void calc_elapsed (struct timeval t1, struct timeval t0, int number_of_runs) {
    printf ("%f,", ((float) (t1.tv_sec-t0.tv_sec) + (float) (t1.tv_usec-t0.tv_usec) / 1000000.0)/ number_of_runs);
}

void timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    printf("%s",asctime( localtime(&ltime) ) );
}

void readMatA (float ** vals, MKL_INT ** cols, MKL_INT ** rows, MKL_INT * m, MKL_INT * nnz); 
void readVecB (float ** vals, MKL_INT ** cols, MKL_INT ** rows, MKL_INT * m, MKL_INT * nnz); 
void set_1based_ind(MKL_INT *rowptr, MKL_INT *colidx, MKL_INT n, MKL_INT nnz);

int main() {


#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = 1;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)
    
    
    //Run on maximum threads available
    MKL_INT max_threads = mkl_get_max_threads();
    printf("Max threads = %lld\n", max_threads);
    mkl_set_num_threads(max_threads);

/* Declaration of values */
    float time = 0;
    int number_of_runs = 0;
    struct timeval t0_beg, t0_end;
    struct timeval t1_beg, t1_end;
    struct timeval t2_beg, t2_end;
    struct timeval t3_beg, t3_end;

    MKL_INT  i, j, status;
    float  *values_A = NULL, *values_B = NULL, *values_C = NULL;
    MKL_INT *columns_A = NULL, *columns_B = NULL, *columns_C = NULL;
    MKL_INT *rowPointer_A = NULL, *rowPointer_B = NULL,*rowPointer_C = NULL;

    //MKL_INT *colA, *rowA;
    //float *valA;
    MKL_INT m, nnz;
    MKL_INT nnz_out;

    //Input variables for MKL function csrmultcsr
    const char trans = 'N';
    MKL_INT sort;
    MKL_INT request1, request2;            
    MKL_INT nzmax = 0;

    MKL_INT job[6] = {2,1,1,0,nnz,0};

    MKL_INT info_coo, info1, info2;

    //sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
    MKL_INT rows, cols, *pointerB_C = NULL, *pointerE_C = NULL;
    ofstream outfile;

//Read Matrix A into COO format
    readMatA(&values_A, &columns_A, &rowPointer_A, &m, &nnz);
    readMatA(&values_B, &columns_B, &rowPointer_B, &m, &nnz);
    printf("N = %lld, NNZ = %lld\n",m, nnz);

/* Allocation of memory */
/*    values_A = (float *)mkl_malloc(sizeof(float) * nnz, ALIGN);
    columns_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * nnz, ALIGN);
    rowPointer_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (m + 1), ALIGN);
*/    
    rowPointer_C = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (m + 1), ALIGN);

/* Set values of the variables*/
    status = 0;                 
    //trans = 'N';                //C = A * B instead of A' * B
    sort = 8;                   //Sort resultant matrix
    request1 = 1;               //Computes only values of the rowPointer array                 
    request2 = 2;               //Computes only values of the column indices and values array given rowPointer array is computed
    nzmax = 0;
    
    //Convert input from COO to CSR
/*    mkl_dcsrcoo(job, &m, values_A, columns_A, rowPointer_A, &nnz, valA, rowA, colA, &info_coo);
    if (info_coo == 1) {
      printf("Error after MKL_DCSRCOO, csrA \n"); fflush(0);                 
      goto memory_free;
    }       
*/
    //mkl_dcsrcsc (job_csc , &m , values_A , columns_A , rowIndex_A , values_csc_A , rows_csc_A , colIndex_csc_A, &info_csc);

    //Compute number of non-zeroes in resultant matrix
    //---------------------------------------------------
    //printf("\nStart of multiplication\n");
    //Timing for Mult phase 1
    gettimeofday (&t0_beg, 0);
    mkl_scsrmultcsr(&trans, &request1, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_B, columns_B, rowPointer_B, values_C, columns_C, rowPointer_C, &nzmax, &info1);
    gettimeofday (&t0_end, 0);

    if (info1 != 0) {
		  printf("Error after MKL_SCSRMULTCSR 1, csrC \n"); fflush(0);
      goto memory_free;
	  }

    gettimeofday (&t2_beg, 0);
    nzmax = rowPointer_C[m] - 1;

    //Allocate memory for Resutant matrix
    values_C = (float *)mkl_malloc(sizeof(float) * nzmax, ALIGN);
    columns_C = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * nzmax, ALIGN);
    cout << "NZMAX for C: " << nzmax << endl;
    gettimeofday (&t2_end, 0);


    //Compute values and column indicies of resultant matrix 
    //------------------------------------------------------
    
    //Timing for Mult phase 2
    gettimeofday (&t1_beg, 0);
    mkl_scsrmultcsr(&trans, &request2, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_B, columns_B, rowPointer_B, values_C, columns_C, rowPointer_C, &nzmax, &info2);
    gettimeofday (&t1_end, 0);
    //printf("\nEnd of multiplication\n");

    if (info2 != 0) {
		  printf("Error after MKL_SCSRMULTCSR 2, csrC \n"); fflush(0);
      mkl_free(rowPointer_C);
      goto memory_free;
	  }

#if (TIME == 1)     
    printf ("Time elapsed for Mult Phase1, Mult Phase2, Malloc time:\n");
    calc_elapsed (t0_end, t0_beg, NUMBER_RUNS);
    calc_elapsed (t1_end, t1_beg, NUMBER_RUNS);
    calc_elapsed (t2_end, t2_beg, 1);
    printf ("\n");
#endif

#if (POWER == 1)
    
    //Get time for number of runs
    gettimeofday (&t0_beg, 0);
    mkl_scsrmultcsr(&trans, &request1, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_A, columns_A, rowPointer_A, values_C, columns_C, rowPointer_C, &nzmax, &info1);
    mkl_scsrmultcsr(&trans, &request2, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_A, columns_A, rowPointer_A, values_C, columns_C, rowPointer_C, &nzmax, &info2);
    gettimeofday (&t0_end, 0);
    
    time = ((float) (t0_end.tv_sec-t0_beg.tv_sec) + (float) (t0_end.tv_usec-t0_beg.tv_usec) / 1000000.0);
    number_of_runs = (int) 50/time; //Make sure power scripts is measuring for time <= 10s
    printf("Number of runs: %d\n", number_of_runs); 

    // system("/home/aporvaa/research/outerspace/gemma_mkl/power_script.sh &");
    system("/home/fengsy/Documents/power_scripts_cpu/power_scripts_vlsi19/power_script.sh &");
    for(int i = 0; i < number_of_runs; i++) {
        mkl_scsrmultcsr(&trans, &request1, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_A, columns_A, rowPointer_A, values_C, columns_C, rowPointer_C, &nzmax, &info1);
        mkl_scsrmultcsr(&trans, &request2, &sort, &m, &m, &m, values_A, columns_A, rowPointer_A, values_A, columns_A, rowPointer_A, values_C, columns_C, rowPointer_C, &nzmax, &info2);
    }
    

#endif



/* Deallocate memory */
memory_free:
 
    //Deallocate arrays for which we allocate memory.
    //mkl_free(valA); mkl_free(colA); mkl_free(rowA);
    mkl_free(values_A); mkl_free(columns_A); mkl_free(rowPointer_A);
    mkl_free(values_C); mkl_free(columns_C); mkl_free(rowPointer_C); 
    return status;
}

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
    set_1based_ind(*rowPtr, *cols, stl_A.size(), *nnz);
}


void set_1based_ind(MKL_INT *rowptr, MKL_INT *colidx, MKL_INT n, MKL_INT nnz)
{
    MKL_INT i;
    for(i=0; i <= n; i++)
        rowptr[i]++;
    for(i=0; i < nnz; i++)
        colidx[i]++;
}


/*    
    ifstream ifile;
    string line, token;
"s!O<Mouse>C"v!O<Mouse>C"y!O
    ifile.open ("./matA.txt");      //matA has matrix in Matrix Market format

    // First, get vals array
    getline (ifile, line);    //comment in file
    getline (ifile, line);    //comment in file
    getline (ifile, line);    //Row_size, col_size, NNZ
    istringstream ss (line);
    int size;
    ss >> size >> *m >> *nnz;
    
    *rows = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * *nnz, ALIGN);
    *cols = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * *nnz, ALIGN);
    *vals = (float *)mkl_malloc(sizeof(float) * *nnz, ALIGN);

    int i = 0;
    while (getline (ifile, line)) {
    MKL_INT a, b;
    istringstream ss0 (line);
        ss0 >> a >> b >> (*vals)[i];
        (*rows)[i] = a;
        (*cols)[i] = b;
        i++;
    }
    ifile.close ();
    */

/* Only for general matrices in mtx format
 
 void readMatA (float ** vals, MKL_INT ** cols, MKL_INT ** rows, MKL_INT * m, MKL_INT * nnz){
    ifstream ifile;
    string line, token;

    ifile.open ("./matA.txt");      //matA has matrix in Matrix Market format

    // First, get vals array
    getline (ifile, line);    //comment in file
    getline (ifile, line);    //comment in file
    getline (ifile, line);    //Row_size, col_size, NNZ
    istringstream ss (line);
    int size;
    ss >> size >> *m >> *nnz;
    
    *rows = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * *nnz, ALIGN);
    *cols = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * *nnz, ALIGN);
    *vals = (float *)mkl_malloc(sizeof(float) * *nnz, ALIGN);

    int i = 0;
    while (getline (ifile, line)) {
    MKL_INT a, b;
    istringstream ss0 (line);
        ss0 >> a >> b >> (*vals)[i];
        (*rows)[i] = a;
        (*cols)[i] = b;
        i++;
    }
    ifile.close ();
}
*/
