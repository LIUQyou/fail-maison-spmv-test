#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "mmio.h"
#include <math.h>
#include <libgen.h>
#include "spmv_csr.h"
#include "spmv_coo.h"
#include "spmv_dia.h"
#include "spmv_ell.h"
#include "conversions.h"
#include "utils.h"
#include "config.h"

#include <sys/ipc.h>
#include <sys/shm.h>

#define MB_1 (1024*1024)
#define GB_1 (1024*MB_1)
#define MB_8 (8*MB_1)

int M, N, entries, anz, i, j, *row, *col, *row_ptr, *colind, *offset, *indices, num_diags, stride, num_cols;
MYTYPE *x, *y, *val, *dia_data, *ell_data, *coo_val;
  
void init_hugetlb_seg(int **vec,int size)
{
  int shmid1;
  shmid1 = shmget(rand(), size*sizeof(int), SHM_HUGETLB
         | IPC_CREAT | SHM_R
         | SHM_W);
  if ( shmid1 < 0 ) {
    perror("shmget");
    exit(1);
  }
  printf("HugeTLB shmid: 0x%x\n", shmid1);
  *vec = shmat(shmid1, 0, 0);
  if (*vec == (int *)-1) {
    printf("Shared memory attach failure");
    shmctl(shmid1, IPC_RMID, NULL);
    exit(2);
  }
  printf("address is %p\n",*vec);
}

void init_hugetlb_seg_mytype(MYTYPE **vec,int size)
{
  int shmid1;
  shmid1 = shmget(rand(), size*sizeof(MYTYPE), SHM_HUGETLB
         | IPC_CREAT | SHM_R
         | SHM_W);
  if ( shmid1 < 0 ) {
    perror("shmget");
    exit(1);
  }
  printf("HugeTLB shmid: 0x%x\n", shmid1);
  *vec = shmat(shmid1, 0, 0);
  if (*vec == (MYTYPE *)-1) {
    printf("Shared memory attach failure");
    shmctl(shmid1, IPC_RMID, NULL);
    exit(2);
  }
  printf("address is %p\n",*vec);
}

// For hugebubbles-00010 on the m400, the MPKI is
// 26.6 for EXE_TIME 20.0 s
// 24.2 for EXE_TIME 10.0 s
// 21.6 for EXE_TIME 5.0 s
// For a trade-off, select 10.0 s for now
// #define EXE_TIME 10.0

int main(int argc, char* argv[])
{
  int ret_code, r, c, k;
  MM_typecode matcode;
  FILE *f, *g;
  double v;
  clock_t start, stop;
  double sum = 0, mean = 0, sd, variance = 0;
  unsigned long inner, inner_max = 1000000;
  double time_span=0;
  
  // A dumb way to indicate when starting the profiling
  // Maybe it can be implemented in another way
  FILE *fptr;
  fptr=fopen("re.txt","w");
  fprintf(fptr,"flag_end_spmv");
  fclose(fptr);

  #ifndef DENSE
  if(argc < 3){
    fprintf(stderr, "Usage: %s [martix-market-filename] [matrix format] [exe time]\n", argv[0]);
      exit(1);
  }
  else{
    if((f = fopen(argv[1], "r")) == NULL)
      exit(1);
  }
  #endif

  #ifndef DENSE
  anz = count_nnz(f);
  rewind(f);

  if (mm_read_banner(f, &matcode) != 0)
  {
    fprintf(stderr, "Could not process Matrix Market banner.\n");
    exit(1);
  }
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &entries)) !=0)
      exit(1);
  
  if(M > N) 
    N = M;
  else
    M = N;
  #endif

  //row = (int*)malloc(anz * sizeof(int));
  init_hugetlb_seg(&row,anz);
  if(row == NULL){
    fprintf(stderr, "couldn't allocate row using malloc\n");
    exit(1);
  }
  
  //col = (int*)malloc(anz * sizeof(int));
  init_hugetlb_seg(&col,anz);
  if(col == NULL){
    fprintf(stderr, "couldn't allocate col using malloc\n");
    exit(1);
  }
  
  //coo_val = (MYTYPE*)malloc(anz * sizeof(MYTYPE));
  init_hugetlb_seg_mytype(&coo_val,anz);
  if(coo_val == NULL){
    fprintf(stderr, "couldn't allocate val using malloc\n");
    exit(1);
  }
  
  row_ptr = (int*)calloc(N+1, sizeof(int));
  //init_hugetlb_seg(&row_ptr,N+1);
  if(row_ptr == NULL){
    fprintf(stderr, "couldn't allocate row_ptr using malloc\n");
    exit(1);
  }
  
  //colind = (int*)malloc(anz * sizeof(int));
  init_hugetlb_seg(&colind,anz);
  if(colind == NULL){
    fprintf(stderr, "couldn't allocate colind using malloc\n");
    exit(1);
  }
  
  //val = (MYTYPE*)malloc(anz * sizeof(MYTYPE));
  init_hugetlb_seg_mytype(&val,anz);
  if(val == NULL){
    fprintf(stderr, "couldn't allocate val using malloc\n");
    exit(1);
  }
  
  //x=(MYTYPE*) malloc( sizeof(MYTYPE)*N );
  init_hugetlb_seg_mytype(&x,N);
  if(x == NULL){
    fprintf(stderr, "couldn't allocate x using malloc\n");
    exit(1);
  }
  init_arr(N,x);
  
  //y=(MYTYPE*) calloc(N, sizeof(MYTYPE));
  init_hugetlb_seg_mytype(&y,N);
  if(y == NULL){
    fprintf(stderr, "couldn't allocate y using calloc\n");
    exit(1);
  }
  
  #ifndef DENSE
  k = 0;
  if(mm_is_symmetric(matcode)){
    if(!mm_is_pattern(matcode)){
      for (i=0; i<entries; i++){
        fscanf(f, "%d %d %lf\n", &r, &c, &v);
        if( v < 0 || v > 0){
          row[k] = r - 1;
          col[k]= c - 1;
          coo_val[k] = v;
          if(fpclassify(coo_val[k]) == FP_NAN){
            fprintf(stderr,"bad value : nan\n");
            exit(1);
          }
          if(fpclassify(coo_val[k]) == FP_INFINITE){
            fprintf(stderr,"bad value : infinite\n");
            exit(1);
          }
          if(fpclassify(coo_val[k]) == FP_SUBNORMAL){
            fprintf(stderr,"bad value : subnormal\n");
            coo_val[k] = 0.0;
          }
          if(r == c){
            k++;
          }
          else{
            row[k+1] = col[k]; 
            col[k+1] = row[k];
            coo_val[k+1] = v;
            if(fpclassify(coo_val[k+1]) == FP_SUBNORMAL){
              fprintf(stderr,"bad value : subnormal\n");
              coo_val[k+1] = 0.0;
            }
            k = k + 2;
          }
        }
      }
    }
    else{
      for (i=0; i<entries; i++){
        fscanf(f, "%d %d\n", &r, &c);
        row[k] = r - 1;
        col[k]= c - 1;
        coo_val[k] = 1.0;
        if(r == c){
          k++;
        }
        else{
          row[k+1] = col[k]; 
          col[k+1] = row[k];
          coo_val[k+1] = 1.0;
          k = k + 2;
        }
      }
    }
  }
  else {
    if(!mm_is_pattern(matcode)){
      for (i=0; i<entries; i++){
        fscanf(f, "%d %d %lf\n", &r, &c, &v);
        if( v < 0 || v > 0){
          row[k] = r - 1;
          col[k]= c - 1;
          coo_val[k] = v;
          if(fpclassify(coo_val[k]) == FP_NAN){
            fprintf(stderr,"bad value : nan\n");
            exit(1);
          }
          if(fpclassify(coo_val[k]) == FP_INFINITE){
            fprintf(stderr,"bad value : infinite\n");
            exit(1);
          }
          if(fpclassify(coo_val[k]) == FP_SUBNORMAL){
            fprintf(stderr,"bad value : subnormal\n");
            coo_val[k] = 0.0;
          }
          k++;
        }
      }
    }
    else {
      for (i=0; i<entries; i++){
        fscanf(f, "%d %d\n", &r, &c);
        row[i] = r - 1;
        col[i]= c - 1;
        coo_val[i] = 1.0;
      }
    }
  }

  quickSort(row, col, coo_val, 0, anz-1);
  #endif

  double EXE_TIME = atof(argv[3]);

  if(!string_compare(argv[2], "coo")){
    // A bit warm up
    zero_arr(N, y);
    spmv_coo(row, col, coo_val, anz, N, x, y);
    // flag to start the profiling
    printf("flag_start_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_start_spmv");
    fclose(fptr);

    // Start running for at least EXE_TIME
    start = clock();
    while (time_span < EXE_TIME){
      zero_arr(N, y);
      spmv_coo(row, col, coo_val, anz, N, x, y);
      stop = clock();
      time_span = ((double)(stop - start)/ CLOCKS_PER_SEC);
    }

    // Reset the flag in the file
    printf("flag_end_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_end_spmv");
    fclose(fptr);
  }

  if(!string_compare(argv[2], "csr")){
    coo_csr(anz, N, row, col, coo_val, row_ptr, colind, val);
    // A bit warm up
    zero_arr(N, y);
    spmv_csr(row_ptr, colind, val, N, x, y);
    // flag to start the profiling
    printf("flag_start_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_start_spmv");
    fclose(fptr);

    // Start running for at least EXE_TIME
    start = clock();
    while (time_span < EXE_TIME){
      zero_arr(N, y);
      spmv_csr(row_ptr, colind, val, N, x, y);
      stop = clock();
      time_span = ((double)(stop - start)/ CLOCKS_PER_SEC);
    }

    // Reset the flag in the file
    printf("flag_end_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_end_spmv");
    fclose(fptr);
  }

// For some matrices, dia may drop the eroor: too large
  if(!string_compare(argv[2], "dia")){
    coo_csr(anz, N, row, col, coo_val, row_ptr, colind, val);
    csr_dia(row_ptr, colind, val, &offset, &dia_data, N, &num_diags, &stride, anz);
    // A bit warm up
    zero_arr(N, y);
    spmv_dia(offset, dia_data, N, num_diags, stride, x, y);
    // flag to start the profiling
    printf("flag_start_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_start_spmv");
    fclose(fptr);

    // Start running for at least EXE_TIME
    start = clock();
    while (time_span < EXE_TIME){
      zero_arr(N, y);
      spmv_dia(offset, dia_data, N, num_diags, stride, x, y);
      stop = clock();
      time_span = ((double)(stop - start)/ CLOCKS_PER_SEC);
    }

    // Reset the flag in the file
    printf("flag_end_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_end_spmv");
    fclose(fptr);
    free(dia_data);
    free(offset);
  }
// For some matrices, ell may drop the eroor: too large
  if(!string_compare(argv[2], "ell")){
    coo_csr(anz, N, row, col, coo_val, row_ptr, colind, val);
    csr_ell(row_ptr, colind, val, &indices, &ell_data, N, &num_cols, anz);
    // A bit warm up
    zero_arr(N, y);
    spmv_ell(indices, ell_data, N, num_cols, x, y);
    // flag to start the profiling
    printf("flag_start_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_start_spmv");
    fclose(fptr);

    // Start running for at least EXE_TIME
    start = clock();
    while (time_span < EXE_TIME){
      zero_arr(N, y);
      spmv_ell(indices, ell_data, N, num_cols, x, y);
      stop = clock();
      time_span = ((double)(stop - start)/ CLOCKS_PER_SEC);
    }

    // Reset the flag in the file
    printf("flag_end_spmv\n");
    fptr=fopen("re.txt","w");
    fprintf(fptr,"flag_end_spmv");
    fclose(fptr);

    free(ell_data);
    free(indices);
  }

  #ifndef DENSE
  if (f !=stdin) 
    fclose(f);
  #endif

  //free(y);
  free(row_ptr);

  return 0;
}

