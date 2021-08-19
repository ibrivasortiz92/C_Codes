#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define printf __mingw_printf
#define eps 1e-12

// PROTOTYPES
void print(int M, long double Matrix[]);
void mult(int M, long double Matrix[], long double X[], long double **B);
void backsub(int M, long double Matrix[], long double **X, long double B[]);
int gelim(int M, long double Matrix[], long double **X, long double B[]);


int main(){

  // INITIALIZATION
  int M = 4;
  long double Matrix[] = {1, 2, 1, 1,
                          5, 3, 1, 2,
                          2, 4, 9, 7,
                          1, 3, 2, 1};
  // long double Matrix[] = {1, 2, 3, 4,
  //                         5, 6, 7, 8,
  //                         9, 10, 11, 12,
  //                         13, 14, 15, 16};
  long double B[] = {1, 1, 1, 1};
  long double *X = NULL;
  long double *IND = NULL;

  // BODY
  if (gelim(M, Matrix, &X, B) == 1){
    mult(M, Matrix, X, &IND);

    printf("SOLUTION:\n");
    for (int i = 0; i < M; i++){
      printf("%.5Lf\n", X[i]);
    }
    printf("\n");
  }

  // FREE MEMORY
  free(X); free(IND);

  return 0;
}


// IMPLEMENTATION

void print(int M, long double Matrix[]){

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      printf("%.5Lf ", Matrix[M * j + i]);
    }
    printf("\n");
  }
  printf("\n");

}


void mult(int M, long double Matrix[], long double X[], long double **B){

  printf("VERIFICATION:\n");

  *B = malloc(sizeof(long double) * M); assert(*B != NULL);
  long double sum;

  for (int j = 0; j < M; j++){
    sum = 0;
    for (int i = 0; i < M; i++){
      sum = sum + Matrix[M * j + i] * X[i];
    }
    (*B)[j] = sum;
    printf("%.5Lf\n", (*B)[j]);
  }
  printf("\n");

}


void backsub(int M, long double Matrix[], long double **X, long double B[]){

  *X = malloc(sizeof(long double) * M); assert(*X != NULL);

  long double sum;

  (*X)[M - 1] = B[M -1] / Matrix[M * M - 1];
  for (int k = M - 2; k >= 0; k--){
    sum = 0;
    for (int i = k; i < M - 1; i++){
      sum = sum + Matrix[M * k + i + 1] * (*X)[i + 1];
    }
    (*X)[k] = (B[k] - sum) / Matrix[M * (k) + k];
  }

}


int gelim(int M, long double Matrix[], long double **X, long double B[]){

  int p;
  long double val, m;

  long double tval, *temp = NULL;
  temp = malloc(sizeof(long double) * M); assert(temp != NULL);

  long double *B2 = NULL, *Matrix2 = NULL;
  B2 = malloc(sizeof(long double) * M); assert(B2 != NULL);
  Matrix2 = malloc(sizeof(long double) * M * M); assert(Matrix2 != NULL);
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      Matrix2[M * j + i] = Matrix[M * j + i];
    }
    B2[j] = B[j];
  }

  for (int j = 0; j < M; j++){

    // PARTIAL PIVOTING
    val = 0;
    for( int k = j; k < M; k++){
      if (val < fabs(Matrix2[M * k + j])){
        val = fabs(Matrix2[M * k + j]);
        p = k;
      }
    }

    for (int i = 0; i < M; i++){
      temp[i] = Matrix2[M * p + i];
      Matrix2[M * p + i] = Matrix2[M * j + i];
      Matrix2[M * j + i] = temp[i];
    }
    tval = B2[p]; B2[p] = B2[j]; B2[j] = tval;

    if (fabs(Matrix2[M * j + j]) <= eps){
      printf("Singular Matrix\n\n");
      free(temp); free(B2); free(Matrix2);
      return 0;
    }

    // TRIANGULIZATION
    for (int jj = j + 1; jj < M; jj++){
      m = Matrix2[M * jj + j] / Matrix2[M * j + j];
      for(int i = j; i < M; i++){
        Matrix2[M * jj + i] = Matrix2[M * jj + i] - m * Matrix2[M * j + i];
      }
      B2[jj] = B2[jj] - m * B2[j];
    }
  }

  // BACK SUBSTITUTION
  backsub(M, Matrix2, &(*X), B2);

  free(temp); free(B2); free(Matrix2);

  return 1;
}
