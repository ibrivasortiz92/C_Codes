#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define printf __mingw_printf
#define eps 1e-12

// PROTOTYPES
void print(int M, long double Matrix[]);
void mult(int M, long double Matrix[], long double X[], long double **B);
void matrix_mult(int M, long double Matrix1[], long double Matrix2[], long double **R);
int inv(int M, long double Matrix[], long double **X);


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
  long double *X = NULL, *IND = NULL, *R = NULL;

  // BODY
  if (inv(M, Matrix, &X)) {
    printf("INVERSE MATRIX:\n");
    print(M, X);

    printf("VERIFICATION:\n");
    matrix_mult(M, Matrix, X, &R);
    print(M, R);
  }
  
  // FREE MEMORY
  free(X); free(IND); free(R);

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

  *B = malloc(sizeof(long double) * M); assert(*B != NULL);
  long double sum;

  for (int j = 0; j < M; j++){
    sum = 0;
    for (int i = 0; i < M; i++){
      sum = sum + Matrix[M * j + i] * X[i];
    }
    (*B)[j] = sum;
  }

}


void matrix_mult(int M, long double Matrix1[], long double Matrix2[], long double **R){

  *R = malloc(sizeof(long double) * M * M); assert(*R != NULL);
  long double sum;

  for (int k = 0; k < M; k++) {
    for (int j = 0; j < M; j++){
      sum = 0;
      for (int i = 0; i < M; i++){
        sum = sum + Matrix1[M * j + i] * Matrix2[M * i + k];
      }
      (*R)[M * j + k] = sum;
    }
  }

}


int inv(int M, long double Matrix[], long double **X){

  int p;
  long double val, m;

  long double *temp = NULL;
  temp = malloc(sizeof(long double) * M); assert(temp != NULL);

  long double *IDEN = NULL, *Matrix2 = NULL;
  IDEN = malloc(sizeof(long double) * M * M); assert(IDEN != NULL);
  Matrix2 = malloc(sizeof(long double) * M * M); assert(Matrix2 != NULL);
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      Matrix2[M * j + i] = Matrix[M * j + i];
      if (i == j) IDEN[M * j + i] = 1.0;
      else IDEN[M * j + i] = 0.0;
    }
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

      temp[i] = IDEN[M * p + i];
      IDEN[M * p + i] = IDEN[M * j + i];
      IDEN[M * j + i] = temp[i];
    }

    if (fabs(Matrix2[M * j + j]) <= eps){
      printf("Singular Matrix\n\n");
      free(temp); free(IDEN); free(Matrix2);
      return 0;
    }

    // GAUSS ELIMINATION
    for (int jj = j + 1; jj < M; jj++){
      m = Matrix2[M * jj + j] / Matrix2[M * j + j];
      for(int i = 0; i < M; i++){
        Matrix2[M * jj + i] = Matrix2[M * jj + i] - m * Matrix2[M * j + i];
        IDEN[M * jj + i] = IDEN[M * jj + i] - m * IDEN[M * j + i];
      }
    }
  }

  // BACK SUBSTITUTION
  *X = malloc(sizeof(long double) * M * M); assert(*X != NULL);
  long double sum;
  for (int k = 0; k < M; k++){
    (*X)[M * (M - 1) + k] = IDEN[M * (M - 1) + k] / Matrix2[M * M - 1];
    for (int j = M - 2; j >= 0; j--){
      sum = 0;
      for (int i = j; i < M - 1; i++){
        sum = sum + Matrix2[M * j + i + 1] * (*X)[M * (i + 1) + k];
      }
      (*X)[M * j + k] = (IDEN[M * j + k] - sum) / Matrix2[M * j + j];
      if (fabs((*X)[M * j + k]) < eps) (*X)[M * j + k] = 0.0;
    }
  }
  
  free(temp); free(IDEN); free(Matrix2);

  return 1;
}
