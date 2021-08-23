#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define printf __mingw_printf
#define eps 1e-12

//////////////////////////////////////////////////////////////////////////////////// PROTOTYPES

////////////////////////////////////////////////////////////////// PRINT VECTOR
void print_vector(int M, long double Vector[]);

////////////////////////////////////////////////////////////////// PRINT MATRIX
void print_matrix(int M, long double Matrix[]);

/////////////////////////////////////////////////////////////////// ZERO MATRIX
long double* zeros(int M);

/////////////////////////////////////////////////////////////// IDENTITY MATRIX
long double* eye(int M);

////////////////////////////////////////////////////////////////// EQUAL MATRIX
long double* equal(int M, long double Matrix[], long double **OUTPUT);

/////////////////////////////////////////////////////////////// NEGATIVE MATRIX
long double* neg(int M, long double Matrix[], long double **OUTPUT);

//////////////////////////////////////////////////////////////////// MATRIX SUM
long double* matrix_sum(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT);

/////////////////////////////////////////////////////// MATRIX MULTIPLICATION 1
long double* matrix_mult1(int M, long double Matrix[], long double X[], long double **OUTPUT);

/////////////////////////////////////////////////////// MATRIX MULTIPLICATION 2
long double* matrix_mult2(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT);

////////////////////////////////////////////////////////////// MATRIX INVERSION
long double* inv(int M, long double Matrix[], long double **OUTPUT);

///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////// MAIN
int main(){

  int M = 4;
  long double *R = eye(M), *RESULT = NULL;
  
  // EQUALITY
  RESULT = equal(M, R, &RESULT);
  print_matrix(M, RESULT);

  // NEGATIVE MATRIX
  RESULT = neg(M, R, &RESULT);
  print_matrix(M, RESULT);

  // MATRIX SUM
  long double *OP1 = NULL, *OP2 = NULL;
  OP1 = equal(M, R, &OP1);
  OP2 = equal(M, R, &OP2);
  RESULT = matrix_sum(M, OP1, OP2, &RESULT);
  print_matrix(M, RESULT);

  // MATRIX MULTIPLICATION 1
  long double Matrix[] = {1, 2, 1, 1,
                          5, 3, 1, 2,
                          2, 4, 9, 7,
                          1, 3, 2, 1};
  long double X[] = {1, 1, 1, 1};
  long double *RVECT = NULL;
  RVECT = matrix_mult1(M, Matrix, X, &RVECT);
  print_vector(M, RVECT);

  // MATRIX INVERSION
  RESULT = inv(M, Matrix, &RESULT);
  if (RESULT == NULL) printf("Singular Matrix\n\n");
  else {
    print_matrix(M, RESULT);

    long double *VER = NULL;
    printf("Verification:\n");
    VER = matrix_mult2(M, RESULT, Matrix, &VER);
    print_matrix(M, VER);

    free(VER);
  }
  
  free(R); free(RESULT); free(OP1); free(OP2);
  free(RVECT);

  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////// IMPLEMENTATION

////////////////////////////////////////////////////////////////// PRINT VECTOR
void print_vector(int M, long double Vector[]){

  for (int i = 0; i < M; i++){
    printf("%.5Lf\n", Vector[i]);
  }
  printf("\n");

}

////////////////////////////////////////////////////////////////// PRINT MATRIX
void print_matrix(int M, long double Matrix[]){

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      printf("%.5Lf ", Matrix[M * j + i]);
    }
    printf("\n");
  }
  printf("\n");

}

/////////////////////////////////////////////////////////////////// ZERO MATRIX
long double* zeros(int M){

  long double *R = NULL;
  R = malloc(sizeof(long double) * M * M); assert(R != NULL);

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      R[M * j + i] = 0.0;
    }
  }

  return R;

}

/////////////////////////////////////////////////////////////// IDENTITY MATRIX
long double* eye(int M){

  long double *R = NULL;
  R = malloc(sizeof(long double) * M * M); assert(R != NULL);

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      if (j == i) R[M * j + i] = 1.0;
      else R[M * j + i] = 0.0;
    }
  }

  return R;

}

////////////////////////////////////////////////////////////////// EQUAL MATRIX
long double* equal(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M); assert(*OUTPUT != NULL);
  }
  
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      (*OUTPUT)[M * j + i] = Matrix[M * j + i];
    }
  }

  return *OUTPUT;

}

/////////////////////////////////////////////////////////////// NEGATIVE MATRIX
long double* neg(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M); assert(*OUTPUT != NULL);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      if (fabs(Matrix[M * j + i]) < eps) (*OUTPUT)[M * j + i] = 0.0;
      else (*OUTPUT)[M * j + i] = - Matrix[M * j + i];
    }
  }

  return *OUTPUT;

}

//////////////////////////////////////////////////////////////////// VECTOR SUM
long double* vector_sum(int M, long double Vector1[], long double Vector2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M); assert(*OUTPUT != NULL);
  }

  for (int i = 0; i < M; i++){
    (*OUTPUT)[i] = Vector1[i] + Vector2[i];
  }

  return *OUTPUT;

}

//////////////////////////////////////////////////////////////////// MATRIX SUM
long double* matrix_sum(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M); assert(*OUTPUT != NULL);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      (*OUTPUT)[M * j + i] = Matrix1[M * j + i] + Matrix2[M * j + i];
    }
  }

  return *OUTPUT;

}

/////////////////////////////////////////////////////// MATRIX MULTIPLICATION 1
long double* matrix_mult1(int M, long double Matrix[], long double X[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M); assert(*OUTPUT != NULL);
  }

  long double sum;

  for (int j = 0; j < M; j++){
    sum = 0.0;
    for (int i = 0; i < M; i++){
      sum = sum + Matrix[M * j + i] * X[i];
    }
    (*OUTPUT)[j] = sum;
  }

  return *OUTPUT;

}

/////////////////////////////////////////////////////// MATRIX MULTIPLICATION 2
long double* matrix_mult2(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M); assert(*OUTPUT != NULL);
  }

  long double sum;

  for (int k = 0; k < M; k++) {
    for (int j = 0; j < M; j++){
      sum = 0;
      for (int i = 0; i < M; i++){
        sum = sum + Matrix1[M * j + i] * Matrix2[M * i + k];
      }
      (*OUTPUT)[M * j + k] = sum;
    }
  }

  return *OUTPUT;

}

////////////////////////////////////////////////////////////// MATRIX INVERSION
long double* inv(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M); assert(*OUTPUT != NULL);
  }

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
      free(temp); free(IDEN); free(Matrix2);
      free(*OUTPUT);
      return NULL;
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
  long double sum;
  for (int k = 0; k < M; k++){
    (*OUTPUT)[M * (M - 1) + k] = IDEN[M * (M - 1) + k] / Matrix2[M * M - 1];
    for (int j = M - 2; j >= 0; j--){
      sum = 0;
      for (int i = j; i < M - 1; i++){
        sum = sum + Matrix2[M * j + i + 1] * (*OUTPUT)[M * (i + 1) + k];
      }
      (*OUTPUT)[M * j + k] = (IDEN[M * j + k] - sum) / Matrix2[M * j + j];
      if (fabs((*OUTPUT)[M * j + k]) < eps) (*OUTPUT)[M * j + k] = 0.0;
    }
  }
  
  free(temp); free(IDEN); free(Matrix2);

  return *OUTPUT;
}

///////////////////////////////////////////////////////////////////////////////////////////////
