//
//  Jacobi_Zeidel.c
//  LUP
//
//  Created by Артем Белов on 20.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//
#define myExp 0.00000001

#include <stdlib.h>
#include <math.h>

#include "Jacobi_Seidel.h"
#include "MatrOperations.h"
#include "VectOperations.h"

double thirdNorm(double **A, int N) {
    int i, j;
    double sum = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            sum += (A[i][j] * A[i][j]);
        }
    }
    sum = sqrt(sum);
    return sum;
}

double vecNorm(double *A, int N) {
    int i;
    double sum = 0;
    for (i = 0; i < N; i++) {
            sum += (A[i] * A[i]);
    }
    sum = sqrt(sum);
    return sum;
}

void runJacobi(double **A, double *x, double *b, int N, int *apost_k1) {
 
   double *dx = (double*) calloc(N,sizeof(double));
   double *y = (double*) calloc(N,sizeof(double));

   double sum;
   apost_k1[0] = 0;
   do {
      apost_k1[0]++;
      sum = 0.0;
      for(int i=0; i<N; ++i) {
         dx[i] = b[i];
          for(int j=0; j<N; ++j) {
              dx[i] -= A[i][j]*x[j];
          }
         dx[i] /= A[i][i];
         y[i] += dx[i];
         sum += ( (dx[i] >= 0.0) ? dx[i] : -dx[i]);
      }
       for(int i = 0; i < N; i++) {
           x[i] = y[i];
       }
   } while (sum >= myExp);
   
}

void runSeidel(double **A, double *x, double *b, int N, int *apost_k2) {
    
    double sum;
    double g;
    apost_k2[0] = 0;
    
    do {
      apost_k2[0]++;
      sum = 0;
      for (int i = 0; i < N; i++) {
                  g = b[i];
                  for (int j = 0; j < N; j++) {
                      g = g + A[i][j] * x[j];
                  }
                  sum += (x[i] - g) * (x[i] - g);
                  x[i] = g;
              }
    } while (sqrt(sum) >= myExp * (1 - thirdNorm(A, N)) / thirdNorm(A, N));
}

void transition(double **A, double **Res, double *x, double *temp_b, double *b, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Res[i][j] = - A[i][j]/A[i][i];
        }
        temp_b[i] = b[i]/A[i][i];
        x[i] = temp_b[i];
        Res[i][i] = 0;
    }
}

int aprioriEst(double **Res, double *temp_b, int N) {
    double alpha = thirdNorm(Res, N);
    double beta = vecNorm(temp_b, N);
    
    if (alpha > 1) {
        printf("CAN'T SOLVE USING THESE METHODS BECAUSE NORM IS BIGGER THAN 1: %lf\n",alpha);
        return 0;
    }
    
    return (int) (log10(myExp) - log10(beta) + log10(1 - alpha))/log10(alpha) + 1;
}



void JacobiSeidel(int N) {
    
    int *apost_k1 = (int*)malloc(sizeof(int) * 1);
    int *apost_k2 = (int*)malloc(sizeof(int) * 1);
    
    double *b = (double*)malloc(sizeof(double) * N);
    double *temp_b = (double*)malloc(sizeof(double) * N);
    
    double *x = (double*)malloc(sizeof(double) * N);
    double *x2 = (double*)malloc(sizeof(double) * N);
    
    double **DiagMatr = (double **)malloc(N*sizeof(double *));
    double **Temp = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        DiagMatr[i] = (double *)malloc(N*sizeof(double));
        Temp[i] = (double *)malloc(N*sizeof(double));
    }
    
    double **NotDiagMatr = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        NotDiagMatr[i] = (double *)malloc(N*sizeof(double));
    }
//  МАТРИЦЫ С ДИАГОНАЛЬНЫМ ПРЕОБЛАДАНИЕМ
    DiagMatr = generateDiagPredomMatr(N);
    printf("Matrix with diagonal predominance:\n");
    showMatr(DiagMatr, N, N);
    
//  ПОЛОЖИТЕЛЬНО ОПРЕДЕЛЕННАЯ БЕЗ ДИАГОНАЛЬНОГО ПРЕОБЛАДАНИЯ
    NotDiagMatr = generatePosDefMatr(N);
    printf("Positive-definitive matrix without diagonal prdominance:\n");
    showMatr(NotDiagMatr, N, N);
    
    int action;
    printf("generate b?\n");
    printf("1)yes!\n");
    printf("2)enter!\n");
    printf("3)exit!\n");
    printf(">>>>>>>>> ");
    scanf("%i", &action);
    switch (action) {
        case 1:
            generateVec(b, N);
            break;
        case 2:
            printf("Введите значения b: ");
            for (int i = 0; i < N; ++i) {
                scanf("%lf", &x[i]);
            }
            break;
        default:
            return;;
    }
    
    
//  МАТРИЦА ПЕРЕХОДА ДЛЯ МАТРИЦЫ С ДИАГОНАЛЬНЫМ ПРЕОБЛАДАНИЕМ
    transition(DiagMatr, Temp, x, temp_b, b, N);
    
//  АПРИОРНАЯ ОЦЕНКА
    int k = aprioriEst(Temp, temp_b, N);
    
    printf("Priori estimation for the matrix with diagonal predominance: %i\n", k);
    runJacobi(DiagMatr, x, b, N, apost_k1);
    runSeidel(Temp, x2, temp_b, N, apost_k2);
//  ПОСТЕРИОРНАЯ ОЦЕНКА
    printf("Posteriori estimation when solving with Jacobi method: %i\n", apost_k1[0]);
    printf("Posteriori estimation when solving with Seidel method: %i\n", apost_k2[0]);
    printf("\n");
    
    printf("Solution using Jacobi method: \n");
    showDVect(x, N);
    check(N, DiagMatr, x, b);
    
    printf("\nSolution using Seidel method\n");
    showDVect(x2, N);
    check(N, DiagMatr, x2, b);
    

printf("\n\nSOLVING POSITIVE DEFINITIVE MATRIX WITHOUT DIAGONAL PREDOMINANCE\n\n");
    
//  МАТРИЦА ПЕРЕХОДА ДЛЯ ПОЛОЖИТЕЛЬНО ОПРЕДЕЛЕННОЙ МАТРИЦЫ БЕЗ ДИАГОНАЛЬНОГО ПРЕОБЛАДАНИЯ
    
    double *temp_b2 = (double*)malloc(sizeof(double) * N);
    double **Temp2 = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        Temp2[i] = (double *)malloc(N*sizeof(double));
    }
    
    transition(NotDiagMatr, Temp2, x, temp_b2, b, N);
    
    showMatr(Temp2, N, N);
    showDVect(temp_b2, N);

   
//  АПРИОРНАЯ ОЦЕНКА
    k = aprioriEst(Temp2, temp_b2, N);
    if (k == 0) {
        return;
    }
    
    printf("Priori estimation for the matrix with diagonal predominance: %i\n", k);
    runJacobi(NotDiagMatr, x, b, N, apost_k1);
    runSeidel(Temp2, x2, temp_b2, N, apost_k2);
    
//  ПОСТЕРИОРНАЯ ОЦЕНКА
    printf("Posteriori estimation when solving with Jacobi method: %i\n", apost_k1[0]);
    printf("Posteriori estimation when solving with Seidel method: %i\n", apost_k2[0]);
    printf("\n");
    
    printf("Solution using Jacobi method: \n");
    showDVect(x, N);
    check(N, NotDiagMatr, x, b);
    
    printf("\nSolution using Seidel method\n");
    showDVect(x2, N);
    check(N, NotDiagMatr, x2, b);
    
}

void check(int N, double **A, double *x, double *b) {
    
    double *u = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j]*x[j];
        }
        u[i] = sum;
    }
    printf("Проверка:\n");
    showDVect(u, N);
}

