//
//  MatrOperations.c
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//
#define myExp 0.00000001
#include "MatrOperations.h"
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

void showMatr(double **A, int N, int M) {
    printf("\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            printf("%lf%s",A[i][j]," ");
        }
        printf("\n");
    }
    printf("\n");
}

void generateMatr(double **A, int N) {
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = rand()%10000;
        }
    }
    printf("Generated matrix:\n");
    showMatr(A, N, N);
}

void generateSingularMatr(double **A, int N) {
    srand(time(0));
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
         A[i][j] = rand()%10000;
       }
    }
    for (int i = 0; i < N; ++i) {
         for (int j = 1; j < N; ++j) {
            A[i][j] = (2-j) * A[i][0] + (j-1) * A[i][2];
       }
    }

}

double **multMatr(double **A, double **B, int N) {
    
    double **R = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        R[i] = (double *)malloc(N*sizeof(double));
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            R[i][j] = 0;
            for (int k = 0; k < N; ++k) {
                R[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return R;
}

void transpose(double **A, int N) {
    
    double **T = (double**)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        T[i] = (double *)malloc(N*sizeof(double));
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            T[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = T[j][i];
}

void swap(double **mat, int row1,int row2, int col) {
    for(int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}

double rankMatr(double** A, int N, int M) {

    double **Temp = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        Temp[i] = (double *)malloc(M*sizeof(double));
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            Temp[i][j] = A[i][j];
        }
    }
    
    double rank = M;
    bool *arr = (bool*)malloc(sizeof(bool) * N);

    for (int i = 0; i < N; ++i) arr[i] = 0;

    for (int i = 0; i < M; ++i) {

        int j;
        for (j = 0; j < N; ++j)
            if (!arr[j] && fabs(Temp[j][i]) > myExp)
                break;
        if (j == N)
            --rank;
        else {
            arr[j] = true;
            for (int p = i + 1; p < M; ++p)
                Temp[j][p] /= Temp[j][i];
            for (int k = 0; k < N; ++k)
                if (k != j && fabs(Temp[k][i]) > myExp)
                    for (int p = i + 1; p < M; ++p)
                        Temp[k][p] -= Temp[j][p] * Temp[k][i];
        }
    }

    return rank;
}

double **numMultMatr(double a, double **A, int N, int M) {

double **Temp = (double **)malloc(N*sizeof(double *));
for(int i = 0; i < N; i++) {
        Temp[i] = (double *)malloc(M*sizeof(double));
}

  for(int i = 0; i < N; ++i) {
    for(int j = 0; j< M; ++j) {
        Temp[i][j] = a*A[i][j];
    }
  }
  return Temp;
}

double **matrPlusMatr(double **A, double **B, int N, int M) {
    double **Temp = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
            Temp[i] = (double *)malloc(M*sizeof(double));
    }
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j< M; ++j) {
          Temp[i][j] = A[i][j]+B[i][j];
      }
    }
    return Temp;
    
}


double **generateDiagPredomMatr(int N) {
    
    double **A = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N*sizeof(double));
    }
    
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                A[i][j] = rand()%1000000;
            } else {
                A[i][j] = rand()%1000;
            }
        }
    }
    return A;
}

double **generatePosDefMatr(int N) {
    
    double **NotDiagMatr = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        NotDiagMatr[i] = (double *)malloc(N*sizeof(double));
    }
    
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
           NotDiagMatr[i][j] = rand()%100;
        }
    }
    double **Temp1 = (double **)malloc(N*sizeof(double *));
       for(int i = 0; i < N; i++) {
           Temp1[i] = (double *)malloc(N*sizeof(double));
       }
    double **Temp2 = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        Temp2[i] = (double *)malloc(N*sizeof(double));
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Temp1[i][j] = NotDiagMatr[i][j];
            Temp2[i][j] = NotDiagMatr[i][j];
        }
    }
    transpose(Temp2, N);
    
    NotDiagMatr = multMatr(Temp1, Temp2, N);
    
    return NotDiagMatr;
}
