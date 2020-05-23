//
//  LUP.c
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//
#define myExp 0.00000001
#include "LUP.h"
#include "MatrOperations.h"
#include "VectOperations.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

void getL(double **A, double **L, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j < i) {
                L[i][j] = A[i][j];
            } else if (j == i) {
                L[i][j] = 1;
            } else {
                L[i][j] = 0;
            }
        }
    }
}
void getU(double **A, double **U, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j >= i) {
                U[i][j] = A[i][j];
            } else {
                U[i][j] = 0;
            }
        }
    }
}

void LUP(double **A, int *p, int N) {

    int maxIndex = -1;
    double maxValue = 0;
        
    for (int i = 0; i < N; ++i) {
        maxValue = 0.0;
        maxIndex = i;

        for (int k = i; k < N; ++k)
            if ((fabs(A[k][i])) > maxValue) {
                maxValue = fabs(A[k][i]);
                maxIndex = k;
            }
            
        if (maxIndex != i) {
            
            int j = p[i];
            p[i] = p[maxIndex];
            p[maxIndex] = j;

            double *r = A[i];
            A[i] = A[maxIndex];
            A[maxIndex] = r;

            p[N]++;
        }

        for (int j = i + 1; j < N; ++j) {
            A[j][i] /= A[i][i];

            for (int k = i + 1; k < N; ++k) {
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }
    printf("LUP-decomposition:\nA:\n");
    showMatr(A, N, N);
    printf("p: ");
    showIVect(p, N);
}


void checkLUP(double **A, int *p, double **L, double **U, int N, double **A1) {
    
    double **T = (double **)malloc(N*sizeof(double *));
    double **B = (double **)malloc(N*sizeof(double *));
    double **C = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        T[i] = (double *)malloc(N*sizeof(double));
        B[i] = (double *)malloc(N*sizeof(double));
        C[i] = (double *)malloc(N*sizeof(double));
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j == p[i]) {
                T[i][j] = 1;
            } else {
                T[i][j] = 0;
            }
        }
    }
    
    getL(A1, L, N);
    getU(A1, U, N);
    
    B = multMatr(L, U, N);
    C = multMatr(T, A, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (fabs(B[i][j]-C[i][j])>myExp) {
                printf("Error!\n");
                return;
            }
        }
    }
    printf("Complete!\n");
}

double determinant(double **A, int *p, int N) {

    double det = 1;

    for (int i = 0; i < N; i++) {
        if (fabs(A[i][i])<myExp) {
            A[i][i] = 0;
        }
         det *= A[i][i];
    }
       
    if (((p[N] - N) % 2 != 0) && (det != 0)) {
         return -det;
    } else {
        return det;
    }
}

double *solve(double **L, double **U, double *b, int *p, int N) {
    
    double *x = (double*)malloc(sizeof(double) * N);
    double *y = (double*)malloc(sizeof(double) * N);
    double *temp = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; ++i) {
        temp[i] = b[p[i]];
    }
    b = temp;
    
    for (int i = 0; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j]*y[j];
        }
        y[i] = b[i] - sum;
    }
    
    showDVect(y, N);
    
    for (int i = N-1; i > -1; --i) {
        
        if (fabs(U[i][i]) < myExp) {
            x[i] = 1;
            printf("%s%i%s%lf%s","!!x_",i," = ",x[i],"\n");
            continue;
        }
        
        double sum = 0;
        for (int j = i; j < N; ++j) {
            sum += U[i][j]*x[j];
        }
        
        x[i] = (y[i] - sum)/U[i][i];
        printf("%s%i%s%lf%s","x_",i," = ",x[i],"\n");
    }
    
    
    return x;
}
void checkSolve(double **A, double *x, double *b, int N) {
    
    double *y = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; ++i) {
        y[i] = 0;
        double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j]*x[j];
        }
        y[i] = b[i] - sum;
        if (fabs(y[i]) > myExp) {
            printf("Error solving SLAE\n");
            return;
        }
    }
    printf("THE TEST IS SUCCESSFUL\n");
}

void invert(double **L, double **U, int *p, int N, double**IA) {
    double *b = (double*)malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        
        for (int j = 0; j < N; ++j) {
            if (j == i) {
                b[j] = 1;
            } else {
                b[j] = 0;
            }
        }

        IA[i] = solve(L, U, b, p, N);
    }
    
    double **T = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        T[i] = (double *)malloc(N*sizeof(double));
    }
    transpose(IA, N);
    printf("INVERSE MATRIX:\n");
    showMatr(IA, N, N);
}
void checkInvert(double **A, double **IA, int N, int *p) {

    printf("A*Aˆ-1:\n");
    showMatr(multMatr(A, IA, N), N, N);
    printf("Aˆ-1*A:\n");
    showMatr(multMatr(IA, A, N), N, N);
}

void condNum(double **A, double **IA, int N, double *cond) {
    
    double normR = 0;
    double normIR = 0;
    for (int i = 0; i < N; ++i) {
        double sum = 0;
        double sumI = 0;
       
        for (int j = 0; j < N; ++j) {
            sum += fabs(A[i][j]);
            sumI += fabs(IA[i][j]);
        }
        if (normR <= sum) {
            normR = sum;
        }
        if (normIR <= sumI) {
            normIR = sumI;
        }
    }
    *cond = normR*normIR;
    printf("CONDITION NUMBER (NORM INF)");
    printf("%f%s", *cond,"\n");
}
