//
//  VectOperations.c
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#include "VectOperations.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void showDVect(double *v, int N) {
    for (int i = 0; i < N; ++i) {
        printf("%f ",v[i]);
    }
    printf("\n");
}

void showIVect(int *v, int N) {
    for (int i = 0; i < N; ++i) {
        printf("%i ",v[i]);
    }
    printf("\n");
}

void generateVec(double *b, int N) {
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        b[i] = rand()%100;
    }
    printf("Generated vector:\n");
    showDVect(b, N);
}

double scalMult(double *a, double *b, int N) {
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res+= a[i]*b[i];
    }
    return res;
}

double *NumMultVect(double a, double *V, int N) {
  double *result = (double*)malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        result[i] = a*V[i];
    }
  return result;
}

double *DifVect(double *a, double *b, int N) {
    double *res = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}

double *SumVect(double *a, double *b, int N) {
    double *res = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

double *orthonormVect(double *v, int N) {
   double *e = (double*)malloc(sizeof(double) * N);
   double c = 1/sqrt(scalMult(v, v, N));
   e = NumMultVect(c, v, N);
   return e;
}

double *proj(double *u, double *v, int N) {
  double *result = (double*)malloc(sizeof(double) * N);
  double c = scalMult(u, v, N)/scalMult(u, u, N);
  result = NumMultVect(c, u, N);
  return result;
}
