//
//  Jacobi_Zeidel.h
//  LUP
//
//  Created by Артем Белов on 20.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#ifndef Jacobi_Seidel_h
#define Jacobi_Seidel_h

#include <stdio.h>

void runJacobi(double **A, double *x, double *b, int N, int *apost_k1);
void runSeidel(double **A, double *x, double *b, int N, int *apost_k2);
void JacobiSeidel(int N);
void check(int N, double **A, double *x, double *b);
void transition(double **A, double **Res, double *x, double *temp_b, double *b, int N);
int aprioriEst(double **Res, double *temp_b, int N);
#endif /* Jacobi_Seidel_h */
