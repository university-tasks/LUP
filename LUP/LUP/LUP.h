//
//  LUP.h
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#ifndef LUP_h
#define LUP_h

#include <stdio.h>

void getL(double **A, double **L, int N);
void getU(double **A, double **U, int N);

void LUP(double **A, int *p, int N);
void checkLUP(double **A, int *p, double **L, double **U, int N, double **A1);

double determinant(double **A, int *p, int N);

double *solve(double **L, double **U, double *b, int *p, int N);
void checkSolve(double **A, double *x, double *b, int N);

void invert(double **L, double **U, int *p, int N, double**IA);
void checkInvert(double **A, double **IA, int N, int *p);

void condNum(double **A, double **IA, int N, double *cond);

#endif /* LUP_h */
