//
//  MatrOperations.h
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#ifndef MatrOperations_h
#define MatrOperations_h

#include <stdio.h>

void showMatr(double **A, int N, int M);
void generateMatr(double **A, int N);
void generateSingularMatr(double **A, int N);
double rankMatr(double** A, int N, int M);
double **multMatr(double **A, double **B, int N);
void transpose(double **A, int N);
void swap(double **mat, int row1,int row2, int col);
double **numMultMatr(double a, double **A, int N, int M);
double **matrPlusMatr(double **A, double **B, int N, int M);

double **generateDiagPredomMatr(int N);
double **generatePosDefMatr(int N);
#endif /* MatrOperations_h */
