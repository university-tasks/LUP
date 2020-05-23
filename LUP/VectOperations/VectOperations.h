//
//  VectOperations.h
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#ifndef VectOperations_h
#define VectOperations_h

#include <stdio.h>

void showDVect(double *v, int N);
void showIVect(int *v, int N);
void generateVec(double *b, int N);
double scalMult(double *a, double *b, int N);
double *NumMultVect(double a, double *V, int N);
double *DifVect(double *a, double *b, int N);
double *SumVect(double *a, double *b, int N);
double *orthonormVect(double *v, int N);
double *proj(double *u, double *v, int N);

#endif /* VectOperations_h */
