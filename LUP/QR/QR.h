//
//  QR.h
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//

#ifndef QR_h
#define QR_h

#include <stdio.h>

double **GrammSchmidtProcess(double **Q, double **A, int N);
double **QR(double **Q, double **R, double **A, int N);
double **cleanMatr(double **Q, int N);
void checkQR(double **A,double **Q, double **R, int N);
double * gauss(double **a, double *y, int N);
double *solveQR(double **Q, double **R, double *b, int N);

#endif /* QR_h */
