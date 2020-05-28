//
//  QR.c
//  LUP
//
//  Created by Артем Белов on 19.05.2020.
//  Copyright © 2020 Артем. All rights reserved.
//
#define myExp 0.00000001
#include "QR.h"
#include "MatrOperations.h"
#include "VectOperations.h"

#include <stdlib.h>
#include <math.h>

double **GrammSchmidtProcess(double **Q, double **A, int N) {
  
  double **T = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        T[i] = (double *)malloc(N*sizeof(double));
    }

  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      T[i][j] = A[i][j];
    }
  }

  transpose(T, N);
  Q[0] = orthonormVect(T[0], N);

  for(int i = 1; i < N; ++i) {

    double *temp = (double*)malloc(sizeof(double) * N);
    for(int k = 0; k < N; ++k) {
       temp[i] = 0;
    }
    for(int j = 0; j < i; ++j) {
       temp = SumVect(temp, proj(Q[j], T[i], N), N);
    }
    Q[i] = DifVect(T[i], temp, N);
    Q[i] = orthonormVect(Q[i], N);
  }
  transpose(Q, N);

  return Q;
}

double **QR(double **Q, double **R, double **A, int N) {

    transpose(Q, N);

    R = multMatr(Q, A, N);
    transpose(Q, N);

    return R;
}

double **cleanMatr(double **Q, int N) {
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        if(isnan(Q[i][j]) != 0) {
          Q[i][j] = 0;
        }
      }
    }
    return Q;
}

void checkQR(double **A,double **Q, double **R, int N) {

    double **tempQ = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
      tempQ[i] = (double *)malloc(N*sizeof(double));
    }
    tempQ = multMatr(Q, R, N);

    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        if(fabs(tempQ[i][j] - A[i][j]) > myExp) {
          printf("ОШИБКА\n");
          return;
        }
      }
    }
    printf("ПРОВЕРКА ПРОЙДЕНА\n");
}
double * gauss(double **a, double *y, int N) {
  double max;
  int k, index;
  double *x = (double*)malloc(sizeof(double) * N);
  k = 0;
  while (k < N) {
    // Поиск строки с максимальным a[i][k]
    max = fabs(a[k][k]);
    index = k;
    for (int i = k + 1; i < N; ++i) {
      if (fabs(a[i][k]) > max) {
        max = fabs(a[i][k]);
        index = i;
      }
    }
    // Перестановка строк
    if (max < myExp) {
      // нет ненулевых диагональных элементов
      printf("Решение получить невозможно из-за нулевого столбца \n");
       printf("%i%s",index," матрицы A\n");
      return 0;
    }
    for (int j = 0; j < N; j++) {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;
    // Нормализация уравнений
    for (int i = k; i < N; i++) {
      double temp = a[i][k];
      if (fabs(temp) < myExp) continue; // для нулевого коэффициента пропустить
      for (int j = 0; j < N; j++)
        a[i][j] = a[i][j] / temp;
      y[i] = y[i] / temp;
      if (i == k)  continue; // уравнение не вычитать само из себя
      for (int j = 0; j < N; j++)
        a[i][j] = a[i][j] - a[k][j];
      y[i] = y[i] - y[k];
    }
    k++;
  }
  // обратная подстановка
  for (k = N - 1; k >= 0; k--) {
    x[k] = y[k];
    for (int i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
}

double *solveQR(double **Q, double **R, double *b, int N) {
  
  double *x = (double*)malloc(sizeof(double) * N);
  double *y = (double*)malloc(sizeof(double) * N);
  y = gauss(Q, b, N);
  showDVect(y, N);

for (int i = N-1; i > -1; --i) {
  double sum = 0;
        for (int j = i; j < N; ++j) {
            sum += R[i][j]*x[j];
        }
        x[i] = (y[i] - sum)/R[i][i];
        printf("%s%i%s%lf%s","x_",i," = ",x[i],"\n");
}

  return x;
}
