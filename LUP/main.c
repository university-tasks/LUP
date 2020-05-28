//
//  main.c
//  LUP
//
//  Created by Артем on 14.04.2020.
//  Copyright © 2020 Артем. All rights reserved.
//
#define myExp 0.00000001
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>


#include "MatrOperations/MatrOperations.h"
#include "VectOperations/VectOperations.h"
#include "LUP/LUP.h"
#include "QR/QR.h"
#include "Jacobi_Seidel/Jacobi_Seidel.h"

//2 4 6
//2 0 2
//6 8 14

//2 0 4

int main() {
    
    int N, M, i, action = 0;
    double det = -100000000;
    int rank;
    int rankExt;
    printf("Size: ");
    scanf("%i",&N);


    M = N + 1;
    
    double *cond = (double*)malloc(sizeof(double) * 1);
    double *b = (double*)malloc(sizeof(double) * N);
    double *x = (double*)malloc(sizeof(double) * N);
    int    *p = (int*)malloc(sizeof(int) * N);
    for (i = 0; i <= N; i++) {
        p[i] = i;
    }

    
    double **A = (double **)malloc(N*sizeof(double *));
    double **IA = (double **)malloc(N*sizeof(double *));
    double **A1 = (double **)malloc(N*sizeof(double *));
    double **AExt = (double **)malloc(N*sizeof(double *));
    
    double **L = (double **)malloc(N*sizeof(double *));
    double **U = (double **)malloc(N*sizeof(double *));
    
    double **Q = (double **)malloc(N*sizeof(double *));
    double **R = (double **)malloc(N*sizeof(double *));
    
    for(int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N*sizeof(double));
        IA[i] = (double *)malloc(N*sizeof(double));
        A1[i] = (double *)malloc(N*sizeof(double));
        AExt[i] = (double *)malloc(M*sizeof(double));
        
        L[i] = (double *)malloc(N*sizeof(double));
        U[i] = (double *)malloc(N*sizeof(double));
        
        Q[i] = (double *)malloc(N*sizeof(double));
        R[i] = (double *)malloc(N*sizeof(double));
    }
  
    int q;
    printf("1)QR\n2)LUP\n3)Jacobi and Seidel iteration methods\n4)exit\n>>>>>>>>> ");
    scanf("%i", &action);
    switch (action) {
        case 1:
            
            generateMatr(A, N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    A1[i][j] = A[i][j];
                    AExt[i][j] = A[i][j];
                }
            }
            printf("Generated matrix:\n");
            showMatr(A, N, N);
            
            Q = GrammSchmidtProcess(Q, A, N);
            Q = cleanMatr(Q, N);
            R = QR(Q, R, A, N);
            R = cleanMatr(R, N);
            
            printf("Q:\n");
            showMatr(Q, N, N);
            
            printf("R:\n");
            showMatr(R, N, N);
            
            checkQR(A, Q, R, N);
            
            printf("1)Generate b\n2)enter\nn)exit\n>>>>>>>  ");
            scanf("%i",&q);
            
            switch (q) {
                case 1:
                    generateVec(b, N);
                    break;
                case 2:
                    printf("Enter coefficients of the vector: ");
                    for (i = 0; i < N; ++i) {
                        scanf("%lf", &b[i]);
                    }
                    break;
                default:
                    return 0;;
            }
            
            
            x = solveQR(Q, R, b, N);

            showDVect(x, N);
            return 0;
        case 2:
            
            generateMatr(A, N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    A1[i][j] = A[i][j];
                    AExt[i][j] = A[i][j];
                }
            }
            
            LUP(A1, p, N);
            checkLUP(A, p, L, U, N, A1);
            break;
        
        case 3:
             JacobiSeidel(N);
        default:
            return 0;
    }
    
    action = -1;
    q = -1;
    while (action!=9) {
        printf("COMMANDS\n");
        printf("1)Find determinant\n");
        printf("2)Solve SLAE\n");
        printf("3)Find inverse matrix\n");
        printf("4)Find condition number\n");
        printf("5)Find rank of a singular matrix and solve if it is joint\n");
        printf("6)Show LUP-decomposition\n");
        printf("7)Show L\n");
        printf("8)Show U\n");
        printf("n)Exit\n");
        printf(">>>> ");
        
        
        
        scanf("%i",&action);
        switch (action) {
            case 1:
                det = determinant(A1, p, N);
                printf("Determinant: %f%s",det,"\n");
                break;
            case 2:
                printf("1)Generate\n2)enter\nn)exit\n>>>>>>>  ");
                scanf("%i",&q);
                
                switch (q) {
                    case 1:
                        generateVec(b, N);
                        break;
                    case 2:
                        printf("Enter coefficients of the vector: ");
                        for (i = 0; i < N; ++i) {
                            scanf("%lf", &b[i]);
                        }
                        break;
                    default:
                        return 0;;
                }
                x = solve(L, U, b, p, N);
                showDVect(x, N);
                checkSolve(A, x, b, N);
                break;
            case 3:
                invert(L, U, p, N, IA);
                checkInvert(A, IA, N, p);
                break;
            case 4:
                condNum(A, IA, N, cond);
                break;
            case 5:
                
                generateSingularMatr(A, N);
                printf("Singular matrix:\n");
                showMatr(A, N, N);

                for(int i = 0; i < N; ++i) {
                  for(int j = 0; j < N; ++j) {
                    A1[i][j] = A[i][j];
                  }
                }
                
                rank = rankMatr(A, N, N);
                printf("Rank of a singular matrix: ");
                printf("%i%s", rank,"\n");
    
                printf("1)Generate vector\n2)enter\nn)exit\n>>>>>>>  ");
                scanf("%i",&q);
                
                switch (q) {
                    case 1:
                        generateVec(b, N);
                        break;
                    case 2:
                        printf("Enter coefficients of the vector: ");
                        for (i = 0; i < N; ++i) {
                            scanf("%lf", &b[i]);
                        }
                        break;
                    default:
                        return 0;;
                }
                
                for (i = 0; i < N; ++i) {
                    AExt[i][N] = b[i];
                }
                
                rankExt = rankMatr(AExt, N, M);
                printf("Rank of an extended singular matrix: ");
                printf("%i%s",rankExt,"\n");
                
                if (rankExt == rank) {

                    x = solve(L, U, b, p, N);
                    printf("Solution:\n");
                    showDVect(x, N);
                    checkSolve(A, x, b, N);
                }
                break;
            case 6:
                showMatr(A1, N, N);
                for (int i = 0; i < N; ++i) {
                    printf("%i ",p[i]);
                }
                printf("\n");
                break;
            case 7:
                showMatr(L, N, N);
                break;
            case 8:
                showMatr(U, N, N);
                break;
            case 9:
                QR(A, Q, R, N);
                break;
            default:
                return 0;;
        }
    }
    
    
    return 0;
}
//-----------------------------------------------------------

