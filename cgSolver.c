#include <stdio.h>
#include <stdlib.h>
#include "Funcoes.h"
#include "sislin.h"


int main(int argc, char *argv[]){
    unsigned int tam;
    unsigned int kDiag;
    int p;
    int it;
    double epsilon;
    char arquivo[MAX_ARQ];
    SistLinear_t *SL;

    
    pegaParametros (argc,argv,&tam,&kDiag,&p,&it,&epsilon,arquivo);
    srand (20222);
    
    SL = lerSisLin (tam,kDiag);
    
	double *x = (double *) malloc (SL->n * (sizeof(double)));

    epsilon = 1.0e-10;
    // it = 5;


    // SL->A[0][0] = 3;
    // SL->A[0][1] = 2;
    // SL->A[1][0] = 2;
    // SL->A[1][1] = 6;

    // SL->b[0] = 2;
    // SL->b[1] = -8;



    // SL->A[0][0] = 4;
    // SL->A[0][1] = 2;
    // SL->A[0][2] = -4;
    
    // SL->A[1][0] = 2;
    // SL->A[1][1] = 10;
    // SL->A[1][2] = 4;

    // SL->A[2][0] = -4;
    // SL->A[2][1] = 4;
    // SL->A[2][2] = 9;

    // SL->b[0] = 6;
    // SL->b[1] = 12;
    // SL->b[2] = 5;



    // ./cgSolver -n 5 -k 3 -p 0 -i 10 -o teste
    // ./cgSolver -n 3 -k 3 -p 0 -i 10 -o teste
    
    prnSisLin (SL);

    gradiente_conjugado (SL, x, p, it, epsilon, arquivo);

    // prnVetor (x, SL->n);

    return 0;
}