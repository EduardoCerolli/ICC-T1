// Eduardo Henrique dos Santos Cerolli - GRR20190397
// Thauan de Souza Tavares da Silva - GRR20171591

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

    gradiente_conjugado (SL, x, p, it, epsilon, arquivo);

    return 0;
}