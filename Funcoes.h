#ifndef __FUNCOES_H__
#define __FUNCOES_H__

#define MAX_ARQ 100

#include "sislin.h"


double normaEuclidiana(double* vetor, unsigned int tam);
int ehPar(unsigned int* k);
void pegaParametros(int argc, char *argv[],unsigned int *tam,unsigned int *kDiag,int *p,int *it,double *epsilon,char *arquivo);


void calcula_residuo (SistLinear_t *SL, double *x, double *r);
double calcula_alpha (double *r, double *d, SistLinear_t *SL);
void atualiza_r (double *r, double alpha, SistLinear_t *SL, double * d);
double calcula_beta (double *r, double *r_antigo, int tam);
int parada (double *x, double *x_antigo, double epsilon, int tam);
void gradiente_conjugado (SistLinear_t *SL, double *x, int p, int it, double epsilon, char *arquivo);


#endif // __FUNCOES_H__
