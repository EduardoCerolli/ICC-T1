#ifndef __FUNCOES_H__
#define __FUNCOES_H__

#define MAX_ARQ 100

#include "sislin.h"


double normaEuclidiana(double* vetor, unsigned int tam);
int ehPar(unsigned int* k);
void pegaParametros(int argc, char *argv[],unsigned int *tam,unsigned int *kDiag,int *p,int *it,double *epsilon,char *arquivo);


void calcula_residuo (SistLinear_t *SL, double *x, double *r);
double calcula_alpha (double *r, double *d, double *z, SistLinear_t *SL);
void atualiza_r (double *r, double alpha, SistLinear_t *SL, double * d);
double calcula_beta (double *r, double *r_antigo, double *z, double *z_antigo, int tam);
int parada (double *x, double *x_antigo, double epsilon, unsigned int tam);
void matriz_identidade (double **M, unsigned int tam);
void matriz_inversa (double **M, SistLinear_t *SL);
void calcula_z (double *z, double **M, double *r, unsigned int tam);
void gradiente_conjugado (SistLinear_t *SL, double *x, int p, int it, double epsilon, char *arquivo);

/*teste*/
void cofactor(SistLinear_t *SL,double **M);
#endif // __FUNCOES_H__
