#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Funcoes.h"
#include "utils.h"



//Função p/ Calcular Norma
double normaEuclidiana(double* vetor, unsigned int tam){
    double soma,pot;
    soma=0.0;
    for (int i = 0; i < tam; i++){
        pot=pow(vetor[i],2);    
        soma+=pot;
    }
    return sqrt(soma);
}

//Verifica se o K passado é impar
int ehPar(unsigned int* k){
    if(*k%2==0){
        return 0;
    }
    return 1;
}


/*FUNÇÃO QUE SEPARA A ENTRADA NO ARGV E ATRIBUI OS VALORES AS VARIAVEIS*/
void pegaParametros(int argc, char *argv[],unsigned int *tam,unsigned int *kDiag,int *p,int *it,double *epsilon,char *arquivo){
    int verifica;
    *tam=atoi(argv[2]);
    //verificar se k é impar
    *kDiag=atoi(argv[4]);
    verifica=ehPar(kDiag);
    if(verifica==0){
        printf("K NÃO É IMPAR \n");
        exit(0);
    }

    *p=atoi(argv[6]);
    *it=atoi(argv[8]);
    if (argc==11){
        *epsilon = -10;
        strcpy(arquivo,argv[10]);
    }else{
        *epsilon=atof(argv[10]);
        strcpy(arquivo,argv[12]);
    }
}



void calcula_residuo (SistLinear_t *SL, double *x, double *r) {

    for (int i = 0; i < SL->n; i++) {
        r[i] = 0;
        for (int j = 0; j < SL->n; j++) {
            r[i] += SL->A[i][j] * x[j];
        }
        r[i] = SL->b[i] - r[i];
    }

    return;
}

double calcula_alpha (double *r, double *d, SistLinear_t *SL) {
    double alpha = 0;
    double aux_2 = 0;
    double aux[SL->n];

    // r.i * r.i
    for (int i = 0; i < SL->n; i++) {
        alpha += r[i] * r[i];
    }

    // (d.i)^T * A
    for (int i = 0; i < SL->n; i++) {
        aux[i] = 0;
        for (int j = 0; j < SL->n; j++) {
            aux[i] += SL->A[i][j] * d[j];
        }
    }

    // (aux)^T * d.i
    for (int i = 0; i < SL->n; i++) {
        aux_2 += aux[i] * d[i];
    }

    return alpha / aux_2;
}

void atualiza_r (double *r, double alpha, SistLinear_t *SL, double *d){

    double aux[SL->n];

    for (int i = 0; i < SL->n; i++) {
        aux[i] = 0;
        for (int j = 0; j < SL->n; j++) {
            aux[i] += (alpha * SL->A[j][i]) * d[j];
        }
    }

    for (int i = 0; i < SL->n; i++) {
        r[i] = r[i] - aux[i];
    }

    return;
}

double calcula_beta (double *r, double *r_antigo, int tam) {
    double beta = 0;
    double aux = 0;

    for (int i = 0; i < tam; i++) {
        beta += r[i] * r[i];
        aux += r_antigo[i] * r_antigo[i];
    }

    return beta / aux;
}

// retorna 0 se tiver que parar
int parada (double *x, double *x_antigo, double epsilon, int tam) {
    double maior = 0;
    double aux;

    if (epsilon == -10) {
        return 1;
    }

    for (int i = 0; i < tam; i++) {
        aux = x[i] - x_antigo[i];
        aux = ABS (aux);
        aux = aux / (ABS (x[i]));

        if (aux > maior) {
            maior = aux;
        }
    }

    if (maior < epsilon) {
        return 0;
    }

    return 1;
}

void gradiente_conjugado (SistLinear_t *SL, double *x, int p, int it, double epsilon, char *arquivo) {
    int cont = 0;
	double *d = (double *) malloc (SL->n * (sizeof(double)));
	double *r = (double *) malloc (SL->n * (sizeof(double)));
	double *r_antigo = (double *) malloc (SL->n * (sizeof(double)));
	double *x_antigo = (double *) malloc (SL->n * (sizeof(double)));
    double alpha;
	double beta;
    

    // 1 x.0
    for (int i = 0; i < SL->n; i++) {
        x[i] = 0;
    }

    // 2 d.0 e r.0
    calcula_residuo (SL, x, r);
    for (int i = 0; i < SL->n; i++) {
        d[i] = r[i];
    }

    do {
        // 4 alpha.i
        alpha = calcula_alpha (r, d, SL);

        // 5 x.i+1
        for (int i = 0; i < SL->n; i++) {
            x_antigo[i] = x[i];
            x[i] = x[i] + (alpha * d[i]);
        }

        // 6 r.i+1
        for (int i = 0; i < SL->n; i++) {
            r_antigo[i] += r[i];
        }
        atualiza_r (r, alpha, SL, d);

        // 7 beta.i+1
        beta = calcula_beta (r, r_antigo, SL->n);

        // 8 d.i+1
        for (int i = 0; i < SL->n; i++) {
            d[i] = r[i] + (beta * d[i]);
        }

        cont++;

    } while ((parada (x, x_antigo, epsilon, SL->n)) && cont < it);
    
    printf ("%d\n", cont);

    return;
}
