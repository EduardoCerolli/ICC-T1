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

double calcula_alpha (double *r, double *d, double *z, SistLinear_t *SL) {
    double alpha = 0;
    double aux_2 = 0;
    double aux[SL->n];

    // r.i * r.i
    for (int i = 0; i < SL->n; i++) {
        alpha += r[i] * z[i];
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

double calcula_beta (double *r, double *r_antigo, double *z, double *z_antigo, int tam) {
    double beta = 0;
    double aux = 0;

    for (int i = 0; i < tam; i++) {
        beta += r[i] * z[i];
        aux += r_antigo[i] * z_antigo[i];
    }

    return beta / aux;
}

// retorna 0 se tiver que parar
int parada (double *x, double *x_antigo, double epsilon, unsigned int tam) {
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

void matriz_identidade (double **M, unsigned int tam) {

    for (int i = 0; i < tam; i++){
        for(int j = 0; j < tam; j++){
            if (i == j) {
                M[i][j] = 1.0;
            }
            else {
                M[i][j] = 0.0;
            }
        }
    }
    
    return;
}

void matriz_inversa(double **M, SistLinear_t *SL) {

    double pivo, aux;
    double matriz_aux[SL->n][SL->n];

    for (int i = 0; i < SL->n; i++) {
        for (int j = 0; j < SL->n; j++) {
            matriz_aux[i][j] = SL->A[i][j];
        }
    }

    for (int j = 0; j < SL->n; j++){
        pivo = matriz_aux[j][j];
        for(int k = j; k < SL->n; k++){
              matriz_aux[j][k] = (matriz_aux[j][k])/(pivo); 
              M[j][k] = (M[j][k])/(pivo); 
        }
        
        for(int i = 0; i < SL->n; i++){
          if(i != j){
               aux = matriz_aux[i][j];
               
               for(int k = 0; k < SL->n; k++){
                  matriz_aux[i][k] = (matriz_aux[i][k]) - (aux * matriz_aux[j][k]); 
                  M[i][k] = (M[i][k]) - (aux * M[j][k]);  
               }
          }
        }  
    }

    for(int linha = 0; linha <  SL->n; linha++){
          for(int coluna = 0; coluna <  SL->n; coluna++){
              printf("%g \t", M[linha][coluna]);
                 
          } 
          
          printf("\n"); 
    }

    return;
}

void calcula_z (double *z, double **M, double *r, unsigned int tam) {

    for (int i = 0; i < tam; i++){
        z[i] = 0.0;
        for(int j=0; j < tam; j++){
            z[i] += M[i][j] * r[j];
        }
    }

    return;
}

void gradiente_conjugado (SistLinear_t *SL, double *x, int p, int it, double epsilon, char *arquivo) {
    int cont = 0;
	double d[SL->n];
	double r[SL->n];
	double z[SL->n];
	double r_antigo[SL->n];
	double x_antigo[SL->n];
	double z_antigo[SL->n];
    double alpha;
	double beta;
    double **M = (double **) malloc (SL->n * (sizeof(double*)));
    for (int i = 0; i < SL->n; i++) {
        M[i] = (double *) malloc (SL->n * (sizeof(double)));
    }

    matriz_identidade (M, SL->n);
    matriz_inversa (M, SL);

    // // 1 x.0
    // for (int i = 0; i < SL->n; i++) {
    //     x[i] = 0;
    // }

    // // 2 d.0 e r.0
    // calcula_residuo (SL, x, r);

    // // z.0
    // calcula_z (z, M, r, SL->n);
    // for (int i = 0; i < SL->n; i++) {
    //     d[i] = z[i];
    // }

    // do {
    //     // 4 alpha.i
    //     alpha = calcula_alpha (r, d, z, SL);

    //     // 5 x.i+1
    //     for (int i = 0; i < SL->n; i++) {
    //         x_antigo[i] = x[i];
    //         x[i] = x[i] + (alpha * d[i]);
    //     }

    //     // 6 r.i+1
    //     for (int i = 0; i < SL->n; i++) {
    //         r_antigo[i] = r[i];
    //     }
    //     atualiza_r (r, alpha, SL, d);

    //     // z.i+1
    //     for (int i = 0; i < SL->n; i++) {
    //         z_antigo[i] = z[i];
    //     }
    //     calcula_z (z, M, r, SL->n);

    //     // 7 beta.i+1
    //     beta = calcula_beta (r, r_antigo, z, z_antigo, SL->n);

    //     // 8 d.i+1
    //     for (int i = 0; i < SL->n; i++) {
    //         d[i] = z[i] + (beta * d[i]);
    //     }

    //     cont++;

    // } while ((parada (x, x_antigo, epsilon, SL->n)) && cont < it);
    
    // printf ("%d\n", cont);

    return;
}
