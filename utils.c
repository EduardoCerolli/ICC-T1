#include <stdio.h>
#include "Funcoes.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>


/******
* Função para medir o Tempo
* Onde o tempo decorrido é medido pela diferença do "timestamp" medidos antes e depois da região de interesse.
**********/
double timestamp(void){
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
}


/***********************
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
 * k: numero de diagonais da matriz A
 * 
 * Antes da primeira chamada da função "generateRandomA()", e somente uma vez em todo código, você deve inicializar a sequência de números aleatóreos chamando a função:
 * srand(20222);
 ***********************/
double generateRandomA( unsigned int i, unsigned int j, unsigned int k ){
  double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

/***********************
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * k: numero de diagonais da matriz A
 ***********************/
double generateRandomB( unsigned int k ){
  double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k<<2) * (double)rand() * invRandMax;
}