#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"

// Alocaçao de matriz em memória. 
SistLinear_t* alocaSisLin (unsigned int n)
{

  
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if ( SL ) {
    
    SL->n = n;
    SL->A = (double **) malloc(n * sizeof(double *));
    SL->b = (double *) malloc(n * sizeof(double));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SL->A[0] = (double *) malloc(n * n * sizeof(double));
    if (!(SL->A[0])) {
      liberaSisLin(SL);
      return NULL;
    }
    
    for (int i=1; i < n; ++i) {
      SL->A[i] = SL->A[i-1]+n;
    }
  }
  
  return (SL);
}

// Liberacao de memória
void liberaSisLin (SistLinear_t *SL)
{
  if (SL) {
    if (SL->A) {
      if (SL->A[0]) free (SL->A[0]);
    free(SL->A);
    }
    
    if (SL->b) free(SL->b);

    free(SL);
  }
}



SistLinear_t *lerSisLin (unsigned int n, unsigned int k)
{

  
  SistLinear_t *SL;
  
  

  SL = alocaSisLin (n);
  
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      SL->A[i][j]=0.0;

  for(int i=0; i < n; ++i)
    SL->b[i]=0.0; 
  
  //gera a matriz com as K-DIAGONAIS
  for ( int l = 0; l < k-1; l++){
    for(int i=0; i < n; ++i){
      for(int j=0; j < n; ++j){      
          if(i==j){
            SL->A[i][j]=generateRandomA(i,j,k);
            if(i+k<=n+1)
              SL->A[i][j+l]=generateRandomA(i,j,k);
            if(i-k>=0) 
              SL->A[i][j-l]=generateRandomA(i,j,k);

          }
      }
    }
  }
  
  for(int i=0; i < n; ++i)
    SL->b[i]=generateRandomB(i+1);
  return SL;
}


void prnSisLin (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n  ");
    for(int j=0; j < n; ++j)
      printf ("%10g", SL->A[i][j]);
    printf ("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (double *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10g ", v[i]);
  printf ("\n\n");

}

