#ifndef __SISLIN_H__
#define __SISLIN_H__

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  double **A; // coeficientes
  double *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;




// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);

// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLin (unsigned int n, unsigned int k);
void prnSisLin (SistLinear_t *SL);
void prnVetor (double *vet, unsigned int n);

#endif // __SISLIN_H__

