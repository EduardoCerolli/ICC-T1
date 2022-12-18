#ifndef __UTILS_H__
#define __UTILS_H__

#define MAX_ARQ 100

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

// real_t: tipo usado para representar valores em ponto flutuante
typedef double real_t;


double timestamp(void);
double generateRandomA( unsigned int i, unsigned int j, unsigned int k );
double generateRandomB( unsigned int k );


#endif // __UTILS_H__
