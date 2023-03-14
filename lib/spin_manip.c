#include "../head/head_&_structures.h"

// FUNZIONE CHE FA IL PRODOTTO S1 * S2 = SUM_i S1[i]*S2[i]
double product(SystemParam_t *Par, double *vec1, double *vec2){
  int i;
  double sum;

  sum = 0;
  for(i=0; i<K; i++){
    sum = sum + vec1[i]*vec2[i];
  }
  return sum;
}

// PROCEDURA CHE FA LA DIFFERENZA TRA DUE SPIN S1-S2
void diff(SystemParam_t *Par, double *vec1, double *vec2, double *diff){
  int i;

  for(i=0;i<K;i++){
    diff[i] = vec1[i] - vec2[i];
  }
}

// PROCEDURA CHE FA LA SOMMA TRA DUE SPIN S1+S2
void two_spin_sum(SystemParam_t *Par, double *vec1, double *vec2, double *somma){
  int i;

  for(i=0;i<K;i++){
    somma[i] = vec1[i] + vec2[i];
  }
}

// PROCEDURA CHE FA LA SOMMA TRA TRE SPIN S1+S2+S3
void three_spin_sum(SystemParam_t *Par, double *vec1, double *vec2, double *vec3, double *somma){
  int i;

  for(i=0;i<K;i++){
    somma[i] = vec1[i] + vec2[i] + vec3[i];
  }
}

// FUNZIONE CHE RINORMALIZZA IL MODULO DEGLI SPIN
// VA RICHIAMATA OGNI TANTO DURANTE GLI UPDATE
// PER RIAGGIUSTARE IL MODULO DEGLI SPIN CHE POTREBBERO
// ASSUMERE VALORI DIVERSI DA UNO PER ERRORI NUMERICI
void renormalize(SystemParam_t *Par, Field_t *Fields){
  int iSite, j;
  double norm;

  #ifdef DEBUG
  resetErr();
  #endif

  for(iSite=0;iSite<(Par->N);iSite++){
    norm = sqrt(product(Par, Fields->spin[iSite], Fields->spin[iSite]));   // CALCOLO LA NORMA
    #ifdef DEBUG
    if(errno == EDOM){
      printf("Funzione renormalize, radice quadrata\n");
      perror("    errno == EDOM");
    }
    if(errno == ERANGE){
      printf("Funzione renormalize, radice quadrata\n");
      perror("    errno == ERANGE");
    }
    if(fetestexcept(FE_INVALID)){
      puts("    FE_INVALID was raised");
      exit(EXIT_FAILURE);
    }
    #endif
    for(j=0; j<K;j++){
      Fields->spin[iSite][j] = Fields->spin[iSite][j] / norm;              // DIVIDO PER LA NORMA, ALMENO TORNA 1
    }
  }
}

// PROCEDURA PER COPIARE UN VETTORE DI SPIN IN UN ALTRO
void spin_copy(SystemParam_t *Par, double *vec_in, double *vec_out){
  int i;
  for(i=0;i<K;i++){
    vec_out[i]=vec_in[i];
  }
}
