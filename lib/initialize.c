#include "../head/head_&_structures.h"

// PROCEDURA GEOMETRY PER MEMORIZZARE I PRIMI VICINI;
void geometry(SystemParam_t *Par){
  int i;

  // SFRUTTO L'ORDINAMENTO LESSICOGRAFICO i_x = x + y*L + d*L^2 con x,y in [0,L-1] e d={0,1}
  // IN QUESTO MODO:
  // x = i_x % L
  // y = (i_x / L) % L
  // d = (i_x / L^2) % L

  // IMPLEMENTAZIONE
  // IL SECONDO INDICE INDICA LA DIREZIONE DEL LEGAME O <-> X, 1 <-> Y, 2 <-> Z
  for(i=0;i<(Par->N);i++){
    // CONDIZIONE SOTTORETICOLO NERO
    if(((i / (int)pow(Par->L,2)) % (Par->L)) == 0){
      // CONDIZIONE x=0,y=/=0, DEVO IMPORRE PBC SU X-BOND
      if(((i % (Par->L)) == 0) && ((i / (Par->L)) % (Par->L) != 0)){
        nn[i][0] = i + (int)pow(Par->L,2) + (Par->L) - 1;
        nn[i][1] = i + (Par->L)*(Par->L - 1);
        nn[i][2] = i + (int)pow(Par->L,2);
      }
      // CONDIZIONE y=0,x=/=0 DEVO IMPORRE PBC SU Y-BOND
      else if(((i / (Par->L)) % (Par->L) == 0) &&  ((i % (Par->L)) != 0)){
        nn[i][1] = i + (Par->L)*(2*(Par->L) - 1);
        nn[i][0] = i + (int)pow(Par->L,2) - 1;
        nn[i][2] = i + (int)pow(Par->L,2);
      }
      // CONDIZIONE x=, y=0 DEVO IMPORRE PBC SU X-BOND E Y-BOND
      else if(((i % (Par->L)) == 0) && ((i / (Par->L)) % (Par->L) == 0)){
        nn[i][0] = i + (int)pow(Par->L,2) + (Par->L) - 1;
        nn[i][1] = i + (Par->L)*(2*(Par->L) - 1);
        nn[i][2] = i + (int)pow(Par->L,2);
      }
      // NON CI SONO DA IMPORRE CONDIZIONI DI BORDO
      else {
        nn[i][0] = i + (int)pow(Par->L,2) - 1;
        nn[i][1] = i + (Par->L)*(Par->L - 1);
        nn[i][2] = i + (int)pow(Par->L,2);
      }
    }
    // CONDIZIONE SOTTORETICOLO BIANCO
    else if(((i / (int)pow(Par->L,2)) % (Par->L)) == 1){
      // CONDIZIONE x=L-1, y=/= L-1 DEVO IMPORRE PBC SU X-BOND
      if(((i % (Par->L)) == (Par->L - 1)) && ((i / (Par->L)) % (Par->L) != (Par->L - 1))){
        nn[i][0] = i - ((int)pow(Par->L,2) + (Par->L) - 1);
        nn[i][1] = i - (Par->L)*(Par->L - 1);
        nn[i][2] = i - (int)pow(Par->L,2);
      }
      // CONDIZIONE y=L-1,x=/=L-1 DEVO IMPORRE PBC SU Y-BOND
      else if(((i / (Par->L)) % (Par->L) == (Par->L - 1)) && ((i % (Par->L)) != (Par->L - 1)) ){
        nn[i][1] = i - (Par->L)*(2*(Par->L) - 1);
        nn[i][0] = i - (int)pow(Par->L,2) + 1;
        nn[i][2] = i - (int)pow(Par->L,2);
      }
      // CONDIZIONE x=L-1,y=L-1 DEVO IMPORRE PBC SU X-BOND E Y-BOND
      else if(((i % (Par->L)) == (Par->L - 1)) && ((i / (Par->L)) % (Par->L) == (Par->L - 1))){
        nn[i][0] = i - ((int)pow(Par->L,2) + (Par->L) - 1);
        nn[i][1] = i - (Par->L)*(2*(Par->L) - 1);
        nn[i][2] = i - (int)pow(Par->L,2);
      }
      // NON CI SONO DA IMPORRE CONDIZIONI DI BORDO
      else {
        nn[i][0] = i - (int)pow(Par->L,2) + 1;
        nn[i][1] = i - (Par->L)*(Par->L - 1);
        nn[i][2] = i - (int)pow(Par->L,2);
      }
    }
  }
}

// PROCEDURA DI ALLOCAZIONE DINAMICA DELLA MEMORIA
void allocation(SystemParam_t *Par , Field_t *Fields, Obs_t *Obs){
  int i;

  nn = (int **)malloc((Par->N)*sizeof(int *));                    // VETTORE PER PRIMI VICINI
  for(i=0;i<(Par->N);i++){
    nn[i] = (int *)malloc((D)*sizeof(int));
  }

  Fields->spin = (double **)malloc((Par->N)*sizeof(double *));    // CONFIGURAZIONE DI SPIN
  for(i=0;i<(Par->N);i++){
    Fields->spin[i] = (double *)malloc((K)*sizeof(double));
  }

  // OSSERVABILI VETTORIALI PER LE DIVERSE FASI
  Obs->N = (double *)malloc((K)*sizeof(double));
  Obs->S = (double *)malloc((K)*sizeof(double));
  Obs->M = (double *)malloc((K)*sizeof(double));
  Obs->Z = (double *)malloc((K)*sizeof(double));
}

// PROCEDURA PER INIZIALIZZARE LA CONFIGURAZIONE DEI CAMPI
// PER IL MOMENTO VENGONO INIZIALIZZATI IN UNA CONFIGURAZIONE ORDINATA
// CON IL CAMPO DI SPIN A COMPONENTI IDENTICHE
void initializeFields(SystemParam_t *Par, Field_t *Fields){
  int i, j, mu;
  double a, b, c;

  if((Par->iStart) == 1){             // COPIO L'ULTIMA CONFIGURAZIONE
    readFields(Par, Fields);
  }
  else if((Par->iStart) == 0){        // INIZIALIZZO IN MODO ORDINATO
    for(i=0;i<(Par->N);i++){
      for(j=0;j<K;j++){
        Fields->spin[i][j] = 1/sqrt(K);
      }
    }
  }
  else if((Par->iStart) == 2){        // INIZIALIZZO IN MODO RANDOM
    srand(time(NULL));
    for(i=0;i<(Par->N);i++){
      for(j=0;j<K;j++){
        Fields->spin[i][j]= (double)rand() / RAND_MAX;   // OGNI COMPONENTE NUMERO RANDOM TRA 0 ED 1
      }
    }
    renormalize(Par, Fields);         // RINORMALIZZO TUTTI GLI SPIN AD 1
  }
  else{
    printf("Errore in iStart\n");
    exit(EXIT_FAILURE);
  }
}

// PROCEDURA DI SUPPORTO PER METTERE A ZERO TUTTE LE OSSERVABILI INIZIALMENTE
// NON SO SE È UTILE MA PER SICUREZZA LO FACCIO
void initializeObs(Obs_t *Obs){
  int i;

  Obs->ene_density = 0.0;
  for(i=0;i<K;i++){
    Obs->N[i]=0.0;
    Obs->S[i]=0.0;
    Obs->M[i]=0.0;
    Obs->Z[i]=0.0;
  }
}

// PROCEDURA DI DEALLOCAZIONE DELLA MEMORIA DINAMICA PRECEDENTEMENTE ALLOCATA
void deallocation(Field_t *Fields, Obs_t *Obs){
  free(nn);
  free(Fields->spin);
  free(Obs->M);
  free(Obs->N);
  free(Obs->Z);
  free(Obs->S);
}

// PROCEDURA PER INIZIALIZZARE I PARAMETRI DI SISTEMA E LE CONFIGURAZIONI
void initializeSystem(SystemParam_t *Par, Field_t *Fields, Obs_t *Obs, char const *finput){

  read_from_input_Param(Par, finput);     // LEGGO DA input.txt I PARAMETRI DEL SISTEMA

  if((Par->iStart) == 1){                 // LEGGO I PARAMETRI DI ACCETTANZA SALVATI
    readEps(Par);
  }

  allocation(Par, Fields, Obs);           // ALLOCO DINAMICAMENTE LA MEMORIA

  geometry(Par);                          // MEMORIZZO I PRIMI VICINI

  initializeFields(Par, Fields);          // INIZIALIZZO I CAMPI

  initializeObs(Obs);                     // INIZIALIZZO LE OSSERVABILI

  #ifdef DEBUG
  resetErr();   // RESETTO LE VARIABILI DI ERROR HANDLING
  #endif
}
