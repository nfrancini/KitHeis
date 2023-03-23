#include "../head/head_&_structures.h"

// FUNZIONE PER CALCOLARE L'ENERGIA DEL SISTEMA
double ene(SystemParam_t *Par, Field_t *Fields){
  double kit_term, heis_term, ene;
  int i;

  kit_term = 0.0;
  heis_term = 0.0;
  // SOMMA SU TUTTI I SITI DEL SOTTORETICOLO NERO
  // CHE EQUIVALE A SOMMA SU TUTTE LE COPPIE
  for(i=0; i<(int)pow(Par->L, 2); i++){
    kit_term = kit_term + (Fields->spin[i][2])*(Fields->spin[nn[i][2]][2]) + (Fields->spin[i][0])*(Fields->spin[nn[i][0]][0]) + (Fields->spin[i][1])*(Fields->spin[nn[i][1]][1]);
    heis_term = heis_term + product(Par, Fields->spin[i], Fields->spin[nn[i][2]]) + product(Par, Fields->spin[i], Fields->spin[nn[i][0]]) + product(Par, Fields->spin[i], Fields->spin[nn[i][1]]);
  }
  ene = 2*sin(Par->phi)*kit_term + cos(Par->phi)*heis_term;
  return ene;
}

// FUNZIONE PER CALCOLARE LA DENSITÀ DI ENERGIA
double ene_dens(SystemParam_t *Par, Field_t *Fields){
  double ene_density;

  #ifdef DEBUG
  resetErr();
  #endif

  // DAL MOMENTO CHE LA SOMMA È SVOLTA SULLE COPPIE DEVO
  // DIVIDERE PER IL NUMERO DI CELLE NON PER IL NUMERO DI SITI
  ene_density = ene(Par, Fields) / ((int)pow(Par->L,2));

  #ifdef DEBUG
    if(errno == EDOM){
      printf("Funzione H_g, pow\n");
      perror("    errno == EDOM");
    }
    if(errno == ERANGE){
      printf("Funzione H_g, pow\n");
      perror("    errno == ERANGE");
      if(fetestexcept(FE_OVERFLOW)){
        puts("    FE_OVERFLOW was raised");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        puts("    FE_UNDERFLOW was raised");
        exit(EXIT_FAILURE);
      }
    }
    if(fetestexcept(FE_INVALID)){
      puts("    FE_INVALID was raised");
      exit(EXIT_FAILURE);
    }
    if(fetestexcept(FE_DIVBYZERO)){
      puts("    FE_DIVBYZERO was raised");
      exit(EXIT_FAILURE);
    }
  #endif

  return ene_density;
}

// PROCEDURA PER CALCOLARE LA MAGNETIZZAZIONE FERROMAGNETICA
void FM_magn(SystemParam_t *Par, Field_t *Fields, double *m){
  int i, j;
  for(j=0;j<K;j++){
    for(i=0;i<(Par->N);i++){
      m[j]= m[j] + (Fields->spin[i][j] / (Par->N));
    }
  }
}

// PROCEDURA PER CALCOLARE LA MAGNETIZZAZIONE DI NEEL
void neel_magn(SystemParam_t *Par, Field_t *Fields, double *n){
  int i, j;
  for(j=0;j<K;j++){
    for(i=0;i<(int)pow(Par->L,2);i++){
      n[j] = n[j] + ((Fields->spin[i][j] - Fields->spin[nn[i][2]][j]) / (Par->N));
    }
  }
}

// PROCEDURA PER CALCOLARE IL PARAMETRO D'ORDINE ZIGZAG
void zigzag_and_stripy_magn(SystemParam_t *Par, Field_t *Fields, double *z, double *s){
  int i, j;
  bool_t flag;
  for(j=0;j<K;j++){
    i=0;                    // PARTO DAL SITO ZERO DEL RETICOLO
    flag=FALSE;
    while(flag == FALSE){   // PER OGNI COMPONENTE USO LA SUPERCELLA FORMATA DA 4 CELLE UNITARIE COSTRUITA RADDOPPIANDO I VETTORI PRIMITIVI
      z[j] = z[j] + ((-Fields->spin[i][j]) + (-Fields->spin[i+Par->L+1][j]) + (Fields->spin[i+(int)pow(Par->L,2)][j]) + (Fields->spin[i+(int)pow(Par->L,2) + Par->L + 1][j]) + (Fields->spin[i+1][j]) + (Fields->spin[i+Par->L][j]) + (-Fields->spin[i+(int)pow(Par->L,2)+1][j]) + (-Fields->spin[i+(int)pow(Par->L,2)+Par->L][j]))/(Par->N);
      s[j] = s[j] + ((-Fields->spin[i][j]) + (-Fields->spin[i+Par->L+1][j]) + (-Fields->spin[i+(int)pow(Par->L,2)][j]) + (-Fields->spin[i+(int)pow(Par->L,2) + Par->L + 1][j]) + (Fields->spin[i+1][j]) + (Fields->spin[i+Par->L][j]) + (Fields->spin[i+(int)pow(Par->L,2)+1][j]) + (Fields->spin[i+(int)pow(Par->L,2)+Par->L][j]))/(Par->N);
      if((i%Par->L)+2<Par->L){
        i=i+2;
      }
      else if(((i%Par->L)+2>=Par->L) && (((i/Par->L)%Par->L)+2<Par->L)){
        i=i-(i%Par->L)+2*Par->L;
      }
      else if(((i%Par->L)+2>=Par->L) && (((i/Par->L)%Par->L)+2>=Par->L)){
        flag=TRUE;
      }
      else{
        printf("Errore in zigzag_and_stripy_magn\n");
        exit(EXIT_FAILURE);
      }
    }
  }
}

// PROCEDURA PER CALCOLARE IL PARAMETRO D'ORDINE DI NEEL
void measures(SystemParam_t *Par, Field_t *Fields, Obs_t *Obs){
  // OGNI VOLTA CHE MISURO INIZIALIZZO TUTTO A ZERO
  initializeObs(Obs);
  Obs->ene_density = ene_dens(Par, Fields);
  FM_magn(Par, Fields, Obs->M);
  neel_magn(Par, Fields, Obs->N);
  zigzag_and_stripy_magn(Par, Fields, Obs->Z, Obs->S);
}
