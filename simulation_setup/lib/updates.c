#include "../head/head_&_structures.h"

// PROCEDURA PER IL CAMPO MEDIO INTORNO ALLO SPIN DEL SITO iSite
// È FONDAMENTALE PER LO STEP METROPOLIS E OVERRELAXATION
void mean_field(SystemParam_t *Par, Field_t *Fields, int site, double *mean_field){
  int i;
  double three_sum[K];

  #ifdef DEBUG
  resetErr();
  #endif

  // SOMMA VETTORIALE DEI TRE SPIN PRIMI VICINI AD site
  three_spin_sum(Par, Fields->spin[nn[site][0]], Fields->spin[nn[site][1]], Fields->spin[nn[site][2]], three_sum);

  // IL CAMPO MEDIO È UN VETTORE CON COMPONENTI f_i = -2*alpha*S{nn(site, i)}_i + (1-alpha)*three_spin_sum_i
  for(i=0;i<K;i++){
    mean_field[i] = 2*(sin(Par->phi))*Fields->spin[nn[site][i]][i] + (cos(Par->phi))*three_sum[i];
  }

  #ifdef DEBUG
    if(fetestexcept(FE_INVALID)){
      puts("    FE_INVALID was raised");
      exit(EXIT_FAILURE);
    }
  #endif
}

// PROCEDURA PER CREARE LO SPIN DI PROVA trial NEL SITO site
void spin_trial(SystemParam_t *Par, Field_t *Fields, int site, double *trial){
  int i, j;
  double a, b;

  #ifdef DEBUG
  resetErr();
  #endif

  spin_copy(Par, Fields->spin[site], trial);    // COPIO IL VETTORE DI SPIN IN TRIAL

  i = floor((K)*rndm());                        // SCELGO IN MODO RANDOM DUE COMPONENTI
  j = floor((K)*rndm());
  while(j==i){                                  // PER EVITARE DI AVERE DUE COMPONENTI UGUALI
    j = floor((K)*rndm());
  }

  a = (2*rndm()-1)*Par->eps1;                   // NUMERO RANDOM TRA [-eps1;eps1]

  b = sqrt(1-pow(a,2));

  #ifdef DEBUG
  if(errno == EDOM){
    printf("Funzione scalar_trial, radice quadrata\n");
    perror("    errno == EDOM");
  }
  if(errno == ERANGE){
    printf("Funzione scalar_trial o pow, radice quadrata\n");
    perror("    errno == ERANGE");
  }
  if(fetestexcept(FE_INVALID)){
    puts("    FE_INVALID was raised");
    exit(EXIT_FAILURE);
  }
  if(fetestexcept(FE_OVERFLOW)){
    puts("    FE_OVERFLOW was raised");
    exit(EXIT_FAILURE);
  }
  #endif

  // RUOTO LE DUE COMPONENTI SCELTE E LE SALVO IN trial
  trial[i] = b*Fields->spin[site][i] - a*Fields->spin[site][j];
  trial[j] = a*Fields->spin[site][i] + b*Fields->spin[site][j];

  // ESTRAGGO DUE NUMERO TRA [-eps1;eps1]
  // SE IL SEGNO È NEGATIVO RIBALTO LE COMPONENTI
  // CON QUESTA PRESCRIZIONE, FISSANDO eps1=0, HO ACCETTANZA 1 (IDENTITÀ)
  if((2*rndm()-1)*Par->eps1<0){
    trial[i]=-trial[i];
  }
  if((2*rndm()-1)*Par->eps1<0){
    trial[j]=-trial[j];
  }
}

// PROCEDURA PER L'UPDATE METROPOLIS
void metro_update(SystemParam_t *Par, Field_t *Fields){
  int iSite, j;
  double r, a, delta, c;
  double prod;
  double trial[K], d[K], mean[K];

  #ifdef DEBUG
  resetErr();
  #endif

  r=0;
  for(iSite=0;iSite<Par->N;iSite++){                      // SPAZZO SU TUTTI I SITI DEL RETICOLO
    spin_trial(Par, Fields, iSite, trial);                // CREO L'ARRAY DI PROVA
    diff(Par, trial, Fields->spin[iSite], d);             // DIFFERENZA TRA SPIN DI PROVA E SPIN ORIGINALE
    mean_field(Par, Fields, iSite, mean);                 // CAMPO MEDIO ATTORNO ALLO SPIN IN iSite
    delta = product(Par, d, mean);                        // DIFFERENZA DI ENERGIA TRA LE DUE CONFIGURAZIONI
    c=-delta/Par->T;                                      // ESPONENTE PER IL CALCOLO DI r

    #ifdef DEBUG
    double ene1, delta_loc;
    double ene2, delta_glob;
    delta_loc = delta;
    ene1 = ene(Par, Fields);
    #endif

    if(c>0){                                              // ACCETTO SICURAMENTE IL PASSO
      spin_copy(Par, trial, Fields->spin[iSite]);
      acc=acc+1.0;
    }
    else{
      a=rndm();
      r=exp(c);

      #ifdef DEBUG
      if(errno == ERANGE){
        printf("Funzione update_metro\n");
        perror("    errno == ERANGE");
      }
      if(fetestexcept(FE_OVERFLOW)){
        printf("    FE_OVERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        printf("    FE_UNDERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_DIVBYZERO)){
        printf("    FE_DIVBYZERO was raised\n");
        exit(EXIT_FAILURE);
      }
      #endif

      if(a<r){                                          // ACCETTO IL PASSO CON PROBABILITÀ r
        spin_copy(Par, trial, Fields->spin[iSite]);     // PASSO ACCETTATO
        acc=acc+1.0;                                    // INCREMENTO L'ACCETTANZA

        #ifdef DEBUG
        ene2 = ene(Par, Fields);
        delta_glob = ene2 - ene1;
        if(fabs(delta_loc - delta_glob)>1.0e-12){
          if((fabs(delta_loc - delta_glob)>1.0e-12)&(fabs(delta_loc - delta_glob)<1.0e-11)){
            err1 = err1+1;
          }
          else if((fabs(delta_loc - delta_glob)>1.0e-11)){
            err2=err2+1;
          }
        }
        #endif
      }
    }
  }
}

// PROCEDURA PER L'UPDATE OVERRELAXATION
void over_update(SystemParam_t *Par, Field_t *Fields){
  int iSite, j;
  double mean[K], trial[K];
  double norm;

  for(iSite=0;iSite<Par->N;iSite++){          // IL PASSO DI OVERRELAXATION È SVOLTO SEQUENZIALMENTE SU TUTTI I SITI
    #ifdef DEBUG
    double ene1;
    ene1 = ene(Par, Fields);
    #endif
    mean_field(Par, Fields, iSite, mean);     // VALUTO IL CAMPO MEDIO ATTORNO AL SITO iSite
    norm = product(Par, mean, mean);          // CALCOLO |f|^2

    if(norm<1.0e-12){                         // SE LA NORMA È TROPPO PICCOLA SOSTITUISCO CON UN VETTORE RANDOM
      for(j=0;j<K;j++){
        trial[j] = 2*rndm()-1.0;
      }
      for(j=0;j<K;j++){
        trial[j] = trial[j]/sqrt(product(Par, trial, trial));     // NORMALIZZO IL VETTORE RANDOM
      }
      spin_copy(Par, trial, Fields->spin[iSite]);
      return;
    }
    else{
      for(j=0;j<K;j++){
        // COSTRUISCO LO SPIN RIFLESSO CON L'OVERRELAXATION
        trial[j] = 2.0*product(Par, Fields->spin[iSite], mean)*mean[j]/norm - Fields->spin[iSite][j];
      }
      spin_copy(Par, trial, Fields->spin[iSite]);

      #ifdef DEBUG
      double ene2;
      ene2 = ene(Par, Fields);
      if(fabs(ene1-ene2)>1e-12){
        if((fabs(ene1-ene2)>1.0e-12)&(fabs(ene1-ene2)<1.0e-11)){
          err1 = err1+1;
        }
        else if((fabs(ene1-ene2)>1.0e-11)){
          err2=err2+1;
        }
      }
      #endif

    }
  }
}

// PROCEDURA DI TERMALIZZAZIONE
void thermalization(SystemParam_t *Par, Field_t *Fields, int count){
  int i, j;
  int cutoff = 20;
  double a1;
  bool_t ctrl_1, ctrl_2;

  // SE PARTO DA UNA CONFIGURAZIONE SALVATA NON DEVO TERMALIZZARE
  if((Par->iStart) == 1){
    return;
  }
  else if(((Par->iStart) == 0) || ((Par->iStart) == 2)){
    if(count>=cutoff){
      return;
    }
    else{
      // PER UN NUMERO iTerm DI VOLTE RIPETO L'UPDATE METROPOLIS DI TUTTI GLI SPIN
      // ED OGNI VOLTA FACCIO iOverr UPDATE DI OVERRELAXATION
      for(i=0;i<(Par->iTerm); i++){
        metro_update(Par, Fields);
        for(j=0;j<(Par->iOverr);j++){
          over_update(Par, Fields);
        }
        // AGGIUSTO A MANO IL MODULO DEL CAMPO SCALARE CHE POTREBBE ESSERE CAMBIATO PER ERRORE NUMERICO
        if(i==((Par->iTerm)/2)){
          renormalize(Par, Fields);
        }
      }

      a1 = acc/((Par->iTerm)*(Par->N));               // ACCETTANZE NELLA PARTE DI TERMALIZZAZIONE
      ctrl_acceptance(a1, &ctrl_1, &ctrl_2);          // CONTROLLO L'ACCETTANZA
      modify_eps(Par, Fields, ctrl_1, ctrl_2, count); // SE LE ACCETTANZE NON PASSANO IL CHECK MODIFICO eps1
      #ifdef DEBUG
      err1=0;
      err2 =0;
      #endif
    }
  }
  else{
    printf("Errore in iStart\n");
  }
}

// PROCEDURA DI UPDATE DELLE CONFIGURAZIONI, PER IDEC VOLTE
void update_configurations(SystemParam_t *Par, Field_t *Fields){
  int i, j;

  for(i=0;i<(Par->iDec); i++){            // RIPETO PER iDec VOLTE GLI UPDATE: 1 METROPOLIS E iOverr OVERRELAXATION
    metro_update(Par, Fields);
    for(j=0;j<(Par->iOverr);j++){
      over_update(Par, Fields);
    }
  }
  renormalize(Par, Fields);               // RINORMALIZZO TUTTI GLI SPIN PER EVENTUALI ERRORI NUMERICI
}
