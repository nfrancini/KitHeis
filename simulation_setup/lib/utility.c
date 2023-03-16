#include "../head/head_&_structures.h"

// DEFINISCO UNA FUNZIONE RNDM CHE RICHIAMA IL DSFMT.
// QUESTO È UTILE PER CAMBIARE GENERATORE SENZA DOVERE OGNI
// VOLTA CORREGGERE IL CODICE
double rndm(){
  double x;

  x = dsfmt_genrand_close_open(&dsfmt);
  return x;
}

// FUNZIONE CHE RESETTA LE VARIABILI DI ERROR HANDLING
// VIENE RICHIAMATA NELLE VARIE FUNZIONI SE DEBUG È DEFINITO
// È IMPORTANTE RESETTARE OGNI VOLTA LE VARIABILI DI ERROR HANDLING ALTRIMENTI
// IN PRESENZA DI UN ERRORE NON SI RIESCE A CAPIRE DA DOVE ARRIVI
void resetErr(){
  errno = 0;
  feclearexcept(FE_ALL_EXCEPT);
}

// FUNZIONE CHE CONTROLLA IL RANGE DELL'ACCETTANZA
// RICHIAMATA NELLA TERMALIZZAZIONE, IN MODO
// DA IMPOSTARE AUTOMATICAMENTE IL PARAMETRO EPS1
// CHE REGOLA L'ACCETTANZA
// IN QUESTO CASO HO PRESO UN LIMITE INFERIORE DI 0.30 E SUPERIORE DI 0.45
void ctrl_acceptance(double ac, bool_t *ctrl_1, bool_t *ctrl_2){
  double a_min = 0.30;
  double a_max = 0.45;

  if((ac>=a_min)&&(ac<=a_max)){       // ACCETTANZA BUONA (TRUE, FALSE)
    *ctrl_1 = TRUE;
    *ctrl_2 = FALSE;
  }
  else if(ac<a_min){                  // ACCETTANZA TROPPO PICCOLA (FALSE, FALSE)
    *ctrl_1 = FALSE;
    *ctrl_2 = FALSE;
  }
  else if(ac>a_max){                  // ACCETTANZA TROPPO GRANDE (TRUE, TRUE)
    *ctrl_1 = TRUE;
    *ctrl_2 = TRUE;
  }
}

// PROCEDURA PER MODIFICA AUTOMATICA DEI PARAMETRI DI ACCETTANZA
// SE QUESTA È TROPPO PICCOLA O TROPPO GRANDE VIENE MODIFICATO
// IL PARAMETRO eps1. VIENE QUINDI CHIAMATA RICORSIVAMENTE
// LA TERMALIZZAZIONE, AUMENTATO IL CONTATORE count
void modify_eps(SystemParam_t *Par, Field_t *Fields, bool_t ctrl_1, bool_t ctrl_2, int count){
  double step=0.05;
  if((Par->eps1+step>=1.0) || (Par->eps1-step)<=0){
    printf("IMPOSSIBILE MODIFICARE IL PARAMETRO DI ACCETTANZA\n");
    exit(EXIT_FAILURE);
  }
  if((ctrl_1==TRUE)&&(ctrl_2==FALSE)){          // ACCETTANZA BUONA
    acc=0;
    return;
  }
  else if((ctrl_1==FALSE)&&(ctrl_2==FALSE)){     // ACCETTANZA TROPPO PICCOLA
    Par->eps1 = Par->eps1 - step;
    acc = 0;
    thermalization(Par, Fields, count+1);
  }
  else if((ctrl_1==TRUE)&&(ctrl_2==TRUE)){     // ACCETTANZA TROPPO GRANDE
    Par->eps1 = Par->eps1 + step;
    acc = 0;
    thermalization(Par, Fields, count+1);
  }
  else{
    printf("ERRORE NELLA TERMALIZZAZIONE\n");
    exit(EXIT_FAILURE);
  }
}
