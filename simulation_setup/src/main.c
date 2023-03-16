// PROGRAMMA DI SIMULAZIONE NUMERICA
// PER IL MODELLO DI KITAEV-HEISENBERG
// SU RETICOLO ESAGONALE
// L'OBIETTIVO È LA RIPRODUZIONE DEI
// RISULTATI NEL PAPER Phys. Rev. B 88, 024410 (2013)

// INCLUDO HEADER DOVE HO DEFINIZIONI DELLE STRUTTURE, MACRO E PROTOTIPI DELLE FUNZIONI
#include "../head/head_&_structures.h"

// STRUTTURE E VARIABILI GLOBALI
dsfmt_t dsfmt;                          // VARIABILE GLOBALE PER IL GENERATORE DI NUMERI RANDOM
int **nn;                               // ARRAY CHE MEMORIZZANO I PRIMI VICINI, nb PER SOTTORETICOLO NERO
                                        // ED nw PER SOTTORETICOLO BIANCO
double acc=0, err1=0, err2=0;           // VARIABILE GLOBALE PER ACCETTANZA ED ERRORE
clock_t start, end;                     // VARIABILI GLOBALI PER IL TEMPO DI ESECUZIONE

// MAIN DEL PROGRAMMA
int main(int argc, char const *argv[]) {
  SystemParam_t Param;                  // STRUTTURA DEI PARAMETRI
  Field_t Fields;                       // CONFIGURAZIONE DEL CAMPO DI SPIN
  Obs_t Obs;                            // STRUTTURA DELLE OSSERVABILI
  int i, count=0;                       // CONTATORI PER TERMALIZZAZIONE E CICLO DI MISURE
  FILE *fptr;                           // PUNTATORE A FILE PER SCRIVERE I DATI

  // INIZIO DELL'ESECUZIONE
  start = clock();

  // INIZIALIZZO IL SISTEMA
  initializeSystem(&Param, &Fields, &Obs, argv[1]);

  // INIZIALIZZO IL GENERATORE NUMERI NATURALI
  dsfmt_init_gen_rand(&dsfmt, Param.dSFMT_seed);

  // TERMALIZZAZIONE
  thermalization(&Param, &Fields, count);

  // SALVO CONFIGURAZIONE
  writeEps(&Param);
  writeFields(&Param, &Fields);

  // APRO IL FILE SU CUI VERRANNO SCRITTE LE MISURE
  fptr = fopen(Param.data_file, "a");
  if (fptr == NULL) {
    perror("Errore in apertura per la scrittura dati");
    exit(EXIT_FAILURE);
  }
  if ((Param.iStart == 0) || (Param.iStart == 2)){
    // SCRIVO I PARAMETRI SU FILE CHE POI RILEGGO NELL'ANALISI (CELLE PER LATO, NUMERO TOTALE DI SITI, NUMERO DI COMPONENTI DELLO SPIN, TEMPERATURA, ALPHA)
    fprintf(fptr, "%d\t%d\t%d\t%lf\t%lf\n", Param.L, Param.N, K, Param.T, Param.alpha);
    // PRIMA LINEA SU FILE DELLE MISURE PER CAPIRE COSA SONO LE COLONNE DI DATI
    fprintf(fptr, "#ene_dens\tmx\tmy\tmz\tnx\tny\tnz\tzx\tzy\tzz\tsx\tsy\tsz\n\n");
  }

  // CICLO DELLE MISURAZIONI VERO E PROPRIO
  for(i=0;i<(Param.iMis);i++){
    update_configurations(&Param, &Fields);             // UPDATE DELLE CONFIGURAZIONI PER iDec VOLTE PRIMA DELLA MISURA
    measures(&Param, &Fields, &Obs);                    // MISURE DELLE VARIE GRANDEZZE SU RETICOLO
    writeObs(fptr, &Obs);                               // SALVATAGGIO DELLE MISURE
    if((i % (Param.iBackup)) == 0){
      writeFields(&Param, &Fields);                     // OGNI iBackup MISURE SALVO UNA CONFIGURAZIONE
    }
  }
  writeFields(&Param, &Fields);                         // SALVO LA CONFIGURAZIONE PRIMA DI CHIUDERE

  // CHIUDO IL FILE DELLE MISURE, DEALLOCO LA MEMORIA E SCRIVO IL FILE DI LOG
  fclose(fptr);
  deallocation(&Fields, &Obs);
  writeLogs(&Param);

  // FINE PROGRAMMA
  return(EXIT_SUCCESS);
}
