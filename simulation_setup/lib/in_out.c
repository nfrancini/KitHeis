#include "../head/head_&_structures.h"

// PROCEDURA PER ELIMINARE GLI SPAZI ED I COMMENTI DAL FILE DI INPUT
void remove_white_line_and_comments(FILE *input){
  int temp_i;

  temp_i=getc(input);
  if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // SCAN PER SPAZI VUOTI E COMMENTI
    {
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\n' || temp_i==' ') // LINEE VUOTE
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i=='\n' || temp_i==' ');
      }
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\043')  // COMMENTI, 043 È IL CODICE ASCII PER #
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i!='\n');
      }
    else
      {
      ungetc(temp_i, input);
      }

    remove_white_line_and_comments(input);
    }
  else
    {
    ungetc(temp_i, input);
    }
  }

// PROCEDURA CHE LEGGE DA FILE I PARAMETRI DI SISTEMA
void read_from_input_Param(SystemParam_t *Par, char const *finput){
  int flag, temp_i;
  int end=1;
  FILE *fInput;
  char str[100];

  fInput = fopen(finput, "r");
  if (fInput == NULL) {
    perror("Errore in apertura in read_from_input_Param");
    exit(EXIT_FAILURE);
  }
  else{
    while(end==1){
      remove_white_line_and_comments(fInput);

      flag = fscanf(fInput, "%s", str);
      if (flag != 1) {
        perror("Errore di lettura della prima stringa");
        exit(EXIT_FAILURE);
      }

      if(strncmp(str, "size", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->L));
        if (flag != 1){
          perror("Errore di lettura nella taglia");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "T", 1) == 0){
        flag = fscanf(fInput, "%lf", &(Par->T));
        if (flag != 1){
          perror("Errore di lettura nella temperatura");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "phi", 3) == 0){
        flag = fscanf(fInput, "%lf", &(Par->phi));
        if (flag != 1){
          perror("Errore di lettura in alpha");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "eps_spin", 8) == 0){
        flag = fscanf(fInput, "%lf", &(Par->eps1));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iTerm", 5) == 0){
        flag = fscanf(fInput, "%d", &(Par->iTerm));
        if (flag != 1){
          perror("Errore di lettura in iTerm");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iDec", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->iDec));
        if (flag != 1){
          perror("Errore di lettura in iDec");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iMis", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->iMis));
        if (flag != 1){
          perror("Errore di lettura in iMis");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iOverr", 6) == 0){
        flag = fscanf(fInput, "%d", &(Par->iOverr));
        if (flag != 1){
          perror("Errore di lettura in iOverr");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iStart", 6) == 0){
        flag = fscanf(fInput, "%d", &(Par->iStart));
        if (flag != 1){
          perror("Errore di lettura in iStart");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iBackup", 7) == 0){
        flag = fscanf(fInput, "%d", &(Par->iBackup));
        if (flag != 1){
          perror("Errore di lettura in iBackup");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "conf_file", 9) == 0){
        flag = fscanf(fInput, "%s", Par->conf_file);
        if (flag != 1){
          perror("Errore di lettura in conf_file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "eps_file", 8) == 0){
        flag = fscanf(fInput, "%s", Par->eps_file);
        if (flag != 1){
          perror("Errore di lettura in eps file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "data_file", 9) == 0){
        flag = fscanf(fInput, "%s", Par->data_file);
        if (flag != 1){
          perror("Errore di lettura in data_file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "log_file", 8) == 0){
        flag = fscanf(fInput, "%s", Par->log_file);
        if (flag != 1){
          perror("Errore di lettura in log_file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "random_seed", 11) == 0){
        flag = fscanf(fInput, "%u", &(Par->dSFMT_seed));
        if (flag != 1){
          perror("Errore di lettura in random_seed");
          exit(EXIT_FAILURE);
        }
      }
      else{
        fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, finput, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }

      remove_white_line_and_comments(fInput);

      // CONTROLLO SE LA LINEA LETTA È L'ULTIMA
      temp_i = getc(fInput);
      if(temp_i == EOF){
        end=0;
      }
      else{
        ungetc(temp_i, fInput);
      }
    }

    fclose(fInput);
    Par->V = (int)pow(Par->L, 2);    // CALCOLO QUA IL "VOLUME" CIOÈ IL NUMERO TOTALE DI CELLE
    Par->N = 2*(int)pow(Par->L, 2);  // CALCOLO IL NUMERO DI SITI TOTALI CIOÈ 2*L*L
  }
}

// PROCEDURA CHE COPIA LOCALMENTE LA STRUTTURA DEGLI SPIN
// QUESTA VIENE POI SALVATA IN BINARIO
void writeFields(SystemParam_t *Par, Field_t *Fields){
  int iSite, j, mu;
  FILE *fptr;

  typedef struct{                           // DEFINISCO UNA STRUTTURA DI CAMPO COPIA CON DIMENSIONI ASSEGNATE
    double spin[Par->N][K];
  } Field_cpy_t;

  Field_cpy_t Fields_cpy;

  for(iSite=0;iSite<(Par->N);iSite++){      // COPIO LA CONFIGURAZIONE DI CAMPO
    for(j=0;j<K;j++){
      Fields_cpy.spin[iSite][j] = Fields->spin[iSite][j];
    }
  }

  fptr = fopen(Par->conf_file, "wb");       // SCRIVO SU FILE LA CONFIGURAZIONE IN BINARIO
  if (fptr == NULL) {
    perror("Errore in apertura");
    exit(1);
  }
  fwrite(&Fields_cpy, sizeof(Field_cpy_t), 1, fptr);
  fclose(fptr);
}

// PROCEDURA CHE LEGGE DA FILE LA CONFIGURAZIONE DEI CAMPI
// PRECEDENTEMENTE SALVATA E LA METTE NELLA CONFIGURAZIONE
// CORRENTE
void readFields(SystemParam_t *Par, Field_t *Fields){
  int iSite, mu, j, t;
  FILE *fptr;

  typedef struct{                           // DEFINISCO UNA STRUTTURA DI CAMPO COPIA CON DIMENSIONI ASSEGNATE
    double spin[Par->N][K];
  } Field_cpy_t;

  Field_cpy_t Fields_cpy;

  fptr = fopen(Par->conf_file, "rb");       // LEGGO LA CONFIGURAZIONE IN BINARIO CHE HO PRECEDENTEMENTE SALVATO
  if (fptr == NULL) {
    perror("Errore in apertura");
    exit(1);
  }

  j = fread(&Fields_cpy, sizeof(Field_cpy_t), 1, fptr);
  if (j!=1){
    printf("ERRORE NELLA LETTURA DA FILE\n");
    exit(EXIT_FAILURE);
  }
  fclose(fptr);

  for(iSite=0;iSite<(Par->N);iSite++){                // LEGGO DA FILE LA CONFIGURAZIONE
    for(t=0;t<K;t++){
      Fields->spin[iSite][t] = Fields_cpy.spin[iSite][t];
    }
  }
}

// PROCEDURA CHE COPIA LOCALMENTE LA STRUTTURA DEI
// PARAMETRI DI ACCETTANZA
void writeEps(SystemParam_t *Par){
  FILE *fptr;

  typedef struct{               // DEFINISCO PER COMODITÀ UNA STRUTTURA COPIA DEI PARAMETRI DI ACCETTANZA
    double eps1;
  } eps_cpy_t;

  eps_cpy_t eps_cpy;

  fptr = fopen(Par->eps_file, "wb");
  if (fptr == NULL) {
    perror("Errore in apertura");
    exit(1);
  }

  eps_cpy.eps1 = Par->eps1;     // COPIO SU FILE I PARAMETRI DI ACCETTANZA

  fwrite(&eps_cpy, sizeof(eps_cpy), 1, fptr);
  fclose(fptr);
}

// PROCEDURA CHE LEGGE DA FILE BINARIO I PARAMETRI DI ACCETTANZA PRECEDENTI
void readEps(SystemParam_t *Par){
  int flag;
  FILE *fptr;

  typedef struct{                 // DEFINISCO PER COMODITÀ UNA STRUTTURA COPIA DEI PARAMETRI DI ACCETTANZA
    double eps1;
  } eps_cpy_t;

  eps_cpy_t eps_cpy;

  fptr = fopen(Par->eps_file, "rb");
  if (fptr == NULL) {
    perror("Errore in apertura");
    exit(1);
  }

  flag = fread(&eps_cpy, sizeof(eps_cpy), 1, fptr);
  if(flag != 1){
    printf("ERRORE NELLA LETTURA DEI FILE\n");
    exit(EXIT_FAILURE);
  }
  fclose(fptr);

  Par->eps1 = eps_cpy.eps1;       // ASSEGNO I VALORI LETTI DA FILE
}

// PROCEDURA CHE SCRIVE LE OSSERVABILI SU FILE DOPO OGNI MISURA
void writeObs(FILE *fptr, Obs_t *Obs){
  fprintf(fptr, "%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", Obs->ene_density, Obs->M[0], Obs->M[1], Obs->M[2], Obs->N[0], Obs->N[1], Obs->N[2], Obs->Z[0], Obs->Z[1], Obs->Z[2], Obs->S[0], Obs->S[1], Obs->S[2]);
}

// PROCEDURA CHE SCRIVE IL FILE DI LOG
void writeLogs(SystemParam_t *Par){
  FILE *fptr;
  double exe_time, seconds;
  int hours, minutes;
  fptr = fopen(Par->log_file, "a");
  if(fptr == NULL){
    perror("Errore in apertura");
    exit(EXIT_FAILURE);
  }

  fprintf(fptr, "+---------------------------------+\n");
  fprintf(fptr, "| Simulation details for kit_heis |\n");
  fprintf(fptr, "+---------------------------------+\n");

  #ifdef DEBUG
  fprintf(fptr, "DEBUG attivo\n\n");
  #endif

  fprintf(fptr, "RETICOLO ESAGONALE CON CELLE UNITARIE CONTENENTI 2 SITI.\n");
  fprintf(fptr, "CELLE PER LATO DEL RETICOLO %d\n", Par->L);
  fprintf(fptr, "NUMERO DI COMPONENTI DELLO SPIN %d\n\n", K);

  fprintf(fptr, "T %.5lf\n", Par->T);
  fprintf(fptr, "phi %.5lf\n\n", Par->phi);

  fprintf(fptr, "eps %.5lf\n\n", Par->eps1);

  fprintf(fptr, "iTerm %d\n", Par->iTerm);
  fprintf(fptr, "iDec %d\n", Par->iDec);
  fprintf(fptr, "iMis %d\n", Par->iMis);
  fprintf(fptr, "iOverr %d\n\n", Par->iOverr);

  fprintf(fptr, "iStart %d (0=ordinata, 1=da file , 2=random)\n", Par->iStart);
  fprintf(fptr, "iBackup %d\n\n", Par->iBackup);

  #ifdef DEBUG
  fprintf(fptr, "1.0e-12<ERRORE<1.0E-11 = %lf\n", err1/((Par->iOverr)*D*(Par->V)*((Par->iDec))));
  fprintf(fptr, "ERRORE>1.0E-11 = %lf\n\n", err2/((Par->iOverr)*D*(Par->V)*((Par->iDec))));
  #endif

  fprintf(fptr, "ACCETTANZA UPDATE METROPOLIS %lf\n\n", acc/((Par->iMis)*(Par->iDec)*(Par->N)));

  // CALCOLO IL TEMPO DI ESECUZIONE
  end = clock();
  exe_time = ((double)end - start)/CLOCKS_PER_SEC;
  minutes = (int)(exe_time/60);
  seconds = exe_time - minutes*60;
  hours = (int)(minutes/60);
  minutes = minutes - hours*60;

  fprintf(fptr, "TEMPO DI ESECUZIONE %d h %d m %lf s\n\n", hours, minutes, seconds);
  fprintf(fptr, "PROGRAMMA ESEGUITO CORRETTAMENTE\n");

  fclose(fptr);
}
