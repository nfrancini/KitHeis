#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "dSFMT.h"
#include <math.h>
#include <string.h>
#include <complex.h>
#include <errno.h>
#include <fenv.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "../config.h"

// MACRO
// #define DEBUG                 // MACRO DI DEBUG
#define D 3                   // NUMERO DI PRIMI VICINI PER OGNI SITO
#define K 3                   // NUMERO DI COMPONENTI DELLO SPIN
#define PI 3.14159265358979   // PI GRECO
#define STD_STRING_LENGTH 50  // LUNGHEZZA STANDARD DELLE STRINGHE DI LUNGHEZZA IGNOTA

typedef struct {  // DEFINISCO IL TIPO SystemParam_t STRUTTURATO NEL SEGUENTE MODO
  int L;          // NUMERO DI CELLE PER LATO
  int V;          // VOLUME DEL SISTEMA, CIOÈ NUMERO TOTALE DI CELLE L*L
  int N;          // NUMERO TOTALE DI SITI (2*V)
  double T;       // TEMPERATURA DEL SISTEMA
  double alpha;   // PARAMETRO CHE REGOLA L'ACCOPPIAMENTO DI KITAEV (-2*alpha) E DI HEISENBERG (1-alpha)
  double eps1;    // PARAMETRO CHE REGOLA L'ACCETTANZA PER LE VARIABILI DI SPIN
  int iTerm;      // NUMERO DI UPDATE DI TERMALIZZAZIONE
  int iDec;       // NUMERO DI UPDATE PER DECORRELARE LA CATENA
  int iMis;       // NUMERO DI MISURAZIONI
  int iOverr;     // NUMERO DI PASSI DI OVERRELAXATION
  int iStart;     // INDICE PER PARTENZA ORDINATA (iStart = 0), DA FILE PRECEDENTE (iStart = 1) O IN MODO RANDOM (iStart = 2)
  int iBackup;    // OGNI iBackup MISURE SALVO UNA CONFIGURAZIONE

  char conf_file[STD_STRING_LENGTH];      // NOME DEL FILE DI SALVATAGGIO DELLE CONFIGURAZIONI
  char eps_file[STD_STRING_LENGTH];       // NOME DEL FILE DI SALVATAGGIO DEI PARAMETRI DI ACCETTANZA
  char data_file[STD_STRING_LENGTH];      // NOME DEL FILE DI SALVATAGGIO DEI DATI
  char log_file[STD_STRING_LENGTH];       // NOME DEL FILE DI SALVATAGGIO DEL LOG

  unsigned int dSFMT_seed;                // SEED PER IL DSFMT
} SystemParam_t;

typedef struct {
  double **spin;                          // CAMPO DI SPIN SPIN[N][K]
} Field_t;

typedef struct {                          // DEFINISCO LA STRUTTURA DELLE OSSERVABILI
  double ene_density;                     // DENSITÀ DI ENERGIA CALCOLATA A PARTIRE DALL'HAMILTONIANA
  double *N;                              // PARAMETRO D'ORDINE PER LA FASE DI NEEL ANTI-FERROMAGNETICA
  double *S;                              // PARAMETRO D'ORDINE PER LA FASE STRIPY ANTIFERROMAGNETICA
  double *M;                              // PARAMETRO D'ORDINE PER LA FASE FERROMAGNETICA
  double *Z;                              // PARAMETRO D'ORDINE PER LA FASE ZIGZAG
} Obs_t;

extern dsfmt_t dsfmt;                     // VARIABILE DI SUPPORTO PER IL DSFMT
extern int **nn;                          // VARIABILE ESTERNA CHE NEL MAIN È GLOBALE
                                          // E CHE MEMORIZZA I PRIMI VICINI
extern double acc, err1, err2;            // VARIABILI GLOBALI PER ACCETTANZA ED ERRORI

typedef enum {FALSE=0, TRUE=1} bool_t;    // DEFINISCO IL TIPO BOOLEANO, UTILE PER ALCUNE PROCEDURE


// FUNZIONI E PROCEDURE AUSILIARIE
double rndm();                            // FUNZIONE CHE GENERA UN NUMERO RANDOM TRA 0 ED 1 CON IL DSFMT
void resetErr();                          // FUNZIONE CHE RESETTA LE VARIABILI DI ERROR HANDLING
void ctrl_acceptance(double ac, bool_t *ctrl_1, bool_t *ctrl_2);        // PROCEDURA CHE CONTROLLA L'ACCETTANZA IN FASE DI TERMALIZZAZIONE
void modify_eps(SystemParam_t *Par, Field_t *Fields, bool_t ctrl_1, bool_t ctrl_2, int count);    // PROCEDURA CHE MODIFICA IL PARAMETRO DI ACCETTANZA IN FASE DI TERMALIZZAZIONE

// FUNZIONI DI INIZIALIZZAZIONE
void geometry(SystemParam_t *Par);        // CREA IL VETTORE CHE MEMORIZZA TUTTI I PRIMI VICINI CON PBC
void allocation(SystemParam_t *Par, Field_t *Fields, Obs_t *Obs);       // ALLOCA LA MEMORIA PER I PARAMETRI DI SISTEMA, IL CAMPO DI SPIN E LE OSSERVABILI
void initializeFields(SystemParam_t *Par, Field_t *Fields);             // INIZIALIZZA I CAMPI IN MODO ORDINATO (iStart=0) O DA FILE (iStart=1)
void initializeObs(Obs_t *Obs);           // INIZIALIZZA A ZERO TUTTE LE OSSERVABILI PER SICUREZZA
void deallocation(Field_t *Fields, Obs_t *Obs);       // DEALLOCA LA MEMORIA ALLOCATA
void initializeSystem(SystemParam_t *Par, Field_t *Fields, Obs_t *Obs, char const *finput);    // METTE INSIEME TUTTO

// FUNZIONI E PROCEDURE CHE REGOLANO IL SALVATAGGIO E LA LETTURA DA FILE DEL SISTEMA
void remove_white_line_and_comments(FILE *input);                         // PROCEDURA CHE RIMUOVE GLI SPAZI, LE LINEE ED I COMMENTI DAL FILE DI INPUT
void read_from_input_Param(SystemParam_t *Par, char const *finput);       // PROCEDURA CHE LEGGE I PARAMETRI DEL SISTEMA DAL FILE DI INPUT
void readFields(SystemParam_t *Par, Field_t *Fields);                     // PROCEDURA CHE LEGGE LE CONFIGURAZIONI DI SPIN PRECEDENTEMENTE SALVATE
void writeFields(SystemParam_t *Par, Field_t *Fields);                    // PROCEDURA CHE SALVA SU FILE LE CONFIGURAZIONI DI SPIN IN BINARIO
void writeEps(SystemParam_t *Par);                                        // PROCEDURA CHE SALVA SU FILE IL PARAMETRO DI ACCETTANZA
void readEps(SystemParam_t *Par);                                         // PROCEDURA CHE LEGGE DA FILE IL PARAMETRO DI ACCETTANZA PRECEDENTEMENTE SALVATO
void writeObs(FILE *fptr, Obs_t *Obs);                                    // PROCEDURA CHE SALVA SU FILE LE MISURE FATTE
void writeLogs(SystemParam_t *Par);                                       // PROCEDURA CHE SCRIVE SU UN FILE DI LOG LE CARATTERISTICHE DELLA SIMULAZIONE

// FUNZIONI E PROCEDURE CHE MANIPOLANO GLI SPIN
double product(SystemParam_t *Par, double *vec1, double *vec2);                                     // PROCEDURA CHE IMPLEMENTA IL PRODOTTO SCALARE TRA DUE SPIN
void diff(SystemParam_t *Par, double *vec1, double *vec2, double *diff);                            // PROCEDURA CHE IMPLEMENTA LA DIFFERENZA TRA DUE SPIN
void two_spin_sum(SystemParam_t *Par, double *vec1, double *vec2, double *somma);                   // PROCEDURA CHE IMPLEMENTA LA SOMMA TRA DUE SPIN
void three_spin_sum(SystemParam_t *Par, double *vec1, double *vec2, double *vec3, double *somma);   // PROCEDURA CHE IMPLEMENTA LA SOMMA TRA TRE SPIN
void renormalize(SystemParam_t *Par, Field_t *Fields);                                              // PROCEDURA CHE RINORMALIZZA A MODULO UNITARIO TUTTI GLI SPIN
void spin_copy(SystemParam_t *Par, double *vec_in, double *vec_out);                                // PROCEDURA CHE COPIA LO SPIN vec_in NELLO SPIN vec_out

// FUNZIONI E PROCEDURE PER L'UPDATE METROPOLIS
void spin_trial(SystemParam_t *Par, Field_t *Fields, int site, double *trial);        // PROCEDURA CHE CREA LO SPIN DI PROVA
void metro_update(SystemParam_t *Par, Field_t *Fields);                               // PROCEDURA CHE IMPLEMENTA L'UPDATE METROPOLIS SU TUTTO IL RETICOLO
void mean_field(SystemParam_t *Par, Field_t *Fields, int site, double *mean_field);   // PROCEDURA CHE CALCOLA IL CAMPO MEDIO INTORNO AD UN DATO SITO
void over_update(SystemParam_t *Par, Field_t *Fields);                                // PROCEDURA CHE IMPLEMENTA L'UPDATE OVERRELAXATION SU TUTTO IL RETICOLO
void thermalization(SystemParam_t *Par, Field_t *Fields, int count);                  // PROCEDURA CHE IMPLEMENTA IL CICLO DI TERMALIZZAZIONE
void update_configurations(SystemParam_t *Par, Field_t *Fields);                      // PROCEDURA CHE IMPLEMENTA L'UPDATE METROPOLI+OVERRELAXATION

// FUNZIONI E PROCEDURE PER IMPLEMENTARE LE MISURE
double ene(SystemParam_t *Par, Field_t *Fields);                                          // CALCOLA L'ENERGIA DEL SISTEMA
double ene_dens(SystemParam_t *Par, Field_t *Fields);                                     // CALCOLA LA DENSITÀ DI ENERGIA
void FM_magn(SystemParam_t *Par, Field_t *Fields, double *m);                             // CALCOLA L'ARRAY DELLA MAGNETIZZ FERROMAGNETICA
void neel_magn(SystemParam_t *Par, Field_t *Fields, double *n);                           // CALCOLA L'ARRAY DEL PARAMETRO D'ORDINE DI NEEL
void zigzag_and_stripy_magn(SystemParam_t *Par, Field_t *Fields, double *z, double *s);   // CALCOLA L'ARRAY DEL PARAMETRO D'ORDINE ZIGZAG E STRIPY
void measures(SystemParam_t *Par, Field_t *Fields, Obs_t *Obs);                           // ESEGUE LE DIVERSE MISURE
