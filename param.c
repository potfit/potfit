
/******************************************************************************
*
*  param.c: Read in parameter files (tag based)
* 
*  $Revision: 1.10 $
*  $Date: 2003/04/04 09:29:13 $
*
******************************************************************************/

#include "potfit.h"

/* the following are needed for gettimeofday */
#include <sys/time.h>
#include <unistd.h>

#if defined(__GNUC__) && defined(__STRICT_ANSI__)
extern char *strdup(char *);
#endif

#ifdef __WATCOMC__
#define strcasecmp strcmpi
#endif

typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_DOUBLE
} PARAMTYPE;

int curline; /* number of current line */

/******************************************************************************
*
*  get parameters from one line
*
******************************************************************************/

/* Parameter:

   param_name ... Parametername (fuer Fehlermeldungen)
   param ........ Adresse der Variable fuer den Parameter
   ptype ........ Parametertyp
                  folgende Werte sind zulaessig:
                  PARAM_STR : String, deklariert als char[]
                  PARAM_STRPTR : String, deklariert als Zeiger auf char*
                  PARAM_INT : Integer-Wert(e)
                  PARAM_DOUBLE : Double-Wert(e)
                  
   pnum_min ..... Minimale Anzahl der einzulesenden Werte
                  (Falls weniger Werte gelesen werden koennen als verlangt,
                  wird ein Fehler gemeldet).
   pnum_max ..... Maximale Anzahl der einzulesenden Werte
                  (Ueberzaehlige Werte werden ignoriert. Falls weniger als
                  pnum_max Werte vorhanden sind, wird das Lesen
                  abgebrochen, es wird aber kein Fehler gemeldet,
                  wenn mindestens pnum_min Werte abgesaettigt wurden.

  Resultat: Die Anzahl der gelesenen Werte wird zurueckgegeben.

*/

int getparam(char *param_name, void *param, PARAMTYPE ptype, 
             int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"Parameter for %s missing in line %u\nstring expected!\n",
              param_name,curline);
      error(errmsg);
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"Parameter for %s missing in line %u\nstring expected!\n",
              param_name,curline);
      error(errmsg);
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\r\n");
      if (str == NULL) {
        sprintf(errmsg,"Parameter for %s missing in line %u!\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),
                "Integer vector of length %u expected!\n",(unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\r\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_DOUBLE) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\r\n");
      if (str == NULL) {
        sprintf(errmsg,"Parameter for %s missing in line %u!\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"Double vector of length %u expected!\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((double*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\r\n")) != NULL) {
        ((double*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  return numread;
} 

/******************************************************************************
*
*  read in parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
******************************************************************************/

void read_paramfile(FILE *pf)
{
  char buffer[1024];
  char *token, *res;

  curline=0;

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\r\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    /* number of atom types */
    if (strcasecmp(token,"ntypes")==0) {
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
    }
    /* file with start potential */
    else if (strcasecmp(token,"startpot")==0) {
      getparam("startpot",startpot,PARAM_STR,1,255);
    }
    /* file for end potential */
    else if (strcasecmp(token,"endpot")==0) {
      getparam("endpot",endpot,PARAM_STR,1,255);
    }
    /* file for IMD potential */
    else if (strcasecmp(token,"imdpot")==0) {
      getparam("imdpot",imdpot,PARAM_STR,1,255);
    }
    /* file for plotting */
    else if (strcasecmp(token,"plotfile")==0) {
	getparam("plotfile",plotfile,PARAM_STR,1,255);
	plot=1;
    }
    /* number of steps in IMD potential */
    else if (strcasecmp(token,"imdpotsteps")==0) {
      getparam("imdpotsteps",&imdpotsteps,PARAM_INT,1,1);
    }
    /* file with atom configuration */
    else if (strcasecmp(token,"config")==0) {
      getparam("config",config,PARAM_STR,1,255);
    }
    /* Optimization flag */
    else if (strcasecmp(token,"opt")==0) {
	getparam("opt",&opt,PARAM_INT,1,1);
    }
    /* break flagfile */
    else if (strcasecmp(token,"flagfile")==0) {
	getparam("flagfile",flagfile,PARAM_STR,1,255);
    }
    /* plotpoint file */
    else if (strcasecmp(token,"plotpointfile")==0) {
	getparam("plotpointfile",plotpointfile,PARAM_STR,1,255);
    }
    /* temporary potential file */
    else if (strcasecmp(token,"tempfile")==0) {
	getparam("tempfile",tempfile,PARAM_STR,1,255);
    }
    /* seed for RNG */
    else if (strcasecmp(token,"seed")==0){
	getparam("seed",&seed,PARAM_INT,1,1);
    }
    /* starting temperature for annealing */
    else if (strcasecmp(token,"anneal_temp")==0) {
      getparam("anneal_temp",&anneal_temp,PARAM_DOUBLE,1,1);
    }
    /* unknown tag */
    else {
      fprintf(stderr,"Unknown tag <%s> ignored!\n",token);
    }
  } while (!feof(pf));
  fclose(pf);

}

/******************************************************************************
*
*  read command line
*
******************************************************************************/

void read_parameters(int argc, char **argv)
{
  char msg[255];
  FILE *pf;

  /* check command line */
  if (argc<2) {
    sprintf(msg, "Usage: %s <paramfile>\n", argv[0]);
    error(msg);
  }

  /* open parameter file, and read it */
  pf = fopen(argv[1], "r");
  if (NULL == pf) {
    fprintf(stderr, "ERROR: Could not open parameter file %s!\n", argv[1] );
    exit(2);
  }
  read_paramfile(pf);
}
