
/******************************************************************************
*
*  IMD -- The ITAP Molecular Dynamics Program
*
*  Copyright 1996-2001 Institute for Theoretical and Applied Physics,
*  University of Stuttgart, D-70550 Stuttgart
*
*  $Revision: 1.1 $
*  $Date: 2004/04/13 14:13:51 $
*
******************************************************************************/

/******************************************************************************
*
*  The utility program pottrans translates potential tables in IMD 
*  formats 1 or 2 into format 3 required by the potfit program. The 
*  input potentials must have the appropriate header for their format.
*
*  Compilation:  make pottrans
*
*  Usage:        pottrans <paramfile>
*
*  pottrans requires the following parameters (here with example values), 
*  which must be given in the parameter file:
*
*  ntypes     2           # number of atom types
*  ncols      3           # number of output potentials 
*  infile     in.pot      # input potential file
*  outfile    out.pot     # output potential file
*  plotfile   plot.pot    # file name for plot version of output potentials
*  nsteps     20 20 20    # number of steps for the output potentials
*
*
*  Pair interactions (default)
*
*  ncols must be equal to (ntypes+1)*ntypes/2. The output contains ncols
*  potential functions. The input functions must vanish at their end. 
*  By default, the output functions are tabulated from the beginning
*  to the end of the corresponding input functions. Optionally, with the
*  parameter r_start modified values for the start of the output functions
*  can be specified. 
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* maximum value in potential file */
#define MAXVAL 40.

/* type of interpolation */
#define POTVAL pot3

/* type of potential */
#define PAIR   1
#define EAM2    2

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define PTR_2D(var,i,j,dim_i,dim_j) (((var) + ((i)*(dim_j)) + (j)))

/* memory allocation increment for potential */
#define PSTEP 50

#define NPLOT 1000

typedef char str255[255];

typedef double real;

/* parameter types */
typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;

/* data structure to store a potential table or a function table */ 
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table (followed by extra zeros) */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  maxsteps;    /* physical length of the table */
  real *table;      /* the actual data */
} pot_table_t;

/* global variables */
str255 infilename="\0", outfilename="\0", plotfilename="\0", 
  plotpointfile="\0", rho_r2_file="\0", f_rho_file="\0", phi_r2_file="\0";
int    ncols=0, ntypes=0, mode=1, *nsteps=NULL, eam;
real   *r_start=NULL, *r_end=NULL;
pot_table_t pt,embed_pt,rho_tab;

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
  exit(2);
}

/*****************************************************************************
*
*  read one parameter (routine from imd_param.c)
*
*****************************************************************************/

int getparam(int curline, char *param_name, void *param, 
             PARAMTYPE ptype, int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected");
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected");
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((int*)param)[i] = ival;
    }
  }
  else if (ptype == PARAM_INTEGER) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int *)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int *)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((int *)param)[i] = (int)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((real*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        rval = atof(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((real*)param)[i] = rval;
    }
  }
  return numread;
} /* getparam */

/*****************************************************************************
*
*  read tag-based parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;
  str255 tmp;
  int  n=0;

  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(tmp,"Cannot open parameter file %s",paramfname);
    error(tmp);
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    n++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"infile")==0) {
      /* file name for input potential */
      getparam(n,"infile",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"atomic_e-density_file")==0) {
      /* file name for input potential */
      getparam(n,"atomic_e-density_file",rho_r2_file,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* file name for input potential */
      getparam(n,"embedding_energy_file",f_rho_file,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* file name for input potential */
      getparam(n,"core_potential_file",phi_r2_file,PARAM_STR,1,255);
    }
   else if (strcasecmp(token,"outfile")==0) {
      /* file name for output potential */
      getparam(n,"outfile",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"plotfile")==0) {
      /* file name for output potential in plot format */
      getparam(n,"plotfile",plotfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam(n,"ntypes",&ntypes,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"ncols")==0) {
      /* number of output columns */
      getparam(n,"ncols",&ncols,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nsteps")==0) {
      if (ncols==0) error("specify ncols before nsteps");
      nsteps = (int *) malloc( ncols * sizeof(int) );
      /* number of output potential values */
      getparam(n,"nsteps",nsteps,PARAM_INT,ncols,ncols);
    }
    else if (strcasecmp(token,"r_start")==0) {
      /* minimal radius */
      if (ncols==0) error("specify ncols before r_start");
      r_start = (real *) malloc( ncols * sizeof(real) );
      if (r_start==NULL) error("allocation of r_start failed");
      getparam(n,"r_start",r_start,PARAM_REAL,ncols,ncols);
    }
    else if (strcasecmp(token,"r_end")==0) {
      /* maximal radius */
      if (ncols==0) error("specify ncols before r_end");
      r_end = (real *) malloc( ncols * sizeof(real) );
      if (r_end==NULL) error("allocation of r_end failed");
      getparam(n,"r_end",r_end,PARAM_REAL,ncols,ncols);
    }
    else if (strcasecmp(token,"plotpoints")==0) {
      /* file for plotpoints */
	getparam(n,"plotpoints",plotpointfile,PARAM_STR,1,255);
    }
  } while (!feof(pf));
  fclose(pf);

} /* getparamfile */

/******************************************************************************
*
*   read the command line parameters, and the parameter file given 
*   on the command line
*
******************************************************************************/

void read_parameters(int argc,char **argv)
{
  str255 progname, paramfilename, msg;

  strcpy(progname,argv[0]);
  if (argc<2) {
    sprintf(msg,"Parameter file missing!\nUsage: %s <paramfile>",progname);
    error(msg);
  }
  strcpy(paramfilename,argv[1]);
  getparamfile(paramfilename);
} 

/*****************************************************************************
*
*  read potential in first format: each line contains
*
*  r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
*
*  N is the number of different atom types
*
*  Note that it is assumed that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table1(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k;
  int tablesize, npot=0;
  real val, numstep, delta;
  real r2, r2_start, r2_step;
  str255 msg;

  /* allocate the function table */
  pt->maxsteps = PSTEP;
  tablesize = ncols * pt->maxsteps;
  pt->table = (real *) malloc(tablesize*sizeof(real));
  if (NULL==pt->table) {
    sprintf(msg,"Cannot allocate memory for function table %s.",filename);
    error(msg);
  }

  /* input loop */
  while (!feof(infile)) {

    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }

    /*  read in potential */
    if ( 1 != fscanf(infile,"%lf",&r2) ) break;
    if (npot==0) r2_start = r2;  /* catch first value */
    for (i=0; i<ncols; ++i) {
      if ( 1 != fscanf(infile,"%lf", &val)) 
        error("Line incomplete in potential file.");
      *PTR_2D(pt->table,npot,i,pt->maxsteps,ncols) = val;
      if (val!=0.0) pt->end[i] = r2; /* catch last non-zero value */
    }
    ++npot;
  }

  r2_step = (r2 - r2_start) / (npot-1);

  if (ncols==ntypes) {
    printf("Read tabulated function %s for %d atoms types.\n",
           filename,ncols);
  } else {
    printf("Read tabulated function %s for %d pairs of atoms types.\n",
           filename,ncols);
  }

  /* fill info block */
  for (i=0; i<ncols; ++i) {
    pt->begin[i]   = r2_start;
    pt->step[i]    = r2_step;
    pt->invstep[i] = 1.0 / r2_step;
    pt->end[i]    += r2_step;
  }

  /* The interpolation uses k+1 and k+2, so we make a few copies 
     of the last value at the end of the table */
  for (k=1; k<=5; ++k) {
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }
    for (i=0; i<ncols; ++i)
      *PTR_2D(pt->table,npot,i,pt->table,ncols) 
          = *PTR_2D(pt->table,npot-1,i,pt->table,ncols);
    ++npot;
  }
}

/*****************************************************************************
*
*  read potential in second format: at the beginning <ncols> times 
*  a line of the form
*
*  r_begin r_end r_step,
*  
*  then the values of the potential (one per line), first those 
*  for atom pair  00, then an empty line (for gnuplot), then 01 and so on.
*  Analogously, if there is only one column per atom type.
*
*  Note that it is assumed that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table2(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k, *len;
  int tablesize;
  real val, numstep;
  str255 msg;

  len = (int  *) malloc(ncols * sizeof(real));
  if (len==NULL) error("allocation failed in read_pot_table");

  /* read the info block of the function table */
  for(i=0; i<ncols; i++) {
    if (3 != fscanf(infile, "%lf %lf %lf",
                  &pt->begin[i], &pt->end[i], &pt->step[i])) {
      sprintf(msg, "Info line in %s corrupt.", filename);
      error(msg);
    }
    pt->invstep[i] = 1.0 / pt->step[i];
    numstep        = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
    len[i]         = (int) (numstep+0.5);  
    pt->maxsteps   = MAX(pt->maxsteps, len[i]);

    /* some security against rounding errors */
    if (fabs(len[i] - numstep) >= 0.1) {
      char msg[255];
      sprintf(msg,"numstep = %f rounded to %d in file %s.",
              numstep, len[i], filename);
      printf(msg);
    }
  }

  /* allocate the function table */
  /* allow some extra values at the end for interpolation */
  tablesize = ncols * (pt->maxsteps+3);
  pt->table = (real *) malloc(tablesize * sizeof(real));
  if (NULL==pt->table) {
    sprintf(msg,"Cannot allocate memory for function table %s.",filename);
    error(msg);
  }

  /* input loop */
  for (i=0; i<ncols; i++) {
    for (k=0; k<len[i]; k++) {
      if (1 != fscanf(infile,"%lf", &val)) {
        sprintf(msg, "wrong format in file %s.", filename);
        error(msg);
      }
      *PTR_2D(pt->table,k,i,pt->maxsteps,ncols) = val;
    }
    /* make some copies of the last value for interpolation */
    for (k=len[i]; k<len[i]+3; k++)
      *PTR_2D(pt->table,k,i,pt->maxsteps,ncols) = val;
  }

  if (ncols==ntypes) {
    printf("Read tabulated function %s for %d atoms types.\n",
           filename,ncols);
  } else {
    printf("Read tabulated function %s for %d pairs of atoms types.\n",
           filename,ncols);
  }
  printf("Maximal length of table is %d.\n",pt->maxsteps);
}

/*****************************************************************************
*
*  read potential table; choose format according to header
*
*****************************************************************************/

void read_pot_table( pot_table_t *pt, char *filename, int ncols )
{
  FILE *infile;
  char buffer[1024], msg[255];
  char *token, *res;
  int  have_header=0, have_format=0, end_header;
  int  i, size=ncols, tablesize, npot=0, format;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read the header */
  do {
    /* read one line */
    res=fgets(buffer,1024,infile);
    if (NULL == res) {
      sprintf(msg,"Unexpected end of file in %s",filename);
      error(msg);
    }

    /* see if it is a header line */
    if (buffer[0]=='#') {
      have_header = 1;
      /* stop after last header line */
      end_header = (buffer[1]=='E');
      /* see if it is the format line */
      if (buffer[1]=='F') {
        /* format complete? */
        if (2!=sscanf( (const char*)(buffer+2), "%d%d", &format, &size )) {
          sprintf(msg,"Corrupted format header line in file %s",filename);
          error(msg);
        }
        /* right number of columns? */
        if (size!=ncols) {
          sprintf(msg,"Wrong number of data columns in file %s",filename);
          error(msg);
        }
        /* recognized format? */
        if ((format!=1) && (format!=2)) {
          sprintf(msg,"Unrecognized format specified for file %s",filename);
          error(msg);
        }
        have_format=1;
      }
    } else if (have_header) { 
      /* header does not end properly */
      sprintf(msg,"Corrupted header in file %s",filename);
      error(msg);
    }
  } while (!end_header);

  /* did we have a format in the header? */
  if (!have_format) {
    sprintf(msg,"Format not specified in header of file %s",filename);
    error(msg);
  }

  /* allocate info block of function table */
  pt->maxsteps = 0;
  pt->begin    = (real *) malloc(size*sizeof(real));
  pt->end      = (real *) malloc(size*sizeof(real));
  pt->step     = (real *) malloc(size*sizeof(real));
  pt->invstep  = (real *) malloc(size*sizeof(real));
  if ((pt->begin   == NULL) || (pt->end == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL)) {
    sprintf(msg,"Cannot allocate info block for function table %s.",filename);
    error(msg);
  }

  /* catch the case where potential is identically zero */
  for (i=0; i<size; ++i) {
    pt->end[i] = 0.0;
  }

  /* read the table */
  if (format==1) read_pot_table1(pt, size, filename, infile);
  if (format==2) read_pot_table2(pt, size, filename, infile);
  fclose(infile);
}

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

real pot2(pot_table_t *pt, int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/*****************************************************************************
*
*  Evaluate potential table with cubic interpolation. 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

real pot3(pot_table_t *pt, int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;
  real fac0, fac1, fac2, fac3, dfac0, dfac1, dfac2, dfac3;
  int  k;

  /* check for distances shorter than minimal distance in table */
  /* we need one extra value at the lower end for interpolation */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a = 0;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);  
  if (k==0) return pot2(pt,col,inc,r2); /* parabolic fit if on left border */
  chi   = (r2a - k * pt->step[col]) * istep; 
  k--;  


  /* intermediate values */ 
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;   /* leftmost value*/
  p1  = *ptr; ptr += inc;   /* next left */
  p2  = *ptr; ptr += inc;   /* first right */
  p3  = *ptr;               /* rightmost value*/

  /* factors for the interpolation */
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);

/* return the potential value */ 
  return fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;
}

/*****************************************************************************
*
*  write potential table (format 3)
*
******************************************************************************/

void write_pot_table_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile, *outfile2;
  char msg[255];
  int  i, j, k, l, col, flag=0;
  real r, r_step;
  if (plotpointfile != "\0") flag=1;
  /* open file */
  outfile = fopen(filename,"w");
  if (NULL == outfile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }
  /* if needed: open file for plotpoints */
  if (flag) {
      outfile2 = fopen(plotpointfile,"w");
      if (NULL == outfile) {
	  sprintf(msg,"Could not open file %s\n",filename);
	  error(msg);
      }
  }

  /* write header */
  fprintf(outfile, "#F 3 %d\n#E\n", ncols ); 

  /* write info block */
  for (i=0; i<ncols; i++) {
    r_step = (r_end[i] - r_start[i]) / (nsteps[i] - 1);
    fprintf( outfile, "%.16e %.16e %d\n", r_start[i], r_end[i], 
	     nsteps[i] );
  }
  fprintf(outfile, "\n");

  /* write data */
  k = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      col    = i * ntypes + j; 
      r      = r_start[k];
      r_step = (r_end[k] - r_start[k]) / (nsteps[k] - 1);
      for (l=0; l<nsteps[k]-1; l++) {
        fprintf(outfile, "%.16e\n", POTVAL(pt, col, ntypes*ntypes, r*r) );
	if (flag) 
	  fprintf(outfile2, "%.6e %.6e\n",r,POTVAL(pt, col, ntypes*ntypes, r*r) );
        r += r_step;
      }
      k++;
      fprintf(outfile, "%.16e\n\n",0.0);
      if (flag) fprintf(outfile2, "%.6e %.6e\n\n\n",r,0.0);
    }
  
  fclose(outfile);
  if (flag) fclose(outfile2);
}

/*****************************************************************************
*
*  write eam potential table (format 3)
*
******************************************************************************/

void write_pot_table_eam(pot_table_t *pt, pot_table_t *embed_pt, pot_table_t *rho_tab, char *filename)
{
  FILE *outfile, *outfile2;
  char msg[255];
  int  i, j, k, l, col, flag=0;
  real r, r_step;
  if (plotpointfile != "\0") flag=1;
  /* open file */
  outfile = fopen(filename,"w");
  if (NULL == outfile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }
  /* if needed: open file for plotpoints */
  if (flag) {
      outfile2 = fopen(plotpointfile,"w");
      if (NULL == outfile) {
	  sprintf(msg,"Could not open file %s\n",filename);
	  error(msg);
      }
  }

  /* write header */
  fprintf(outfile, "#F 3 %d\n#E\n", ncols ); 

  /* write info block */
  for (i=0; i<ncols; i++) {
    r_step = (r_end[i] - r_start[i]) / (nsteps[i] - 1);
    fprintf( outfile, "%.16e %.16e %d\n", r_start[i], r_end[i], 
	     nsteps[i] );
  }
  fprintf(outfile, "\n");

  /* write data */
  k = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      col    = i * ntypes + j; 
      r      = r_start[k];
      r_step = (r_end[k] - r_start[k]) / (nsteps[k] - 1);
      for (l=0; l<nsteps[k]-1; l++) {
        fprintf(outfile, "%.16e\n", POTVAL(pt, col, ntypes*ntypes, r*r) );
	if (flag) 
	  fprintf(outfile2, "%.6e %.6e\n",r,POTVAL(pt, col, ntypes*ntypes, r*r) );
        r += r_step;
      }
      k++;
      fprintf(outfile, "%.16e\n\n",0.0);
      if (flag) fprintf(outfile2, "%.6e %.6e\n\n\n",r,0.0);
    }
  for (i=0; i<ntypes; i++) {
    r      =  r_start[k];
    r_step =  (r_end[k] - r_start[k]) / (nsteps[k] - 1);
    for (l=0; l<nsteps[k]-1; l++) {
      fprintf(outfile, "%.16e\n", POTVAL(rho_tab, i, ntypes*ntypes, r*r) );
      if (flag) 
	fprintf(outfile2, "%.6e %.6e\n",r,POTVAL(rho_tab, i, ntypes*ntypes, r*r) );
      r += r_step;
    }
    k++;
    fprintf(outfile, "%.16e\n\n",0.0);
    if (flag) fprintf(outfile2, "%.6e %.6e\n\n\n",r,0.0);
  }

  for (i=0; i<ntypes; i++) {
    r      =  r_start[k];
    r_step =  (r_end[k] - r_start[k]) / (nsteps[k] - 1);
    for (l=0; l<nsteps[k]; l++) {
      fprintf(outfile, "%.16e\n", POTVAL(embed_pt, i, ntypes, r) );
      if (flag) 
	fprintf(outfile2, "%.6e %.6e\n",r,POTVAL(embed_pt, i, ntypes, r) );
      r += r_step;
    }
    k++;
    fprintf(outfile, "\n");
    if (flag) fprintf(outfile2, "\n\n");
  }

  fclose(outfile);
  if (flag) fclose(outfile2);
}

/*****************************************************************************
*
*  write plot version of potential table
*
******************************************************************************/

void write_plotpot_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  char msg[255];
  int  i, j, k, l, col;
  real r, r_step;

  /* open file */
  outfile = fopen(filename,"w");
  if (NULL == outfile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* write data */
  k = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      col    = i * ntypes + j; 
      r      = r_start[k];
      r_step = (r_end[k] - r_start[k]) / (NPLOT - 1);
      for (l=0; l<NPLOT; l++) {
        fprintf( outfile, "%e %e\n", r, POTVAL(pt, col, ntypes*ntypes, r*r) );
        r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
  if (eam) {
    for (i=0; i<ntypes; i++) {
      r=r_start[k];
      r_step=(r_end[k] - r_start[k]) / (NPLOT - 1);
      for (l=0;l<NPLOT;l++) {
	fprintf(outfile, "%e %e\n", r, POTVAL(&rho_tab, i, ntypes*ntypes, r*r) );
        r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
    for (i=0; i<ntypes; i++) {
      r=r_start[k];
      r_step=(r_end[k] - r_start[k]) / (NPLOT - 1);
      for (l=0;l<NPLOT;l++) {
	fprintf(outfile, "%e %e\n", r, POTVAL(&embed_pt, i, ntypes, r) );
        r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
  }

  fclose(outfile);
}

/*****************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv)
{
  int i, j, k;

  read_parameters(argc, argv);

  /* pair interactions */
  if (mode==PAIR) {

    /* read potential table */
    if (2*ncols==(ntypes+1)*ntypes) {
      printf("Reading Pair Potential\n");
      eam=0;
    } else if (2*ncols==(ntypes+5)*ntypes){
      printf("Reading EAM Potentials\n");
      eam=1;
    } else {
      error("ntypes and ncols are incompatible");
    }


    if (eam) {
      read_pot_table(&pt,phi_r2_file,ntypes*ntypes);
      read_pot_table(&embed_pt,f_rho_file,ntypes);
      read_pot_table(&rho_tab,rho_r2_file,ntypes*ntypes);
    } else {
      read_pot_table(&pt,infilename,ntypes*ntypes);
    }
    /* always set r_end from potential table */
    if (r_end==NULL) {
      r_end = (real *) malloc( ncols * sizeof(real) );
      if (NULL==r_end) error("allocation of r_end failed");
    }
    k=0;
    for (i=0; i<ntypes; i++)
      for (j=i; j<ntypes; j++)
        r_end[k++] = sqrt(pt.end[i*ntypes+j]);
    if (eam) {
      for (i=0; i<ntypes; i++)
	r_end[k++] = sqrt(rho_tab.end[i]);
      for (i=0; i<ntypes; i++)
	r_end[k++] = embed_pt.end[i];
    }
    /* set r_start, if not read in */
    if (r_start==NULL) {
      r_start = (real *) malloc( ncols * sizeof(real) );
      if (r_start==NULL) error("allocation of r_start failed");
    }

    k=0;
    for (i=0; i<ntypes; i++)
      for (j=i; j<ntypes; j++) {
	  r_start[k++] = sqrt(pt.begin[i*ntypes+j]);
	  while(POTVAL(&pt,i*ntypes+j,ntypes*ntypes,
		       r_start[k-1]*r_start[k-1])>MAXVAL){
	      r_start[k-1]+=(r_end[k-1]-r_start[k-1])/1000.;
	  }
      }
    if (eam) {
      for (i=0; i<ntypes; i++)
	r_start[k++] = sqrt(rho_tab.begin[i]);
      for (i=0; i<ntypes; i++)
	r_start[k++] = embed_pt.begin[i];
    }
    
    /* write potential table */
    if (eam)
      write_pot_table_eam(&pt,&embed_pt,&rho_tab,outfilename);
    else
      write_pot_table_pair(&pt,outfilename);
    write_plotpot_pair(&pt,plotfilename);
  }

  return 0;
}


