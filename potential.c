#define NPLOT 1000
#include "potfit.h"

/******************************************************************************
*
* read potential table
*
******************************************************************************/

void read_pot_table( pot_table_t *pt, char *filename, int ncols )
{
  FILE *infile;
  char buffer[1024], msg[255], *res;
  int  have_format=0, end_header=0;
  int  format, size, i, j, k, l, *nvals;
  real *val;

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
    /* check if it is a header line */
    if (buffer[0]!='#') {
      sprintf(msg,"Header corrupt in file %s",filename);
      error(msg);
    }
    /* stop after last header line */
    if (buffer[1]=='E') {
      end_header = 1;
    }
    /* see if it is the format line */
    else if (buffer[1]=='F') {
      /* format complete? */
      if (2!=sscanf( (const char*)(buffer+2), "%d %d", &format, &size )) {
        sprintf(msg,"Corrupt format header line in file %s",filename);
        error(msg);
      }
      /* right number of columns? */
      if (size!=ncols) {
        sprintf(msg,"Wrong number of data columns in file %s",filename);
        error(msg);
      }
      /* recognized format? */
      if (format!=3) {
        sprintf(msg,"Unrecognized format specified for file %s",filename);
        error(msg);
      }
      have_format = 1;
    }
  } while (!end_header);

  /* did we have a format in the header? */
  if (!have_format) {
    sprintf(msg,"Format not specified in header of file %s",filename);
    error(msg);
  }

  /* allocate info block of function table */
  pt->len     = 0;
  pt->ncols   = size;
  pt->begin   = (real *) malloc(size*sizeof(real));
  pt->end     = (real *) malloc(size*sizeof(real));
  pt->step    = (real *) malloc(size*sizeof(real));
  pt->invstep = (real *) malloc(size*sizeof(real));
  pt->first   = (int  *) malloc(size*sizeof(int ));
  pt->last    = (int  *) malloc(size*sizeof(int ));
  nvals       = (int  *) malloc(size*sizeof(int ));
  if ((pt->begin   == NULL) || (pt->end   == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL) || 
      (nvals       == NULL)) {
    sprintf(msg,"Cannot allocate info block for potential table %s",filename);
    error(msg);
  }

  /* read the info block of the function table */
  for(i=0; i<ncols; i++) {
    if (3>fscanf(infile,"%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i])) {
        sprintf(msg, "Premature end of potential file %s", filename);
        error(msg);
    }
    pt->step[i] = (pt->end[i] - pt->begin[i]) / (nvals[i]-1);
    pt->invstep[i] = 1.0 / pt->step[i];
    if (i==0) pt->first[i] = 0; else pt->first[i] = pt->last[i-1] + 1;
    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }

  /* allocate the function table */
  pt->table = (real *) malloc(pt->len * sizeof(real));
  pt->d2tab = (real *) malloc(pt->len * sizeof(real));
  pt->idx   = (int  *) malloc(pt->len * sizeof(int ));
  if ((NULL==pt->table) || (NULL==pt->idx) || (NULL==pt->d2tab)) {
    error("Cannot allocate memory for potential table");
  }

  /* input loop */
  val = pt->table;
  k=0; l=0;
  for (i=0; i<ncols; i++) {
    for (j=0; j<nvals[i]; j++) {
      if (1>fscanf(infile, "%lf\n", val)) {
        sprintf(msg, "Premature end of potential file %s", filename);
        error(msg);
      } else val++;
      if (j<nvals[i]-1) pt->idx[k++] = l++;
      else l++;
    }
  }
  pt->idxlen = k;

  fclose(infile);

  /* compute rcut */
  rcut = (real *) malloc( ntypes * ntypes * sizeof(real) );
  if (NULL==rcut) error("Cannot allocate rcut");
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1))/2) 
	           : j * ntypes + i - ((j * (j + 1))/2);
      rcut[i * ntypes + j] = pt->end[k];
    }

}


/*****************************************************************************
*
*  Evaluate derivative of potential with quadratic interpolation. 
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real grad2(pot_table_t *pt, real *xi, int col, real r)
{
  real rr, istep, chi, p0, p1, p2, dv, d2v;
  int  k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  chi   = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];

  /* intermediate values */
  p0  = (k<=pt->last[col]) ? xi[k++] : 0.0;
  p1  = (k<=pt->last[col]) ? xi[k++] : 0.0;
  p2  = (k<=pt->last[col]) ? xi[k++] : 0.0;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the derivative */
  return istep * (dv + (chi - 0.5) * d2v);
}

/*****************************************************************************
*
*  Evaluate derivative of potential with cubic interpolation. 
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real grad3(pot_table_t *pt, real *xi, int col, real r)
{
  real rr, istep, chi, p0, p1, p2, p3;
  real dfac0, dfac1, dfac2, dfac3;
  int  k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  if (k==0) return grad2(pt,xi,col,r); /* parabolic fit if on left border */
  chi   = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];
  k--;  

  /* intermediate values */
  if (k<=pt->last[col]) p0 = xi[k++]; else return 0.0;
  if (k<=pt->last[col]) p1 = xi[k++]; else return 0.0;
  if (k<=pt->last[col]) p2 = xi[k++]; else return 0.0;
  if (k<=pt->last[col]) p3 = xi[k  ]; else { /* p2 = 0.0;  dp2 = 0.0; */
      dfac0 = -0.25 * (3.0 * chi - 1.0) * (chi - 1.0);
      dfac1 =         (3.0 * chi + 1.0) * (chi - 1.0);
      /* dfac2 = -0.25 * (9.0 * chi + 5.0) * (chi - 1.0);*/
      /* dfac3 =  0.5  * (3.0 * (chi*chi - 1));*/
      return istep * (dfac0 * p0 + dfac1 * p1);
  }

  /* factors for the interpolation of the 1. derivative */
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);

  /* return the derivative */
  return istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);
}


/*****************************************************************************
*
*  Evaluate potential with quadratic interpolation. 
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real pot2(pot_table_t *pt, int col, real r)
{
  real rr, istep, chi, p0, p1, p2, dv, d2v;
  int  k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  chi   = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];

  /* intermediate values */
  p0  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  p1  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  p2  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}


/*****************************************************************************
*
*  Evaluate potential with cubic interpolation. 
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real pot3(pot_table_t *pt, int col, real r)
{
  real rr, istep, chi, p0, p1, p2, p3;
  real fac0, fac1, fac2, fac3;
  int  k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  if (k==0) return pot2(pt,col,r); /* parabolic fit if on left border */
  chi   = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];
  k--;

  /* intermediate values */
  if (k<=pt->last[col]) p0 = pt->table[k++]; else return 0.0;
  if (k<=pt->last[col]) p1 = pt->table[k++]; else return 0.0;
  if (k<=pt->last[col]) p2 = pt->table[k++]; else return 0.0;
  if (k<=pt->last[col]) p3 = pt->table[k  ]; else {/* p2 = 0.0; dp2 = 0.0 */
      fac0 = -0.25 * chi * SQR(chi-1.0);
      fac1 =         (chi*chi - 1) * (chi - 1);
      /* fac2 = -0.25 * chi * (chi + 1) * (3.0*chi - 5.0); */
      /* fac3 = -0.5  * (chi*chi - 1) * chi;           */ 
      return fac0 * p0 + fac1 * p1;
  }  /* go smoothly: interpolate with f=f'=0 at chi=1. */ 

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

void write_pot_table(pot_table_t *pt, char *filename)
{
  FILE *outfile, *outfile2;
  char msg[255];
  int  i, j, flag=0;
  real r;

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
  fprintf(outfile, "#F 3 %d\n#E\n", pt->ncols ); 

  /* write info block */
  for (i=0; i<pt->ncols; i++) {
    fprintf(outfile, "%.16e %.16e %d\n", 
            pt->begin[i], pt->end[i], pt->last[i] - pt->first[i] + 1);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i=0; i<pt->ncols; i++) {
    r = pt->begin[i];
    for (j=pt->first[i]; j<=pt->last[i]; j++) {
      fprintf(outfile, "%.16e\n", pt->table[j] );
      if (flag) 
	fprintf(outfile2, "%.6e %.6e\n",r,pt->table[j] );
      r += pt->step[i];
    }
    fprintf(outfile, "\n");
    if (flag) fprintf(outfile2,"\n\n");
  }
  fclose(outfile);
  if (flag) fclose(outfile2);
}


/*****************************************************************************
*
*  write potential table for IMD (format 2)
*
******************************************************************************/

void write_pot_table_imd(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  char msg[255];
  real *r2begin, *r2end, *r2step, r2;
  int  i, j, k, m, m2, col1, col2;

  /* allocate memory */
  r2begin = (real *) malloc( ntypes * ntypes *sizeof(real) );
  r2end   = (real *) malloc( ntypes * ntypes *sizeof(real) );
  r2step  = (real *) malloc( ntypes * ntypes *sizeof(real) );
  if ((r2begin==NULL) || (r2end==NULL) || (r2step==NULL)) 
    error("Cannot allocate memory in  write_pot_table_imd");

  /* open file */
  outfile = fopen(filename,"w");
  if (NULL == outfile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes ); 

  /* write info block */
  m=0;
  for (i=0; i<ntypes; i++){ 
      m+=i;
      m2=0;
      for (j=0; j<ntypes; j++) {
	  m2+=j;
	  col1 = i<j ? i * ntypes + j - m : j * ntypes + i - m2;
	  col2 = i * ntypes + j;
	  r2begin[col2] = SQR(pt->begin[col1]);
	  r2end  [col2] = SQR(pt->end[col1]);
	  r2step [col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
	  fprintf(outfile, "%.16e %.16e %.16e\n", 
		  r2begin[col2], r2end[col2], r2step[col2]);
      }
    }
  fprintf(outfile, "\n");

  /* write data */
  m=0;
  for (i=0; i<ntypes; i++) {
      m+=i;
      m2=0;
      for (j=0; j<ntypes; j++) {
	  m2+=j;
	  col1 = i<j ? i * ntypes + j - m : j * ntypes + i - m2;
	  col2 = i * ntypes + j;
	  r2 = r2begin[col2];
	  for (k=0; k<=imdpotsteps; k++) { 
	      fprintf(outfile, "%.16e\n", splint_ed(pt, col1, sqrt(r2) ));
	      r2 += r2step[col2];
	  }
	  fprintf(outfile, "\n");
      }
    }
  fclose(outfile);
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
      r      = pt->begin[k];
      r_step = (pt->end[k] - pt->begin[k]) / (NPLOT-1);
      for (l=0; l<NPLOT; l++) {
        fprintf( outfile, "%e %e\n", r, splint_ed(pt, k, r) );
        r += r_step;
      }
      fprintf( outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
  fclose(outfile);
}
