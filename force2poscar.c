
/******************************************************************************
*
*  IMD -- The ITAP Molecular Dynamics Program
*
*  Copyright 1996-2001 Institute for Theoretical and Applied Physics,
*  University of Stuttgart, D-70550 Stuttgart
*
*  $Revision: 1.1 $
*  $Date: 2003/01/07 19:30:07 $
*
*  Convert an IMD force file to a VASP POSCAR file
*
*  Usage: force2poscar <infile>
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXTYP 10
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef double real;
typedef struct { real x, y, z; } vector;

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
*  main
*
******************************************************************************/

int main(int argc, char **argv)
{
  FILE *infile, *outfile;
  vector *pos, box_x, box_y, box_z, tmp;
  int *typ, num[MAXTYP], max=0, n, i, j;

  /* open input file */
  infile = fopen(argv[1], "r" );
  if (infile==NULL) error("cannot open input file");

  /* read number of atoms, allocate memory */
  fscanf(infile, "%d\n", &n);
  pos = (vector *) malloc( n * sizeof(vector) );
  typ = (int    *) malloc( n * sizeof(int   ) );
  if ((NULL==pos) || (NULL==typ)) error("cannot allocate memory");

  /* read box vectors */
  fscanf(infile, "%lf %lf %lf\n", &box_x.x, &box_x.y, &box_x.z ); 
  fscanf(infile, "%lf %lf %lf\n", &box_y.x, &box_y.y, &box_y.z ); 
  fscanf(infile, "%lf %lf %lf\n", &box_z.x, &box_z.y, &box_z.z ); 

  /* read atoms */
  for (i=0; i<MAXTYP; i++) num[i]=0;
  for (i=0; i<n; i++) {
    fscanf(infile, "%d %lf %lf %lf %lf %lf %lf\n",
      typ+i, &(pos+i)->x, &(pos+i)->y, &(pos+i)->z, &tmp.x, &tmp.y, &tmp.z );
    if (typ[i]<MAXTYP) {
      num[typ[i]]++;
      max = MAX( max, typ[i]+1 );
    } else error("not enough types - increase MAXTYP");
  }

  /* close input file */
  fclose(infile);

  /* open POSCAR file */
  outfile = fopen("POSCAR", "w");
  if (outfile==NULL) error("cannot open POSCAR file");

  /* write header */
  fprintf(outfile, "Al-Co-Ni structure from %s\n", argv[1] );
  fprintf(outfile, " %f\n", (real) 1.0 );
  fprintf(outfile, "%f %f %f\n", box_x.x, box_x.y, box_x.z );
  fprintf(outfile, "%f %f %f\n", box_y.x, box_y.y, box_y.z );
  fprintf(outfile, "%f %f %f\n", box_z.x, box_z.y, box_z.z );
  for (i=0; i<max; i++) fprintf(outfile, " %d", num[i]);
  fprintf(outfile, "\nCartesian\n");

  /* write atoms */
  for (i=0; i<max; i++)
    for (j=0; j<n; j++)
      if (typ[j]==i) 
        fprintf(outfile, "%f %f %f\n", pos[j].x, pos[j].y, pos[j].z ); 

  /* close POSCAR file */
  fclose(outfile);

  return 0;
}
