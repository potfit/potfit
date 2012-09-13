
/******************************************************************************
*
*  IMD -- The ITAP Molecular Dynamics Program
*
*  Copyright 1996-2001 Institute for Theoretical and Applied Physics,
*  University of Stuttgart, D-70550 Stuttgart
*
*  $Revision: 1.3 $
*  $Date: 2007/10/04 14:55:35 $
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
typedef struct {
  real  x, y, z;
} vector;

static real swaptemp;
#define SWAP(a,b) swaptemp=(a) ; (a)=(b) ; (b)=swaptemp

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr, "Error: %s\n", msg);
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
  char *res;
  char  buffer[1024];
#ifdef STRESS
  vector tmp2;
#endif
  int  *typ, num[MAXTYP], max = 0, n, i, j, tag_format, h_boxx, h_boxy, h_boxz;
  real  triple_prod;

  /* open input file */
  infile = fopen(argv[1], "r");
  if (infile == NULL)
    error("cannot open input file");

  /* old format or new format? */
  res = fgets(buffer, 1024, infile);
  if (NULL == res)
    error("Unexpected end of file\n");
  if (res[0] == '#') {		/* new file type */
    tag_format = 1;
    if (res[1] == 'N') {	/* Atom number line */
      if (sscanf(res + 3, "%d %d", &n, &i) < 2)
	error("Error in atom number specification\n");
    } else
      error("Error - number of atoms missing\n");
  } else {
    /* number of atoms in this configuration */
    tag_format = 0;
    if (1 > sscanf(buffer, "%d", &n))
      error("Unexpected end of file");
  }

  /* read number of atoms read, now allocate memory */
  pos = (vector *)malloc(n * sizeof(vector));
  typ = (int *)malloc(n * sizeof(int));
  if ((NULL == pos) || (NULL == typ))
    error("cannot allocate memory");

  if (tag_format) {
    do {
      res = fgets(buffer, 1024, infile);
      /* read the box vectors */
      if (res[1] == 'X') {
	if (sscanf(res + 3, "%lf %lf %lf\n", &box_x.x, &box_x.y, &box_x.z) == 3)
	  h_boxx++;
	else
	  error("Error in box_x vector\n");
      } else if (res[1] == 'Y') {
	if (sscanf(res + 3, "%lf %lf %lf\n", &box_y.x, &box_y.y, &box_y.z) == 3)
	  h_boxy++;
	else
	  error("Error in box_y vector\n");
      } else if (res[1] == 'Z') {
	if (sscanf(res + 3, "%lf %lf %lf\n", &box_z.x, &box_z.y, &box_z.z) == 3)
	  h_boxz++;
	else
	  error("Error in box_z vector\n");
      }
    } while (res[1] != 'F');
    /* read energy and stress -> discard */
    if (!(h_boxx && h_boxy && h_boxz))
      error("Incomplete force file!");
  } else {
    /* read box vectors */
    fscanf(infile, "%lf %lf %lf\n", &box_x.x, &box_x.y, &box_x.z);
    fscanf(infile, "%lf %lf %lf\n", &box_y.x, &box_y.y, &box_y.z);
    fscanf(infile, "%lf %lf %lf\n", &box_z.x, &box_z.y, &box_z.z);

    /* read energy and stress -> discard */
    fscanf(infile, "%lf\n", &tmp.x);
    //printf("%f\n",tmp.x);
#ifdef STRESS
    fscanf(infile, "%lf %lf %lf %lf %lf %lf\n", &tmp.x, &tmp.y, &tmp.z, &tmp2.x,
      &tmp2.y, &tmp2.z);
//  printf("%f %f\n",tmp.x,tmp2.z);
#endif
  }
  /* read atoms */
  for (i = 0; i < MAXTYP; i++)
    num[i] = 0;
  for (i = 0; i < n; i++) {
    fscanf(infile, "%d %lf %lf %lf %lf %lf %lf\n", typ + i, &(pos + i)->x,
      &(pos + i)->y, &(pos + i)->z, &tmp.x, &tmp.y, &tmp.z);
    if (typ[i] < MAXTYP) {
      num[typ[i]]++;
      max = MAX(max, typ[i] + 1);
//      printf("Atom %d, Type %d\n",i,typ[i]);
    } else
      error("not enough types - increase MAXTYP");
  }

  /* close input file */
  fclose(infile);

  /* if triple product of box vectors is negative: swap vectors 1 and 2 */
  triple_prod = ((box_y.y * box_z.z) - (box_y.z * box_z.y)) * box_x.x;
  triple_prod += ((box_y.z * box_z.x) - (box_y.x * box_z.z)) * box_x.y;
  triple_prod += ((box_y.x * box_z.y) - (box_y.y * box_z.x)) * box_x.z;
  if (triple_prod < 0) {
    SWAP(box_x.x, box_y.x);
    SWAP(box_x.y, box_y.y);
    SWAP(box_x.z, box_y.z);
  }

  /* open POSCAR file */
  outfile = fopen("POSCAR", "w");
  if (outfile == NULL)
    error("cannot open POSCAR file");

  /* write header */
  fprintf(outfile, "Structure from %s\n", argv[1]);
  fprintf(outfile, " %f\n", (real)1.0);
  fprintf(outfile, "%f %f %f\n", box_x.x, box_x.y, box_x.z);
  fprintf(outfile, "%f %f %f\n", box_y.x, box_y.y, box_y.z);
  fprintf(outfile, "%f %f %f\n", box_z.x, box_z.y, box_z.z);
  for (i = 0; i < max; i++)
    fprintf(outfile, " %d", num[i]);
  fprintf(outfile, "\nCartesian\n");

  /* write atoms */
  for (i = 0; i < max; i++)
    for (j = 0; j < n; j++)
      if (typ[j] == i)
	fprintf(outfile, "%f %f %f\n", pos[j].x, pos[j].y, pos[j].z);

  /* close POSCAR file */
  fclose(outfile);

  return 0;
}
