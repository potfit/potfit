/****************************************************************
* 
* neightab.c (obsolete): Used to create static neighbour tables.
*      As those are generated on the fly in current versions, this
*      file is no longer used.
* 
*****************************************************************/

/****************************************************************
* $Revision: 1.1 $
* $Date: 2004/04/13 14:13:50 $
*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAXNEIGH 200
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))

typedef double real;
typedef struct { real x; real y; real z; } vektor;

typedef struct {
  int    typ;
  real   r;
  vektor dist;
} neigh_t;

typedef struct {
  int    typ;
  int    n_neigh;
  vektor pos;
  vektor force;
  neigh_t *neigh;
} atom_t;

vektor box_x,  box_y,  box_z;
vektor tbox_x, tbox_y, tbox_z;
atom_t *atoms;
real   rcut;
int    natoms;

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

/* vector product */ 
vektor vec_prod(vektor u, vektor v)
{
  vektor w;
  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;
  return w;
}

/******************************************************************************
*
*  compute box transformation matrix
*
******************************************************************************/

void make_box( void )
{
  real volume;

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  /* first unnormalized */
  tbox_x = vec_prod( box_y, box_z );
  tbox_y = vec_prod( box_z, box_x );
  tbox_z = vec_prod( box_x, box_y );

  /* volume */
  volume = SPROD( box_x, tbox_x );
  if (0==volume) error("Box edges are parallel");

  /* normalization */
  tbox_x.x /= volume;  tbox_x.y /= volume;  tbox_x.z /= volume;
  tbox_y.x /= volume;  tbox_y.y /= volume;  tbox_y.z /= volume;
  tbox_z.x /= volume;  tbox_z.y /= volume;  tbox_z.z /= volume;
}

/******************************************************************************
*
*  read the configuration
*
******************************************************************************/

atom_t *read_config(char *filename)
{
  FILE *infile;
  char buffer[1024], msg[255], *res;
  int  i, end_header;
  atom_t *atom;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read the header */
  end_header=0;
  while (end_header==0) {
    /* read one line */
    res=fgets(buffer,1024,infile);
    if (NULL == res) {
      sprintf(msg, "Unexpected end of file %s", filename);
      error(msg);
    }
    /* check if it is a header line */
    if (buffer[0]!='#') {
      sprintf(msg, "Corrupt header in file %s", filename);
      error(msg);
    }
    /* stop after last header line */
    if (buffer[1]=='E') {
      end_header = 1;
    }
    /* number of atoms */
    else if (buffer[1]=='N') {
      if (1!=sscanf( (const char*)(buffer+2), "%d", &natoms )) {
        sprintf(msg, "Corrupt header in file %s", filename);
        error(msg);
      }
    }
    /* box_x */
    else if (buffer[1]=='X') {
      if (3!=sscanf( (const char*)(buffer+2), "%lf %lf %lf", 
                     &box_x.x, &box_x.y, &box_x.z )) {
        sprintf(msg, "Corrupt header in file %s", filename);
        error(msg);
      }
    }
    /* box_y */
    else if (buffer[1]=='Y') {
      if (3!=sscanf( (const char*)(buffer+2), "%lf %lf %lf", 
                     &box_y.x, &box_y.y, &box_y.z )) {
        sprintf(msg, "Corrupt header in file %s", filename);
        error(msg);
      }
    }
    /* box_z */
    else if (buffer[1]=='Z') {
      if (3!=sscanf( (const char*)(buffer+2), "%lf %lf %lf", 
                     &box_z.x, &box_z.y, &box_z.z )) {
        sprintf(msg, "Corrupt header in file %s", filename);
        error(msg);
      }
    }
  }

  make_box();

  /* allocate atoms */
  atoms = (atom_t *) malloc( natoms * sizeof(atom_t) );
  if (atoms==NULL) error("Cannot allocate atoms list");

  /* allocate neighbor tables */
  for (i=0; i<natoms; i++) {
    atoms[i].neigh = (neigh_t *) malloc( MAXNEIGH * sizeof(neigh_t) );
    if (atoms[i].neigh==NULL) error("Cannot allocate neighbor list");
    atoms[i].n_neigh=0;
  }

  /* read atoms */
  for (i=0; i<natoms; i++) {
    atom = atoms+i;
    fscanf( infile, "%d %lf %lf %lf %lf %lf %lf\n",
            &atom->typ, &atom->pos.x, &atom->pos.y, &atom->pos.z, 
            &atom->force.x, &atom->force.y, &atom->force.z ); 
  }
  fclose(infile);
  return atoms;
}

/******************************************************************************
*
*  make the neighbor table
*
******************************************************************************/

void make_neightab(void)
{
  int i, j, k, ix, iy, iz;
  vektor d, dd;
  real r;

  for (i=0; i<natoms; i++)
    for (j=0; j<natoms; j++) {
      d.x  = atoms[j].pos.x - atoms[i].pos.x; 
      d.y  = atoms[j].pos.y - atoms[i].pos.y;
      d.z  = atoms[j].pos.z - atoms[i].pos.z;
      for (ix=-1; ix<=1; ix++)
        for (iy=-1; iy<=1; iy++)
          for (iz=-1; iz<=1; iz++) {
            if ((i==j) && (ix==0) && (iy==0) && (iz==0)) continue;
            dd.x = d.x + ix * box_x.x + iy * box_y.x + iz * box_z.x;
            dd.y = d.y + ix * box_x.y + iy * box_y.y + iz * box_z.y;
            dd.z = d.z + ix * box_x.z + iy * box_y.z + iz * box_z.z;
            r = sqrt(SPROD(dd,dd));
            if (r <= rcut) {
              if (atoms[i].n_neigh==MAXNEIGH) 
                error("Neighbor table is too small");
              dd.x /= r;
              dd.y /= r;
              dd.z /= r;
              k = atoms[i].n_neigh;
              atoms[i].neigh[k].typ  = atoms[j].typ;
              atoms[i].neigh[k].r    = r;
              atoms[i].neigh[k].dist = dd;
              atoms[i].n_neigh++;
	    }
	  }
    }
}

void write_neightab(char *filename)
{
  FILE *outfile;
  char msg[255];
  int  i, j;
  atom_t  *atom;
  neigh_t *neigh;

  /* open file */
  outfile = fopen(filename,"w");
  if (NULL == outfile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* write the data */
  fprintf(outfile, "%d\n", natoms );
  for (i=0; i<natoms; i++) {
    atom = atoms + i;
    fprintf( outfile, "%d %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
             atom->typ, atom->pos.x, atom->pos.y, atom->pos.z, 
             atom->force.x, atom->force.y, atom->force.z,
             atom->n_neigh );
    for (j=0; j<atom->n_neigh; j++) {
      neigh = atom->neigh + j;
      fprintf( outfile, "  %d %.16e %.16e %.16e %.16e\n", 
               neigh->typ, neigh->r, 
               neigh->dist.x, neigh->dist.y, neigh->dist.z ); 
    }
  }
  fclose(outfile);
}

int main(int argc, char **argv)
{
  char msg[255];

  if (argc<3) {
    printf("Usage: %s <rcut> <infile> <outfile>\n", argv[0] );
    return 1;
  }
  rcut  = atof(argv[1]);
  atoms = read_config( argv[2] );
  make_neightab();
  write_neightab( argv[3] );
  return 0; 
}
