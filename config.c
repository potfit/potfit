
#include "potfit.h"

/*****************************************************************************
*
*  read the configurations
*
******************************************************************************/

void read_config(char *filename)
{
  FILE    *infile;
  char    msg[255];
  int     nconf=0, count, i, j, k;
  atom_t  *atom;
  neigh_t *neigh;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read configurations until the end of the file */
  do {

    /* number of atoms in this configuration */
    if (1>fscanf(infile, "%d", &count)) error("Unexpected end of file");

    /* increase memory for this many additional atoms */
    atoms = (atom_t *) realloc(atoms, (natoms+count) * sizeof(atom_t));
    if (NULL==atoms)   error("Cannot allocate memory for atoms");
    force_0 = (real *) realloc(force_0, 3 * (natoms+count) * sizeof(real));
    if (NULL==force_0) error("Cannot allocate memory for forces");

    /* read the atoms */
    for (i=0; i<count; i++) {

      k    = 3 * (natoms + i);
      atom = atoms + natoms + i;
      if (8>fscanf(infile,"%d %lf %lf %lf %lf %lf %lf %d\n", &(atom->typ), 
                   &(atom->pos.x), &(atom->pos.y), &(atom->pos.z), 
                   force_0+k, force_0+k+1, force_0+k+2, &(atom->n_neigh)))
        error("Corrupt configuration file");

      /* allocate memory for neighbors */
      atom->neigh = (neigh_t *) malloc(atom->n_neigh * sizeof(neigh_t));
      if (NULL==atom->neigh) error("Cannot allocate memory for neighbors");

      /* read the neighbors */
      for (j=0; j<atom->n_neigh; j++) {
	neigh = atom->neigh + j;
        if (5>fscanf(infile,"%d %lf %lf %lf %lf\n", 
                     &(neigh->typ), &(neigh->r ), 
                     &(neigh->dist.x), &(neigh->dist.y), &(neigh->dist.z)))
          error("Corrupt configuration file");
      }
    }

    /* increment natoms and configuration number */
    natoms += count;
    nconf++;

  } while (!feof(infile));

  /* print diagnostic message and close file */
  printf("Read %d configurations with a total of %d atoms\n", nconf, natoms);
  fclose(infile);
}


