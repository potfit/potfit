/****************************************************************
* 
*  config.c: Reads atomic configurations and forces.
* 
*****************************************************************/
/****************************************************************
* $Revision: 1.16 $
* $Date: 2003/11/20 08:38:44 $
*****************************************************************/

#include "potfit.h"

/*****************************************************************************
*
*  read the configurations
*
******************************************************************************/

void read_config_old(char *filename)
{
  FILE    *infile;
  char    msg[255];
  int     count, i, j, k, maxneigh=0;
  atom_t  *atom;
  neigh_t *neigh;
  
  nconf=0;

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

      /* check if n_neigh < MAXNEIGH */
      maxneigh = MAX( maxneigh, atom->n_neigh );
      if (maxneigh >= MAXNEIGH) {
        sprintf(msg, "%d neighbors, MAXNEIGH is %d!", maxneigh, MAXNEIGH );
        error(msg);
      }

      /* read the neighbors */
      for (j=0; j<atom->n_neigh; j++) {
	neigh = atom->neigh + j;
        if (5>fscanf(infile,"%d %lf %lf %lf %lf\n", 
                     &(neigh->typ), &(neigh->r ), 
                     &(neigh->dist.x), &(neigh->dist.y), &(neigh->dist.z)))
          error("Corrupt configuration file");
      }
    }
    mdim=k+3;			/* mdim is dimension of force vector */

    /* increment natoms and configuration number */
    natoms += count;
    nconf++;

  } while (!feof(infile));

  /* print diagnostic message and close file */
  printf("Maximal number of neighbors is %d, MAXNEIGH is %d\n",
         maxneigh, MAXNEIGH );
  printf("Read %d configurations with a total of %d atoms\n", nconf, natoms);
  fclose(infile);
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

/*****************************************************************************
*
*  read the configurations
*
******************************************************************************/

void read_config(char *filename)
{
  int     maxneigh=0, count;
  int     i, j, k, ix, iy, iz;
  FILE    *infile;
  char    msg[255];
  atom_t  *atom;
  neigh_t *neigh;
  vektor  d, dd;
  real    r;
  
  nconf = 0;

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
    coheng = (real *) realloc(coheng, (nconf+1) * sizeof (real));
    if (NULL==coheng)  error("Cannot allocate memory for cohesive energy");
    inconf = (int *) realloc(inconf, (nconf+1) * sizeof (int));
    if (NULL==inconf)  error("Cannot allocate memory for atoms in conf");
    cnfstart = (int *) realloc(cnfstart, (nconf+1) * sizeof (int));
    if (NULL==cnfstart)  error("Cannot allocate memory for start of conf");
 
    inconf[nconf]=count;
    cnfstart[nconf]=natoms;

    /* read the box vectors */
    fscanf( infile, "%lf %lf %lf\n", &box_x.x, &box_x.y, &box_x.z );
    fscanf( infile, "%lf %lf %lf\n", &box_y.x, &box_y.y, &box_y.z );
    fscanf( infile, "%lf %lf %lf\n", &box_z.x, &box_z.y, &box_z.z );
    make_box();

    /* read cohesive energy */
    if (1!=fscanf(infile, "%lf\n", &(coheng[nconf])))
	error("Configuration file without cohesive energy -- old format!");

    /* read the atoms */
    for (i=0; i<count; i++) {
      k    = 3 * (natoms + i);
      atom = atoms + natoms + i;
      if (7>fscanf( infile, "%d %lf %lf %lf %lf %lf %lf\n", &(atom->typ), 
                    &(atom->pos.x), &(atom->pos.y), &(atom->pos.z), 
                    &(atom->force.x), &(atom->force.y), &(atom->force.z)))
        error("Corrupt configuration file");
      atom->conf = nconf;
    }

    /* compute the neighbor table */
    for (i=natoms; i<natoms+count; i++) {
      atoms[i].n_neigh = 0;
      for (j=natoms; j<natoms+count; j++) {
        d.x = atoms[j].pos.x - atoms[i].pos.x; 
        d.y = atoms[j].pos.y - atoms[i].pos.y;
        d.z = atoms[j].pos.z - atoms[i].pos.z;
        for (ix=-1; ix<=1; ix++)
          for (iy=-1; iy<=1; iy++)
            for (iz=-1; iz<=1; iz++) {
              if ((i==j) && (ix==0) && (iy==0) && (iz==0)) continue;
              dd.x = d.x + ix * box_x.x + iy * box_y.x + iz * box_z.x;
              dd.y = d.y + ix * box_x.y + iy * box_y.y + iz * box_z.y;
              dd.z = d.z + ix * box_x.z + iy * box_y.z + iz * box_z.z;
              r = sqrt(SPROD(dd,dd));
              if (r <= rcut[ atoms[i].typ * ntypes + atoms[j].typ ]) {
		if (r <= rmin[atoms[i].typ * ntypes + atoms[j].typ ]){
		  sprintf(msg,"Distance too short between atom %d and %d in conf %d",
			  i,j,nconf);
		  error(msg);
		}
                if (atoms[i].n_neigh==MAXNEIGH) 
                  error("Neighbor table is too small");
                dd.x /= r;
                dd.y /= r;
                dd.z /= r;
                k = atoms[i].n_neigh;
                atoms[i].neigh[k].typ  = atoms[j].typ;
		atoms[i].neigh[k].nr   = j;
                atoms[i].neigh[k].r    = r;
                atoms[i].neigh[k].dist = dd;
                atoms[i].n_neigh++;
	      }
	    }
      }
      maxneigh = MAX( maxneigh, atoms[i].n_neigh );
    }

    /* increment natoms and configuration number */
    natoms += count;
    nconf++;

  } while (!feof(infile));
  fclose(infile);

  mdim=3*natoms+nconf;       /* mdim is dimension of force vector 
				3*natoms are real forces, 
				nconf cohesive energies, */ 
#ifdef EAM  
  if (eam) mdim+=1+2*ntypes;          /* 1+2*ntypes dummy constraints */
#ifdef LIMIT
  if (eam) mdim+=nconf;		/* nconf limiting constraints */
#endif LIMIT
#endif EAM
  /* copy forces into single vector */
  if (NULL==(force_0 = (real *) malloc( mdim * sizeof(real) ) ) )
    error("Cannot allocate forces");
  k = 0;
  for (i=0; i<natoms; i++) {
    force_0[k++] = atoms[i].force.x;
    force_0[k++] = atoms[i].force.y;
    force_0[k++] = atoms[i].force.z;
  }
  for (i=0; i<nconf; i++) {
    force_0[k++] = coheng[i] * (real) ENG_WEIGHT; }
#ifdef EAM
  if (eam) {
#ifdef LIMIT 
   for(i=0; i<nconf; i++) force_0[k++]=0.; /* punishment rho out of bounds */
#endif
    force_0[k++]=DUMMY_WEIGHT * dummy_rho;		/* dummy constraints */
#ifdef LIMIT
    force_0[k-1]=0.; 		/* ignore dummy_rho if LIMIT is used */
#endif
    for (i=0; i<ntypes; i++) { 	/* constraints on phi */
      force_0[k++]=DUMMY_WEIGHT * dummy_phi[i];}
    for (i=0; i<ntypes;i++) {  /* constraint on U(n=0):=0 */
      force_0[k++]=0.;}
  }

#endif
  /* print diagnostic message and close file */
  printf("Maximal number of neighbors is %d, MAXNEIGH is %d\n",
         maxneigh, MAXNEIGH );
  printf("Read %d configurations with a total of %d atoms\n", nconf, natoms);
  return;
}
