/****************************************************************
* 
*  config.c: Reads atomic configurations and forces.
* 
*****************************************************************/
/****************************************************************
* $Revision: 1.28 $
* $Date: 2004/12/03 17:33:25 $
*****************************************************************/

#include "potfit.h"
#include "nrutil_r.h"

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

real make_box( void )
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
  return volume;
}

/*****************************************************************************
*
*  read the configurations
*
******************************************************************************/

void read_config(char *filename)
{
  int     maxneigh=0, count;
  int     i, j, k, ix, iy, iz,typ1,typ2,col,slot,klo,khi;
  FILE    *infile;
  char    msg[255];
  atom_t  *atom;
  neigh_t *neigh;
  stens   *stresses;
  vektor  d, dd;
  real    r,rr,istep,shift,step;
#ifdef DEBUG
  real    *mindist;
  mindist = (real *) malloc(ntypes*ntypes*sizeof(real));
  for (i=0;i<ntypes*ntypes;i++) mindist[i]=rcut[i];
#endif
  
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
    volumen = (real *) realloc(volumen, (nconf+1) * sizeof (real));
    if (NULL==volumen)  error("Cannot allocate memory for volume");
    stress = (stens *) realloc(stress, (nconf+1) * sizeof (stens));
    if (NULL==stress)  error("Cannot allocate memory for stress");
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
    volumen[nconf]=make_box();

    /* read cohesive energy */
    if (1!=fscanf(infile, "%lf\n", &(coheng[nconf])))
	error("Configuration file without cohesive energy -- old format!");

    /* read stress tensor */
    stresses=stress+nconf;
    if (6!=fscanf(infile, "%lf %lf %lf %lf %lf %lf\n", &(stresses->xx),
		  &(stresses->yy), &(stresses->zz), &(stresses->xy), 
		  &(stresses->yz), &(stresses->zx))) 
      error("No stresses given -- old format");


    /* read the atoms */
    for (i=0; i<count; i++) {
      k    = 3 * (natoms + i);
      atom = atoms + natoms + i;
      if (7>fscanf( infile, "%d %lf %lf %lf %lf %lf %lf\n", &(atom->typ), 
                    &(atom->pos.x), &(atom->pos.y), &(atom->pos.z), 
                    &(atom->force.x), &(atom->force.y), &(atom->force.z)))
        error("Corrupt configuration file");
      atom->absforce = sqrt(SQR(atom->force.x)+
			    SQR(atom->force.y)+
			    SQR(atom->force.z));
      /* ++++++++++++++ */
//      printf("Atom %d, x %f, y %f, z %f, abs %f\n", natoms+i, atom->force.x, atom->force.y, atom->force.z, atom->absforce);
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
	      typ1=atoms[i].typ;
	      typ2=atoms[j].typ;
              if (r <= rcut[ typ1 * ntypes + typ2 ]) {
		if (r <= rmin[typ1 * ntypes + typ2 ]){
		  printf("%d: %f %f %f\n", i-natoms, atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
		  printf("%d: %f %f %f\n", j-natoms, dd.x, dd.y, dd.z);
		  sprintf(msg,
 "Distance %f too short between atom %d (type %d) and %d (type %d) in conf %d",
			  r,i-natoms, typ1, j-natoms, typ2,nconf);
		  error(msg);
		}
                if (atoms[i].n_neigh==MAXNEIGH) 
                  error("Neighbor table is too small");
                dd.x /= r;
                dd.y /= r;
                dd.z /= r;
                k = atoms[i].n_neigh;
                atoms[i].neigh[k].typ  = typ2;
		atoms[i].neigh[k].nr   = j;
                atoms[i].neigh[k].r    = r;
                atoms[i].neigh[k].dist = dd;
                atoms[i].n_neigh++;
#ifdef DEBUG
		/* Minimal distance check */
		mindist[ntypes*typ1 +typ2]=
		  MIN(mindist[ntypes*typ1 +typ2],r);
#endif
		
		/* pre-compute index and shift into potential table */
		/* pair potential */
		col = (typ1 <= typ2) ? 
		   typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
		   : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);
		if (format==3) {
		  rr    = r - pair_pot.begin[col];
		  if (rr < 0) {
		    printf("%f %f %d %d %d\n",r,pair_pot.begin[col],col,nconf,i-natoms);
//		  printf("%f %f %f %f %f %f\n", d.x,d.y,d.z,coheng[nconf],stresses->xx,stresses->yz);
			 
		    fflush(stdout);
		    error("short distance in config.c!");
		  }		
		  istep = pair_pot.invstep[col];
		  slot  = (int)( rr * istep);
		  shift = (rr - slot * pair_pot.step[col]) * istep;
		  slot  += pair_pot.first[col];		
		  step  = pair_pot.step[col];

		} else {	/* format == 4 ! */

		  klo=pair_pot.first[col];
		  khi=pair_pot.last[col];
		  /* bisection */
		  while (khi-klo > 1) {
		    slot=(khi+klo) >> 1;
		    if ( pair_pot.xcoord[slot] > r ) khi=slot;
		    else klo=slot;
		  }
		  slot=klo;
		/* Check if we are at the last index - we should be lower */
                /* should be impossible anyway */
		/*  if (slot>=pair_pot.last[col]) {
		    klo--;khi--;
		    } */
		  step=pair_pot.xcoord[khi]-pair_pot.xcoord[klo];
		  shift=(r-pair_pot.xcoord[klo])/step;

		}
		/* independent of format - we should be left of last index */
		if (slot>=pair_pot.last[col]) {
		  slot--;shift+=1.0;
		}
		atoms[i].neigh[k].shift[0]  = shift;
		atoms[i].neigh[k].slot[0]   = slot;
		atoms[i].neigh[k].step[0]   = step;
#ifdef EAM
		/* EAM part */
		col=paircol+typ2;
		if (format==3) {
		  rr    = r - pair_pot.begin[col];
		  if (rr < 0) {
		    printf("%f %f %d %d %d\n",r,pair_pot.begin[col],col,typ1,typ2);
		    fflush(stdout);
		    error("short distance in config.c!");
		  }		
		  istep = pair_pot.invstep[col];
		  slot  = (int)( rr * istep);
		  shift = (rr - slot * pair_pot.step[col]) * istep;
		  slot  += pair_pot.first[col];
		  step  = pair_pot.step[col];
		} else {	/* format == 4 ! */
		  klo=pair_pot.first[col];
		  khi=pair_pot.last[col];
		  /* bisection */
		  while (khi-klo > 1) {
		    slot=(khi+klo) >> 1;
		    if ( pair_pot.xcoord[slot] > r ) khi=slot;
		    else klo=slot;
		  }
		  slot=klo;
		/* Check if we are at the last index - we should be lower */
                /* should be impossible anyway */
		/*   if (slot>=pair_pot.last[col]) {  */
 		/*    klo--;khi--; */
 		/*  } */
		  step=pair_pot.xcoord[khi]-pair_pot.xcoord[klo];
		  shift=(r-pair_pot.xcoord[klo])/step;
		}
		  
                /* Check if we are at the last index */
		if (slot>=pair_pot.last[col]) {
		  slot--;shift+=1.0;
		}
		atoms[i].neigh[k].shift[1]  = shift;
		atoms[i].neigh[k].slot[1]   = slot;
		atoms[i].neigh[k].step[1]   = step;
#endif		
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

  mdim=3*natoms+7*nconf;       /* mdim is dimension of force vector 
				3*natoms are real forces, 
				nconf cohesive energies,
			        6*nconf stress tensor components*/ 
#ifdef EAM  
  mdim+=2*ntypes;          /* ntypes dummy constraints */
  mdim+=nconf;		/* nconf limiting constraints */
#endif /* EAM */
  /* copy forces into single vector */
  if (NULL==(force_0 = (real *) malloc( mdim * sizeof(real) ) ) )
    error("Cannot allocate forces");
  k = 0;
  for (i=0; i<natoms; i++) {	/* first forces */
    force_0[k++] = atoms[i].force.x;
    force_0[k++] = atoms[i].force.y;
    force_0[k++] = atoms[i].force.z;
  }
  for (i=0; i<nconf; i++) {	/* then cohesive energies */
    force_0[k++] = coheng[i] * eweight; }
#ifdef STRESS
  for (i=0; i<nconf; i++) {	/* then stresses */
    force_0[k++] = stress[i].xx * sweight;
    force_0[k++] = stress[i].yy * sweight;
    force_0[k++] = stress[i].zz * sweight;
    force_0[k++] = stress[i].xy * sweight;
    force_0[k++] = stress[i].yz * sweight;
    force_0[k++] = stress[i].zx * sweight;
  }
#else                           
  for (i=0; i<6*nconf; i++)
    force_0[k++] = 0.;
#endif /* STRESS */
#ifdef EAM
  for(i=0; i<nconf; i++) force_0[k++]=0.; /* punishment rho out of bounds */
  for (i=0; i<2*ntypes;i++) {  /* constraint on U(n=0):=0 */
                               /* XXX and U'(n_mean)=0  */
    force_0[k++]=0.;}
#endif

#ifdef DEBUG
  printf("Minimal Distances Matrix \n");
  printf("Atom\t"); 
  for (i=0;i<ntypes;i++) printf("%8d\t",i);
  printf("with\n"); 
  for (i=0;i<ntypes;i++) { 
    printf("%d\t",i);
    for (j=0;j<ntypes;j++)
      printf("%f\t",mindist[ntypes*i+j]);
    printf("\n");
  }
#endif

  /* print diagnostic message and close file */
  printf("Maximal number of neighbors is %d, MAXNEIGH is %d\n",
         maxneigh, MAXNEIGH );
  printf("Read %d configurations with a total of %d atoms\nfrom file %s\n", nconf, natoms,filename);
  return;
}
