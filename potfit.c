/****************************************************************
* 
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/

/****************************************************************
* $Revision: 1.25 $
* $Date: 2004/04/13 14:10:11 $
*****************************************************************/


#define MAIN

#include "potfit.h"

/******************************************************************************
*
*  error -- complain and abort
*
******************************************************************************/

void error(char *msg)
{
  real* force;
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
#ifdef MPI
  calc_forces(pair_pot.table,force,1); /* go wake up other threads */
  shutdown_mpi();
#endif
  exit(2);
}

/******************************************************************************
 *
 *  warning -- just complain, don't abort
 *
 *****************************************************************************/

void warning(char *msg)
{
    fprintf(stderr,"Warning: %s\n",msg);
    fflush(stderr);
    return;
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv)
{
  real *force;
  real tot, min, max, sqr;
  int  i, diff;
  pi = 4.0 * atan(1.);
#ifdef MPI
  init_mpi(&argc,argv);
#endif
  srandom(seed+myid);random();random();random();random();
  if (myid==0) {
    read_parameters(argc, argv);
#ifdef EAM
    /* check whether dummy_phi was specified... */
    if (NULL==dummy_phi) 
      error("dummy_phi not specified in parameter file, aborting...");
#endif
    read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
    read_config(config);
#ifdef EAM
    printf("Dummy constraints are:\n");
    printf("rho[%f]\t = %f \t for atom type %d\n", 
	     dummy_r,dummy_rho,DUMMY_COL_RHO);
    for (i=0;i<ntypes;i++) {
      printf("phi[%f]\t = %f \t between atoms of type %d\n",
	     dummy_r,dummy_phi[i],i);
    }
#ifndef NORESCALE
    rescale(&pair_pot,1.,1); 	/* rescale now... */
    embed_shift(&pair_pot);	/* and shift */
#endif NORESCALE
#endif EAM
  }
#ifdef MPI
  broadcast_params(); 	     /* let the others know what's going on */
#endif MPI
 /*   mdim=3*natoms+nconf; */
  ndim=pair_pot.idxlen;
  ndimtot=pair_pot.len;
  paircol=(ntypes*(ntypes+1))/2;
  idx=pair_pot.idx;
  calc_forces = calc_forces_pair;
  force = (real *) malloc( (mdim) * sizeof(real) );

  if (myid > 0) 		/* all but root go to calc_forces */
    calc_forces(pair_pot.table,force,0);
  else {			/* root thread does minimization */
    if (opt) {
      anneal(pair_pot.table);
      powell_lsq(pair_pot.table);
    }
/*  for (i=0; i<pair_pot.ncols; i++) 
      spline_ed(pair_pot.step[i],pair_pot.table+pair_pot.first[i],
      pair_pot.last[i]-pair_pot.first[i]+1,
      1e30,0,pair_pot.d2tab+pair_pot.first[i]);*/

//    rescale(&pair_pot,1.);
    tot = calc_forces(pair_pot.table,force,0);
    write_pot_table( &pair_pot, endpot );
    printf("Potential written to file %s\n",endpot);
    printf("Plotpoint file written to file %s\n", plotpointfile);
    
    write_pot_table_imd( &pair_pot, imdpot );
    if (plot) write_plotpot_pair(&pair_pot, plotfile);
#ifdef PDIST
#ifndef MPI 			/* will not work with MPI */
    write_pairdist(&pair_pot,distfile);
#endif
#endif
#ifdef EAM
    printf("Local electron density rho\n");
    for (i=0; i<natoms;i++) printf("%d %d %f\n",i,atoms[i].typ,atoms[i].rho);
#endif
    
    max = 0.0;
    min = 100000.0;
    
    for (i=0; i<3*natoms; i++) {
      sqr = SQR(force[i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n",i/3,sqr,
	     force[i]+force_0[i],force_0[i],force[i]/force_0[i]);
    }
    printf("Cohesive Energies\n");
    for (i=0; i<nconf; i++){
      sqr = SQR(force[3*natoms+i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", i, sqr, force[3*natoms+i]+force_0[3*natoms+i],
	     force_0[3*natoms+i],force[3*natoms+i]/force_0[3*natoms+i]);
    }
#ifdef STRESS
    printf("Stresses on unit cell\n");
    for (i=3*natoms+nconf; i<3*natoms+7*nconf; i++) {
      sqr = SQR(force[i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", (i-(3*natoms+nconf))/6, sqr, 
	     force[i]+force_0[i], force_0[i],force[i]/force_0[i]);
    }    
#endif 
#ifdef EAM
#ifdef LIMIT
    printf("Punishment Constraints\n");
#ifdef STRESS
    diff = 6*nconf;
#else
    diff = 0;
#endif
    for (i=3*natoms+nconf+diff; i<3*natoms+2*nconf+diff; i++){
      sqr = SQR(force[i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", i-(3*natoms+nconf+diff), sqr, 
	     force[i]+force_0[i],
	     force_0[i], 
	     force[i]/force_0[i]);
    }
#endif    
    printf("Dummy Constraints\n");
    sqr = SQR(force[mdim-(2*ntypes+1)]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    printf("%d %f %f %f %f\n", 0, sqr, 
	   force[mdim-(2*ntypes+1)]+force_0[mdim-(2*ntypes+1)],
	   force_0[mdim-(2*ntypes+1)],
	   force[mdim-(2*ntypes+1)]/force_0[mdim-(2*ntypes+1)]);
    for (i=2*ntypes; i>0; i--){
      sqr = SQR(force[mdim-i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", 2*ntypes+1-i, sqr, 
	     force[mdim-i]+force_0[mdim-i],
	     force_0[mdim-i],force[mdim-i]/force_0[mdim-i]);
    }
#endif
    printf("av %e, min %e, max %e\n", tot/mdim, min, max);
    printf("Used %d function evaluations.\n",fcalls);
    calc_forces(pair_pot.table,force,1); /* go wake up other threads */
  }
#ifdef MPI
  /* kill MPI */
  shutdown_mpi();
#endif
  return 0;
}
