/****************************************************************
* 
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/

/****************************************************************
* $Revision: 1.19 $
* $Date: 2003/12/11 12:44:47 $
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
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
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
  int  i;
  pi = 4.0 * atan(1.);
  read_parameters(argc, argv);
  srandom(seed);random();random();random();random();
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  read_config(config);
#ifdef EAM
  if (eam) {
    printf("Dummy constraints are:\n");
    printf("rho[%f]\t = %f \t for atom type %d\n", 
	   dummy_r,dummy_rho,DUMMY_COL_RHO);
    for (i=0;i<ntypes;i++) {
      printf("phi[%f]\t = %f \t between atoms of type %d\n",
	     dummy_r,dummy_phi[i],i);
    }
  }
#endif
 /*   mdim=3*natoms+nconf; */
  ndim=pair_pot.idxlen;
  ndimtot=pair_pot.len;
  idx=pair_pot.idx;

  calc_forces = calc_forces_pair;

  if (opt) {
      anneal(pair_pot.table);
      powell_lsq(pair_pot.table);
  }
/*  for (i=0; i<pair_pot.ncols; i++) 
      spline_ed(pair_pot.step[i],pair_pot.table+pair_pot.first[i],
		pair_pot.last[i]-pair_pot.first[i]+1,
		1e30,0,pair_pot.d2tab+pair_pot.first[i]);*/

  force = (real *) malloc( (mdim) * sizeof(real) );
  tot = calc_forces(pair_pot.table,force);
  write_pot_table( &pair_pot, endpot );
  printf("Potential written to file %s\n",endpot);
  printf("Plotpoint file written to file %s\n", plotpointfile);

  write_pot_table_imd( &pair_pot, imdpot );
  if (plot) write_plotpot_pair(&pair_pot, plotfile);

#ifdef EAM
  if (eam) {
    printf("Local electron density rho\n");
    for (i=0; i<natoms;i++) printf("%d %d %f\n",i,atoms[i].typ,atoms[i].rho);
  }
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
    printf("%d %f %f %f %f\n", i-(3*natoms+nconf), sqr, force[i]+force_0[i],
	   force_0[i],force[i]/force_0[i]);
  }    
#endif 
#ifdef EAM
  if (eam) {
#ifdef LIMIT
    printf("Punishment Constraints\n");
    for (i=0; i<nconf; i++){
      sqr = SQR(force[3*natoms+nconf+i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", i, sqr, 
	     force[3*natoms+nconf+i]+force_0[3*natoms+nconf+i],
	     force_0[3*natoms+nconf+i], 
	     force[3*natoms+nconf+i]/force_0[3*natoms+nconf+i]);
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
  }
#endif
  printf("av %e, min %e, max %e\n", tot/mdim, min, max);
  return 0;
}
