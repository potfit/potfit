/****************************************************************
* 
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/

/****************************************************************
* $Revision: 1.14 $
* $Date: 2003/05/16 12:17:01 $
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
  real *force, tot, min, max, sqr;
  int  i;
  pi = 4.0 * atan(1.);
  read_parameters(argc, argv);
  srandom(seed);random();random();random();random();
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  read_config(config);

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

  force = (real *) malloc( (3 * natoms + nconf) * sizeof(real) );
  tot = calc_forces(pair_pot.table,force);
  write_pot_table( &pair_pot, endpot );
  write_pot_table_imd( &pair_pot, imdpot );
  if (plot) write_plotpot_pair(&pair_pot, plotfile);

  printf("Local electron density rho\n");
  for (i=0; i<natoms;i++) printf("%d %f\n",i,atoms[i].rho);

  max = 0.0;
  min = 100000.0;
  
  for (i=0; i<3*natoms; i++) {
    sqr = SQR(force[i]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    printf("%d %f %f %f\n",i/3,sqr,force[i]+force_0[i],force_0[i]);
  }
  printf("Cohesive Energies\n");
  for (i=0; i<nconf; i++){
    sqr = SQR(force[3*natoms+i]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    printf("%d %f %f %f\n", i, sqr, force[3*natoms+i]+force_0[3*natoms+i],
	   force_0[3*natoms+i]);
  }
#ifdef EAM
  printf("Punishment Constraints\n");
  for (i=0; i<nconf; i++){
   sqr = SQR(force[3*natoms+nconf+i]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    printf("%d %f %f %f\n", i, sqr, 
	   force[3*natoms+nconf+i]+force_0[3*natoms+nconf+i],
	   force_0[3*natoms+nconf+i]);
  }
    
  printf("Dummy Constraints\n");
  for (i=2; i>=1; i--){
    sqr = SQR(force[mdim-i]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    printf("%d %f %f %f\n", 2-i, sqr, force[mdim-i]+force_0[mdim-i],
	   force_0[mdim-i]);
  }
#endif
  printf("av %e, min %e, max %e\n", tot/mdim, min, max);
  return 0;
}
