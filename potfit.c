
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

  mdim=3*natoms+nconf;
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

  force = (real *) malloc( 3 * natoms * sizeof(real) );
  tot = calc_forces(pair_pot.table,force);
  write_pot_table( &pair_pot, endpot );
  write_pot_table_imd( &pair_pot, imdpot );
  if (plot) write_plotpot_pair(&pair_pot, plotfile);


  max = 0.0;
  min = 100000.0;
  for (i=0; i<3*natoms; i++) {
    sqr = SQR(force[i]);
    max = MAX( max, sqr );
    min = MIN( min, sqr );
    
    printf("%d %f %f %f\n",i/3,sqr,force[i]+force_0[i],force_0[i]);
    
  }
  printf("av %e, min %e, max %e\n", tot/(3*natoms), min, max);
  return 0;
}
