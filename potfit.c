
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

  read_parameters(argc, argv);
  read_config(config);
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  ndim=pair_pot.len;  
  if (opt) powell_lsq(pair_pot.table,calc_forces_pair);
  write_pot_table( &pair_pot, endpot );
  write_pot_table_imd( &pair_pot, imdpot );

  force = (real *) malloc( 3 * natoms * sizeof(real) );
  tot = calc_forces_pair(pair_pot.table,force);

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
