
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
*  main
*
******************************************************************************/

int main(int argc, char **argv)
{
  real *force, tot, max;
  int  i;

  read_parameters(argc, argv);
  read_config(config);
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  write_pot_table( &pair_pot, endpot );
  write_pot_table_imd( &pair_pot, imdpot );

  force = (real *) malloc( 3 * natoms * sizeof(real) );
  tot = calc_forces_pair(force);

  max = 0.0;
  for (i=0; i<3*natoms; i++) max = MAX( max, SQR(force[i]-force_0[i]) );
  printf("av %e, max %e\n", tot/(3*natoms), max);

}
