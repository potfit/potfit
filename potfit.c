/****************************************************************
* 
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/
/*
*   Copyright 2002-2008 Peter Brommer, Franz G"ahler
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*  
*   This file is part of potfit.
*
*   potfit is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   potfit is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor, 
*   Boston, MA  02110-1301  USA
*/
/****************************************************************
* $Revision: 1.39 $
* $Date: 2008/04/02 15:05:39 $
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
  real* force=NULL;
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
#ifdef MPI
  calc_forces(calc_pot.table,force,1); /* go wake up other threads */
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
  real tot, min, max, sqr,*totdens;
  int  i, diff,*ntyp;
  pi = 4.0 * atan(1.);
#ifdef MPI
  init_mpi(&argc,argv);
#endif
  srandom(seed+myid);random();random();random();random();
  calc_forces = calc_forces_pair;
  if (myid==0) {
    read_parameters(argc, argv);
    read_pot_table( &opt_pot, startpot, ntypes*(ntypes+1)/2 );
    read_config(config);
    printf("Energy weight: %f\n",eweight);
#ifdef STRESS
    printf("Stress weight: %f\n",sweight);
#endif
  /* Select correct spline interpolation and other functions */
    if (format==3) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table= write_pot_table3;
#ifdef PARABEL
      parab = parab_ed;
      parab_comb = parab_comb_ed;
      parab_grad = parab_grad_ed;
#endif /* PARABEL */
    } else { /*format == 4 ! */
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
      write_pot_table=write_pot_table4;
#ifdef PARABEL
      parab = parab_ne;
      parab_comb = parab_comb_ne;
      parab_grad = parab_grad_ne;
#endif /* PARABEL */
      
    }
    /* set spline density corrections to 0 */
    lambda = (real *) malloc(ntypes * sizeof(real));
    totdens = (real *) malloc(ntypes *sizeof(real));
    ntyp=(int*) malloc(ntypes*sizeof(int));
    for (i=0;i<ntypes;i++) lambda[i]=0.;

#ifdef EAM
#ifndef NORESCALE
    rescale(&opt_pot,1.,1); 	/* rescale now... */
#ifdef WZERO
//    embed_shift(&opt_pot);	/* and shift  - diabolical */
#endif /* WZERO */
#endif /* NORESCALE */
#endif /* EAM */
  }
#ifdef MPI
  broadcast_params(); 	     /* let the others know what's going on */
#else  /* MPI */
  /* Identify subset of atoms/volumes belonging to individual process
     with complete set of atoms/volumes */
  conf_atoms = atoms;
  conf_vol = volumen;
  conf_uf = useforce;
  conf_us = usestress;
#endif /* MPI */
 /*   mdim=3*natoms+nconf; */
  ndim=opt_pot.idxlen;
  ndimtot=opt_pot.len;
  paircol=(ntypes*(ntypes+1))/2;
  idx=opt_pot.idx;

  force = (real *) malloc( (mdim) * sizeof(real) );

  if (myid > 0) {
  /* Select correct spline interpolation and other functions */
    /* Root process has done this earlier */
    if (format==3) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table= write_pot_table3;
#ifdef PARABEL
      parab = parab_ed;
      parab_comb = parab_comb_ed;
      parab_grad = parab_grad_ed;
#endif /* PARABEL */
    } else { /*format == 4 ! */
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
      write_pot_table=write_pot_table4;
#ifdef PARABEL
      parab = parab_ne;
      parab_comb = parab_comb_ne;
      parab_grad = parab_grad_ne;
#endif /* PARABEL */
    }
  /* all but root go to calc_forces */
    calc_forces(calc_pot.table,force,0);
  }  else {			/* root thread does minimization */
    if (opt) {
      anneal(opt_pot.table);
      powell_lsq(opt_pot.table);
    }
/*  for (i=0; i<opt_pot.ncols; i++) 
      spline_ed(opt_pot.step[i],opt_pot.table+opt_pot.first[i],
      opt_pot.last[i]-opt_pot.first[i]+1,
      1e30,0,opt_pot.d2tab+opt_pot.first[i]);*/

//    rescale(&opt_pot,1.);
    tot = calc_forces(calc_pot.table,force,0);
    write_pot_table( &opt_pot, endpot );
    printf("Potential in format %d written to file %s\n",format,endpot);
    printf("Plotpoint file written to file %s\n", plotpointfile);
    write_pot_table_imd( &opt_pot, imdpot );
    if (plot) write_plotpot_pair(&opt_pot, plotfile);


#ifdef PDIST
#ifndef MPI 			/* will not work with MPI */
    write_pairdist(&opt_pot,distfile);
#endif
#endif
    if (format == 3) {		/* then we can also write format 4 */
      sprintf(endpot,"%s_4",endpot);
      write_pot_table4(&opt_pot,endpot);
      printf("Potential in format 4 written to file %s\n",endpot);
    }
#ifdef EAM 
#ifndef MPI /* Not much sense in printing rho when not communicated... */
    printf("Local electron density rho\n");
    for (i=0;i<ntypes;i++) {totdens[i]=0.; ntyp[i]=0;}
    for (i=0; i<natoms;i++) {
      printf("%d %d %f\n",i,atoms[i].typ,atoms[i].rho);
      totdens[atoms[i].typ]+=atoms[i].rho;
      ntyp[atoms[i].typ]++;
    }
    for(i=0;i<ntypes;i++) {
      totdens[i] /= (real) ntyp[i];
      printf("Average local electron density at atom sites type %d: %f\n",
	     i,totdens[i]);
    }
#ifdef NEWSCALE
    for (i=0;i<ntypes;i++) {
      lambda[i]=splint_grad(&opt_pot,opt_pot.table,
			    paircol+ntypes+i,totdens[i]);
      printf("lambda[%d] = %f \n", i, lambda[i]) ;
	}
    sprintf(plotfile,"%s_new",plotfile);
    sprintf(imdpot,"%s_new",imdpot);
    /* write new potential plotting table */
    if (plot) write_altplot_pair(&opt_pot, plotfile);
    /* write NEW imd potentials */
    write_pot_table_imd( &opt_pot, imdpot );
#endif /* NEWSCALE */


#endif /* MPI */
#endif /* EAM */

    max = 0.0;
    min = 100000.0;
//    printf("conf-atom df^2 f f0 df/f0 |f|\n ");
    for (i=0; i<3*natoms; i++) {
      sqr = SQR(force[i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
#ifdef FWEIGHT
      printf("%d-%d %f %f %f %f %f\n",atoms[i/3].conf,i/3,sqr,
	     force[i]*(FORCE_EPS+atoms[i/3].absforce)+force_0[i],
	     force_0[i],
	     (force[i]*(FORCE_EPS+atoms[i/3].absforce))/force_0[i],
	     atoms[i/3].absforce);
#else  /* FWEIGHT */
      printf("%d-%d %f %f %f %f\n",atoms[i/3].conf,i/3,sqr,
      force[i]+force_0[i],force_0[i],force[i]/force_0[i]);
#endif /* FWEIGHT */
    }
    printf("Cohesive Energies\n");
    printf("conf w*de^2 e e0 de/e0\n");
    for (i=0; i<nconf; i++){
      sqr = SQR(force[3*natoms+i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", i, sqr, force[3*natoms+i]+force_0[3*natoms+i],
	     force_0[3*natoms+i],force[3*natoms+i]/force_0[3*natoms+i]);
    }
#ifdef STRESS
    printf("Stresses on unit cell\n");
    printf("conf w*ds^2 s s0 ds/s0\n");
    for (i=3*natoms+nconf; i<3*natoms+7*nconf; i++) {
      sqr = SQR(force[i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", (i-(3*natoms+nconf))/6, sqr, 
	     force[i]+force_0[i], force_0[i],force[i]/force_0[i]);
    }    
#endif 
#ifdef EAM
    printf("Punishment Constraints\n");
//    printf("conf dp p p0 dp/p0");
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
    printf("Dummy Constraints\n");
    for (i=2*ntypes; i>0; i--){
      sqr = SQR(force[mdim-i]);
      max = MAX( max, sqr );
      min = MIN( min, sqr );
      printf("%d %f %f %f %f\n", ntypes-i, sqr, 
	     force[mdim-i]+force_0[mdim-i],
	     force_0[mdim-i],force[mdim-i]/force_0[mdim-i]);
    }
#endif
    printf("av %e, min %e, max %e\n", tot/mdim, min, max);
    printf("Sum %f, count %d\n", tot, mdim);
    printf("Used %d function evaluations.\n",fcalls);
#ifdef MPI
    calc_forces(calc_pot.table,force,1); /* go wake up other threads */
#endif /* MPI */
    free(ntyp);
    free(totdens);
    free(lambda);
  }
  free(force);
#ifdef MPI
  /* kill MPI */
  shutdown_mpi();
#endif
  return 0;
}
