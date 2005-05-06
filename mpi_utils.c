/********************************************************
 *
 *  mpi_utils.c: Contains utilities to be used with MPI
 *
 *******************************************************/
/*
*   Copyright 2004-2005 Peter Brommer, Franz G"ahler
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
* $Revision: 1.10 $
* $Date: 2005/05/06 13:38:31 $
*****************************************************************/

#include "potfit.h"
#define N 20
#ifdef MPI

/******************************************************************************
*
* set up mpi
*
******************************************************************************/

void init_mpi(int *argc_pointer, char **argv)
{
  /* Initialize MPI */
  MPI_Init(argc_pointer, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if (0 == myid) {
    printf("%s\n", argv[0]);
    printf("Starting up MPI with %d processes.\n", num_cpus);
  }
}


/******************************************************************************
*
* shut down mpi
*
******************************************************************************/

void shutdown_mpi(void)
{
  MPI_Barrier(MPI_COMM_WORLD);   /* Wait for all processes to arrive */
  MPI_Finalize();                /* Shutdown */
}

void dbb(int i) {
/*   int j,k; */
/*   MPI_Status status; */
/*   if (myid==0)  */
/*     for (k=1;k<num_cpus;k++) { */
/*       MPI_Send(&i,1,MPI_INT,k,12,MPI_COMM_WORLD); */
/*       printf("Node 0, sent %d to node %d\n",i,k); */
/*     } */
/*   else {  */
/*     MPI_Recv(&j,1,MPI_INT,0,12,MPI_COMM_WORLD,&status); */
/*     printf("Node %d, received %d at point %d.\n",myid,j,i); */
/*   } */
  MPI_Barrier(MPI_COMM_WORLD);
  printf("I am node %d. Still running at point %d.\n",myid,i);
  fflush(stdout);
}

/***************************************************************************
 *
 * broadcast_param: Broadcast parameters etc to other nodes
 *
 **************************************************************************/
  /* 9: number of entries in struct atom_t  */
#define MAX_MPI_COMPONENTS 9

void broadcast_params() {
  int ierr,blklens[MAX_MPI_COMPONENTS]; 
  MPI_Aint displs[MAX_MPI_COMPONENTS]; 
  MPI_Datatype typen[MAX_MPI_COMPONENTS];
  neigh_t testneigh;
  atom_t testatom;
  int size,i,each,odd,nodeatoms=0;


  /* Define Structures */
  /* first the easy ones: */
  /* MPI_VEKTOR */
  MPI_Type_contiguous(3,REAL,&MPI_VEKTOR);
//MPI_Type_commit(&MPI_VEKTOR);
  /* MPI_STENS */
  MPI_Type_contiguous(6,REAL,&MPI_STENS);
//  MPI_Type_commit(&MPI_STENS);

  /* MPI_NEIGH */
  blklens[0]=1;         typen[0]=MPI_INT;  /* typ */
  blklens[1]=1;         typen[1]=MPI_INT;  /* nr */
  blklens[2]=1;         typen[2]=REAL;     /* r */
  blklens[3]=1;         typen[3]=MPI_VEKTOR; /* dist */
  blklens[4]=2;         typen[4]=MPI_INT; /* slot */
  blklens[5]=2;         typen[5]=REAL; /* shift */
  blklens[6]=2;         typen[6]=REAL; /* step */
  MPI_Address(&testneigh.typ,displs);
  MPI_Address(&testneigh.nr,&displs[1]);
  MPI_Address(&testneigh.r,&displs[2]);
  MPI_Address(&testneigh.dist,&displs[3]);
  MPI_Address(testneigh.slot,&displs[4]);
  MPI_Address(testneigh.shift,&displs[5]);
  MPI_Address(testneigh.step,&displs[6]);
  
  for(i=1;i<7;i++) {
      displs[i]-=displs[0];
  }
  displs[0]=0; /* set displacements */
  MPI_Type_struct(7,blklens,displs,typen,&MPI_NEIGH );
  MPI_Type_commit(&MPI_NEIGH);

  /* MPI_ATOM */
  blklens[0]=1;         typen[0]=MPI_INT; /* typ */
  blklens[1]=1;         typen[1]=MPI_INT; /* n_neigh */
  blklens[2]=1;         typen[2]=MPI_VEKTOR;    /* pos */
  blklens[3]=1;         typen[3]=MPI_VEKTOR;    /* force */
  blklens[4]=1;         typen[4]=REAL; /* absforce */
  blklens[5]=MAXNEIGH;  typen[5]=MPI_NEIGH; /* neigh */
  blklens[6]=1;         typen[6]=MPI_INT; /* conf */
  size=7;
#ifdef EAM
  blklens[7]=1;         typen[7]=REAL; /* rho */
  blklens[8]=1;         typen[8]=REAL; /* gradF */
  size+=2;
#endif
  MPI_Address(&testatom.typ,&displs[0]);
  MPI_Address(&testatom.n_neigh,&displs[1]);
  MPI_Address(&testatom.pos,&displs[2]);
  MPI_Address(&testatom.force,&displs[3]);
  MPI_Address(&testatom.absforce,&displs[4]);
  MPI_Address(testatom.neigh,&displs[5]);
  MPI_Address(&testatom.conf,&displs[6]);
#ifdef EAM
  MPI_Address(&testatom.rho,&displs[7]);
  MPI_Address(&testatom.gradF,&displs[8]);
#endif
  for(i=1;i<size;i++) {
    displs[i]-=displs[0];
  }
  displs[0]=0; /* set displacements */

  MPI_Type_struct(size,blklens,displs,typen,&MPI_ATOM);
  MPI_Type_commit(&MPI_ATOM);

  /* Distribute fundamental parameters */
  MPI_Bcast( &mdim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &anneal_temp, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &opt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid>0) {
    inconf=(int *) malloc(nconf*sizeof(int));
    cnfstart=(int *) malloc(nconf*sizeof(int));
    force_0=(real *) malloc(mdim*sizeof(real));
  }
  MPI_Bcast( inconf, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( cnfstart, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( force_0, mdim, REAL, 0, MPI_COMM_WORLD);

  /* Broadcast weights... */
  MPI_Bcast ( &eweight, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast ( &sweight, 1, REAL, 0, MPI_COMM_WORLD);
  /* Broadcast the potential... */
  MPI_Bcast(&format,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&pair_pot.len,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&pair_pot.idxlen,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&pair_pot.ncols,1,MPI_INT,0,MPI_COMM_WORLD);
  size=pair_pot.ncols;
  ndimtot=pair_pot.len;
  ndim=pair_pot.idxlen;
  if (myid>0) {
    pair_pot.begin = (real *) malloc(size * sizeof(real));
    pair_pot.end = (real *) malloc(size * sizeof(real));
    pair_pot.step = (real *) malloc(size * sizeof(real));
    pair_pot.invstep = (real *) malloc(size * sizeof(real));
    pair_pot.first = (int *) malloc(size * sizeof(int));
    pair_pot.last = (int *) malloc(size * sizeof(int));
    pair_pot.table = (real *) malloc(ndimtot * sizeof(real)); 
    pair_pot.xcoord = (real *) malloc(ndimtot * sizeof(real));
    pair_pot.d2tab = (real *) malloc(ndimtot * sizeof(real)); 
    pair_pot.idx = (int *) malloc(ndim * sizeof(int));
  }
  MPI_Bcast(pair_pot.begin,size,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.end,size,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.step,size,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.invstep,size,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.first,size,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.last,size,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.table,ndimtot,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.d2tab,ndimtot,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.xcoord,ndimtot,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.idx,ndim,MPI_INT,0,MPI_COMM_WORLD);

  /* Distribute configurations */
  /* Each node: nconf/num_cpus configurations. 
     Last nconf%num_cpus nodes: 1 additional config */
  each=(nconf/num_cpus);
  odd=(nconf%num_cpus)-num_cpus;
  if (myid==0){
    atom_len = (int *) malloc(num_cpus * sizeof(int));
    atom_dist  = (int *) malloc(num_cpus * sizeof(int));
    conf_len = (int *) malloc(num_cpus * sizeof(int));
    conf_dist  = (int *) malloc(num_cpus * sizeof(int));
    for(i=0;i<num_cpus;i++) 
      conf_dist[i]=i*each+(((i+odd)>0)?(i+odd):0);
    for(i=0;i<num_cpus-1;i++) conf_len[i]=conf_dist[i+1]-conf_dist[i];
    conf_len[num_cpus-1]=nconf-conf_dist[num_cpus-1];
    for(i=0;i<num_cpus;i++) 
      atom_dist[i]=cnfstart[i*each+(((i+odd)>0)?(i+odd):0)];
    for(i=0;i<num_cpus-1;i++) atom_len[i]=atom_dist[i+1]-atom_dist[i];
    atom_len[num_cpus-1]=natoms-atom_dist[num_cpus-1];
  }
  MPI_Scatter(atom_len,1,MPI_INT,&myatoms,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(atom_dist,1,MPI_INT,&firstatom,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(conf_len,1,MPI_INT,&myconf,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(conf_dist,1,MPI_INT,&firstconf,1,MPI_INT,0,MPI_COMM_WORLD);
  conf_atoms = (atom_t *) malloc(myatoms*sizeof(atom_t));
  MPI_Scatterv(atoms,atom_len,atom_dist,MPI_ATOM,
		 conf_atoms,myatoms,MPI_ATOM,
		 0,MPI_COMM_WORLD);
  conf_vol = (real *) malloc(myconf*sizeof(real));
  MPI_Scatterv(volumen,conf_len,conf_dist,REAL,
	       conf_vol,myconf,REAL,0,MPI_COMM_WORLD);
}

/***************************************************************************
 *
 * potsync: Broadcast parameters etc to other nodes
 *
 **************************************************************************/
#ifdef EAM

void potsync() {
  int firstcol,firstval,nvals;
  firstcol=paircol+ntypes;
  /* Memory is allocated - just bcast that changed potential... */
  /* bcast begin/end/step/invstep of embedding energy  */
  MPI_Bcast(pair_pot.begin+firstcol,ntypes,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.end+firstcol,ntypes,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.step+firstcol,ntypes,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.invstep+firstcol,ntypes,REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pair_pot.first+firstcol,ntypes,MPI_INT,0,MPI_COMM_WORLD);
  /* bcast table values of transfer fn. and embedding energy */
  firstval=pair_pot.first[paircol];
  nvals=ndimtot-firstval;
  MPI_Bcast(pair_pot.table+firstval,nvals,REAL,0,MPI_COMM_WORLD);
  /* das war's auch schon... */
}
#endif
#endif

