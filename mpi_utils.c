/********************************************************
 *
 *  mpi_utils.c: Contains utilities to be used with MPI
 *
 *******************************************************/

/****************************************************************
* $Revision: 1.1 $
* $Date: 2004/02/25 16:35:49 $
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

void broadcast_params() {
  int ierr,blklens[7];   /* 7: number of entries in struct atom_t  */
  MPI_Aint displs[7]; 
  MPI_Datatype typen[7];
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
  MPI_Address(&testneigh.typ,displs);
  MPI_Address(&testneigh.nr,&displs[1]);
  MPI_Address(&testneigh.r,&displs[2]);
  MPI_Address(&testneigh.dist,&displs[3]);
  MPI_Address(testneigh.slot,&displs[4]);
  MPI_Address(testneigh.shift,&displs[5]);
  
  for(i=1;i<6;i++) {
      displs[i]-=displs[0];
  }
  displs[0]=0; /* set displacements */
  MPI_Type_struct(6,blklens,displs,typen,&MPI_NEIGH );
  MPI_Type_commit(&MPI_NEIGH);

  /* MPI_ATOM */
  blklens[0]=1;         typen[0]=MPI_INT; /* typ */
  blklens[1]=1;         typen[1]=MPI_INT; /* n_neigh */
  blklens[2]=1;         typen[2]=MPI_VEKTOR;    /* pos */
  blklens[3]=1;         typen[3]=MPI_VEKTOR;    /* force */
  blklens[4]=MAXNEIGH;  typen[4]=MPI_NEIGH; /* neigh */
  blklens[5]=1;         typen[5]=MPI_INT; /* conf */
  size=6;
#ifdef EAM
  blklens[6]=1;         typen[6]=REAL; /* rho */
  size+=1;
#endif
    MPI_Address(&testatom.typ,&displs[0]);
    MPI_Address(&testatom.n_neigh,&displs[1]);
    MPI_Address(&testatom.pos,&displs[2]);
    MPI_Address(&testatom.force,&displs[3]);
    MPI_Address(testatom.neigh,&displs[4]);
    MPI_Address(&testatom.conf,&displs[5]);
#ifdef EAM
    MPI_Address(&testatom.rho,&displs[6]);
#endif
    for(i=1;i<size;i++) {
      displs[i]-=displs[0];
   }
    displs[0]=0; /* set displacements */
/*  printf("node %d: blklens: %d %d %d %d, displs: %d %d %d %d \n",
	 myid, blklens[0], blklens[1], blklens[2], blklens[3],
	 displs[0], displs[1], displs[2], displs[3]); */

  MPI_Type_struct(size,blklens,displs,typen,&MPI_ATOM);
  MPI_Type_commit(&MPI_ATOM);

  /* Distribute fundamental parameters */
  MPI_Bcast( &mdim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
  /* Broadcast the potential... */
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
/*  if (myid==0) {
    printf("Sent information... eg:\n");
    printf("Atom 0: x %f y %f z %f typ %d n_neigh %d \n", atoms[0].pos.x, atoms[0].pos.y, atoms[0].pos.z,atoms[0].typ,atoms[0].n_neigh);
    printf("Atom 0, Neigh 102: dist.x %f y %f z %f typ %d\n", atoms[0].neigh[102].dist.x, atoms[0].neigh[102].dist.y, atoms[0].neigh[102].dist.z, atoms[0].neigh[102].typ);
    printf("Atom %d: x %f y %f z %f typ %d n_neigh %d\n",natoms-1, atoms[natoms-1].pos.x, atoms[natoms-1].pos.y, atoms[natoms-1].pos.z,atoms[natoms-1].typ,atoms[natoms-1].n_neigh);
    printf("Atom %d, Neigh 102: dist.x %f y %f z %f typ %d\n", natoms-1, atoms[natoms-1].neigh[102].dist.x, atoms[natoms-1].neigh[102].dist.y, atoms[natoms-1].neigh[102].dist.z, atoms[natoms-1].neigh[102].typ);
  }
  printf("Received information... eg:\n");
  printf("Node %d, Atom 0: x %f y %f z %f typ %d n_neigh %d \n", myid, conf_atoms[0].pos.x, conf_atoms[0].pos.y, conf_atoms[0].pos.z,conf_atoms[0].typ,conf_atoms[0].n_neigh);
  printf("Node %d, Atom 0, Neigh 102: dist.x %f y %f z %f typ %d\n", myid, conf_atoms[0].neigh[102].dist.x, conf_atoms[0].neigh[102].dist.y, conf_atoms[0].neigh[102].dist.z, conf_atoms[0].neigh[102].typ);
  printf("Node %d, Atom %d: x %f y %f z %f typ %d n_neigh %d\n", myid, myatoms-1, conf_atoms[myatoms-1].pos.x, conf_atoms[myatoms-1].pos.y, conf_atoms[myatoms-1].pos.z,conf_atoms[myatoms-1].typ,conf_atoms[myatoms-1].n_neigh);
  printf("Node %d, Atom %d, Neigh 102: dist.x %f y %f z %f typ %d\n", myid, myatoms-1, conf_atoms[myatoms-1].neigh[102].dist.x, conf_atoms[myatoms-1].neigh[102].dist.y, conf_atoms[myatoms-1].neigh[102].dist.z, conf_atoms[myatoms-1].neigh[102].typ); */
  /*distribute additional information: energy, volume, stress..*/
/*  conf_eng = (real *) malloc(myconf*sizeof(real)); */
  conf_vol = (real *) malloc(myconf*sizeof(real));
/*  conf_stress = (stens *) malloc(myconf*sizeof(stens));*/
//  MPI_Scatterv(coheng,conf_len,conf_dist,REAL,
//	       conf_eng,myconf,REAL,0,MPI_COMM_WORLD);
  MPI_Scatterv(volumen,conf_len,conf_dist,REAL,
	       conf_vol,myconf,REAL,0,MPI_COMM_WORLD);
#ifdef STRESS
//  MPI_Scatterv(stress,conf_len,conf_dist,MPI_STENS,
//	       conf_stress,myconf,MPI_STENS,0,MPI_COMM_WORLD);
#endif STRESS
/*   printf("I now have %d atoms in %d configurations. I am node %d.\n", */
/* 	   myatoms,myconf,myid); */
}


#endif
