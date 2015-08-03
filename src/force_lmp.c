// Created by Kubo 20120518 based on "force_pair.c"

/****************************************************************
 *
 * force_pair.c: Routines used for calculating pair forces/energies
 *
 ****************************************************************
 *
 * Copyright 2002-2011
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
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
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#ifdef LMP
#include "potfit.h"

#include "functions.h"
#include "potential_input.h"
#include "splines.h"
#include "utils.h"
#include "time.h"

double lammps( // RETURN: Total Energy
int me,
int   natm, //-------- IN: Num of Atoms
int   ntps, //-------- IN: Num of Atomic Types (=ntypes)
double *lvec, //-------- IN: Lattice Vector Length (x,y,z)
double *avec, //-------- IN: Lattice Vector Angle (xy,xz,yz)
double *coordtmp, //---- IN: Coords[3N] (rx1,ry1,rz1,rx2,ry2,rz2,...)
int  *asp, //--------- IN: Atomic Species[N] (a1,a2,...)
double *forcetmp, //--- OUT: Forces[3N] (fx1,fy1,fz1,fx2,fy2,fz2,...)
double *potrearr,
void* lammpsOBJ //------- IN: Potaential Param[npot+tags]
);

double calc_forces(double *xi_opt, double *forces, int flag)
{
  void write_ff(int n,double *p);
  void rearrange_ff(int n,double *p,double *q);
  void show_ff(int n,double *p);

  int   first, col, i;
  double  tmpsum = 0., sum = 0.;
  double *xi = NULL;
  apot_table_t *apt = &g_pot.apot_table;

  // create communicator for each node
  //MPI_Comm new_comm;

  //MPI_Comm_split(MPI_COMM_WORLD, myid, myid, &new_comm);
  //MPI_Group orig_group, new_group;
  //MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
 // int ranks[1];
  //ranks[0] = myid;
  //MPI_Group_incl(orig_group, 1, ranks, &new_group);
  //printf("Entered force calculation: %d. Will create lammps object now.\r\n", myid);
  //fflush(stdout);  
//MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
//  clock_t t1, t2;

  //t1 = clock();
  //LAMMPS* lammpsOBJ = new LAMMPS(0,NULL, MPI_COMM_SELF); 
  //t2 = clock(); 
  //float diff = (((float)t2 - (float)t1) / 1000000.0F ) * 1000;   
  //printf("MyID: %d, LAMMPS Init.: %f ms\n", myid, diff);
  // Edited by Kubo 20120528 ===================================//
// From ------------------------------------------------------//
/*
  switch (format) {
      case 0:
        xi = calc_pot.table;
        break;
      case 3:			// fall through
      case 4:
        xi = xi_opt;		// calc-table is opt-table
        break;
      case 5:
        xi = calc_pot.table;	// we need to update the calc-table
  }
  */

// To --------------------------------------------------------//
  //if (flag != 1)
     xi = xi_opt;
   
  //if (myid > 0)
 //show_ff(4, xi);
 //return 0;
// End of Edition ============================================//
 
  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.;		/* sum of squares of local process */

#if defined APOT && !defined MPI
    if (format == 0) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif /* APOT && !MPI */

   #ifdef MPI
#ifndef APOT
    /* exchange potential and flag value */
    MPI_Bcast(xi,calc_pot.len,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif /* APOT */
    //if (flag==1) {
    //printf("OUTPUT: %d %d", myid, flag);
   
    //return 0;
    //break; }
    
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //printf("CLOSE: %d %d\n", myid, flag);
    //fflush(stdout);

    if(flag==1){break;}   /* Exception: flag 1 means clean up */

#ifdef APOT
    if(g_mpi.myid==0){apot_check_params(xi_opt);}
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else /* APOT */
    /* if flag==2 then the potential parameters have changed -> sync */
    if(flag==2){potsync();}
#endif /* APOT */
#endif /* MPI */

//return 0;

    /* init second derivatives for splines */

// Edited by Kubo 20120528 ===================================//
// Comm-Out
/*
    // pair potentials
    for(col=0;col<paircol;col++){
      first = calc_pot.first[col];
      if (format == 0 || format == 3)
        spline_ed(calc_pot.step[col], xi + first,
                  calc_pot.last[col] - first + 1,
                  *(xi + first - 2), 0.0,
                  calc_pot.d2tab + first);
      else			// format >= 4 ! 
        spline_ne(calc_pot.xcoord + first, xi + first,
                  calc_pot.last[col] - first + 1,
                  *(xi + first - 2), 0.0,
                  calc_pot.d2tab + first);
    }
*/
// End of Edition ============================================//

#ifndef MPI
    myconf = nconf;
#endif /* MPI */
  //printf("CHKPNT 1: %d\n", myid);
  //fflush(stdout);  
    /* region containing loop over configurations,
       also OMP-parallelized region */
    {
      atom_t *atom;
      int   h, j, k, l;
      int   self, uf;
#ifdef STRESS
      int   us, stresses;
#endif /* STRESS */

      neigh_t *neigh;

      // create communicator for each node
      //MPI_Comm new_comm;
      //MPI_Comm_split(MPI_COMM_WORLD, myid, myid, &new_comm);
      //LAMMPS* lammpsOBJ = new LAMMPS(0,NULL, new_comm /*MPI_COMM_WORLD/*); 

      /* pair variables */
      double  phi_val, phi_grad;
      vector tmp_force;

      /* loop over configurations */
      for(h=g_mpi.firstconf;h<g_mpi.firstconf+g_mpi.myconf;h++){
        uf = g_config.conf_uf[h-g_mpi.firstconf];
#ifdef STRESS
        us = g_config.conf_us[h-g_mpi.firstconf];
#endif /* STRESS */
	/* reset energies and stresses */
        forces[g_calc.energy_p+h] = 0.;

#ifdef STRESS
        for(i=0;i<6;i++){forces[g_calc.stress_p+6*h+i] = 0.;}
#endif /* STRESS */

// Edited by Kubo 20120522 ===//
// From ----------------------//
//#ifdef APOT
//          if (enable_cp)
//            forces[energy_p + h] +=
//                chemical_potential(ntypes, na_type[h], xi_opt + cp_start);
//#endif /* APOT */
// To ------------------------//
          forces[g_calc.energy_p+h] = g_config.force_0[g_calc.energy_p+h];
// End of Edition ============//

	/* first loop over atoms: reset forces, densities */
          for (i = 0; i < g_config.inconf[h]; i++) {
            if (uf) {
              k = 3 * (g_config.cnfstart[h] + i);
              forces[k] = -g_config.force_0[k];
              forces[k + 1] = -g_config.force_0[k + 1];
              forces[k + 2] = -g_config.force_0[k + 2];
            } else {
              k = 3 * (g_config.cnfstart[h] + i);
              forces[k] = 0.;
              forces[k + 1] = 0.;
              forces[k + 2] = 0.;
            }
          }
	/* end first loop */


//====== This Loop is to be Modified...! =====================//
  int   natm = g_config.inconf[h]; //----- Num of Atoms in h-th Config
  int   nprm = apt->total_par; // Num of Params in ReaxFF Function
  int   asp[natm]; //------------ Atomic Species List
  double *pottmp=xi_opt; //-------- Potential Param
  double  coordtmp[3*natm]; //----- 
  double  forcetmp[3*natm]; //----- Force Calced by Subroutine
  double  forcetmp_sep[3*natm]; // force calculated by subroutine in the case if separation is needed
  double  en_sep;
  double potrearr[rxp_num];
 
  for(i=0;i<natm;i++){
    atom = g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
    asp[i] = atom->type+1;
    coordtmp[3*i+0] = atom->pos.x;
    coordtmp[3*i+1] = atom->pos.y;
    coordtmp[3*i+2] = atom->pos.z;
    forcetmp[3*i+0] = 0.0;
    forcetmp[3*i+1] = 0.0;
    forcetmp[3*i+2] = 0.0; 
  }

  double lvec[3],avec[3];
  double lat_a = sqrt(g_config.lattice[h].xx*g_config.lattice[h].xx+g_config.lattice[h].xy*g_config.lattice[h].xy+g_config.lattice[h].xz*g_config.lattice[h].xz);
  double lat_b = sqrt(g_config.lattice[h].yx*g_config.lattice[h].yx+g_config.lattice[h].yy*g_config.lattice[h].yy+g_config.lattice[h].yz*g_config.lattice[h].yz);
  double lat_c = sqrt(g_config.lattice[h].zx*g_config.lattice[h].zx+g_config.lattice[h].zy*g_config.lattice[h].zy+g_config.lattice[h].zz*g_config.lattice[h].zz);
  double lat_cosalpha = (g_config.lattice[h].xx*g_config.lattice[h].yx+g_config.lattice[h].xy*g_config.lattice[h].yy+g_config.lattice[h].xz*g_config.lattice[h].yz)/(lat_a*lat_b);
  double lat_cosbeta  = (g_config.lattice[h].xx*g_config.lattice[h].zx+g_config.lattice[h].xy*g_config.lattice[h].zy+g_config.lattice[h].xz*g_config.lattice[h].zz)/(lat_a*lat_c);
  double lat_cosgamma = (g_config.lattice[h].yx*g_config.lattice[h].zx+g_config.lattice[h].yy*g_config.lattice[h].zy+g_config.lattice[h].yz*g_config.lattice[h].zz)/(lat_b*lat_c);

  // artificial increase of cell size along z direction to be able to move Ni Slab
  double d = 0;

  avec[0] = lat_b*lat_cosalpha;
  lvec[0] = lat_a;
  avec[1] = lat_c*lat_cosbeta;
  lvec[1] = sqrt(lat_b*lat_b-avec[0]*avec[0]);
  avec[2] = (lat_b*lat_c*lat_cosgamma-avec[0]*avec[1])/lvec[1];
  lvec[2] = sqrt(lat_c*lat_c-avec[1]*avec[1]-avec[2]*avec[2]) + 2 * d;
  
  rearrange_ff(g_param.ntypes,xi,potrearr);
  forces[g_calc.energy_p+h] = lammps(g_mpi.myid, g_config.inconf[h], g_param.ntypes, lvec, avec, coordtmp, asp, forcetmp, potrearr, g_pot.lammpsObj); 

  
  // shift of ni
  /*for(i=0;i<natm;i++){
    if (asp[i] == nickelId + 1)
	coordtmp[3*i+2] +=d;
  }*/
  
  // subtract basis B2 nizr state
  lvec[0] = 3.22;
  lvec[1] = 3.22;
  lvec[2] = 3.22;

  coordtmp[0] = 0;
  coordtmp[1] = 0;
  coordtmp[2] = 0;
  coordtmp[3] = 1.61;
  coordtmp[4] = 1.61;
  coordtmp[5] = 1.61;
  // end of subtract B2 Nizr
  
  forces[g_calc.energy_p+h] -= lammps(g_mpi.myid, g_config.inconf[h], g_param.ntypes, lvec, avec, coordtmp, asp, forcetmp_sep, potrearr, g_pot.lammpsObj);
  
  if(uf){
    for(i=0;i<g_config.inconf[h];i++){
      k = 3*(g_config.cnfstart[h]+i);
      forces[k]   += forcetmp[3*i+0] - forcetmp_sep[3*i+0];
      forces[k+1] += forcetmp[3*i+1] - forcetmp_sep[3*i+1];
      forces[k+2] += forcetmp[3*i+2] - forcetmp_sep[3*i+2];

     //printf("%d %d %f %f %f  %f  %f  %f\n", asp[i], i+1, coordtmp[3*i+0], coordtmp[3*i+1], coordtmp[3*i+2],  forcetmp[3*i+0]*23.045126829,forcetmp[3*i+1]*23.045126829, forcetmp[3*i+2]*23.045126829);
      //printf("OR:  %f  %f  %f\n",forces[k],forces[k+1],forces[k+2]);
    }
  } 

 // printf("CHKPNT 2: %d\n", myid);
  //fflush(stdout);  
 
  //return 7;
  // 2nd loop: calculate pair forces and energies
	for(i=0;i<g_config.inconf[h];i++){
	  atom = g_config.conf_atoms+i+g_config.cnfstart[h]-g_mpi.firstatom;
	  k = 3*(g_config.cnfstart[h]+i);
	//then we can calculate contribution of forces right away
	  if(uf){
#ifdef FWEIGHT
	    // Weigh by absolute value of force
	    forces[k]     /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif // FWEIGHT
	    // sum up forces
	    tmpsum += g_config.conf_weight[h]
             *(dsquare(forces[k])
              +dsquare(forces[k+1])
              +dsquare(forces[k+2]));
	  }			// second loop over atoms
	}
//====== End of Loop to be Modified ==========================//
   //printf("CHKPOINT 2.1: %d\n", myid);
   //fflush(stdout);
	/* energy contributions */
	forces[g_calc.energy_p+h] /= (double)g_config.inconf[h]; // Total Ene to Ene Per-Atom
	forces[g_calc.energy_p+h] -= g_config.force_0[g_calc.energy_p+h];
#ifdef COMPAT
	tmpsum += g_config.conf_weight[h]*dsquare(g_param.eweight*forces[g_calc.energy_p+h]);
#else
	tmpsum += g_config.conf_weight[h]*g_param.eweight*dsquare(forces[g_calc.energy_p+h]);
#endif /* COMPAT */
#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stress_p+6*h+i] /= conf_vol[h-firstconf];
	    forces[stress_p+6*h+i] -= g_config.force_0[stress_p+6*h+i];
	    tmpsum +=
#ifdef COMPAT
	      conf_weight[h]*dsquare(sweight*forces[stress_p+6*h+i]);
#else
	      conf_weight[h]*sweight*dsquare(forces[stress_p+6*h+i]);
#endif /* COMPAT */
	  }
	}
#endif /* STRESS */
	/* limiting constraints per configuration */
      }				/* loop over configurations */
   }
   //printf("CHKPNT prefinal: %d %10.4f\n", myid, tmpsum);
//   fflush(stdout);  
 
 

		
    /* parallel region */
    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if(g_mpi.myid == 0){tmpsum += apot_punish(xi_opt,forces);}
   //printf("CHKPNT final: %d %10.4f\n", myid, tmpsum);

    
#endif /* APOT */
    //printf("\n%f %f\n\n", tmpsum, apot_punish(xi_opt,forces));
    sum = tmpsum;		/* global sum = local sum  */
#ifdef MPI
    //return 0;  
    /* reduce global sum */
    sum = 0.;
    MPI_Reduce(&tmpsum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   
   //printf("CHKPNT 4 id %d: %10.4f %10.4f\n", myid, tmpsum, sum);
//   fflush(stdout); 
   /* gather forces, energies, stresses */
    /* forces */
   MPI_Gatherv(forces+g_mpi.firstatom*3,g_mpi.myatoms, g_mpi.MPI_VECTOR,
        forces,g_mpi.atom_len,g_mpi.atom_dist, g_mpi.MPI_VECTOR,0,MPI_COMM_WORLD);
    /* energies */
    MPI_Gatherv(forces+g_config.natoms*3+g_mpi.firstconf,g_mpi.myconf,MPI_DOUBLE,
        forces+g_config.natoms*3,g_mpi.conf_len,g_mpi.conf_dist,MPI_DOUBLE,0,MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces+g_config.natoms*3+g_config.nconf+6*g_mpi.firstconf,g_mpi.myconf, g_mpi.MPI_STENS,
        forces+g_config.natoms*3+g_config.nconf,g_mpi.conf_len,g_mpi.conf_dist, g_mpi.MPI_STENS,0,MPI_COMM_WORLD);
#endif /* MPI */

    //printf("CHKPNT final: %d\n", myid);
  //fflush(stdout);  
 

    /* root process exits this function now */
    if(g_mpi.myid==0){
      g_calc.fcalls++;			/* Increase function call counter */
      if(isnan(sum)){
   #ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
   #endif /* DEBUG */
        //delete lammpsOBJ;
	return 10e10;
      } else{
	 //delete lammpsOBJ;
	return sum; }
    }

  }

  //delete lammpsOBJ;
  /* once a non-root process arrives here, all is done. */
  return -1.0;
}

// Rearrange Force-Field Array ===============================//
void rearrange_ff(
    int n, //--------- IN: Num of Atomic Species
    double *p, //----- IN: Potential Table to Arrange
    double *q //----- OUT: Arranged Table
){
  int i,j;
  int i1,i2,i3,i4;
  int sr0,sr1,sr2,sr3,sr4;
  int srO,srH;
  int comb0,comb1,comb2,comb3,comb4;
  int combO,combH;
  int offset;
  int iq; // Index for q

  iq = 0;

  sr0 = 39; comb0 = 1;
  sr1 = 32; comb1 = n;
  sr2 = 16; comb2 = n*(n+1)/2;
  sr3 =  7; comb3 = n*n*(n+1)/2;
  sr4 =  7; comb4 = n*n*(n*n+1)/2;
  srO =  6; combO = n*(n+1)/2;
  srH =  4; combH = n*n*n;

// 0-Body --------------------//
  q[iq] = sr0; iq++;
  for(i=0;i<sr0;i++){q[iq] = p[i]; iq++;}

// 1-Body --------------------//
  offset = sr0;
  q[iq] = comb1; iq++;

  j = 0;
  for(i1=0;i1<n;i1++){
    for(i= 0;i<32;i++){q[iq] = p[offset+i*comb1+j]; iq++;}
    j++;
  }

// 2-Body --------------------//
  offset += comb1*sr1;
  q[iq] = comb2; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    for(i= 0;i<16;i++){q[iq] = p[offset+i*comb2+j]; iq++;}
    j++;
  }
  }

// Off-Diag
  offset += comb2*sr2;
  q[iq] = combO; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    for(i= 0;i<srO;i++){q[iq] = p[offset+i*combO+j]; iq++;}
    j++;
  }
  }

// 3-Body --------------------//
  offset += combO*srO;
  q[iq] = comb3; iq++;

  j = 0;
  for(i2= 0;i2<n;i2++){ // <- Attention!! Order of Loop is 2-1-3
    for(i1= 0;i1<n;i1++){
    for(i3=i1;i3<n;i3++){ // i1 and i3 are Commutative. 
      q[iq] = i1+1; iq++;
      q[iq] = i2+1; iq++;
      q[iq] = i3+1; iq++;
      for(i=0;i<sr3;i++){q[iq] = p[offset+i*comb3+j]; iq++;}
      j++;
    }
    }
  }

// 4-Body --------------------//
  offset += comb3*sr3;
  q[iq] = comb4; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=0 ;i2<n;i2++){
  for(i3=0 ;i3<n;i3++){
  for(i4=i1;i4<n;i4++){
    if((i1==i4)&&(i2>i3)){continue;}
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    q[iq] = i3+1; iq++;
    q[iq] = i4+1; iq++;
    for(i=0;i<sr4;i++){q[iq] = p[offset+i*comb4+j]; iq++;}
    j++;
  }
  }
  }
  }

// H-Bond
  offset += comb4*sr4;
  q[iq] = combH; iq++;

  j = 0;
  for(i1=0;i1<n;i1++){
  for(i2=0;i2<n;i2++){
  for(i3=0;i3<n;i3++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    q[iq] = i3+1; iq++;
    for(i=0;i<srH;i++){q[iq] = p[offset+i*combH+j]; iq++;}
    j++;
  }
  }
  }

}
/*
void rearrange_ff(
    int n, //--------- IN: Num of Atomic Species
    double *p, //----- IN: Potential Table to Arrange
    double *q //----- OUT: Arranged Table
){
  int i,j;
  int i1,i2,i3,i4;
  int sr0,sr1,sr2,sr3,sr4;
  int sr2b,sr2o,sr2h;
  int comb1,comb2,comb3,comb4;
  int offset;
  int iq; // Index for q

  iq = 0;

  sr0 = 39;
  sr1 = 32;
  sr2 = 26; sr2b = 16; sr2o = 6; sr2h = 4;
  sr3 =  7;
  sr4 =  7;

// 0-Body --------------------//
  q[iq] = sr0; iq++;
  for(i=0;i<sr0;i++){
    q[iq] = p[i]; iq++;
  }

// 1-Body --------------------//
  offset = sr0;
  comb1 = n;
  q[iq] = comb1; iq++;

  j = 0;
  for(i1=0;i1<n;i1++){
    for(i= 0;i< 8;i++){q[iq] = p[offset+i*comb1+j]; iq++;}
    for(i= 8;i<16;i++){q[iq] = p[offset+i*comb1+j]; iq++;}
    for(i=16;i<24;i++){q[iq] = p[offset+i*comb1+j]; iq++;}
    for(i=24;i<32;i++){q[iq] = p[offset+i*comb1+j]; iq++;}
    j++;
  }

// 2-Body --------------------//
  offset = sr0 + comb1*sr1;
  comb2 = n*(n+1)/2;
  q[iq] = comb2; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    for(i= 0;i< 8;i++){q[iq] = p[offset+i*comb2+j]; iq++;}
    for(i= 8;i<16;i++){q[iq] = p[offset+i*comb2+j]; iq++;}
    j++;
  }
  }

// Off-Diag
  offset = sr0 + comb1*sr1 + comb2*sr2b;
  q[iq] = comb2; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    q[iq] = i1+1; iq++;
    q[iq] = i1+2; iq++;
    for(i= 0;i<sr2o;i++){q[iq] = p[offset+i*comb2+j]; iq++;}
    j++;
  }
  }

// 3-Body --------------------//
  offset = sr0 + comb1*sr1 + comb2*sr2;
  comb3 = n*(n+1)*(n+2)/6;
  q[iq] = comb3; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
  for(i3=i2;i3<n;i3++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    q[iq] = i3+1; iq++;
    for(i=0;i<sr3;i++){q[iq] = p[offset+i*comb3+j]; iq++;}
    j++;
  }
  }
  }

// 4-Body --------------------//
  offset = sr0 + comb1*sr1 + comb2*sr2 + comb3*sr3;
  comb4 = n*(n+1)*(n+2)*(n+3)/24;
  q[iq] = comb4; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
  for(i3=i2;i3<n;i3++){
  for(i4=i3;i4<n;i4++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    q[iq] = i3+1; iq++;
    q[iq] = i4+1; iq++;
    for(i=0;i<sr4;i++){q[iq] = p[offset+i*comb4+j]; iq++;}
    j++;
  }
  }
  }
  }

// H-Bond
  offset = sr0 + comb1*sr1 + comb2*(sr2b+sr2o);
  q[iq] = comb2; iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    q[iq] = i1+1; iq++;
    q[iq] = i2+1; iq++;
    for(i=0;i<sr2h;i++){q[iq] = p[offset+i*comb2+j]; iq++;}
    j++;
  }
  }

}
*/

/*
// Write Force-Field File ====================================//
void write_ff(
    int n, //----- Num of Atomic Species
    double *p //-- Potential Parameter Table
){
  int i,j;
  int i1,i2,i3,i4;
  int sr0,sr1,sr2,sr3,sr4;
  int sr2b,sr2o,sr2h;
  int comb1,comb2,comb3,comb4;
  int offset;
  FILE *fp;

  sr0 = 39;
  sr1 = 32;
  sr2 = 26; sr2b = 16; sr2o = 6; sr2h = 4;
  sr3 =  7;
  sr4 =  7;

  fp = fopen("ffield.reax","w");
  fprintf(fp,"Reactive MD-force field\n");

// 0-Body
  fprintf(fp,"%3d",sr0);
  fprintf(fp,"       ! Number of general parameters\n");
  for(i=0;i<sr0;i++){fprintf(fp,"%10.4f\n",p[i]);}

// 1-Body
  offset = sr0;
  comb1 = n;
  fprintf(fp,"%3d",comb1);
  fprintf(fp,"    ! Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;#\n");
  fprintf(fp,"            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n");
  fprintf(fp,"            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.\n");
  fprintf(fp,"            ov/un;val1;n.u.;val3,vval4\n");

  j = 0;
  for(i1=0;i1<n;i1++){
//    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3s",elements[i1]);
    for(i= 0;i< 8;i++){fprintf(fp,"%10.4f",p[offset+i*comb1+j]);}
    fprintf(fp,"\n");
    fprintf(fp,"   ");
    for(i= 8;i<16;i++){fprintf(fp,"%10.4f",p[offset+i*comb1+j]);}
    fprintf(fp,"\n");
    fprintf(fp,"   ");
    for(i=16;i<24;i++){fprintf(fp,"%10.4f",p[offset+i*comb1+j]);}
    fprintf(fp,"\n");
    fprintf(fp,"   ");
    for(i=24;i<32;i++){fprintf(fp,"%10.4f",p[offset+i*comb1+j]);}
    fprintf(fp,"\n");
    j++;
  }

// 2-Body
  offset = sr0 + comb1*sr1;
  comb2 = n*(n+1)/2;
  fprintf(fp,"%3d",comb2);
  fprintf(fp,"      ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n");
  fprintf(fp,"                         pbe2;pbo3;pbo4;Etrip;pbo1;pbo2;ovcorr\n");

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3d",i2+1);
    for(i= 0;i< 8;i++){fprintf(fp,"%10.4f",p[offset+i*comb2+j]);}
    fprintf(fp,"\n");
    fprintf(fp,"   ");
    fprintf(fp,"   ");
    for(i= 8;i<16;i++){fprintf(fp,"%10.4f",p[offset+i*comb2+j]);}
    fprintf(fp,"\n");
    j++;
  }
  }

// Off-Diag
  offset = sr0 + comb1*sr1 + comb2*sr2b;
  fprintf(fp,"%3d",comb2);
  fprintf(fp,"    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2\n");

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3d",i2+1);
    for(i= 0;i<sr2o;i++){fprintf(fp,"%10.4f",p[offset+i*comb2+j]);}
    fprintf(fp,"\n");
    j++;
  }
  }

// 3-Body
  offset = sr0 + comb1*sr1 + comb2*sr2;
  comb3 = n*(n+1)*(n+2)/6;
  fprintf(fp,"%3d",comb3);
  fprintf(fp,"    ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2\n");

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
  for(i3=i2;i3<n;i3++){
    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3d",i2+1);
    fprintf(fp,"%3d",i3+1);
    for(i=0;i<sr3;i++){fprintf(fp,"%10.4f",p[offset+i*comb3+j]);}
    fprintf(fp,"\n");
    j++;
  }
  }
  }

// 4-Body
  offset = sr0 + comb1*sr1 + comb2*sr2 + comb3*sr3;
  comb4 = n*(n+1)*(n+2)*(n+3)/24;
  fprintf(fp,"%3d",comb4);
  fprintf(fp,"    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n\n");

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
  for(i3=i2;i3<n;i3++){
  for(i4=i3;i4<n;i4++){
    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3d",i2+1);
    fprintf(fp,"%3d",i3+1);
    fprintf(fp,"%3d",i4+1);
    for(i=0;i<sr4;i++){fprintf(fp,"%10.4f",p[offset+i*comb4+j]);}
    fprintf(fp,"\n");
    j++;
  }
  }
  }
  }

// H-Bond
  offset = sr0 + comb1*sr1 + comb2*(sr2b+sr2o);
  fprintf(fp,"%3d",comb2);
  fprintf(fp,"    ! Nr of hydrogen bonds;at1;at2;at3;Rhb;Dehb;vhb1\n");

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    fprintf(fp,"%3d",i1+1);
    fprintf(fp,"%3d",i2+1);
    for(i=0;i<sr2h;i++){fprintf(fp,"%10.4f",p[offset+i*comb2+j]);}
    fprintf(fp,"\n");
    j++;
  }
  }

// End Process
  fclose(fp);
}
*/

// Show Rearranged FF Array ----------------------------------//
void show_ff(
    int n, //--------- IN: Num of Atomic Species
    double *p //------ IN: Arranged Table
){
  int i,j;
  int i1,i2,i3,i4;
  int sr0,sr1,sr2,sr3,sr4;
  int srO,srH;
  int comb0,comb1,comb2,comb3,comb4;
  int combO,combH;
  int offset;
  int iq; // Index for q

  iq = 0;

  sr0 = 39; comb0 = 1;
  sr1 = 32; comb1 = n;
  sr2 = 16; comb2 = n*(n+1)/2;
  sr3 =  7; comb3 = n*n*(n+1)/2;
  sr4 =  7; comb4 = n*n*(n*n+1)/2;
  srO =  6; combO = n*(n+1)/2;
  srH =  4; combH = n*n*n;

// 0-Body --------------------//
  printf("%4d\n",(int)p[iq]); iq++;
  for(i=0;i<sr0;i++){printf("%10.4f\n",p[iq]); iq++;}

// 1-Body --------------------//
  offset = sr0;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0;i1<n;i1++){
    for(i= 0;i< 8;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    for(i= 8;i<16;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    for(i=16;i<24;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    for(i=24;i<32;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }

// 2-Body --------------------//
  offset += comb1*sr1;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    for(i= 0;i< 8;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    printf("    ");
    printf("    ");
    for(i= 8;i<16;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }
  }

// Off-Diag
  offset += comb2*sr2;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=i1;i2<n;i2++){
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    for(i= 0;i<srO;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }
  }

// 3-Body --------------------//
  offset += combO*srO;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=0 ;i2<n;i2++){
  for(i3=i1;i3<n;i3++){
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    for(i=0;i<sr3;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }
  }
  }

// 4-Body --------------------//
  offset += comb3*sr3;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0 ;i1<n;i1++){
  for(i2=0 ;i2<n;i2++){
  for(i3=0 ;i3<n;i3++){
  for(i4=i1;i4<n;i4++){
    if((i1==i4)&&(i2>i3)){continue;}
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    for(i=0;i<sr4;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }
  }
  }
  }

// H-Bond
  offset = comb4*sr4;
  printf("%4d\n",(int)p[iq]); iq++;

  j = 0;
  for(i1=0;i1<n;i1++){
  for(i2=0;i2<n;i2++){
  for(i3=0;i3<n;i3++){
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    printf("%4d",(int)p[iq]); iq++;
    for(i=0;i<srH;i++){printf("%10.4f",p[iq]); iq++;} printf("\n");
    j++;
  }
  }
  }

}

#endif /* LMP */
