/****************************************************************
* 
*  potscale.c: Potfit rescaling program
*
*****************************************************************/

/****************************************************************
* $Revision: 1.5 $
* $Date: 2004/03/19 16:08:31 $
*****************************************************************/


#define MAIN

#include "potscale.h"

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
  real tot, min, max, sqr;
  int  i;
  int rescale=0.;
  pi = 4.0 * atan(1.);
  read_parameters(argc, argv);
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  if (dummy_phi!=NULL) rescale=1;
 /*   mdim=3*natoms+nconf; */
  ndim=pair_pot.idxlen;
  ndimtot=pair_pot.len;
  idx=pair_pot.idx;
  init_splines(&pair_pot);
  print_scaling(&pair_pot);

  if (rescale) {
    rescale_pot(&pair_pot, dummy_r, dummy_rho, dummy_phi);
    printf("Rescaled Potential:\n");
    print_scaling(&pair_pot);
  }
  write_pot_table( &pair_pot, endpot );
  write_pot_table_imd( &pair_pot, imdpot );
  if (plot) write_plotpot_pair(&pair_pot, plotfile);

  return 0;
}


/**************************************************************************
 *
 *  print the current scaling to the screen
 *
 *************************************************************************/

void print_scaling(pot_table_t *pt)
{
  int i,col1;
  paircol=(ntypes*(ntypes+1))/2;
  printf("Constraint on rho:\n");
  printf("rho[%f]\t = %f \t for atom type %d\n", 
	 dummy_r,
	 splint_ed(pt,pt->table,paircol+DUMMY_COL_RHO,dummy_r),
	 DUMMY_COL_RHO);
  printf("Constraints on phi:\n");
  col1=0;
  for (i=0;i<ntypes;i++) {
    printf("phi[%f]\t = %f \t between atoms of type %d\n",
	   dummy_r,splint_ed(pt,pt->table,col1,dummy_r),i);
    col1+=ntypes-i;
  }
  printf("Constraints on F:\n");
  col1=paircol+ntypes;
  for (i=0;i<ntypes;i++) {
#ifdef PARABEL
    printf("F[%f]\t = %f \t for atoms of type %d\n",
	   0.,parab_ed(pt,pt->table,col1,0.),i);
#else
    printf("F[%f]\t = %f \t for atoms of type %d\n",
	   0.,splint_ed(pt,pt->table,col1,0.),i);
#endif
    col1++;
  }
  
/*   printf("FYI: \n"); */
/*   for (i=0;i<ntypes;i++) { */
/*     printf("rho[%f]\t = %f \t for atom type %d\n", */
/* 	   dummy_r, */
/* 	   splint_ed(pt,pt->table,paircol+i,dummy_r), */
/* 	   i); */
/*   } */

  printf("\n");
  return;
}


void init_splines(pot_table_t *pt)
{
  int col1, first;
  int paircol=(ntypes*(ntypes+1))/2;
  real grad0, y0,y1,x0,x1;
 /* init second derivatives for splines */
  for (col1=0; col1<paircol; col1++){  /* just pair potentials */
    first=pt->first[col1];
    x0=pt->begin[col1];
    x1=x0+pt->step[col1];
    y0=pt->table[first];
    y1=pt->table[first+1];
    if (y0*y1>0)
      grad0=(y0*log(y0/y1)) /  (x0*log(x0/x1));
    else
      grad0=1e30;
    if (!((grad0>-1e10) && (grad0<1e10))) grad0=1e30;
    if (  (grad0>-1e-20)&& (grad0<1e-20)) grad0=0.;

    spline_ed(pt->step[col1], pt->table+first,
	      pt->last[col1]-first+1,
	      grad0, 0.0, pt->d2tab+first);
  }
  for (col1=paircol; col1<paircol+ntypes; col1++) { /* rho */
    first=pair_pot.first[col1];
    spline_ed(pt->step[col1], pt->table+first, 
	      pt->last[col1]-first+1,
	      1e30,0.0,pt->d2tab+first);
  }

#ifndef PARABEL
  for  (col1=paircol+ntypes; col1<paircol+2*ntypes; col1++) { /* F */
    first=pt->first[col1];
    spline_ed(pt->step[col1], pt->table+first, 
	      pt->last[col1]-first+1,
	      0.,1.e30,pt->d2tab+first);
  }

  return;
#endif
}

/****************************************************************************
 *
 * This one does the rescaling...
 *
 ***************************************************************************/

void rescale_pot(pot_table_t *pt,real dummy_r, real dummy_rho, real *dummy_phi)
{
  real a,r,step,temp;
  int i,j,k,first,col1;
  int paircol=(ntypes*(ntypes+1))/2;
  real *b;
  b=(real *) malloc( ntypes * sizeof(real)); 
  /* Calculate value of rho at dummy_r, compare with dummy_rho */
  a=dummy_rho/splint_ed(pt,pt->table,paircol+DUMMY_COL_RHO,dummy_r);
  /*rescale all rho by a */
  for (i=paircol;i<paircol+ntypes;i++){
/*    first =pt->first[i];*/
    for (j=pt->first[i];j<=pt->last[i];j++){
      pt->table[j]*=a;
    }
  }
  /* rescale all embed. by a */
  for (i=paircol+ntypes;i<paircol+2*ntypes;i++){
    pt->begin[i]*=a;
    pt->end[i]*=a;
    pt->step[i]*=a;
    pt->invstep[i]/=a;
  }
  printf("Gauge constant a:\t%f\n",a);
  init_splines(pt);
  col1=0;
  for (i=0;i<ntypes;i++){
    b[i]=(dummy_phi[i]-splint_ed(pt,pt->table,col1,dummy_r))/
      (2*splint_ed(pt,pt->table,paircol+i,dummy_r));
    col1+=ntypes-i;
    printf("Gauge constant b[%d]:\t%f\n",i,b[i]);    
  }
  col1=0;
  for (i=0;i<ntypes;i++){
    for(j=i;j<ntypes;j++) {
      r=pt->begin[col1];
      step=pt->step[col1];
      for (k=pt->first[col1];k<pt->last[col1];k++){
	if (r<pt->end[paircol+i])
	  pt->table[k]+=b[j]*splint_ed(pt,pt->table,paircol+i,r);
	if (r<pt->end[paircol+j])
	  pt->table[k]+= b[i]*splint_ed(pt,pt->table,paircol+j,r);
	r+=step;
      }
      col1++;
    }
  }
  for (i=paircol+ntypes;i<paircol+2*ntypes;i++){
    r=pt->begin[i];
    step=pt->step[i];
    for (k=pt->first[i];k<=pt->last[i];k++){
      pt->table[k]-=b[i-(paircol+ntypes)]*r;
      r+=step;
    }
  }
  init_splines(pt);
  /* translate F(0) to 0 */
  col1=paircol+ntypes;
  for (i=0;i<ntypes;i++) {
#ifdef PARABEL
    temp= parab_ed(pt,pt->table,col1,0.);
#else
    temp= splint_ed(pt,pt->table,col1,0.);
#endif
    for (k=pt->first[col1];k<=pt->last[col1];k++){
      pt->table[k]-=temp;
    }
    col1++;
  }
  init_splines(pt);
  return;
}
