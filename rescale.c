/****************************************************************
* 
* rescale.c: Routines used to automatically rescale 
*     EAM potential.
*
*****************************************************************/
/****************************************************************
* $Revision: 1.3 $
* $Date: 2004/08/16 13:02:50 $
*****************************************************************/

#include "potfit.h"

/* Doesn't make much sense without EAM  */

#ifdef EAM


/****************************************************************
* 
* rescale: Routine used to automatically rescale 
*     EAM potential. Flag indicates whether to force update...
*     upper is upper limit of electron density.
*
*****************************************************************/

real rescale(pot_table_t *pt, real upper, int flag)
{
  int mincol,maxcol,col,col2,first,vals,h,i,j,typ1,typ2,sign,dimneuxi;
  real *xi,*neuxi,*neuord,*neustep,*maxrho,*minrho,*left,*right;
  atom_t *atom;
  neigh_t *neigh;
  real fnval,pos,grad,a;
  real min=1e100,max=-1e100;

  xi=pt->table;
  dimneuxi=pt->last[paircol+2*ntypes-1]-pt->last[paircol+ntypes-1];
  neuxi=(real *) malloc(dimneuxi*sizeof(real));
  neuord=(real *) malloc(dimneuxi*sizeof(real));
  neustep=(real *) malloc(ntypes*sizeof(real));
  maxrho=(real *) malloc(ntypes*sizeof(real));
  minrho=(real *) malloc(ntypes*sizeof(real));
  left=(real *) malloc(ntypes*sizeof(real));
  right=(real *) malloc(ntypes*sizeof(real));
  for (i=0;i<ntypes;i++){
    maxrho[i]=-1e100;
    minrho[i]=1e100;
  }
  /* Max/Min rho finden */
  /* Splines initialisieren - sicher ist sicher*/
  for (col=paircol; col<paircol+ntypes; col++) { /* rho */
    first=pt->first[col];
    if (format == 3)
      spline_ed(pt->step[col], xi+first, 
		pt->last[col]-first+1,
		1e30,0.0,pt->d2tab+first);
    else			/* format == 4 ! */
      spline_ne(pt->xcoord+first, xi+first, pt->last[col]-first+1,
		1e30,0.0,pt->d2tab+first);
  }
  for  (col=paircol+ntypes; col<paircol+2*ntypes; col++) { /* F */
    first=pt->first[col];
    /* Steigung 0 am rechten Rand */
    if (format == 3) 
      spline_ed(pt->step[col], xi+first, 
		pt->last[col]-first+1,
		0,1e30,pt->d2tab+first);
    else 			/* format == 4 */
      spline_ne(pt->xcoord+first, xi+first, 
		pt->last[col]-first+1,
		0,1e30,pt->d2tab+first);
  }
  
  /* atom_rho berechnen (eigentlich Verschwendung...) */
  for (h=0; h<nconf; h++) {
    for (i=0; i<inconf[h]; i++) 
      atoms[cnfstart[h]+i].rho=0.0;
    for (i=0; i<inconf[h]; i++) {
      atom = atoms + i + cnfstart[h];
      typ1 = atom->typ;
      for (j=0; j<atom->n_neigh; j++) {
	neigh = atom->neigh+j;
	if (neigh->nr > i+cnfstart[h]) {
	  typ2  = neigh->typ;
	  col2 = paircol+typ2;
	  if (typ2==typ1) {
	    if (neigh->r < pt->end[col2]) {
	      fnval = splint_dir(pt,xi,col2,neigh->slot[1],
				 neigh->shift[1],neigh->step[1]);
	      atom->rho += fnval;
	      atoms[neigh->nr].rho += fnval;
	    }
	  } else {
	    col = paircol+typ1;
	    if (neigh->r < pt->end[col2]) {
	      atom->rho += splint_dir(pt,xi,col2,neigh->slot[1],
				      neigh->shift[1],neigh->step[1]);
	    }
	    if (neigh->r < pt->end[col])
	      atoms[neigh->nr].rho += 
		splint(pt,xi,col,neigh->r);
	  }
	}
      }
      maxrho[typ1]=MAX(maxrho[typ1],atom->rho);
      minrho[typ1]=MIN(minrho[typ1],atom->rho);
    }
  }
  for (i=0;i<ntypes;i++) {
    //printf("maxrho[%d]=%f\tminrho[%d]=%f\n",i,maxrho[i],i,minrho[i]);
    if (maxrho[i]> max) { max=maxrho[i];maxcol=i;}
    if (minrho[i]< min) { min=minrho[i];mincol=i;}
  }
  /* dominante Seite bestimmen */
  sign = ( max>=-min ) ? 1 : -1;
  
  /* Neue linke und rechte Grenze ermitteln, 60 Prozent zugeben... */
  
  for (i=0;i<ntypes;i++) {
    j=paircol+ntypes+i;
    left[i]=minrho[i]-0.6*pt->step[j];
    right[i]=maxrho[i]+0.6*pt->step[j];
    /* Ist Erweiterung nötig? */
    if ( flag ||
	 minrho[i]-pt->begin[j] < 0. ||
	 minrho[i]-pt->begin[j] > .95*pt->step[j] ||
	 maxrho[i]-pt->end[j]   > 0 ||
	 maxrho[i]-pt->end[j]   < -.95*pt->step[j])
      flag=1;
  }
  
  /* Skalierungsfaktor bestimmen */
  a= (sign==1) ? upper/right[maxcol] : upper/left[mincol];
  
  if (flag || fabs(a) > 1.05 || fabs(a) < 0.95) flag=1;
  
  /* Haben wir Grund zum Updaten? */
  
  if (!flag) return 0.;		/* offensichtlich nicht */

  /* Also machen wir ein Update... */

  /* Null muss dabei sein ???  -- nein! */
/*   if (sign==1) { */
/*     for (i=0;i<ntypes;i++) */
/*       left[i]=MIN(0.,left[i]); */
/*   } else { */
/*     for (i=0;i<ntypes;i++) */
/*       right[i]=MAX(0.,right[i]); */
/*   } */

  /* Potenzial erweitern */
  h=0;
  for(i=0;i<ntypes;i++) {
    col=paircol+ntypes+i;
    vals=pt->last[col]-pt->first[col];
    neustep[i]=(right[i]-left[i])/(double) vals;
    pos=left[i];
    for (j=0;j<=vals;j++) {
      if (pt->begin[col]>pos) { /* extrapolating... */
	grad=splint_grad(pt,xi,col,pt->begin[col]); /* Steigung */
	neuxi[h]=grad*(pt->begin[col]-pos)+xi[pt->first[col]];
      } else if (pt->end[col]<pos) { /* also extrapolating */
	grad=splint_grad(pt,xi,col,pt->end[col]); /* Steigung */
	neuxi[h]=grad*(pos-pt->end[col])+xi[pt->last[col]];
      } else 			/* direct calculation */
	neuxi[h]=splint(pt,xi,col,pos);
      
      neuord[h]=pos;
      h++;
      pos+=neustep[i];
    }
  }
  
  /* Werte zurückschreiben */
  col=pt->first[paircol+ntypes]; /* erster zu ändernder Wert */
  for (i=0;i<dimneuxi;i++){
    xi[i+col]=neuxi[i];
    pt->xcoord[i+col]=neuord[i];
  }
//#ifndef DEBUG
  printf("Skalierungsfaktor %f\n",a);
//#endif
  /* skalieren */
  for (i=paircol;i<paircol+ntypes;i++){
/*    first =pt->first[i];*/
    for (j=pt->first[i];j<=pt->last[i];j++){
      pt->table[j]*=a;
    }
  }
  
  /* rescale all embed. by a */
  if (sign==1){
    j=0;
    for (i=paircol+ntypes;i<paircol+2*ntypes;i++){
      pt->begin[i]=a*left[j];
      pt->end[i]=a*right[j];
      pt->step[i]=a*neustep[j];
      pt->invstep[i]=1.0/pt->step[i];
      j++;
    }
  } else { 			/* wir drehen um - a negativ */
    j=0;
    for (i=paircol+ntypes;i<paircol+2*ntypes;i++){
      pt->begin[i]=a*right[j];
      pt->end[i]=a*left[j];
      pt->step[i]=-a*neustep[j];
      pt->invstep[i]=1.0/pt->step[i];
      j++;
    }
    h=0;
    for (i=0;i<ntypes;i++) { 	/* Werte in umgekehrter Reihenfolge */
      col=paircol+ntypes+i;
      for (j=pt->last[col];j>=pt->first[col];j--){
	neuxi[h]=xi[j];
	neuord[h]=pt->xcoord[j];
	h++;
      }
    }
    col=pt->first[paircol+ntypes];	/* und wieder zurückschreiben */
    for (i=0;i<dimneuxi;i++){
      xi[i+col]=neuxi[i];
      pt->xcoord[i+col]=neuord[h];
    }
  }
  free(neuxi);
  free(neustep);
  free(maxrho);
  free(minrho);
  free(left);
  free(right);
    

  /* Faktor zurückgeben */
  return a;
}

/**********************************************************************
 *
 * embed_shift: Shift embedding function to U(0)=0
 *                   
 **********************************************************************/

void embed_shift(pot_table_t *pt){
  real shift;
  real *xi;
  int i,j,first;
  xi=pt->table;
  for (i=paircol+ntypes;i<paircol+2*ntypes;i++){
    first=pt->first[i];
    /* init splines - sicher ist sicher */
    /* Steigung 0 am linken Rand */
/********* GEFAHR HIER ****************/
/** NICHT IDIOTENSICHER ***************/
    if (pt->begin[i]<=0) {	/* 0 in domain of U(n) */
      if (format == 3) 
	spline_ed(pt->step[i], xi+first,pt->last[i]-first+1,0.,
		  1e30,pt->d2tab+first);
      else 			/* format == 4 ! */
	spline_ne(pt->xcoord+first,xi+first,pt->last[i]-first+1,
		  0.,1e30,pt->d2tab+first);
      shift=splint(pt,xi,i,0.);
#ifdef DEBUG
      printf("shifting by %f\n",shift);
#endif /* DEBUG */
    } else
      shift=xi[first];
    for (j=first;j<=pt->last[i];j++)
      xi[j]-=shift;
  }
}
#endif /* EAM */
