/****************************************************************
* 
*  potscale.c: Contains program to rescale potential for IMD or plotting.
*
*****************************************************************/

/****************************************************************
* $Revision: 1.1 $
* $Date: 2004/12/03 17:31:41 $
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
  real* force;
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
  int col1,first;
  real x0, y0, x1, y1, grad0, *xi;
  read_parameters(argc, argv);
  read_pot_table( &pair_pot, startpot, ntypes*(ntypes+1)/2 );
  if (format==3) {
    splint = splint_ed;
    splint_comb = splint_comb_ed;
    splint_grad = splint_grad_ed;
  } else { /*format == 4 ! */
    splint = splint_ne;
    splint_comb = splint_comb_ne;
    splint_grad = splint_grad_ne;
  }
  xi=pair_pot.table;
  /* init splines */
  for (col1=0; col1<paircol; col1++){  /* just pair potentials */
    first=pair_pot.first[col1];
    x0=pair_pot.begin[col1];
    x1=pair_pot.xcoord[pair_pot.first[col1]+1];
    y0=xi[first];
    y1=xi[first+1];
    /* use power law for inclination at left border */
    if (y0*y1>0)
      grad0=(y0*log(y0/y1)) /  (x0*log(x0/x1));
    else
      grad0=1e30; 		/* natural spline: curvature 0 */
    if (!((grad0>-1e10) && (grad0<1e10))) grad0=1e30;
    if (  (grad0>-1e-20)&& (grad0<1e-20)) grad0=0.;
    if (format==3) 
      spline_ed(pair_pot.step[col1], xi+first,
		pair_pot.last[col1]-first+1,
		grad0, 0.0, pair_pot.d2tab+first);
    else 			/* format == 4 ! */
      spline_ne(pair_pot.xcoord+first,xi+first,
		pair_pot.last[col1]-first+1,
		grad0, 0.0, pair_pot.d2tab+first);
  }
  
  for (col1=paircol; col1<paircol+ntypes; col1++) { /* rho */
    first=pair_pot.first[col1];
    if (format==3)
      spline_ed(pair_pot.step[col1], xi+first, 
		pair_pot.last[col1]-first+1,
		1e30,0.0,pair_pot.d2tab+first);
    else                   /* format == 4 ! */
      spline_ne(pair_pot.xcoord+first,xi+first,
		pair_pot.last[col1]-first+1,
		1e30,0.0,pair_pot.d2tab+first);		  
  }
  for  (col1=paircol+ntypes; col1<paircol+2*ntypes; col1++) { /* F */
    first=pair_pot.first[col1];
    /* Steigung  am linken Rand passt an Wurzelfunktion 
       falls 0 nicht in domain(F), sonst natürlicher Spline */
    if (format==3)
      spline_ed(pair_pot.step[col1], xi+first, 
		pair_pot.last[col1]-first+1,
#ifdef WZERO
		((pair_pot.begin[col1]<=0.) ? 1.e30
		 : .5/xi[first]),
		((pair_pot.end[col1]>=0.) ? 1e30
		 : -.5/xi[pair_pot.last[col1]]),
#else  /* WZERO: F is natural spline in any case */ 
		1.e30,1.e30,
#endif /* WZERO */
		pair_pot.d2tab+first); /* XXX */
    else                   /* format == 4 ! */
      spline_ne(pair_pot.xcoord+first,xi+first,
		pair_pot.last[col1]-first+1,
#ifdef WZERO
		(pair_pot.begin[col1]<=0.?1.e30
		 : .5/xi[first]),
		(pair_pot.end[col1]>=0.?1e30
		 : -.5/xi[pair_pot.last[col1]]),
#else  /* WZERO */
		1.e30,1.e30,
#endif /* WZERO */
		pair_pot.d2tab+first);
  }
  
  write_altplot_pair(&pair_pot, plotfile);
  /* write NEW imd potentials */
  write_pot_table_imd( &pair_pot, imdpot );
  return 0;
}
