#include <math.h>
#include "potfit.h"
#include "nrutil_r.h"
#define RAND_MAX 2147483647
#define EPS 0.1
#define NEPS 4
#define NSTEP 20
#define STEPVAR 2.0
#define NTEMP (5*ndim)
#define TEMPVAR 0.85
/* #define TSTART 6 */
#define KMAX 1000
#define GAUSS(a) (1.0/sqrt(2*pi)*(exp(-(DSQR(a))/2.)))

/*int RAND_MAX=(1<<31)-1;*/


real normdist(void)
{
    static int have = 0;
    static real nd2;
    real x1,x2,sqr,cnst;
    if (!(have)) {
	do {
	    /* x1=2.0*drand48()-1.0;
	       x2=2.0*drand48()-1.0;*/
	    x1=2.0*random()/(RAND_MAX+1.0)-1.0;
	    x2=2.0*random()/(RAND_MAX+1.0)-1.0;
	    sqr=x1*x1+x2*x2;
	} while (!(sqr <= 1.0 && sqr>0));
	/* Box Muller Transformation */
	cnst=sqrt(-2.0*log(sqr)/sqr);
	nd2 = x2 * cnst;
	have = 1;
	return x1* cnst;}
    else {
	have = 0;
	return nd2;
    }
}

void makebump(real *x,real width,real height,int center)
{
    int i,j=0;
    /* find pot to which center belongs */
    while (pair_pot.last[j]<center) j++;
    for (i=0;i<=4.*width;i++){
	if (center+i<pair_pot.last[j]){
	    x[center+i]+=GAUSS((double) i / width)*height;
	}
    }
    for (i=1;i<=4.*width;i++){
	if (center-i>=pair_pot.first[j]){
	    x[center-i]+=GAUSS((double) i / width) * height;
	}
    }
    return;
}


void anneal(real *xi)
{
    int i=0,j=0,m=0,k=0,n,h=0;			/* counters */
    int nstep=NSTEP,ntemp=NTEMP;
    int loopagain;		/* loop flag */
    real c=STEPVAR;
    real p;			/* Probability */
    real T;	   	        /* Temperature */
    real F, Fopt, F2;		/* Fn value */
    real *Fvar;			/* backlog of Fn vals */
    real *v; 			/* step vector */
    real *xopt,*xi2;		/* optimal value */
    real *fxi1;			/* two latest force vectors */
    real width, height;		/* gaussian bump size */
    FILE *ff;			/* exit flagfile */
    int *naccept;		/* number of accepted changes in dir */
    /* init starting temperature for annealing process */
    T=anneal_temp;
    if (T==0.) return; 		/* don't anneal if starttemp equal zero */
    Fvar=dvector(-NEPS+1,KMAX+5);
    v=dvector(0,ndim-1);
    xopt=dvector(0,ndimtot-1);
    xi2=dvector(0,ndimtot-1);
    fxi1=dvector(0,mdim-1);
    naccept=ivector(0,ndim-1);
    /* init step vector and optimum vector */
    for (n=0;n<ndim;n++) {
	v[n]=1.;
	naccept[n]=0;
    }
    for (n=0;n<ndimtot;n++) {xi2[n]=xi[n];xopt[n]=xi[n];}
    F=(*calc_forces)(xi,fxi1);
    Fopt=F;
    printf("k\tT        \tm\tF          \tFopt\n");   
    printf("%d\t%f\t%d\t%f\t%f\n", 0, T,0, F, Fopt);
    for (n=0;n<NEPS;n++) Fvar[-n]=F;
    do {
	for (m=0;m<ntemp;m++){
	    for (j=0;j<nstep;j++) {
		for (h=0;h<ndim;h++){
		    /* Step #1 */
		    /* Create a gaussian bump, 
		       width & hight distributed normally */
		    
		    width=fabs(normdist()); height=normdist()*v[h];
		    for (n=0;n<ndimtot;n++) xi2[n]=xi[n];
		    makebump(xi2,width,height,idx[h]);
		    F2=(*calc_forces)(xi2,fxi1);
		    if (F2<=F) {		/* accept new point */
			for (n=0;n<ndimtot;n++) xi[n]=xi2[n];
			F=F2;
			i++;naccept[h]++;
			if(F<Fopt) {
			    for (n=0;n<ndimtot;n++) xopt[n]=xi2[n];
			    Fopt=F2;
			}
		    }
		    else if (random()/(RAND_MAX+1.0)<exp((F-F2)/T)) {
			for (n=0;n<ndimtot;n++) xi[n]=xi2[n];
			F=F2;
			i++;naccept[h]++;
		    }
		}
	    }
	    /* Step adjustment */
	    for (n=0;n<ndim;n++) {
		if (naccept[n]>0.6*nstep) 
		    v[n]*=(1+c*((double) naccept[n]/nstep - 0.6)/0.4);
		else if (naccept[n]<0.4*nstep)
		    v[n]/=(1+c*(0.4-(double) naccept[n]/nstep)/0.4);
		naccept[n]=0;
	    }
	    printf("%d\t%f\t%d\t%f\t%f\n", k, T,m+1, F, Fopt);
	    /* End fit if break flagfile exists */
	    ff = fopen(flagfile,"r");
	    if (NULL != ff) {
	      printf("Annealing terminated in presence of break flagfile %s!\n", flagfile);
	      printf("Temperature was %f, returning optimum configuration\n", T);	    
	      for (n=0;n<ndimtot;n++) xi[n]=xopt[n];
	      F=Fopt;
	      k=KMAX+1;
	      break;
	    }
	
	} 
	/*Temp adjustment */
	T*=TEMPVAR;
	k++;
	Fvar[k]=F;
	loopagain = 0; 
	for (n=1;n<=NEPS;n++) { 
	  /* printf("%d %f %f",n,fabs(F-Fvar[k-n]),EPS); */
	  if (fabs(F-Fvar[k-n])>EPS) loopagain=1; }
	/* if ((F-Fopt)>EPS) loopagain=1; */
	if (loopagain) {
	    i++;
	}
	else if ((F-Fopt)>EPS){
	    for (n=0;n<ndimtot;n++) xi[n]=xopt[n];
	    F=Fopt;
	    loopagain=1;
	}
    } while (k<KMAX && loopagain);
    free_dvector(Fvar,-NEPS+1,KMAX+5);
    free_dvector(v,0,ndim-1);
    free_dvector(xopt,0,ndimtot-1);
    free_dvector(xi2,0,ndimtot-1);
    free_dvector(fxi1,0,mdim-1);
    free_ivector(naccept,0,ndim-1);
    return;
}
/*int main(void)
{
    double sum = 0;
    int hist[BINS+2];
    int i,j;
    for (i=0;i<BINS+2;i++) hist[i]=0;
    srandom(12345);
    drand48();drand48();drand48();
    random();random();random();
    for (i=0;i<1e8;i++) {
	j=(int) (normdist()*BINS/(2*5.0) + BINS*0.5);
	if (j<0) j=0;
	if (j>=BINS) j=BINS-1;
	hist[j]++;
    }
    for (i=0;i<BINS+2;i++) printf("%f %f %f\n", (2*5.0)/(double)BINS*((i+.5)-BINS*.5),(double)hist[i]/1.e8, sum+=(double) hist[i]/1.e8);
    return 1;
}*/
