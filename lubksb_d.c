/**** rewritten for double precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ***/

void lubksb_d(double **a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j; /* zo */
	double sum;

	for (i=0;i<n;i++) { /* zo */
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>=0) /* zo */
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j]; /* zo */
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) { /* zo */
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j]; /* zo */
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
