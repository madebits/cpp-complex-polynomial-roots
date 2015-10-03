/* Iterative methods for complex coefficients/roots polynomials demo.
 * Vasian CEPA 1996.
 */

/* nsh.c */

#include "nsh.h"

/* calculates f(z) and its derivate f'(z) for
 * a polynomial f with coefficients in a[] and order n
 * using Horner's factorization method
 * the result is returned in f[] where
 * f[0] is f(z) and
 * f[1] is f'(z)
 */
void hornerfd(cplx a[],int n,cplx z,cplx f[])
{
	cplx b1,c1,b0=a[n];
	cplx c0=b0;

	int i=n-1;
	for(;i>=1;i--)
	{
		b1=simple(a[i],simple(z,b0,M),A);
		c1=simple(b1,simple(z,c0,M),A);
		b0=b1;
		c0=c1;
	}

	f[0]=simple(a[0],simple(z,b0,M),A);   /*f(z) */
	f[1]=c1;                              /*f'(z)*/
	return;
}

/* calculates f(z) for a polynomial f denoted by
 * its coefficients in a[]
 * a light-weight version of hornerfd
 */
cplx horner(cplx a[],int n,cplx z)
{
	cplx b1,b0=a[n];
	int i=n-1;
	for(;i>=0;i--)
	{
		b1=simple(a[i],simple(z,b0,M),A);
		b0=b1;
	}
	return b0;
}

/* Newton iterative method for finding a root of an polynomial equation
 * a[] - the coefficients of the polynomial
 * n - the polynomial order (size of a[])
 * e - epsilon: the acceptable margin of error
 * count - the max number of iteration allowed, in end contains how iterations were used
 */
cplx npr(cplx a[],int n,cplx r,double e,int *count)
{
	int icount=0;
	cplx f[2],z1;
	cplx z0=r;
	cplx dl=init(e+100.0,0.0,ALG); /* make sure we do not exit in beginning */

	while(fabs(take_part(dl,MOD))>e && icount<=*count)
	{
		hornerfd(a,n,z0,f);
		dl=simple(f[0],f[1],D);
		z1=simple(z0,dl,S);
		z0=z1;
		if(iszero(horner(a,n,z1))) break;
		icount++;
	}
	*count= (icount>=*count) ? -1:icount;
	return z1;
}

/* Secant iterative method
 * See npr() for argument explanation.
 */
cplx secf(cplx a[],int n,cplx r,double e,int *count)
{
	int icount=0;
	cplx z0=r;
	cplx dl=init(-0.1,0.0,ALG);
	cplx z1=simple(z0,dl,S);
	cplx f0=horner(a,n,z0);
	cplx t,f1,z2;

	while(fabs(take_part(dl,MOD))>e && icount<=*count)
	{
		f1=horner(a,n,z1);

		if(fabs(take_part(f1,MOD))>fabs(take_part(f0,MOD)))
		{
			t=z1,z1=z0,z0=t;
			t=f1,f1=f0,f0=t;
			dl=simple(dl,init(-1.0,0.0,ALG),M);
		}

		t=simple(f1,f0,D);
		dl=simple(simple(dl,t,M),simple(init(1.0,0.0,ALG),t,S),D);
		z2=simple(z1,dl,S);
		z0=z1;
		z1=z2;
		f0=f1;

		if(iszero(horner(a,n,z1))) break;
		icount++;
	}

	*count= (icount>=*count) ? -1:icount;
	return z2;
}

/* end of nsh.c */

