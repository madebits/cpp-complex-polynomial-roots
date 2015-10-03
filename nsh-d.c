/* Iterative methods for complex coefficient/roots polynomials demo.
 * Vasian CEPA 1996.
 */

/* nsh-d.c */

#include "nsh.h"

#define maxiter 50 /* no more than that iterations allowed */
#define eps 0.0001 /* allowed error margin in approximating the result */

/* helper function to print a polynomial in screen */
void printPoly(cplx a[], int n)
{
	char *f = "x^%d + ";
	int i = n;
	printf("\nYou are using this polynomial for testing:\n");
	for(;i > -1; i--){
		if(i == 0) f = "x^%d = 0";
		show(a[i]); printf(f, i);
	}
	printf("\nThis is a polynomial of order: %d", n);
}

/* some test code for the library */
int main()
{
	/* We will use this polynomial for demonstration:
	 * 2*x^2 - 8 = 0 this is the function we are using for demo
	 * It is a polynomial of order n=2
	 *
	 * Note the first coefficient goes last. This is a matter of implementation,
	 * not of the methods used.
	 */

	cplx a[]={{-8.0,0.0},{0.0,0.0},{2.0,0.0}};

	/* cplx a[]={{0.0,0.0},{1.0,0.0},{1.0,0.0},{1.0,0.0}}; */
	/*
	 * x^3 + x^2 + x = 0 - alternative test of complex roots
	 * It is a polynomial of order n=3
	 */

	int i = 0;
	int n = sizeof(a)/sizeof(a[0]) - 1; /* this the order of polynomial */
	cplx t = init(-0.4,-0.2,ALG); /* a starting point, somewhere near the root */
	cplx g, z[2]; /* some variables for holding results */

	/* nice-print the polynomial */
	printPoly(a,n); 

	/* newton - finds a approximate root in the complex plane */

	i = maxiter;
	g = npr(a,n,t,eps,&i);
	printf("\n\nnewt root="); show(g);
	printf("\niterations=%d",i);

	/* secant - another algorithm to find the root in the complex plane */

	i = maxiter;
	g = secf(a,n,t,eps,&i);
	printf("\n\nsec root="); show(g);
	printf("\niterations=%d",i);

	/* horner - finds a polynomial value and its first derivate in a given point */

	hornerfd(a,n,init(5.0,0.0,ALG),z);
	printf("\n\nhorner calculation:");
	printf("\nf((5.0,0.0)) = "); show(z[0]);
	printf("\nf'((5.0,0.0)) = "); show(z[1]);

	return 0;
}

/* end of nsh-d.c */

