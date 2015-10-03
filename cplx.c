/* Iterative methods for complex coefficient/roots polynomials demo.
 * Vasian CEPA 1996.
 */

/* cplx.c */

#include "cplx.h"

/* internally used function prototype */
static double kfangle(double rang);

/* returns the 'part' of a complex z noted by constant c */
double take_part(cplx z,cpart c)
{
	switch(c)
	{
	case RE:  return z.x;
	case IM:  return z.y;
	case MOD: return sqrt(z.x*z.x+z.y*z.y);
	case ARG: return (z.y!=0)? atan2(z.x,z.y) : 0.0;
	default:  return 0.0;
	}
}

/* return the n-th root of a complex z
 * the result is a double array 'result' where
 * result[0] is (equal) the module of the roots
 * result[1] - result[n] are the arguments of n-roots
 * NOTE: The user is responsible to -free- the result[] array.
 */
double *nrootc(cplx z,int root)
{
	double take_part(cplx,cpart);
	double arg=take_part(z,ARG);
	double mod=take_part(z,MOD);
	double* buffer;
	double t=2*PI;
	int i;

	if(root<0) return NULL;
	buffer=(double *)malloc((root+1)*sizeof(double));

	if(buffer==NULL) return NULL;
	else { *buffer=pow(mod,1.0/root);
    	for(i=0;i<root;i++) *(buffer+i+1)=(arg+t*i)/root;
	}
	return buffer;
}

/* for an angle greater than 2*PI returns the corresponding [0,2PI] angle */
static double kfangle(double rang)
{
	double t=2*PI;
	if(rang>0) while(rang>t) rang-=t;
	else       while(rang<-t) rang+=t;
	return rang;
}

/* constructor for a complex number z based on it representation c
 * if c = ALG then (rm,ia) are the (real, imaginary) parts of z
 * if c = TRG, EXP (rm,ia) are the (module, angle) parts of z
 * At any time we use algebraic representation of a complex inside.
 */
cplx init(double rm,double ia,ctype c)
{
	cplx t;
	switch(c)
	{
		case ALG: t.x=rm,t.y=ia; break;
		case TRG:
		case EXP: ia=kfangle(ia);
		  t.x=rm*cos(ia),t.y=rm*sin(ia); break;
		default:  t.x=t.y=0.0;
	}
	return t;
}

/* performs simple arithmetic operations
 * on two complex numbers z and v based on operation type p
 * which can be:
 * A - add, S - sub, M - mul, D - div
 */
cplx simple(cplx z,cplx v,optype p)
{
	cplx t; double m;
	switch(p)
	{
		case A: t.x=z.x+v.x,t.y=z.y+v.y; break;
		case S: t.x=z.x-v.x,t.y=z.y-v.y; break;
		case M: t.x=z.x*v.x-z.y*v.y,t.y=z.x*v.y+z.y*v.x; break;
		case D: m=v.x*v.x+v.y*v.y;
			if(m!=0){t.x=(z.x*v.x+z.y*v.y)/m,t.y=(z.y*v.x-z.x*v.y)/m;}
			else t.x=t.y=0.0;        break;
		default: t.x=t.y=0.0;
	}
	return t;
}

/* end of cplx.c */

