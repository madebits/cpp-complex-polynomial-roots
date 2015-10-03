/* Iterative methods for complex coefficient/roots polynomials demo.
 * Vasian CEPA 1996.
 */

/* cplx.h */

#include <stdio.h>
#include <math.h>
#include <malloc.h> /* alloc.h */

/* some macros to make calculations faster:
 * PI - the number constant PI
 * dtor(n) - degree to radians converter of n
 * isodd - true if k is odd, false if k is even
 * show(c) - nicely print a complex c in screen
 */

#define PI (acos(-1.0))
#define dtor(n) ((n)*PI/180.0)
#define isodd(k)  ((k)&01)
#define iszero(d) ((d.x==d.y) ? ((d.x==0.0) ? 1:0):0)
#define show(c) printf("(%lf%s%lfj)",c.x,((c.y>=0.0) ? "+" : ""),c.y)

/* we represent the complex number as a structure of real and imaginary parts
 * optype: A - add, S - sub, M - mul, D - div, simple shorthand arithmetic operations
 * cpart: (RE - real, IM - imaginary), (MOD - modulus, ARG - argument)
 */

typedef struct{ double x,y;}cplx;
typedef enum{A,S,M,D}optype;
typedef enum{ALG,TRG,EXP}ctype;
typedef enum{RE,IM,MOD,ARG}cpart;

/* module functions:
 * simple() - carries out the optype operations on two complex numbers
 * take_part() - return the required cpart from a complex number
 * nroot() - raises a complex in an integer power
 *  (if power is negative, it is called root, hence its name)
 */

cplx simple(cplx ,cplx ,optype );
double take_part(cplx ,cpart );
cplx init(double ,double ,ctype );
double* nrootc(cplx ,int );

/* end of cplx.h */

