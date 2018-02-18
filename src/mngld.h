#ifndef __MNGLD__
#define __MNGLD__

#ifdef __cplusplus
extern "C" {
#endif

/* Routine from Numerical Recipes in C */
void  nrerror(char error_text[]); /* from nrutil.h */

//---------------------------------------------------------------------//
  
/* Routine for initially bracketing a minimum */
void   mnbrak(double *ax, double *bx, double *cx,
	      double *fa, double *fb, double *fc, 
	      double (*func)(double, double*), ... );
/* Routine for Golden Section Search */
double golden(double ax, double bx, double cx,
	      double (*f)(double,double*), double tol, double *xmin, ... );
/* Routine for Brent's method */
double  brent(double ax, double bx, double cx, double (*f)(double, double*),
	      double tol, double *xmin, ... );
/* Routine for a one-demensional search with first derivatives */
double dbrent(double ax, double bx, double cx, double (*f)(double, double*),
	      double (*df)(double, double*), double tol, double *xmin, ...);
#ifdef __cplusplus
}
#endif

#endif

