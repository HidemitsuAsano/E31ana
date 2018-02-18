#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "mngld.h"

#define GOLD   1.618034 	/* ratio for "golden sencion" */
#define GLIMIT 100.0		/* limit for a parabolic fit step */
#define TINY   1.0e-20
/* Utility routines */
#define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define SIGN(a,b) ( (b) >= 0.0 ? fabs(a) : -fabs(a) )

#define R 0.61803399		/* The golden ratio */
#define C (1.0 - R)
#define SHIFT2(a,b,c) (a)=(b);(b)=(c);

#define ITMAX 100 		/* maximum number of iteration */
#define CGOLD 0.3819680		/* golden ratio */
#define ZEPS  (1.0e-10)		/* fractional accuracy */

#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);

void nrerror( char error_text[] )
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

void mnbrak( double *ax, double *bx, double *cx,
	     double *fa, double *fb, double *fc,
	     double (*func)(double, double*), ... )
/* Given a function func, and given distinct ititial interval points ax and  */
/* bx, this routine searches in the downhill direction (defined by the       */
/* function as evaluated at the initial points) and returns new points ax,   */
/* bx, cx that bracket a minimum of the function. Also returned are the      */
/* function values at the three points, fa, fb, fc. ... is used as the       */
/* second argument of the function func.                                    */
{
  va_list list;
  va_start( list, func );
  
  double ulim,u,r,q,fu,dum;
  double *par;
  par = va_arg( list, double* );

  *fa = (*func)(*ax, par);
  *fb = (*func)(*bx, par);
  if( *fb > *fa ){		/* Switch roles of a and b so that we can go */
    SHIFT(dum,*ax,*bx,dum);	/* downhill in the direction from a to b.    */
    SHIFT(dum,*fb,*fa,dum);
  }

  *cx = (*bx) + GOLD*(*bx - *ax); /* First guess for c */
  *fc = (*func)(*cx, par);
  
  while( *fb > *fc ){		/* Keep returning here until we braket */
    /* Compute u by parabolic extrapolation from a, b, c. TINY is used to  */
    /* prevent any possible division by zero.                              */
    r = (*bx - *ax)*(*fb - *fc); 
    q = (*bx - *cx)*(*fb - *fa); 
    u = (*bx) - ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX( fabs(q-r), TINY),q-r));
    ulim = (*bx) + GLIMIT*(*cx - *bx);

    /* We won't go farther than this. Test various possibilities. */
    if( (*bx - u)*(u - *cx) > 0.0 ){ /* Parabolic u is between b and c: Try it */
      fu = (*func)(u, par);
      if( fu < *fc ){ 		/* Got a minimum between b and c. */
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	va_end( list );
	return;
      }
      else if( fu > *fb ){	/* Got a minimum between a and u. */
	*cx = u;
	*fc = fu;
	va_end( list );
	return;
      }
      
      u = (*cx) + GOLD*(*cx - *bx); /* Parabolic fit is no use.   */
      fu = (*func)(u, par);	    /* Use default magnification. */
    }      
    else if( (*cx - u)*(u - ulim) > 0.0 ){ /* Parabolic fit is between c and its allowed limit. */
      fu = (*func)(u, par);
      if( fu < *fc ){
	SHIFT(*bx,*cx,u,*cx+GOLD*(*cx- *bx));
	SHIFT(*fb,*fc,fu,(*func)(u, par));
      }
    }
    else if( (u - ulim)*(ulim - *cx) >= 0.0 ){ /* Limit parabolic u to maximum allowd value. */
      u = ulim;
      fu = (*func)(u, par);
    }
    else {			/* Reject parabolic u, Use default magnification. */
      u = (*cx) + GOLD*(*cx - *bx);
      fu = (*func)(u, par);
    }
    
    SHIFT(*ax,*bx,*cx,u);	/* Eliminate oldest point and continue. */
    SHIFT(*fa,*fb,*fc,fu);
  }
  va_end( list );
}


double golden( double ax, double bx, double cx, double (*f)(double, double*),
	       double tol, double *xmin, ... )
/* Given a function f, and given a braketing triplet of abscissas ax, bx, cx */
/* (such that bx is between ax and cx, and f(bx) is less than both f(ax) and */
/* f(cx)), this routine performs a golden section search for the minimum,    */
/* isolating it ti a fractional precision of about tol. The abscissas of the */
/* minimum is returend as xmin, and the minimum function value is returned   */
/* as golden, the returned function value. ... is used as the second         */
/* argument of the function f.                                               */
{
  va_list list;
  va_start( list, xmin );

  double f1,f2,x0,x1,x2,x3;
  double *par;
  par = va_arg( list, double* );

  x0 = ax;			/* At any given time, we will keep track of */
  x3 = cx;			/* four points, x0, x1, x2, x3.             */
  if( fabs(cx - bx) > fabs(bx - ax) ){ /* Make x0 to x1 the smaller segment. */
    x1 = bx;		               /* and fill in the new point to be    */
    x2 = bx + C*(cx - bx);             /* tried.                             */
  }
  else{
    x2 = bx;
    x1 = bx - C*(bx - ax);
  }
  f1 = (*f)(x1,par);   	/* The initial function evaluation. Note that we     */
  f2 = (*f)(x2,par);	/* never need to evaluate the function at the        */
                        /* endpoints.                                        */
  while( fabs(x3 - x0) > tol*( fabs(x1) + fabs(x2) ) ){
    if( f1 > f2 ){		      /* One possible outcome,          */
      SHIFT(x0, x1, x2, R*x1+C*x3);   /* its housekeeping,              */
      SHIFT2(f1, f2, (*f)(x2,par) );  /* and a new function evaluation. */
    }
    else{			     /* The other outcome, */
      SHIFT(x3, x2, x1, R*x2+C*x0);
      SHIFT2(f2, f1, (*f)(x1,par));  /* and its new function evaluation. */
    }
  }                             /* Back to see if we are done.             */
  if( f1 < f2 ){		/* We are done. Output the best of the two */
    *xmin = x1;			/* current values.                         */
    va_end( list );
    return f1;
  }
  else{
    *xmin = x2;
    va_end( list );
    return f2;
  }
  va_end( list );
}


double brent( double ax, double bx, double cx, double (*f)(double, double*),
	      double tol, double *xmin, ... )
/* Given a function f, and given a bracketing triplet of abscissas ax, bx,   */
/* cx (such that bx is between ax and cx, and f(bx) is less than both f(ax)  */
/* and f(cx)), this routine isolates the minimum to a fraction precision of  */
/* about tol using Brent's method. The abscissa of the minimum is returned   */
/* as xmin, and the minimum function value is returned as brent, the         */
/* returned function value. ... is used as the second argument of the        */
/* function f.                                                               */
{
  va_list list;
  va_start( list, xmin );

  int iter;
  double *par;
  double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
  double e=0.0;			/* This will be the distance moved on the    */
                                /* step before last.                         */
  par = va_arg( list, double* );
  //  d=0.0;			/* d is initialized. This statement is       */
                                /* originally not written.                   */
  a = ( ax < cx ? ax : cx );	/* a and b must be in ascending order, but   */
  b = ( ax > cx ? ax : cx );	/* the input abscissas need not be .         */
  x=w=v=bx;			/* Initialization ...                        */
  fw=fv=fx=(*f)(x,par);
  for( iter=1; iter<=ITMAX; iter++ ){ /* Main program loop. */
    xm = 0.5*(a+b);
    tol2 = 2.0*(tol1=tol*fabs(x)+ZEPS);
    if( fabs(x-xm) <= (tol2-0.5*(b-a)) ){ /* Test for done here. */
      *xmin = x;
      va_end( list );
      return fx;
    }
    if( fabs(e) > tol1 ){	/* Construct a trial parabolic fit. */
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if( q > 0.0 ) p = -p;
      q = fabs(q);
      etemp=e;
      e = d;  
      if( fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) )
	d = CGOLD*( e=( x >= xm ? a-x : b-x ) );
      /* The above conditions determine the acceptability of the parabolic   */
      /* fit. Here we take the golden section step into the two larger       */
      /* segments.                                                           */
      else{
	d = p/q;		/* Take the parabolic step. */
	u = x+d;
	if( u-a < tol2 || b-u < tol2 )
	  d = SIGN(tol1,xm-x);
      }
    }
    else{
      d = CGOLD*( e=( x >= xm ? a-x : b-x ) );
    }
    u = ( fabs(d) > tol1 ? x+d : x+SIGN(tol1,d) );
    fu = (*f)(u,par);  /* This is the one function evaluation per iteration. */
    if( fu <= fx ){   /* Now decide what to do with our function evaluation. */
      if( u >= x ) a=x; else b=x;
      SHIFT(v,w,x,u); 		/* Housekeeping follows. */
      SHIFT(fv,fw,fx,fu);
    }
    else{
      if( u < x ) a=u; else b=u;
      if( fu <= fw || w == x ){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if( fu <= fv || v == x || v == w ){
	v=u;
	fv=fu;
      }
    }
  } /* Done with housekeeping, Back for another iteration. */
  nrerror("Too many interations in brent");
  *xmin=x;
  va_end( list );
  //  return fx;
  return -999.;			/* for our applications. */
}


double dbrent( double ax, double bx, double cx,
	       double (*f)(double,double*), double (*df)(double,double*),
	       double tol, double *xmin, ... )
/* Given a function f and its derivative function df, and given a bracketing */
/* triplet of abscissas ax, bx, cx [ such a that bx is between ax and cx,    */
/* and f(bx) is less than both f(ax) and f(cx)], this routine isolates the   */
/* minimum to a fractional precision of about tol using a modification of    */
/* Brent's method that uses derivatives. The abscissa of the minimum is      */
/* returned as xmin, and the minimum function value is returned as dbrent,   */
/* the returned function value. par is used as the second argument of the    */
/* function f.                                                               */
{
  va_list list;
  va_start( list, xmin );
  
  int iter,ok1,ok2;                     /* Will be used as flags for whether */
  double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;	/* proposed steps are acceptable or  */
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm; /* not.                */
  double *par;
  par = va_arg( list, double* );

  /* Comments following will point out only difference from the routine      */
  /* brent, Read that routine first.                                         */
  a = ( ax < cx ? ax : cx );
  b = ( ax > cx ? ax : cx );
  x=w=v=bx;
  fw=fv=fx=(*f)(x,par);
  dw=dv=dx=(*df)(x,par);
  for( iter=1; iter<=ITMAX; iter++ ){ /* All our housekeeping chores are     */
    xm = 0.5*(a+b);		      /* doubled by the necessity of moving  */
    tol1 = tol*fabs(x)+ZEPS;	      /* derivative values around as well as */
    tol2 = 2.0*tol1;		      /* function values.                    */
    if( fabs(x-xm) <= (tol2-0.5*(b-a)) ){
      *xmin=x;
      va_end( list );
      return fx;
    }
    if( fabs(e) > tol1 ){
      d1 = 2.0*(b-a);		/* Initialize these d's to an out-of-bracket */
      d2=d1;			/* value.                                    */
      if( dw != dx ) d1 = (w-x)*dx/(dx-dw); /* Secant method with one point. */
      if( dv != dx ) d2 = (v-x)*dx/(dx-dv); /* And the other.                */
      /* Which of these two estimates of d shall we take?. We will insist    */
      /* that they be within the bracket, and on the side pointted to by the */
      /* derivative at x.                                                    */
      u1 = x+d1;
      u2 = x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e;			/* Movement on the step before last. */
      e=d;
      if( ok1 || ok2 ){	       /* Take only an acceptable d, and if both are */
	if( ok1 && ok2 )       /* acceptable, then take the smallest one.    */
	  d = (fabs(d1) < fabs(d2) ? d1 : d2);
	else if( ok1 )
	  d=d1;
	else
	  d=d2;
	if( fabs(d) <= fabs(0.5*olde) ){
	  u=x+d;
	  if( u-a < tol2 || b-u < tol2 )
	    d = SIGN(tol1, xm-x);
	}
	else{ 			/* Bisect, not golden section. */
	  d = 0.5*( e=(dx >= 0.0 ? a-x : b-x) );
	  /* Decide which segment by the sign of the derivative. */
	}
      }
      else{
	d = 0.5*( e=(dx >= 0.0 ? a-x : b-x) );
      }
    }
    else{
      d = 0.5*( e=(dx >= 0.0 ? a-x : b-x) );
    }
    if( fabs(d) >= tol1 ){
      u=x+d;
      fu=(*f)(u,par);
    }
    else{
      u=x+SIGN(tol1,d);
      fu=(*f)(u,par);
      if( fu > fx ){  /* If the minimum step in the downhill direction takes */
	*xmin = x;    /* us uphill, then we are done.                        */
	va_end( list );
	return fx;
      }
    }
    du=(*df)(u,par);
    if( fu <= fx ){
      if( u >= x ) a=x; else b=x;
      MOV3(v,fv,dv,w,fw,dw);
      MOV3(w,fw,dw,x,fx,dx);
      MOV3(x,fx,dx,u,fu,du);
    }
    else{
      if( u < x ) a=u; else b=u;
      if( fu <= fw || w == x ){
	MOV3(v,fv,dv,w,fw,dw);
	MOV3(w,fw,dw,u,fu,du);
      } 
      else if( fu < fv || v == x || v == w ){
	MOV3(v,fv,dv,u,fu,du);
      }
    }
  }
  nrerror("Too many iterations in routine dbrent");
  va_end( list );
  //  return 0.0; 			/* Never get here. */
  return -999.; 		/* for our applications. */
}
