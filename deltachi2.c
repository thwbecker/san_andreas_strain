#include "deltachi2.h"
/* 


   comput delta chi2 as a function of p and degrees of freedom


   numerical recipes p 697
*/



/* to find delta chi2 at 1e-6 accuracy and nu parameters of freedom
   for probability p

delta_chi2_rtbis(1e-6,nu,p)


 */
/* 

   gammq(nu/2,delta/2) = 1-p

*/
CPREC rfunc(CPREC delta, CPREC nu, CPREC p)
{
  return gammq(nu/2,delta/2)+p-1.;

}


#define JMAX 100
/* bisection */
CPREC delta_chi2_rtbis(CPREC nu,CPREC p)
{
  int j;
  CPREC dx,f,fmid,xmid,rtb,x1,x2,xacc;
  x1=0;
  x2=1e4;
  xacc = EPS;
  f=rfunc(x1,nu,p);
  fmid=rfunc(x2,nu,p);
  if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=JMAX;j++) {
    fmid=rfunc(xmid=rtb+(dx *= 0.5),nu,p);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  fprintf(stderr,"%g %g\n",nu, p);
  nrerror("Too many bisections in rtbis");
  return 0.0;
}
#undef JMAX


CPREC gammq(CPREC a,CPREC x)
{
  CPREC gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) {
    fprintf(stderr,"%g %g\n",a,x);
    nrerror("Invalid arguments in routine gammq");
  }
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}





#define ITMAX 100


void gser(CPREC *gamser,CPREC a,CPREC x,CPREC *gln)
{
  int n;
  CPREC sum,del,ap;

  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}
#undef ITMAX


#define ITMAX 100
#define FPMIN 1.0e-30

void gcf(CPREC *gammcf,CPREC a,CPREC x,CPREC *gln)
{
  int i;
  CPREC an,b,c,d,del,h;

  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef FPMIN



CPREC gammln(CPREC xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

