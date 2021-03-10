#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef SINGLE_PREC
#define FMT_LINE1 "%f %f %f %f %f %f %f %s"
#define FMT_1 "%f"
#define EPS 1.2e-7
#define CPREC float
#else
#define FMT_LINE1 "%lf %lf %lf %lf %lf %lf %lf %s"
#define FMT_1 "%lf"
#define EPS 4e-15
#define CPREC double
#endif
#define CLEN 200
#define boolean unsigned short
#define TRUE 1
#define FALSE 0
#define MEMERROR(x) {fprintf(stderr,"%s: out of memory\n",x);exit(-1);}
#define RADIUS_KM 6371.
#define TWOPI 6.28318530717958647
#define ONEEIGHTYOVERPI   57.295779513082321 
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.1415926535897932384626
#endif
#define PIOVERONEEIGHTY 0.0174532925199433
#define RAD2DEG(x) ( (x) * ONEEIGHTYOVERPI )
#define DEG2RAD(x) ( (x) / ONEEIGHTYOVERPI )
#define LON2PHI(x) (DEG2RAD(x))
#define LAT2THETA(x) (DEG2RAD(90.0-x))
#define PHI2LON(x) ( RAD2DEG(x) )
#define THETA2LAT(x) ( (90.0-RAD2DEG(x)))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define RAD2KM(x) (((x)*RADIUS_KM))
#define KM2RAD(x) (((x)/RADIUS_KM))



CPREC delta_dist(CPREC , CPREC , CPREC , CPREC );

struct data{
  CPREC lon,lat,phi,theta;
  CPREC v[2],sv[2],err;
  CPREC corr;
  char code[100];

};

/* 

   read in GPS data and condense colocated stations

*/
int main(int argc, char **argv)
{

  /* max velocity error bound */
  CPREC max_err = 12;
  /* max velocity bound */
  CPREC max_vel = 150;
  /* distance */
  CPREC eps_dist_km = 2;
  /* geographic range */
  const CPREC xrange[2] = {0,360}, yrange[2]={-90,90};
  /* 

   */
  FILE *in;
  int i,ndata,nread;
  char datafile[CLEN];
  boolean use;
  struct data *d,tmpd;
  CPREC min_dist_rad,dist,dist_min;
  
  if(argc < 2){
    fprintf(stderr,"%s datafile [eps_km, %g]\n\n",argv[0],eps_dist_km);
    fprintf(stderr,"\twill select data with vel <= %g, errors <= %g, and prune them for colocations of eps_km\n",
	    max_vel, max_err);
    fprintf(stderr,"\tby selecting data with smaller error, geographic range -R%g/%g/%g/%g\n\n",
	    xrange[0],xrange[1],yrange[0],yrange[1]);
    exit(-1);
  }
  /* options */
  sprintf(datafile,"%s",argv[1]); /* filename */
  if(argc>2)			  /* different range */
    sscanf(argv[2],FMT_1,&eps_dist_km);

  
  min_dist_rad = KM2RAD(eps_dist_km);
  /* 
     
     read in data
     
  */
  in = fopen(datafile,"r");
  if(!in){fprintf(stderr,"%s: error opening %s\n",argv[0],datafile);exit(-1);}
  d = (struct data *)malloc(sizeof(struct data));
  
  ndata = nread = 0;
  while(fscanf(in,FMT_LINE1,&(tmpd.lon),&(tmpd.lat),&(tmpd.v[0]),&(tmpd.v[1]),
	       &(tmpd.sv[0]),&(tmpd.sv[1]),&(tmpd.corr),(tmpd.code)) == 8){
    /* process this datum */
    while(tmpd.lon < 0)
      tmpd.lon += 360;
    while(tmpd.lon > 360)
      tmpd.lon -= 360;
    
    tmpd.phi = LON2PHI(tmpd.lon);
    tmpd.theta = LAT2THETA(tmpd.lat);
    tmpd.err = hypot(tmpd.sv[0],tmpd.sv[1]);

    if((hypot(tmpd.v[0],tmpd.v[1]) <= max_vel) && (tmpd.err <= max_err) && ( tmpd.lon <= xrange[1]) &&
       (tmpd.lon >= xrange[0]) && (tmpd.lat <= yrange[1]) && (tmpd.lat >= yrange[0])){
      use = TRUE;
      dist_min = 1e20;
      for(i=0;i < ndata;i++){	/* check for existing data */
	if((dist=delta_dist(tmpd.theta, tmpd.phi,d[i].theta, d[i].phi)) <= min_dist_rad){
	  /* 
	     close to previous location
	  */
	  if(tmpd.err < d[i].err){	       /* overwrite previous
						  if this datum
						  appears better
						  constrained */
	    fprintf(stderr,"%s: colocated: replacing old x: %8.3f,%8.3f v: %8.3f,%8.3f (s: %5.3f), with         ",
		    argv[0],d[i].lon,d[i].lat,d[i].v[0],d[i].v[1],d[i].err);
	    /* don't add this one */
	    fprintf(stderr,"x: %8.3f,%8.3f v: %8.3f,%8.3f (s: %5.3f)\n",
		    tmpd.lon,tmpd.lat,tmpd.v[0],tmpd.v[1],tmpd.err);
	    memcpy((d+i),&tmpd,sizeof(struct data));
	  }else{
	    fprintf(stderr,"%s: colocated: keeping   old x: %8.3f,%8.3f v: %8.3f,%8.3f (s: %5.3f), disregarding ",
		    argv[0],d[i].lon,d[i].lat,d[i].v[0],d[i].v[1],d[i].err);
	    /* don't add this one */
	    fprintf(stderr,"x: %8.3f,%8.3f v: %8.3f,%8.3f (s: %5.3f)\n",
		    tmpd.lon,tmpd.lat,tmpd.v[0],tmpd.v[1],tmpd.err);
	  }
	  use = FALSE;		/* either way, don't use */
	  break;		/* exit search loop */
	}
	if(dist<dist_min)dist_min=dist;
      }
      fprintf(stderr,"%s: datum %4i (%7.2f, %7.2f) was at least %6.1f km from any other station\n",
	      argv[0],ndata+1,tmpd.lon,tmpd.lat,RAD2KM(dist_min));
      if(use){
	memcpy((d+ndata),&tmpd,sizeof(struct data));
	ndata++;
	d = (struct data *)realloc(d,sizeof(struct data)*(ndata+1));
      }
    }else{
      fprintf(stderr,"%s: disregarding:            x: %8.3f,%8.3f v: %8.3f,%8.3f (s: %5.3f)\n",
	      argv[0],tmpd.lon,tmpd.lat,tmpd.v[0],tmpd.v[1],tmpd.err);

    }
    nread++;
  }
  fclose(in);
  /* 
     throw out removed data and compute min/max 
  */
  fprintf(stderr,"%s: selected %i data out of %i from %s, eps %g km\n",argv[0],ndata,nread,datafile, RAD2KM(min_dist_rad));
  for(i=0;i<ndata;i++)
    printf("%10.6f %10.6f\t%10.4f %10.4f\t%6.2f %6.2f\t%6.3f\t%s\n",
	   d[i].lon,d[i].lat,d[i].v[0],d[i].v[1],d[i].sv[0],d[i].sv[1],d[i].corr,d[i].code);
  return 0;
}
/* 

   compute epicenteral distance on sphere given two locations given as 

   theta1,phi1  and theta2, phi2

   in radians

*/

CPREC delta_dist(CPREC theta1, CPREC phi1, CPREC theta2, CPREC phi2)
{
  double dphi,dtheta,tmp1,tmp2,tmp3;
  dphi =   ((double)phi1 - (double)phi2)/2.;
  dtheta = ((double)theta1 - (double)theta2)/2.;
  tmp1 = sin(dtheta);tmp1 *= tmp1;
  tmp2 = sin(dphi);tmp2 *= tmp2;
  tmp3 = tmp1 + sin(theta1) * sin(theta2) * tmp2;
  tmp3 = sqrt(tmp3);
  return (CPREC)(2.0*asin(tmp3));
}
