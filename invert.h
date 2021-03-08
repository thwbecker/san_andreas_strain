//#define DEBUG
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "precision.h"
#include "deltachi2.h"

#define CLEN 200



#define boolean unsigned short
#define TRUE 1
#define FALSE 0
#define MEMERROR(x) {fprintf(stderr,"%s: out of memory\n",x);exit(-1);}
#define NRERROR(x) {fprintf(stderr,"NRERROR: %s\n",x);exit(-1);}



#include "nrutil.h"

#define NPLIM 100

/* samples per degree of freedom */
//#define NRANDOM_SAMPLE  700000
//#define NRANDOM_SAMPLE 2000
#define NRANDOM_SAMPLE 500

/* 

   sample for monte carlo iteration over GPS uncertainties

*/
#define MC_SAMPLE 150000
//#define MC_SAMPLE 25000
//#define MC_SAMPLE 5000
//#define MC_SAMPLE 1

/* 
   exit value for solutions in inner probability range
*/
#define MC_SAMPLE_INNER_MIN 20000


/* 
   for LEVMAR routines
*/
#define LEVMAR_INIT -1
#define LEVMAR_RUN 0
#define LEVMAR_FINISH 1

/* 

solver modes

*/
#define RANDOM_SOLVER 0 
#define LEVMAR_SOLVER 1
#define SIMPLEX_SOLVER 2


/* 

model parameters


*/


struct data{			/* data structure */
  CPREC x, y;
  CPREC sy;
  boolean use ;
};

struct flt{			/* fault structure */
  CPREC xoff,w,u;
  boolean xoff_fixed, w_fixed;
};
struct sol{			/* solution structure */
  float *par;			/* parameters */
  float chi2;
};
struct mdl{			/* model structure */
  int ndata,npar,dof,nfault,nfree_par,nsave_par,nfixed_par;
  struct data *d,*orig_d;
  struct flt *fault;
  struct sol *solution;
  struct sol best_solution;
  CPREC xmin,xmax,xrange,ymean,yoff,wmin,wrange,umin,wmax;
  char datafile[CLEN];
  long seed;
};

/* deltachi2.c */
CPREC rfunc(CPREC, CPREC, CPREC);
CPREC delta_chi2_rtbis(CPREC, CPREC);
CPREC gammq(CPREC, CPREC);
void gser(CPREC *, CPREC, CPREC, CPREC *);
void gcf(CPREC *, CPREC, CPREC, CPREC *);
CPREC gammln(CPREC);
/* invert.c */
void init_levmar_para(int *, CPREC *, struct mdl *);
void solve_for_fault_slip(struct mdl *);
int sort_by_xoff(const void *, const void *);
int sort_by_chi2(const void *, const void *);
void print_fitting_function(struct mdl, FILE *, int);
void consolidate_data(struct mdl *);
void generate_random_data_realization(struct mdl *);
void print_best_result(struct mdl *, FILE *, int, unsigned short, unsigned short, CPREC **);
CPREC fitting_function(struct flt *, int, CPREC, CPREC);
void fault_yoff_to_para(struct flt *, int, CPREC, CPREC **, int *);
CPREC misfit(struct mdl *);
CPREC scaled_atan(CPREC, CPREC, CPREC, CPREC);
void copy_fault_par_to_sol(struct mdl *, CPREC *, unsigned short);
void copy_sol_to_fault_par(CPREC *, struct mdl *, unsigned short);
void read_fault_data(struct flt *, int, FILE *, char *);
CPREC urandom(CPREC, CPREC, struct mdl *);
CPREC ran2(long *);
CPREC gasdev(long *);
void nr_fitting_function(CPREC, CPREC *, CPREC *, CPREC *, int);
int mrqmin(CPREC *, CPREC *, CPREC *, int, CPREC *, int *, int, CPREC **, CPREC **, CPREC *, CPREC *, int);
void covsrt(CPREC **, int, int *, int);
int gaussj(CPREC **, int, CPREC **, int);
void mrqcof(CPREC *, CPREC *, CPREC *, int, CPREC *, int *, int, CPREC **, CPREC *, CPREC *);
int lu_solve(CPREC **, int, CPREC **, int);
void ludcmp(CPREC **, int, int *, CPREC *);
void lubsksb(CPREC **, int, int *, CPREC *);
boolean valid_solution(struct mdl *);

/* nrutil.c */
void nrerror(char []);
float *vector(long, long);
int *ivector(long, long);
unsigned char *cvector(long, long);
unsigned long *lvector(long, long);
double *dvector(long, long);
float **matrix(long, long, long, long);
double **dmatrix(long, long, long, long);
int **imatrix(long, long, long, long);
float **submatrix(float **, long, long, long, long, long, long);
float **convert_matrix(float *, long, long, long, long);
float ***f3tensor(long, long, long, long, long, long);
void free_vector(float *, long, long);
void free_ivector(int *, long, long);
void free_cvector(unsigned char *, long, long);
void free_lvector(unsigned long *, long, long);
void free_dvector(double *, long, long);
void free_matrix(float **, long, long, long, long);
void free_dmatrix(double **, long, long, long, long);
void free_imatrix(int **, long, long, long, long);
void free_submatrix(float **, long, long, long, long);
void free_convert_matrix(float **, long, long, long, long);
void free_f3tensor(float ***, long, long, long, long, long, long);

/* prune_data.c */
CPREC delta_dist(CPREC, CPREC, CPREC, CPREC);

CPREC amotry(CPREC **,CPREC *,CPREC *,int ,int ,CPREC ,struct mdl *);
CPREC amo_cost_func(CPREC *, struct mdl *,boolean *);
int amoeba(CPREC **,CPREC *,int , CPREC ,int *, struct mdl *);
