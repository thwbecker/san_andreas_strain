#include "invert.h"

void init_levmar_para(int *nr_ia, CPREC *nr_a, struct mdl *model){

  int i,j;

  j = 1;
  /* set initial values */
  nr_ia[j]  = 1;		/* always fit offset */
  nr_a[j] = 0;
  j = 2;
  for(i=0;i < model->nfault; i++){
    /* offset */
    if(model->fault[i].xoff_fixed){
      nr_ia[j] = 0;		/* fixed */
      nr_a[j]  = model->fault[i].xoff;
    }else{
      nr_ia[j] = 1;		/* solve */
      nr_a[j]  = model->xmin + 
	(model->xrange/(CPREC)(model->nfault+1)) * (CPREC)(i+1); /* default location */
      nr_a[j] += urandom(-model->xrange/50,model->xrange/50,model); /* add some noise */
      model->fault[i].xoff = nr_a[j];
    }
    j++;
    /* locking depth */
    if(model->fault[i].w_fixed){ /* fixed */
      nr_ia[j] = 0;
      nr_a[j]  = model->fault[i].w;
    }else{			/* solve */
      nr_ia[j] = 1;
      nr_a[j]  = 7.5;
      nr_a[j] += urandom(-2.5,2.5,model); /* noise */
      model->fault[i].w = nr_a[j];
    }
    j++;
    /* slip */
    nr_ia[j] = 1;
    nr_a[j] = 50/(CPREC)model->nfault;
    nr_a[j] += urandom(-2,2,model);
    model->fault[i].u = nr_a[j];
    j++;
  }
  if(j-1 != model->npar){fprintf(stderr,"init_lev_mar_para: logic error lev mar %i %i\n",
				 j-1,model->npar);
    exit(-1);
  }

}


/* 

   for the given data set and fault geometry, solve for fault slip and
   best fit offset


*/

void solve_for_fault_slip(struct mdl *model)
{

  /* 
     
     given assigned data vector, compute new A matrix based on fault
     geometry, solve, and assign slip values to fault structure
     
     
  */
  int i,j,imode;  
  CPREC rnorm, *a,*u,*b;

  a = (CPREC *)malloc(sizeof(CPREC)*(model->nfault+1) * model->ndata);
  b = (CPREC *)malloc(sizeof(CPREC) * model->ndata);
  u = (CPREC *)malloc(sizeof(CPREC) * (model->nfault+1));

  
  /* data */
  for(i=0;i < model->ndata;i++){
    b[i] = model->d[i].y / model->d[i].sy;
    /* first row is for mean offset */
    a[i] = 1.0 /  model->d[i].sy;
  }

  /* 
     
     assemble design matrix 

  */
  model->yoff = 0.0;
  for(j=0;j < model->nfault;j++){	/* loop through basis functions */
    for(i=0;i < model->ndata;i++){	/* loop through data */
      a[(j+1) * model->ndata + i] = 
	scaled_atan((model->fault+j)->xoff,(model->fault+j)->w,1.0,model->d[i].x) /  
	model->d[i].sy;
    }
  }
  /* 
     solve for the non-negative least squares solution 
  */
  nnls_driver(model->ndata,(model->nfault+1),1,a,b,u,&rnorm,&imode);
  if(imode){
    fprintf(stderr,"solver: error in NNLS, %i\n",imode);
    exit(-1);
  }
  /* 
     assign solution 
  */
  model->yoff = u[0];
  for(i=0;i < model->nfault;i++){
    //fprintf(stderr,"%i  %g %g %g\n",i, model->fault[i].xoff,  model->fault[i].w, u[i]);
    model->fault[i].u = u[i+1];
  }
  free(a);free(u);free(b);
}

int sort_by_xoff(const void *x, const void *y)
{

  struct flt *one,*two;
  one = (struct flt *)x;
  two = (struct flt *)y;

  if(one->xoff < two->xoff)
    return -1;
  else if(one->xoff == two->xoff)
    return 0;
  else
    return 1;
}

int sort_by_chi2(const void *x, const void *y)
{
  struct sol *one, *two;
  one = (struct sol *)x;
  two = (struct sol *)y;

  if(one->chi2 < two->chi2)
    return -1;
  else if(one->chi2 == two->chi2)
    return 0;
  else
    return 1;
}

void print_fitting_function(struct mdl model, FILE *out, int mode)
{
  /* 
     mode == 1:  add all up 
     mode < 0:  only show individual function and add header
  */
  static int nstep = 1000;
  CPREC x,dx,x1,xoff,yoff;
  dx = (model.xmax-model.xmin)/(nstep-1);
  x1 = model.xmax + 1e-7;
  if(mode == 1){		/* whole function */
    for(x = model.xmin;x <= x1;x += dx)
      fprintf(out,"%g %g\n",x,fitting_function(model.fault,model.nfault,x,model.yoff));
  }else{
    /* 
       individual, shifted by the y offset of the total function at
       the offset point
    */
    if((mode > -1)||(mode < -(int)model.nfault)){
      fprintf(stderr,"print_fitting_function: mode %i out of bounds\n",mode);
      exit(-1);
    }
    mode = -mode-1;
    xoff =  model.fault[mode].xoff; /* offset of this function */
    yoff =  fitting_function(model.fault,model.nfault,xoff,model.yoff);
    /* header line */
    fprintf(out,"# func %i out of %i off %g W %g slip %g\n",
	    mode+1,model.nfault,model.fault[mode].xoff,model.fault[mode].w,model.fault[mode].u);
    for(x = model.xmin;x <= x1;x += dx)
      fprintf(out,"%g %g\n",x, 
	      scaled_atan((model.fault+mode)->xoff,
			  (model.fault+mode)->w,
			  (model.fault+mode)->u,x)+yoff);
  }
}

/* 

   go through all data and throw out those that were not used and do
   some ymean xmin xmax stats

*/
void consolidate_data(struct mdl *model)
{
  struct data *tmp;
  int n,i;
  tmp = (struct data *)malloc(sizeof(struct data));
  model->ymean = 0.0;
  model->xmin = 1e20;model->xmax = -1e20;
  n = 0;
  for(i=0;i < model->ndata;i++){
    if(model->d[i].use){
      memcpy((tmp+n),(model->d+i),sizeof(struct data));
      /* do some stats */
      model->ymean += tmp[n].y;
      if(tmp[n].x < model->xmin)
	model->xmin = tmp[n].x;
      if(tmp[n].x > model->xmax)
	model->xmax = tmp[n].x;
      n++;
      tmp = (struct data *)realloc(tmp,sizeof(struct data)*(n+1));
    }
  }
  fprintf(stderr,"consolidate_data: using %i out of %i data\n",n,model->ndata);
  model->ndata = n;
  model->ymean /= (CPREC)model->ndata;
  /* copy and shrink */
  memcpy(model->d,tmp,sizeof(struct data)*model->ndata);
  model->d = (struct data *)realloc(model->d,sizeof(struct data)*model->ndata);
  free(tmp);
  model->xrange = model->xmax - model->xmin;
  /* 
     store to original data copy 
  */
  model->orig_d = (struct data *)realloc(model->orig_d,sizeof(struct data)*model->ndata);
  memcpy(model->orig_d,model->d,sizeof(struct data)*model->ndata);
  
}

void generate_random_data_realization(struct mdl *model)
{
  int i;
  CPREC dnorm;
  /* restore */
  memcpy(model->d,model->orig_d,sizeof(struct data)*model->ndata);
  /* add noise */
  for(i=0;i < model->ndata;i++){
    model->d[i].y += gasdev(&model->seed) * model->d[i].sy;
  }
#ifdef SUPER_DEBUG  
  dnorm = 0.0;
  for(i=0;i < model->ndata;i++)
    dnorm += (model->d[i].y)*(model->d[i].y);
  fprintf(stderr,"generate_random_data_realization: data norm: %20.7e\n",sqrt(dnorm));
#endif
}


void print_best_result(struct mdl *model, FILE *out, int mode,
		       boolean store_slip, boolean use_covar, CPREC **nr_covar)
{
  int i,j;
  FILE *fileout;
  char filename[CLEN];
  /*  */
  copy_sol_to_fault_par(model->solution[0].par,model,store_slip);

  switch(mode){
  case 0:
    /* 
       best-fit parameters 
    */
    fprintf(out,"ndata= %i ymean= %.4e nfunc= %i npar= %i dof= %i chi2(WSSR)= %.5e sqrt(red_chi2)= %.5f uleft: %.5f uright: %.5f\n",
	    model->ndata,model->ymean,model->nfault,model->npar,model->dof,
	    model->solution[0].chi2,sqrt(model->solution[0].chi2/(CPREC)model->dof),
	    fitting_function(model->fault,model->nfault,-1e10,model->yoff),
	    fitting_function(model->fault,model->nfault, 1e10,model->yoff));
    if(use_covar){
      fprintf(out,"best yoff: %11g (%8.3e)\n",model->yoff,sqrt(nr_covar[1][1]));
    }else
      fprintf(out,"best yoff: %11g\n",model->yoff);
    for(j=2,i=0;i < model->nfault;i++,j+=3){
      fprintf(out,"func%03i: ",i+1);
      if(use_covar){
	fprintf(out,"par_%02i= %11g ( %8.3e ) par_%02i= %11g ( %8.3e ) par_%02i= %11g ( %8.3e )",
		1,model->fault[i].xoff,sqrt(nr_covar[j][j]),
		2,model->fault[i].w,sqrt(nr_covar[j+1][j+1]),
		3,model->fault[i].u,sqrt(nr_covar[j+2][j+2]));
      }else{
	fprintf(out,"par_%02i= %11g par_%02i= %11g par_%02i= %11g",
		1,model->fault[i].xoff,
		2,model->fault[i].w,
		3,model->fault[i].u);
      }
      fprintf(out,"\n");
    }
    break;
  case 1:
    /* 
       fitting functions to file 
    */
    print_fitting_function(*model,out,1);
    /* 
       individual 
    */
    for(i=0;i < model->nfault;i++){
      sprintf(filename,"invert.fit.%i.out",i+1);
      fileout = fopen(filename,"w"); /* total */
      print_fitting_function(*model,fileout,-(i+1));
      fclose(fileout);
    }
    break;
  case 2:			/* original data file name and data
				   actually used for fitting */
    fprintf(out,"# %s\n",model->datafile);
    for(i=0;i<model->ndata;i++)
      if(model->d[i].use)
	fprintf(out,"%.6e %.6e %.6e\n",model->d[i].x,model->d[i].y,model->d[i].sy);
    break;
  default:
    fprintf(stderr,"mode %i undefined in print_best_result\n",mode);
    exit(-1);
  }

}
/* 

   evaluate fitting function at location x, with nfaults and
   parameters in fault, offset in yoff

   generalized wrapper

*/
CPREC fitting_function(struct flt *fault,int nfault,CPREC x, CPREC yoff)
{
  CPREC *a=NULL,*dyda, y=0;
  int npara;
  /* convert from fault to list of parameters, which are automatically allocates */
  fault_yoff_to_para(fault,nfault,yoff,&a,&npara);
  /* make room for derivatives (not needed) */
  dyda = (CPREC *)malloc(sizeof(CPREC)*npara); 
  /* evaluate function and discard derivatives, call numerical recipes
     style */
  nr_fitting_function(x,(a-1),&y,(dyda-1),npara);
  /* free storage */
  free(a);free(dyda);
  return y;
}

/* 

   assigns parameters, allocates space for a, and assign to a

   pass a initialized


*/
void fault_yoff_to_para(struct flt *fault,int nfault,CPREC yoff, CPREC **a, int *npara)
{
  int i,j;
  *npara = 1 + nfault * 3;
  *a = (CPREC *)realloc(*a, sizeof(CPREC)* (*npara));
  *(*a+0) = yoff;
  for(i=j=0;i < nfault;i++,j+=3){
    *(*a+j+1) = fault[i].xoff;
    *(*a+j+2) = fault[i].w;
    *(*a+j+3) = fault[i].u;
  }
}

/* 

   compute chi2 misfit for current yoffset and fault geometry and slip

*/
CPREC misfit(struct mdl *model)
{
  int i;
  CPREC chi2, tmp, ymod;
  chi2 = 0.0;
  for(i=0;i < model->ndata;i++){	/* loop through data */
    if(model->d[i].use){
      ymod = fitting_function(model->fault,model->nfault,model->d[i].x,model->yoff);
      tmp = ((model->d[i].y) - ymod)/(model->d[i].sy);
      //fprintf(stderr,"o: %i %11g %11g %11g %11g\n",i+1,d[i].x,d[i].y,d[i].sy,ymod);
      chi2 += tmp*tmp;
    }
  }
  /*  */
  return chi2;
}

/* actual basis function */
CPREC scaled_atan(CPREC xoff, CPREC w, CPREC u, CPREC x)
{
  return atan(-(x - xoff)/w)/M_PI * u;
}



void copy_fault_par_to_sol(struct mdl *model, CPREC *par,
			   boolean store_slip)
{
  int i,pc;
  pc = 0;
  if(store_slip){
    par[pc] = model->yoff;
    pc++;
  }
  for(i=0;i < model->nfault;i++){
    if(!model->fault[i].xoff_fixed){
      par[pc] = model->fault[i].xoff;
      pc++;
    }
    if(!model->fault[i].w_fixed){
      par[pc] = model->fault[i].w;
      pc++;
    }
    if(store_slip){
      par[pc] = model->fault[i].u;
      pc++;
    }
  }
#ifdef DEBUG
  if(pc != model->nsave_par){
    fprintf(stderr,"copy_fault_par_to_sol: logic error: %i vs %i\n",
	    pc,model->nsave_par);
    exit(-1);
  }
#endif
}

void copy_sol_to_fault_par(CPREC *par, struct mdl *model,boolean store_slip)
{
  int i,pc;
  pc = 0;
  if(store_slip){			/* offset */
    model->yoff = par[pc];
    pc++;
  }
  for(i=0;i < model->nfault;i++){
    if(!model->fault[i].xoff_fixed){
      model->fault[i].xoff = par[pc];
      pc++;
    }
    if(!model->fault[i].w_fixed){
      model->fault[i].w = par[pc];
      pc++;
    }
    if(store_slip){		/* store the slip values */
      model->fault[i].u =  par[pc];
      pc++;
    }
  }
  if(!store_slip)			/* solve NNLS */
    solve_for_fault_slip(model);
#ifdef DEBUG
  if(pc != model->nsave_par){
    fprintf(stderr,"copy_sol_to_fault_par: logic error: %i vs %i\n",
	    pc,model->nsave_par);
    exit(-1);
  }
#endif
}

/* read information about the fault parameters */
void read_fault_data(struct flt *fault,int nbasis, FILE *in, char *progname)
{
  int i,yes;
  /* 
     read in fault parameters 
  */
  for(i=0;i < nbasis;i++){
    /* location */
    fprintf(stderr,"%s: fault %i, fixed location (1/0)?\n",progname,i+1);
    fscanf(in,"%i",&yes);
    if(yes){
      fault[i].xoff_fixed = TRUE;
      fprintf(stderr,"%s: fault %i, provide location (distance along profile)\n",progname,i+1);
      fscanf(in,FMT_1,&(fault[i].xoff));
    }else{
      fault[i].xoff_fixed = FALSE;
    }
    /* locking depth */
    fprintf(stderr,"%s: fault %i, fixed locking depth (1/0)?\n",progname,i+1);
    fscanf(in,"%i",&yes);
    if(yes){
      fault[i].w_fixed = TRUE;
      fprintf(stderr,"%s: fault %i, provide locking depth\n",progname,i+1);
      fscanf(in,FMT_1,&(fault[i].w));
    }else{
      fault[i].w_fixed = FALSE;
    }
  }
  /* sort faults by x */
  qsort(fault,nbasis,sizeof(struct flt),sort_by_xoff);
  /* 
     report fault settings 
  */
  for(i=0;i < nbasis;i++){
    fprintf(stderr,"%s: fault %3i: ",progname,i+1);
    if(fault[i].xoff_fixed)
      fprintf(stderr,"location fixed at %11g ",fault[i].xoff);
    else
      fprintf(stderr,"location free                 ");
    if(fault[i].w_fixed)
      fprintf(stderr,"locking depth fixed at %11g\n",fault[i].w);
    else
      fprintf(stderr,"locking depth free\n");
  }


}
/* return a uniformly distributed random number between x0 and x1 */
CPREC urandom(CPREC x0, CPREC x1, struct mdl *model)
{
  return x0 + ran2(&(model->seed)) * (x1-x0);
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-EPS)

CPREC ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  CPREC temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


CPREC gasdev(long *idum)
{
  static int iset=0;
  static CPREC gset;
  CPREC fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran2(idum)-1.0;
      v2=2.0*ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

/* 
   actual fitting function, expects NR style (1..n) vectors

   x: location
   a: parameters in 1,2,...na format
   y: return value 
   dyda[1,....,na]: derivatives wrt. to para 
   na: number of parameters

*/
void nr_fitting_function(CPREC x, CPREC *a, CPREC *y, CPREC *dyda, int na)
{
  int i,i3,j,nfault;
  CPREC tmp_atan,tmp_atand,u,w,xoff,dx;

  nfault = (na - 1)/3;
  j= 1;
  *y        = a[j];			/* offset */
  dyda[j++] = 1;			/* derivative */

  for(i=i3=0;i < nfault;i++,i3+=3){
    xoff = a[2+i3];
    w =    a[3+i3];
    u =    a[4+i3];

    dx = x - xoff;
    tmp_atan =  atan(-dx/w)/M_PI; 
    /* function value */
    *y += tmp_atan * u;	/* u/pi * atan(-(x-xoff)/w)  */
    //fprintf(stderr,"%5i  - %11g - %11g %11g %11g %11g - %11g\n",i+1,a[1],xoff,dx,w,u,*y);
    /* derivative of atan(x) = 1/(x^2+1) */
    tmp_atand = dx/w;
    tmp_atand = 1./(tmp_atand * tmp_atand + 1.0)/M_PI*u/w; /* 1/(((x-xoff)/w)^2 + 1) * u/pi/w a*/
    /* dy/dx_off */
    dyda[j++] = tmp_atand;
    /* dy/dw  */
    dyda[j++] = tmp_atand/w * dx;
    /* dy/du */
    dyda[j++] = tmp_atan;
  }



}
/* 

   mode = -1: init
   mode =  1, close, free arrays 
   else mode = 0


   return < 0 on error, zero else
*/
/* pass numerical recipes style */

int mrqmin(CPREC *x,CPREC *y,CPREC *sig,int ndata,CPREC *a,
	   int *ia,int ma,CPREC **covar,CPREC **alpha,CPREC *chisq,
	   CPREC *alamda,int mode)
{

  int j,k,l,m,err;
  static int mfit;
  static CPREC ochisq,*atry,*beta,*da,**oneda;
  if (mode == LEVMAR_INIT) {		/* init */
    atry = nr_vector(1,ma);
    beta = nr_vector(1,ma);
    da=nr_vector(1,ma);
    for (mfit=0,j=1;j<=ma;j++)
      if (ia[j]) 
	mfit++;			/* count the parameters to be fit */
    oneda=nr_matrix(1,mfit,1,1);
    *alamda=0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq);
    ochisq=(*chisq);
    for (j=1;j<=ma;j++) 
      atry[j]=a[j];
  }
  if(mode == LEVMAR_FINISH)	/* last call */
    *alamda = 0.0;		/* finish */
  /*  */
  for (j=1;j <= mfit;j++){
    for (k=1;k <= mfit;k++) 
      covar[j][k] = alpha[j][k];

    covar[j][j] = alpha[j][j]*(1.0+(*alamda));
    oneda[j][1] = beta[j];
  }
  /* 

     solver covar * x = oneda
     and return covar^-1

  */
  if((err = gaussj(covar,mfit,oneda,1)) < 0){
    //if((err = lu_solve(covar,mfit,oneda,1)) < 0){
    /* 
       free pointers on error exit
    */
    free_nr_matrix(oneda,1,mfit,1,1);
    free_nr_vector(da,1,ma);
    free_nr_vector(beta,1,ma);
    free_nr_vector(atry,1,ma);
    return err;
  }
  for (j=1;j<=mfit;j++) 
    da[j]=oneda[j][1];
  if (mode == LEVMAR_FINISH) {
    /* finalize */
    *alamda = 0.0;
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    free_nr_matrix(oneda,1,mfit,1,1);
    free_nr_vector(da,1,ma);
    free_nr_vector(beta,1,ma);
    free_nr_vector(atry,1,ma);
    return 0;
  }
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) 
      atry[l] = a[l] + da[++j];
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq);
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    for(j=1;j<=mfit;j++){
      for (k=1; k <= mfit;k++) 
	alpha[j][k]=covar[j][k];
      beta[j]=da[j];
    }
    for(l=1;l<=ma;l++)
      a[l] = atry[l];
  } else {
    *alamda *= 10.0;
    *chisq=ochisq;
  }
  return 0;
}




#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/* pass numerical recipes style */

void covsrt(CPREC **covar,int ma,int *ia,int mfit)
{
  int i,j,k;
  CPREC swap;
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) {
      covar[i][j]=covar[j][i]=0.0;
    }
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) {
	SWAP(covar[i][k],covar[i][j]);
      }
      for (i=1;i<=ma;i++) {
	SWAP(covar[k][i],covar[j][i]);
      }
      k--;
    }
  }
}
#undef SWAP

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
/* pass numerical recipes style */

int gaussj(CPREC **a,int n,CPREC **b,int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  CPREC big,dum,pivinv,temp;
#ifdef SUPER_DEBUG
  fprintf(stderr,"\nA\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      fprintf(stderr,"%25.17e ",a[i][j]);
    fprintf(stderr,"\n");
  }
#endif

  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) {
	    fprintf(stderr,"gaussj: Singular Matrix-1\n");
	    return -1;
	  }
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			   for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
						}
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) {
      fprintf(stderr,"gaussj: Singular Matrix-2\n");
      return -2;
    }
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
#ifdef SUPER_DEBUG
  fprintf(stderr,"\nA^{-1}\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      fprintf(stderr,"%25.15e ",a[i][j]);
    fprintf(stderr,"\n");
  }
#endif

  return 0;
}
#undef SWAP
/* pass numerical recipes style */

/* solve by LU */
int lu_solve(CPREC **a,int n,CPREC **b,int m)
{
  int *indx,i,j;
  CPREC d,*col,**ainv,*bloc;
#ifdef SUPER_DEBUG
  fprintf(stderr,"\nA\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      fprintf(stderr,"%25.15e ",a[i][j]);
    fprintf(stderr,"\n");
  }
#endif
  col = nr_vector(1,n);
  bloc = nr_vector(1,n);
  indx = ivector(1,n);
  ainv = nr_matrix(1,n,1,n);
  /*  */
  ludcmp(a,n,indx,&d);		/* LU decompose */
  for(i=1;i<=m;i++){
    for(j=1;j<=n;j++)
      bloc[j] = b[j][i];
    lubsksb(a,n,indx,bloc);		/* solve for all B(i) */
    for(j=1;j<=n;j++)
      b[j][i] = bloc[j];
  }
  /* compute inverse of A */
  for(j=1;j <= n;j++){
    for(i=1;i <= n;i++)
      col[i] = 0.0;
    col[j] = 1.0;
    lubsksb(a,n,indx,col);
    for(i=1;i<=n;i++)
      ainv[i][j] = col[i];
  }
  /* copy to A for return */
  for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
      a[i][j] = ainv[i][j];
  /*  */
  free_nr_matrix(ainv,1,n,1,n);
  free_ivector(indx,1,n);
  free_nr_vector(col,1,n);
  free_nr_vector(bloc,1,n);
#ifdef SUPER_DEBUG
  fprintf(stderr,"\nA^{-1}\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      fprintf(stderr,"%25.15e ",a[i][j]);
    fprintf(stderr,"\n");
  }
#endif
  return 0;
}

#define TINY 1.0e-20;
//#define TINY 1.0e-25;
/* pass numerical recipes style */

void ludcmp(CPREC **a,int n,int *indx,CPREC *d)
{
  int i,imax,j,k;
  float big,dum,sum,temp;
  float *vv;

  vv=vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,1,n);
}
#undef TINY
/* pass numerical recipes style */
void lubsksb(CPREC **a,int n,int *indx,CPREC *b)
{
  int i,ii=0,ip,j;
  float sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/* pass numerical recipes style */
void mrqcof(CPREC *x,CPREC *y,CPREC *sig,int ndata,CPREC *a,
	    int *ia,int ma,CPREC **alpha,CPREC *beta,CPREC *chisq)
{
  int i,j,k,l,m,mfit=0,nf;
  CPREC ymod,wt,sig2i,dy,*dyda;
  const int add_apriori_bounds = 0; /* add physical limits (not a good
				       idea!!!) */

  dyda=nr_vector(1,ma);
  for (j=1;j<=ma;j++)
    if (ia[j]) 
      mfit++;
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=j;k++) 
      alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  *chisq=0.0;
  for (i=1;i <= ndata;i++) {
    nr_fitting_function(x[i],a,&ymod,dyda,ma);
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	wt=dyda[l]*sig2i;
	for (j++,k=0,m=1;m<=l;m++)
	  if (ia[m]) 
	    alpha[j][++k] += wt*dyda[m];
	beta[j] += dy*wt;
      }
    }
    //fprintf(stderr,"n: %i %11g %11g %11g %11g\n",i,x[i],y[i],sig[i],ymod);
    *chisq += dy*dy*sig2i;
  }
  //for(i=1;i<=mfit;i++)fprintf(stderr,"%20.5e ",alpha[i][i]);fprintf(stderr,"\n");
  if(add_apriori_bounds){

    /* check for negative slip and large W */
    nf = (ma -1)/3;		/* number of faults */
    for(i=0;i < nf;i++){
      if(a[1+i*3+3] < 0)		/* slip */
	*chisq += 1e5;
      if(a[1+i*3+2] > 15)
	*chisq += 1e5;
    }
  }
  
  
  //fprintf(stderr,"n: chi2: %11g\n",*chisq);
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++) 
      alpha[k][j]=alpha[j][k];
  
  free_nr_vector(dyda,1,ma);
}



#define NMAXA 15000
#define GET_PSUM					\
  for (j=1;j<=ndim;j++) {				\
    for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];	\
    psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/* pass numerical receipes style */
int amoeba(CPREC **p,CPREC *y,int ndim, CPREC ftol,int *nfunk, struct mdl *model)
{
  int i,ihi,ilo,inhi,j,mpts=ndim+1;
  CPREC rtol,sum,swap,ysave,ytry,*psum;
  boolean valid;
  psum=nr_vector(1,ndim);
  *nfunk=0;
  GET_PSUM;
  for (;;) {
    ilo=1;
    ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
    for (i=1;i<=mpts;i++) {
      if (y[i] <= y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi] && i != ihi) inhi=i;
    }
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) {
      SWAP(y[1],y[ilo]);
      for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]);
      break;
    }
    if (*nfunk >= NMAXA) 
      return -1;
    *nfunk += 2;
    ytry=amotry(p,y,psum,ndim,ihi,-1.0,model);
    if (ytry <= y[ilo])
      ytry=amotry(p,y,psum,ndim,ihi,2.0,model);
    else if (ytry >= y[inhi]) {
      ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,ihi,0.5,model);
      if (ytry >= ysave) {
	for (i=1;i<=mpts;i++) {
	  if (i != ilo) {
	    for (j=1;j<=ndim;j++)
	      p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	    y[i]=amo_cost_func((psum+1),model,&valid);
	  }
	}
	*nfunk += ndim;
	GET_PSUM;
      }
    } else --(*nfunk);
  }
  free_nr_vector(psum,1,ndim);
  return 0;
}
#undef SWAP
#undef GET_PSUM
#undef NMAXA
/* pass numerical recipes styles */
CPREC amotry(CPREC **p,CPREC *y,CPREC *psum,int ndim,int ihi,CPREC fac,struct mdl *model)
{
  int j;
  CPREC fac1,fac2,ytry,*ptry;
  boolean valid;
  ptry=nr_vector(1,ndim);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=1;j<=ndim;j++) 
    ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=amo_cost_func((ptry+1),model,&valid);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free_nr_vector(ptry,1,ndim);
  return ytry;
}

/* 
   evaluate the cost function for a set of parameters 

   pass para in C style, not numerical recipes style

*/
CPREC amo_cost_func(CPREC *para, struct mdl *model,boolean *valid)
{				
  int i,j=0;
  CPREC cost;
  /* copy parameters over and solve */
  copy_sol_to_fault_par(para,model,FALSE);
  /* chi2 */
  cost = misfit(model);
  /* check if solution is out of bounds */
  *valid = valid_solution(model);
  if(!(*valid))
    cost += 1e6;
  return cost;
}
/* 

check if this solution is out of bounds

*/
boolean valid_solution(struct mdl *model)
{
  int i;
  /* bounds */
  for(i=0;i < model->nfault;i++){
    if(model->fault[i].u < model->umin)
      return FALSE;
    if(model->fault[i].w < model->wmin)
      return FALSE;
    if(model->fault[i].w > model->wmax)
      return FALSE;
    if(!model->fault[i].xoff_fixed){
      if(model->fault[i].xoff < model->xmin)
	return FALSE;
      if(model->fault[i].xoff > model->xmax)
	return FALSE;
    }
  }
  return TRUE;
}
