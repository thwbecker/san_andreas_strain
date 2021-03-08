#include "invert.h"
/* 


   Levenberg-Marquardt or random inversion for best fitting fault
   locations and locking depths, fitting atan profiles

   An outer Monte Carlo approach attempts to parameter robustness
   output and uses non-negative least squares for slip rate inferences

   USAGE:


   invert data_file nbasis [xmin, -800] [xmax, 800] [add_random, 0] [wrange {100,60}] [umin, 0]


   the program will then, for each nbasis fault, if the fault locatio
   is fixed, and if the locking depth is fixed. if so, the respective
   next line from stdin has to specify the value

   if nbasis > 0, will use random
   if nbasis < 0, will use Levenberg-Marquardt
   if nbasis > 100, will use simplex solver

   
   - datafile is expected to have 

   x v_normal sigma_v_normal in column 1, 2, and 3 of the datafile, according to 

   #define FMT_LINE1 "%f %f %f"

   - xmin and xmax will limit the data range
   
   - if add_random > 0, will add random faults for constrained inversion

   - if wrange is set, will limit locking depth to [wmin, wmin +
   wrange] where wmin = 1 km by default, and wrange 100 and 60 km for
   single or multiple faults, respectively

   - umin is the minimum slip allowed


   $Id: invert.c,v 1.2 2009/07/23 19:02:02 becker Exp becker $
   
*/
int main(int argc, char **argv)
{

  int nbasis,add_random,nbig;

  /* 

     bounds on parameters

  */
  const int iter_max_umin = 1000;     /* max trials to achieve minimum slip */
  

  /*  */
  const CPREC wrange0[2] = {100, 60};  /* maximum locking depth for
					  single and multiple faults */

  /* minimum spacing between faults */
  const CPREC min_xoff  = 10;
  /* range within xmin ... xmax to use for faults */
  const CPREC xuse_range = 0.9;
  /* 
     
     data selection parameters

  */
  /* max velocity error bound */
  const CPREC max_err = 4;

  /* min distance between two data points */
  const CPREC eps_dist_km = .5;
  //const CPREC eps_dist_km = -1;

  /* replace or average close locations? */
  const boolean replace_close = FALSE;

  const boolean use_orig_data_for_chi2 = TRUE;
  /* 

     should we iterate to throw out outliers?

  */
  const int big_sample = 2;	/* 2: yes, iterate 1: nope, just one go */
  const CPREC du_outlier = 6; /* treat all data further than
				   du_outlier as outliers */

  /* maximum distance for data, xmin <= x <= xmax */
  CPREC xmin = -800, xmax = 800;

  /* 
     
     for delta chi2 statistic
     for inner circle (determines iteration break)
  */
  //const CPREC prob_del_chi2_inner = 0.954;/* 95.4% */
  const CPREC prob_del_chi2_inner = 0.99;

  //const CPREC prob_del_chi2_inner = 0.68;/* 68% */
  CPREC del_chi2,chi2_lim;	
  
  /* 

   */
  FILE *in,*out;
  int i,j,ran_sample,mc_sample,npara,istart,i_ran,ii_ran,iii_ran,iiii_ran,iiiii_ran,
    i_mc,i_big,i_levmar,ii_levmar,ierror,np,simplex_sample,
    nchi2_close;
  CPREC  umin_loc,off_min,off_range,*apara=NULL,yp,dist,last_chi2,dchi2,simplex_eps;
  boolean ok,inner_init,inner_exit,store_slip,modified,valid;
  int solver_mode,amoeba_return;
  /*  */
  CPREC *nr_x, *nr_y, **nr_alpha, *nr_sig, **nr_covar, **best_nr_covar,*nr_a, nr_chi2, nr_alamda,min_chi2,avg_chi2;
  int *nr_ia;
  struct sol solution;
  
  /* master structure */
  struct mdl model;
  /* 

     defaults 
  */
  nbasis = 1;
  add_random = 0;
  model.wmin = 1;		/* min locking depth */
  model.umin = 0;		/* minimum slip-rate allowed */

  /* general init */
  model.seed = -1;
  /* init random */
  ran2(&(model.seed));


  if(argc < 2){
    fprintf(stderr,"%s datafile [nbasis, %i] [xmin, %g] [xmax, %g] [add_random, %i] [wrange, {%g, %g}] [umin, %g]\n",
	    argv[0],nbasis,xmin,xmax,add_random,wrange0[0],wrange0[1],model.umin);
    exit(-1);
  }
  /* options */
  sprintf(model.datafile,"%s",argv[1]);
  if(argc > 2)			/* number of faults */
    sscanf(argv[2],"%i",&nbasis);
   /* 
     
     decide on solver

  */
  if(nbasis < 0){
    solver_mode = LEVMAR_SOLVER;
    nbasis = -nbasis;
    store_slip = TRUE;
  }else if(nbasis > 100){
    solver_mode = SIMPLEX_SOLVER;
    nbasis -= 100;
    store_slip = FALSE;
  }else{
    solver_mode = RANDOM_SOLVER;
    store_slip = FALSE;
  }
  if(argc > 3)			/* xrange */
    sscanf(argv[3],FMT_1,&xmin);
  if(argc > 4)
    sscanf(argv[4],FMT_1,&xmax);
  if(argc > 5)			/* random additional faults? */
    sscanf(argv[5],"%i",&add_random);
  /* range of locking depths */
  model.wrange = (model.nfault == 1)?(wrange0[0]):(wrange0[1]);
  /*  */
  if(argc > 6)			/* adjust max locking depth? */
    sscanf(argv[6],FMT_1,&(model.wrange));
  if(argc > 7)			
    sscanf(argv[7],FMT_1,&(model.umin));
  /*  */
  model.wmax = model.wmin + model.wrange;
  fprintf(stderr,"%s: nbasis: %i xmin: %g  xmax: %g add_random: %i w: %g-%g umin: %g solver type: %i\n",
	  argv[0],nbasis,xmin,xmax,add_random,model.wmin,model.wmin+model.wrange,model.umin,solver_mode);

  /* 

     init faults

  */
  model.nfault = nbasis + add_random;
  model.yoff = 0;
  model.fault = (struct flt *)calloc(model.nfault,sizeof(struct flt));
  if(!model.fault)MEMERROR(argv[0]);
  fprintf(stderr,"%s: using %i regular faults, plus an additional %i random ones\n",
	  argv[0],nbasis,add_random);
  /* input of nbasis fault parameters */
  read_fault_data(model.fault,nbasis,stdin,argv[0]);

  /*  */
  model.npar = 1 + model.nfault * 3; /* total parameters, offset plus faults */
  /* 
     
     count free degrees of freedom 

  */
  for(model.nfree_par = model.nfixed_par = i = 0;
      i < model.nfault;i++){
    if(model.fault[i].xoff_fixed) /* offset */
      model.nfixed_par++;
    else
      model.nfree_par++;
    if(model.fault[i].w_fixed)	/* locking depth */
      model.nfixed_par++;
    else
      model.nfree_par++;
    model.nfree_par++;		/* slip */
  }
  /* do not store any prescribed values */
  if(store_slip){
    model.nsave_par = model.npar - model.nfixed_par;
  }else{
    model.nsave_par = model.nfault * 2 - model.nfixed_par; /* for random and simplex,
								 only save xoff and W */
  }

  fprintf(stderr,"%s: number of faults: %i, parameters: %i (free: %i, fixed %i, saved: %i)\n",
	  argv[0],model.nfault,model.npar,model.nfree_par,model.nfixed_par,model.nsave_par);
 


 /*  
      random para exploration
      
      don't simulate for slip
  */
  ran_sample = NRANDOM_SAMPLE * (model.nfree_par-model.nfault) + 1;
  if(solver_mode == RANDOM_SOLVER)
    fprintf(stderr,"%s: counted %i free parameters for random, using %i iterations\n",
	    argv[0],model.nfree_par,ran_sample);

  /* 
     monte carlo realization
  */
  
  /* 
     
     read in data
     
  */
  model.d = (struct data *)malloc(sizeof(struct data));
  model.orig_d = (struct data *)malloc(sizeof(struct data));
  model.ndata = 0;
  in = fopen(model.datafile,"r");
  if(!in){fprintf(stderr,"%s: error opening %s\n",argv[0],model.datafile);exit(-1);}
  while(fscanf(in,FMT_LINE1,&(model.d[model.ndata].x),
	       &(model.d[model.ndata].y),&(model.d[model.ndata].sy)) == 3){
    if((model.d[model.ndata].sy < max_err)&&
       (model.d[model.ndata].x <= xmax)&&(model.d[model.ndata].x >= xmin)){
      if(model.d[model.ndata].sy <= 0){
	fprintf(stderr,"%s: error line %i, cannot deal with non-positive errors\n",
		argv[0],model.ndata+1);
	exit(-1);
      }
      model.d[model.ndata].use = TRUE;
      if(model.ndata){			
	/*
	  more than one data point 
	*/
	if(model.d[model.ndata].x < model.d[model.ndata-1].x){ 
	  /* decrease in x is no good */
	  fprintf(stderr,"%s: error, x locations have to be monotonously increasing, please sort\n",argv[0]);
	  fprintf(stderr,"%s: %i x %g %i x %g\n",argv[0],model.ndata,model.d[model.ndata-1].x,
		  model.ndata+1,model.d[model.ndata].x);
	  exit(-1);
	}
	if(fabs(model.d[model.ndata].x - model.d[model.ndata-1].x) < eps_dist_km){ 
	  /* 
	     close to previous location
	  */
	  if(replace_close){
	    if(model.d[model.ndata].sy < model.d[model.ndata-1].sy){	       /* overwrite
										  previous if
										  better
										  constrained */
#ifdef SUPER_DEBUG
	      fprintf(stderr,"%s: replacing %11g %11g %11g with             ",
		      argv[0],model.d[model.ndata-1].x,model.d[model.ndata-1].y,model.d[model.ndata-1].sy);
	      fprintf(stderr,"%11g %11g %11g\n",model.d[model.ndata].x,model.d[model.ndata].y,model.d[model.ndata].sy);
#endif
	      model.d[model.ndata-1].use = FALSE; /* kill the last one */
	    }else{
#ifdef SUPER_DEBUG
	      fprintf(stderr,"%s: keeping   %11g %11g %11g and disregarding ",argv[0],
		      model.d[model.ndata-1].x,model.d[model.ndata-1].y,model.d[model.ndata-1].sy);
	      fprintf(stderr,"%11g %11g %11g\n",model.d[model.ndata].x,model.d[model.ndata].y,model.d[model.ndata].sy);
#endif
	      model.d[model.ndata].use = FALSE;
	    }
	  }else{
	    /* average */
#ifdef SUPER_DEBUG
	    fprintf(stderr,"%s: averaging %11g %11g %11g and ",argv[0],
		    model.d[model.ndata-1].x,model.d[model.ndata-1].y,model.d[model.ndata-1].sy);
	    fprintf(stderr,"%11g %11g %11g to ",model.d[model.ndata].x,model.d[model.ndata].y,model.d[model.ndata].sy);
#endif
	    model.d[model.ndata-1].x = (model.d[model.ndata-1].x + model.d[model.ndata].x)/2.0;
	    model.d[model.ndata-1].y = (model.d[model.ndata-1].y / model.d[model.ndata-1].sy + 
					model.d[model.ndata].y / model.d[model.ndata].sy)/
	      (1/model.d[model.ndata-1].sy + 1/model.d[model.ndata].sy);
	    model.d[model.ndata-1].sy = (model.d[model.ndata-1].sy + model.d[model.ndata].sy)/2.0;
#ifdef SUPER_DEBUG
	    fprintf(stderr,"%11g %11g %11g\n",model.d[model.ndata-1].x,model.d[model.ndata-1].y,model.d[model.ndata-1].sy);
#endif
	    model.d[model.ndata].use = FALSE;
	  }
	}
      }
      model.ndata++;
      model.d = (struct data *)realloc(model.d,sizeof(struct data)*(model.ndata+1));
    }
  }
  fclose(in);
  /* throw out removed data and compute min/max  */
  consolidate_data(&model);
  fprintf(stderr,"%s: selected %i original data from %s, xmin/xmax: %g/%g (imposed limits: %g/%g),\n%s: mean y %g\n",
	  argv[0],model.ndata,model.datafile,model.xmin,model.xmax,
	  xmin,xmax,argv[0],model.ymean);
  /* 

     end original data  input
  
  */
  i_big = 0;
  do{				/* major loop for data outlier iteration */
    if(i_big > 0){
      /* 
	 throw out all data outside sigma_outlier if a previous
	 iteration has found a best-fit model
      */
      fprintf(stderr,"%s: removing outliers outside du %g mm/yr of best solution\n",
	      argv[0],du_outlier);
      copy_sol_to_fault_par(model.solution[0].par,&model,store_slip);
      for(i = 0;i < model.ndata;i++){
	yp = fitting_function(model.fault,model.nfault,model.d[i].x,model.yoff);
	dist = fabs(yp - model.d[i].y); /* use absolute distance */
	if(dist > du_outlier){
	  fprintf(stderr,"%s: rejecting %g, %g, %g based on distance %g from best-fit\n",
		  argv[0],model.d[i].x,model.d[i].y,model.d[i].sy,dist);
	  model.d[i].use = FALSE;
	}
      }
      consolidate_data(&model);
      fprintf(stderr,"%s: consolidation done\n",argv[0]);
      /* realloc */
      model.solution = (struct sol *)realloc(model.solution,MC_SAMPLE * sizeof(struct sol));
      if(!model.solution)
	MEMERROR(argv[0]);
      for(i=0;i < MC_SAMPLE;i++){
	model.solution[i].par = (float *)realloc(model.solution[i].par,
						 sizeof(float)*model.nsave_par);
	if(!model.solution[i].par)
	  MEMERROR(argv[0]);
      }
    }else{
      /* 

	 first time around, common, initial steps

      */
      /* solutions to save */
      model.solution = (struct sol *)malloc(MC_SAMPLE * sizeof(struct sol));
      if(!model.solution)
	MEMERROR(argv[0]);
      for(i=0;i < MC_SAMPLE;i++){
	model.solution[i].par = (float *)malloc(sizeof(float)*model.nsave_par);
	if(!model.solution[i].par)
	  MEMERROR(argv[0]);
      }
    }
    /* degrees of freedom */
    model.dof = model.ndata - model.npar;
    /* 
       range for offsets 
    */
    off_min   = model.xmin + (1-xuse_range)/2. * model.xrange; /* from 0.05 to 0.95 of range */
    off_range = model.xrange * xuse_range;

    /* 

       allocate storage for all solutions

    */
    nchi2_close = 0;
    avg_chi2 = 0.0;
    if(i_big == big_sample - 1)
      mc_sample = MC_SAMPLE;	/* only last big iter will MC sample */
    else
      mc_sample = 1;
    for(i_mc = 0;(i_mc < mc_sample) && (nchi2_close <  MC_SAMPLE_INNER_MIN);i_mc++){
      /* 


	 outer monte carlo iteration for GPS uncertainties

      */
      inner_init = (i_mc==0)?(TRUE):(FALSE);
      inner_exit = (i_mc==(mc_sample-1))?(TRUE):(FALSE);
      if(!inner_init){		
#ifdef SUPER_DEBUG
	fprintf(stderr,"%s: working on random MC realization %i\n",argv[0],i_mc);
#endif
	generate_random_data_realization(&model);
      }
       /* 
	 decide on inner least squares solver
      */
      switch(solver_mode){
      case LEVMAR_SOLVER:
	/* 

	   levenberg marquardt part, optimize chi2 for current
	   (perhaps statistical dataset)
	   
	*/
	
	/* assign data */
	if(inner_init){
	  nr_x = nr_vector(1,model.ndata);
	  nr_y = nr_vector(1,model.ndata);
	  nr_sig = nr_vector(1,model.ndata);
	  for(i=0,j=1;i < model.ndata;i++,j++){
	    if(!model.d[i].use){
	      fprintf(stderr,"%s: error, Lev Mar assumes all data to be used, at least one use flag false\n",
		      argv[0]);
	      exit(-1);
	    }
	    nr_x[j]   = model.d[i].x; /* y is assigned later */
	    nr_sig[j] = model.d[i].sy;
	  }
	  /* add random faults? */
	  if(add_random > 0){
	    fprintf(stderr,"%s: adding random faults not implemented yet for Lev-Mar\n",argv[0]);
	    exit(-1);
	  }
	  /* initial values for parameters */
	  nr_ia    = ivector(1,model.npar); /* should we invert for this parameter? */
	  nr_a     = nr_vector(1,model.npar);
	  nr_covar = nr_matrix(1,model.npar,1,model.npar);
	  best_nr_covar = nr_matrix(1,model.npar,1,model.npar);
	  nr_alpha = nr_matrix(1,model.npar,1,model.npar);
	}
	
	for(i=0,j=1;i < model.ndata;i++,j++) /* data assign to NR
						array (redo every
						time) */
	  nr_y[j]   = model.d[i].y;
	
	/* 
	   init the lev mar algorithm
	 */
	init_levmar_para(nr_ia,nr_a,&model);

	/* init Lev Mar */
	ok = FALSE;
	/* reset min chi2 */
	min_chi2 = 1e20;
  
	i_levmar = ii_levmar = 0;
	do{			
	  /* 
	     inner Lev Mar iteration 
	  */
	  ierror = mrqmin(nr_x,nr_y,nr_sig,model.ndata,nr_a,nr_ia,
		       model.npar,nr_covar,nr_alpha,&nr_chi2,&nr_alamda,
		       (i_levmar == 0)?(LEVMAR_INIT):(LEVMAR_RUN));
	  if(ierror < 0){
	    fprintf(stderr,"%s: matrix should not be singular\n",argv[0]);
	    exit(-1);
	  }
	  i_levmar++;
	  if(nr_chi2 > min_chi2)
	    ok = FALSE;
	  else{
	    if((min_chi2 - nr_chi2)/min_chi2 < 1e-3)
	      ii_levmar++;		/* exit criterion */
	    if(ii_levmar > 2)
	      ok = TRUE;	/* exit criterion fulfilled a number of times */
	    min_chi2 = nr_chi2;
	  }
#ifdef DEBUG
	  fprintf(stderr,"%s: Lev Mar i %5i ii %5i alamda %10.2e chi2' %11g best chi2 %11g\n",
		  argv[0],i_levmar,ii_levmar,nr_alamda,nr_chi2,min_chi2);
#endif
	}while(!ok);
	/* 
	   converged, compute errors 
	*/
	ierror = mrqmin(nr_x,nr_y,nr_sig,model.ndata,nr_a,nr_ia,model.npar,
		     nr_covar,nr_alpha,&nr_chi2,&nr_alamda,LEVMAR_FINISH);
	if(ierror < 0){
	  fprintf(stderr,"%s: finalizing mrqmin should not return error\n",argv[0]);
	  exit(-1);
	}
	/* assign solution to fault structure */
	model.yoff = nr_a[1];
	for(j=2,i=0;i < model.nfault;i++,j+=3){
	  model.fault[i].xoff = nr_a[j];
	  model.fault[i].w    = nr_a[j+1];
	  model.fault[i].u    = nr_a[j+2];
	}
	/* sort faults by location */
	qsort(model.fault,model.nfault,sizeof(struct flt),sort_by_xoff);

	/* 

	   compute chi2 
	
	*/
	/* restore original data */
	if(use_orig_data_for_chi2)
	  memcpy(model.d,model.orig_d,sizeof(struct data)*model.ndata);
	/* assign to solution */
	model.solution[i_mc].chi2 = misfit(&model);
	/* store geometry and slip parameters */
	copy_fault_par_to_sol(&model,model.solution[i_mc].par,store_slip);
	/* 

	   end lev mar fit part

	 */
	if(i_mc == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: Lev Mar done: chi2: %11g \n",
		  argv[0],i_big,big_sample,i_mc,model.solution[i_mc].chi2);
	else if(i_mc % 500 == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: Lev Mar done: avg chi2: %11g ndchi2: %i\n",
		  argv[0],i_big,big_sample,i_mc,avg_chi2/(CPREC)i_mc,nchi2_close);
	/*  */
	if(inner_exit){
	  free_nr_vector(nr_x,1,model.ndata);free_nr_vector(nr_y,1,model.ndata);
	  free_nr_vector(nr_sig,1,model.ndata);free_ivector(nr_ia,1,model.npar);
	  free_nr_vector(nr_a,1,model.npar);
	  free_nr_matrix(nr_alpha,1,model.npar,1,model.npar);
	  free_nr_matrix(nr_covar,1,model.npar,1,model.npar);
	}
	break;
      case RANDOM_SOLVER:
	/* 
	 
	   random fit part 
	   
	*/
	last_chi2 = 0;
	/* reset min chi2 */
	min_chi2 = 1e20;
  	i_ran = 0;
	while(i_ran < ran_sample){
	
	  ii_ran = 0;		/* for positivity */
	  do{				/* iteration to have all u[] > umin */
	    /* 
	       assign all random xoff locations 
	    */
	    /* 
	       pick random fault locations
	       
	    */
	    for(modified=FALSE,i=0;i < model.nfault;i++){
	      if(!model.fault[i].xoff_fixed){
		ok = FALSE;
		iii_ran = 0;
		while((!ok) && (iii_ran < NPLIM)){
		  model.fault[i].xoff = off_min + off_range * ran2(&model.seed);
		  ok = TRUE;
		  for(j=0;j < model.nfault;j++)
		    if((j != i) && (fabs(model.fault[i].xoff - model.fault[j].xoff) < min_xoff)){
		      ok = FALSE;
		      break;
		    }
		  iii_ran++;
		}
		if(iii_ran > NPLIM){
		  fprintf(stderr,"%s: error, xoff loop exceeds limit of %i\n",
			  argv[0],NPLIM);
		  exit(-1);
		}
		modified = TRUE;
	      }
	    }
	    if(modified){
	      /* sort the faults by location */
	      qsort(model.fault,model.nfault,sizeof(struct flt),sort_by_xoff);
	    }
	    /* 
	       
	       assign all random locking  depths
	       
	    */
	    for(i=0;i < model.nfault;i++){	/* loop through number of tans */
	      if(!model.fault[i].w_fixed){
		/* 
		   random locking depth 
		*/
		model.fault[i].w = model.wmin + ran2(&model.seed) * model.wrange;	/*  */
	      }
	    }
	    /* 
	       
	       solve for the best fitting, non-negative slip values 
	       
	    */
	    solve_for_fault_slip(&model);
	    /* 

	     */
	    umin_loc = 1e20;
	    for(i=0;i < model.nfault;i++)
	      if(model.fault[i].u < umin_loc)
		umin_loc = model.fault[i].u;
	    ii_ran++;
	  }while((umin_loc < model.umin) && (ii_ran < iter_max_umin)); /* inner loop to check for 
								    umin solution */
	  if(ii_ran == iter_max_umin){
	    fprintf(stderr,"%s: ERROR: could not achieve umin %g in %i iterations\n",
		    argv[0],model.umin,iter_max_umin);
	    exit(-1);
	  }
	  /* 
	     
	     evaluate misfit and store the non-constrained part of the solution
	     
	  */

	  nr_chi2 = misfit(&model);
	  if(nr_chi2 < min_chi2){
	    /* store solution */
	    min_chi2 = nr_chi2;
	    model.solution[i_mc].chi2 = nr_chi2;
	    copy_fault_par_to_sol(&model,model.solution[i_mc].par,store_slip);
	  }
	  i_ran++;
#ifdef DEBUG
	  if(i_ran % 100000 == 0){
	    dchi2 = (min_chi2 - last_chi2)/min_chi2;
	    fprintf(stderr,"%s: iteration %10i out of %10i %5.1f%% done, min chi2: %12.5e dchi2: %12.5e\r",
		    argv[0],i_ran,ran_sample,(CPREC)i_ran/(CPREC)ran_sample*100,
		  min_chi2,dchi2);
	    last_chi2 = min_chi2;
	  }
#endif
	} /* end random loop */
	if(i_mc == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: random   done: chi2 %11g\n",
		  argv[0],i_big,big_sample,i_mc,model.solution[i_mc].chi2);
	else if(i_mc % 500 == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: random   done: avg chi2 %11g ndchi2: %5i\n",
		  argv[0],i_big,big_sample,i_mc,avg_chi2/(CPREC)i_mc,nchi2_close);
	if(use_orig_data_for_chi2){
	  /* compute misfit wrt. to original data */
	  memcpy(model.d,model.orig_d,sizeof(struct data)*model.ndata);
	}
	model.solution[i_mc].chi2 = misfit(&model);
	/* end random branch */
	break;
      case SIMPLEX_SOLVER:
	/* 

	   downhill simplex of nelder and meade
	   
	*/
	if(inner_init){
	  nr_alpha = nr_matrix(1,model.nsave_par+1,1,model.nsave_par);
	  nr_y = nr_vector(1,model.nsave_par+1);
	  nr_ia = ivector(1,model.nsave_par);
	}
	i_ran = 0;
	do{
	  /* initial guess */
	  for(i=0,j=0;i < model.nfault;i++){
	    if(!model.fault[i].xoff_fixed){
	      j++;
	      nr_alpha[1][j] = model.xmin + (model.xrange/(CPREC)(model.nfault+1)) * (CPREC)(i+1) + urandom(-20,20,&model);
	      nr_ia[j] = 1;	/* code for xoff */
	    }
	    if(!model.fault[i].w_fixed){
	      j++;
	      nr_alpha[1][j] = 7.5 + urandom(-5,5,&model);
	      nr_ia[j] = 0;	/* code for W */
	    }
	  }
	  if(j != model.nsave_par){
	    fprintf(stderr,"%s: amoeba logic error %i %i\n",argv[0],j,model.nsave_par);
	    exit(-1);
	  }
	  /* check if valid guess */
	  amo_cost_func(&(nr_alpha[1][1]),&model,&valid);
	  i_ran++;
	}while((!valid) && (i_ran < 100));
	if(i_ran == 100){
	  fprintf(stderr,"%s: could not find valid initial guess\n",argv[0]);
	  exit(-1);
	}

	if(i_mc == 0){
	  simplex_sample = 5 * model.nfault;	/* number of restarts */
	  simplex_eps = 3e-6;
	}else{
	  simplex_sample = 2 + model.nfault;
	  simplex_eps = 5e-6;
	}
	  
	i_ran = iii_ran = 0;
	do{			/* run a simplex a few times */
	  iiii_ran = 0;
	  do{
	    iiii_ran++;
	    iiiii_ran=0;
	    do{
	      /* perturbations from P0 */
	      for(i=1;i <= model.nsave_par;i++)
		for(j=1;j <= model.nsave_par;j++){
		  dist = ((nr_ia[j]==1)?(10):(2.5));
		  nr_alpha[i+1][j] = nr_alpha[1][j] + ((i==j)?(urandom(-dist,dist,&model)):(0));
		  if(nr_ia[j]){
		    if(nr_alpha[i+1][j] < model.xmin)
		      nr_alpha[i+1][j] = model.xmin + 10;
		    if(nr_alpha[i+1][j] > model.xmax)
		      nr_alpha[i+1][j] = model.xmax - 10;
		  }else{
		    if(nr_alpha[i+1][j] < model.wmin)
		      nr_alpha[i+1][j] = model.wmin + 0.5;
		    if(nr_alpha[i+1][j] > model.wmax)
		      nr_alpha[i+1][j] = model.wmax - 0.5;
		  }
		}
	      /* misfits */
	      for(i=1;i <= model.nsave_par+1;i++){
		nr_y[i] = amo_cost_func(&(nr_alpha[i][1]),&model,&valid);
		if(!valid)
		  break;
	      }
	      iiiii_ran++;
	    }while((!valid) && (iiiii_ran <100));
	    if(iiiii_ran == 100){
	      fprintf(stderr,"%s: errror, could not find random starting sequences\n",argv[0]);
	      exit(-1);
	    }
	    /* make sure amoeba converges and yields a valid
	       solution*/
	    amoeba_return = amoeba(nr_alpha,nr_y,model.nsave_par,simplex_eps,&ii_ran,&model);
	    valid = valid_solution(&model);
	    //fprintf(stderr,"iiii %3i amoeba %i valid %i\n",iiii_ran,amoeba_return,valid);
	  }while((iiii_ran < 500) && ((amoeba_return != 0)||(!valid)));

	  if(iiii_ran == 500){
	    fprintf(stderr,"%s: amoeba inner loop problem, no proper solution found within 500 steps\n",argv[0]);
	    exit(-1);
	  }
	  nr_chi2 = misfit(&model);
#ifdef DEBUG
	  fprintf(stderr,"%s: amoeba run %3i, avg %5i iterations, chi2 %11g\n",
		  argv[0],i_ran,ii_ran,nr_chi2);
#endif
	  iii_ran += ii_ran;
	  i_ran++;
	}while(i_ran < simplex_sample);
	iii_ran = ((CPREC)iii_ran/(CPREC)simplex_sample + 0.5);
	/* save solution */
	qsort(model.fault,model.nfault,sizeof(struct flt),sort_by_xoff);
	if(use_orig_data_for_chi2){
	  /* restore original data */
	  memcpy(model.d,model.orig_d,sizeof(struct data)*model.ndata);
	}
	model.solution[i_mc].chi2 = misfit(&model);
	copy_fault_par_to_sol(&model,model.solution[i_mc].par,store_slip);
	if(i_mc == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: simplex avg iter: %5i chi2 %11g\n",
		  argv[0],i_big,big_sample,i_mc,iii_ran,model.solution[i_mc].chi2);
	else if(i_mc % 500 == 0)
	  fprintf(stderr,"%s: i_big: %i/%i i_mc: %5i: simplex avg iter: %5i avg chi2 %11g ndchi2: %5i\n",
		  argv[0],i_big,big_sample,i_mc,iii_ran,avg_chi2/(CPREC)i_mc,nchi2_close);

	if(inner_exit){
	  free_nr_vector(nr_y,1,model.nsave_par+1);
	  free_nr_matrix(nr_alpha,1,model.nsave_par+1,1,model.nsave_par);
	}
	break;
	/* end simplex */
      default:
	fprintf(stderr,"%s: solver mode %i undefined\n",argv[0],solver_mode);
	exit(-1);
      }


      if(i_mc == 0){
	/* original solution */
	if(solver_mode == LEVMAR_SOLVER){
	  /* store estimates of error */
	  for(i=1;i<= model.npar;i++)
	    for(j=1; j <= model.npar;j++)
	      best_nr_covar[i][j] = nr_covar[i][j];
	}
	/* 
	   compute del chi2 stats
	*/
	del_chi2 = delta_chi2_rtbis((CPREC)model.nfree_par,prob_del_chi2_inner);
	fprintf(stderr,"%s: min chi2: %11g del npara(dof): %4i prob: %.2f%% del chi2: %11g\n",
		argv[0],model.solution[0].chi2,model.nfree_par,prob_del_chi2_inner*100,
		del_chi2);
	chi2_lim = model.solution[0].chi2 + del_chi2;
      }
      /* avg chi2 */
      avg_chi2 += model.solution[i_mc].chi2;
      if(model.solution[i_mc].chi2 < chi2_lim)
	nchi2_close++;
    } /* end MC iteration */
    
    /* restore data to original data  */
    memcpy(model.d,model.orig_d,sizeof(struct data)*model.ndata);
    i_big++;
  }while(i_big < big_sample); /* big outer loop, if we want to redo
				 the solution after outliers */

  mc_sample = i_mc;		/* in case of early exit */
  
  /* 
     
     print the orirginal solution set of parameters to stdout

  */
  print_best_result(&model,stdout,0,store_slip,(solver_mode == LEVMAR_SOLVER),best_nr_covar);
  /* 

     total and individual fitting functions 

  */
  out = fopen("invert.fit.out","w"); /* total */
  print_best_result(&model,out,1,store_slip,(solver_mode == LEVMAR_SOLVER),best_nr_covar);	/* fitting
												   functions
												   (for
												   whole
												   model
												   and
												   individual)
												   to
												   file  */
  fclose(out);
  if(solver_mode == LEVMAR_SOLVER){
    /*
      
      output of error bounds and solution 
      
    */
    out = fopen("invert.lmcov.out","w"); /* covariance */
    for(i=1;i<= model.npar;i++){
      for(j=1; j <= model.npar;j++)
	fprintf(out,"%11g ",best_nr_covar[i][j]);
      fprintf(out,"\n");
    }
    fclose(out);
  }
  /* 

     actually used data, could be subset of all

  */
  out = fopen("invert.dat.out","w");
  print_best_result(&model,out,2,store_slip,(solver_mode == LEVMAR_SOLVER),best_nr_covar);	/* data used for fitting to file */
  fclose(out);

  /* 

     print MC solutions
     
     format:
     
     1     2          3          4         5          6           7        8 
     chi2 yoff fault_1_xoff fault_1_D fault_1_u fault_2_xoff fault_2_D fault_2_u  ....
     
     
  */
  out = fopen("invert.bestsol.dat","w");
  for(i=0,np=0;i < mc_sample;i++){
    fprintf(out,"%16.7e\t",model.solution[i].chi2);
    copy_sol_to_fault_par(model.solution[i].par,&model,store_slip);
    fault_yoff_to_para(model.fault,model.nfault,model.yoff,&apara,&npara);
    for(j=0;j < npara;j++)
      fprintf(out,"%16.7e ",apara[j]);
    fprintf(out,"\n");
  }
  fprintf(stderr,"%s: printed %i MC simulation parameters\n",argv[0],mc_sample);
  fclose(out);
  
  if(solver_mode == LEVMAR_SOLVER)
    free_nr_matrix(best_nr_covar,1,model.npar,1,model.npar);
  free(apara);


  fprintf(stderr,"%s: done\n",argv[0]);
  return 0;
}

