#ifndef __PREC_READ__
#include "slatec_src.h"

#ifdef SINGLE_PREC
#define nr_vector vector
#define nr_matrix matrix
#define free_nr_vector free_vector
#define free_nr_matrix free_matrix

#define FMT_LINE1 "%f %f %f"
#define FMT_1 "%f"
#define EPS 1.2e-7

#else
#define nr_vector dvector
#define nr_matrix dmatrix
#define free_nr_vector free_dvector
#define free_nr_matrix free_dmatrix

#define FMT_LINE1 "%lf %lf %lf"
#define FMT_1 "%lf"
#define EPS 4e-15

#endif
#define CPREC COMP_PRECISION	/* from slatec */

#define __PREC_READ__
#endif
