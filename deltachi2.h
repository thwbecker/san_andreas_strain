#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "precision.h"
#include "nrutil.h"

CPREC delta_chi2_rtbis(CPREC, CPREC);
CPREC rfunc(CPREC, CPREC, CPREC);
CPREC gammq(CPREC, CPREC);

void gser(CPREC *, CPREC, CPREC, CPREC *);

void gcf(CPREC *, CPREC, CPREC, CPREC *);
CPREC gammln(CPREC);

