#
# compile the fault slip inversion code used for Platt and Becker (2010)
#

#F77 = /usr/geosys/x86_64/intel/fce/10.1.015/bin/ifort
#CC = /usr/geosys/x86_64/intel/cce/10.1.015/bin/icc
#FLAGS = -O3 -ipo -no-prec-div -static -DSINGLE_PREC
#F77 = ifort
#CC = icc
CC = gcc
#
# double precision
#FLAGS = -O2 
#LDFLAGS = -L$(HOME)/progs/src/slatec/$(ARCH)/ -lslatec_dbl $(MATHLIBS) -lm
#
# single precision
FLAGS = -O2 -DSINGLE_PREC
LDFLAGS = -L$(HOME)/progs/src/slatec/$(ARCH)/ -lslatec_sgl $(MATHLIBS) -lm

# DEBUG
#
#FLAGS = -g -ftrapuv -inline-debug-info -debug all  -DSINGLE_PREC -DDEBUG
#LDFLAGS = -L$(HOME)/progs/src/slatec/$(ARCH)/ -lslatec_sgl $(MATHLIBS) -lm


FFLAGS = $(FLAGS) -nofor_main
CFLAGS = $(FLAGS) 

ODIR = objects/
BDIR = bin/


INCLUDES = -I$(HOME)/progs/src/slatec/ 
HDR_FLS = invert.h precision.h deltachi2.h nrutil.h 

all:	$(BDIR)/invert $(BDIR)/prune_data $(BDIR)/test_minv

clean:
	rm $(BDIR)/invert;\
	rm $(ODIR)/*.o


$(BDIR)/invert: $(ODIR)/invert.o $(ODIR)/nrutil.o $(ODIR)/invert_util.o \
	$(ODIR)/deltachi2.o \
	$(HOME)/progs/src/slatec/$(ARCH)/libslatec_sgl.a
	$(F77) $(FFLAGS) $(ODIR)/nrutil.o  $(ODIR)/deltachi2.o \
	$(ODIR)/invert.o  $(ODIR)/invert_util.o -o $(BDIR)/invert $(LDFLAGS)

$(BDIR)/prune_data: prune_data.c 
	$(CC) $(CFLAGS) prune_data.c  -o $(BDIR)/prune_data $(LDFLAGS)

$(BDIR)/test_minv: $(ODIR)/test_minv.o $(ODIR)/nrutil.o $(ODIR)/invert_util.o \
	$(ODIR)/deltachi2.o \
	$(HOME)/progs/src/slatec/$(ARCH)/libslatec_sgl.a
	$(F77) $(FFLAGS) $(ODIR)/nrutil.o  $(ODIR)/deltachi2.o \
	$(ODIR)/test_minv.o  $(ODIR)/invert_util.o -o $(BDIR)/test_minv $(LDFLAGS)





$(ODIR)/%.o: %.F   $(HDR_FLS)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.o


$(ODIR)/%.o: %.f90  $(HDR_FLS)
	$(F90) $(FFLAGS) -c $< -o $(ODIR)/$*.o


$(ODIR)/%.o: %.f  $(HDR_FLS)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.o

$(ODIR)/invert.o: invert.c  $(HDR_FLS)
	$(CC) $(CFLAGS) $(INCLUDES) -c invert.c -o $(ODIR)/invert.o



$(ODIR)/%.o: %.c  $(HDR_FLS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

