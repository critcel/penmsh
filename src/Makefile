# penmshxp Makefile
# prerequisite:
#   DISLIN graphic library  

PROGRAM = penmshxp

SRCS = subcode0.f90 mytecio.f90  dislin.f90 subcode1.f90 \
       subcode2.f90 subcode3.f90 subcode4.f90 subcode5.f90 \
       subcode6.f90 subcode7.f90 subcode8.f90 supplement01.f90 \
       supplement02.f90 inpred.f90 
OBJS = $(SRCS:.f90=.o)


FC = gfortran
FCFLAG = -O3 

# gfortran option to accept longer lines
ifeq ($(FC), gfortran)
FCFLAG += -ffree-line-length-0
endif

# DISLIN library
.PHONY: check-dislin

LDFLAG = -L$(DISLIN) -ldislin

all: $(PROGRAM)
check-dislin: 
ifndef DISLIN
	$(error DISLIN is requied to build this version of $(PROGRAM) )
endif

$(PROGRAM):  check-dislin $(OBJS)
	$(FC) $(FCFLAG) $(OBJS) -o $(PROGRAM)  $(LDFLAG) 
%.o: %.f90 
	$(FC) -c $(FCFLAG) $< -o $@
clean:
	rm -f $(OBJS)
	rm -f *.mod
clear: clean
	rm -f $(PROGRAM)
