#OPTION FOR Makefile defined here
FC      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -Wall -fbacktrace -finit-real=nan
FC      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -Wall -fbacktrace -finit-real=nan

FC      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -fbacktrace -finit-real=nan
FC      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -fbacktrace 

FL = gfortran

#OPTIONS COMPILE

OPTF =  -O3   
OPTC =  -O3  
OPTL =  -O3   

OPTFF = ${OPTF}

########################## implicit rules ###################


.SUFFIXES: .f .o
.f.o:
	$(F77) $(OPTF) $(INC) -c $*.f
	 
.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) $(INC) -c $*.f90
	 
.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -c $*.c

########################## archivage #########################
AREXTRA= vx
ARFUNCT= cruv
AR= ar
RANLIB= ranlib
