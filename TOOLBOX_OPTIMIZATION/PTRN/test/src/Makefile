include ../../../Makefile.inc

CXX = $(FC)
FLAGS = $(OPTF)

##SOURCES 

COMMON_SRC = ../../../COMMON/test/
CXX_SOURCES = $(COMMON_SRC)rosenbrock.f90 test_PTRN.f90 
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

INCLUDE = -I ../../../COMMON/include
TARGET = ../bin/test_PTRN.bin 

#LIB =  -L ../../lib -lPTRN
LIB =  -L ../../../lib -lSEISCOPE_OPTIM


$(TARGET) : $(CXX_OBJECTS) 
	$(CXX) $(FLAGS) $(CXX_OBJECTS) -o $(TARGET) $(LIB)

%.o : %.f90
	$(CXX) $(FLAGS) $(INCLUDE) -c -o $@ $<

clean :
	rm -f *.o *~ *backup *.mod *genmod*
	rm -f $(COMMON_SRC)*.o $(COMMON_SRC)*~ $(COMMON_SRC)*backup $(COMMON_SRC)*.mod $(COMMON_SRC)*genmod*
