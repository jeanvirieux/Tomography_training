include ../../Makefile.inc

CXX = $(FC)
FLAGS = $(OPTF)


##SOURCES 
COM = /src/

CXX_SOURCES_common = $(COM)normL2.f90 $(COM)scalL2.f90 $(COM)std_linesearch.f90 $(COM)std_test_conv.f90 $(COM)print_info.f90 $(COM)project.f90

CXX_SOURCES = $(CXX_SOURCES_common) 
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

INCLUDE = -I ../COMMON/include

%.o : %.f90
	$(CXX) $(FLAGS) $(INCLUDE) -c -o $@ $<

donothing : # for giving the information when typing make
	@echo please type -make lib- for compiling the library

lib : $(CXX_OBJECTS) 

#$(AR) $(ARFUNCT) $(LIB)libPSTD.a $(CXX_OBJECTS)
#cp $(LIB)libPSTD.a ../lib

clean :	
	rm -f $(COM)*.o $(COM)*~ $(COM)*backup $(COM)*.mod $(COM)*genmod* 
	rm -f *.mod *genmod* 