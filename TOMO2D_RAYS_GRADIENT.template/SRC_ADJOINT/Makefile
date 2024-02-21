include Makefile.inc

all: tomo2D_gradient

INCLUDE = -I../INCLUDE -I ../../TOOLBOX_OPTIMIZATION/COMMON/include

LIB =  -L ../../TOOLBOX_OPTIMIZATION/lib -lSEISCOPE_OPTIM

SUB = read_src_rec.o read_data.o read_model.o ray_over_data.o ray_sngl_data_src_rec.o ray_segment.o 

%.o : %.f90
	$(FC) $(OPTF) $(INCLUDE) -c -o $@ $<

clean :
	rm -f *.o *~ *backup *.mod *genmod*

tomo2D_gradient : $(SUB) tomo2D_gradient.o
	$(FL) $(SUB) tomo2D_gradient.o -o ../BIN/tomo2D_gradient $(LIB)
