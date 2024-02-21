include Makefile.inc

all: tomo2D_Frechet

INCLUDE = -I../INCLUDE 

SUB1 = lsmrDataModule.o lsmrReverse.o 

SUB = lsmrblas.o read_src_rec.o read_data.o read_model.o ray_over_data.o ray_sngl_data_src_rec.o ray_segment.o 

LIB = 


%.o : %.f90
	$(FC) $(OPTF) $(INCLUDE) -c -o $@ $<

clean :
	rm -f *.o *~ *backup *.mod *genmod*

tomo2D_Frechet : $(SUB1) $(SUB) tomo2D_Frechet.o
	$(FL) $(SUB1) $(SUB) tomo2D_Frechet.o -o ../BIN/tomo2D_Frechet $(LIB)
