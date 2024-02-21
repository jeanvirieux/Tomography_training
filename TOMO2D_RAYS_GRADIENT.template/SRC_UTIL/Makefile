include ../SRC_ADJOINT/Makefile.inc


all: line_array model_homo temps_homo model_hetero

INCLUDE = -I../INCLUDE


%.o : %.f90
	$(FC) $(OPTF) $(INCLUDE) -c -o $@ $<

clean :
	rm -f *.o *~ *backup *.mod *genmod*

line_array : line_array.o
	$(FL) $(SUB) line_array.o -o ../BIN/line_array

model_homo : model_homo.o
	$(FL) $(SUB) model_homo.o -o ../BIN/model_homo

model_hetero : model_hetero.o
	$(FL) $(SUB) model_hetero.o -o ../BIN/model_hetero

temps_homo : temps_homo.o ../SRC_ADJOINT/read_src_rec.o ../SRC_ADJOINT/read_model.o
	$(FL) $(SUB) ../SRC_ADJOINT/read_src_rec.o ../SRC_ADJOINT/read_model.o temps_homo.o -o ../BIN/temps_homo
