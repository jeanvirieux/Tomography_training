#MASTER Makefile for the OPTIMIZATION
#include ../Makefile.inc
include Makefile.inc

LIB = ./lib/
COM = ./COMMON/src/
KER_PSTD = ./PSTD/kernel/src/
KER_PNLCG = ./PNLCG/kernel/src/
KER_LBFGS = ./LBFGS/kernel/src/
KER_PLBFGS = ./PLBFGS/kernel/src/
KER_TRN = ./TRN/kernel/src/
KER_PTRN = ./PTRN/kernel/src/

TEST_PSTD = ./PSTD/test/
TEST_PNLCG = ./PNLCG/test/
TEST_LBFGS = ./LBFGS/test/
TEST_PLBFGS = ./PLBFGS/test/
TEST_TRN = ./TRN/test/
TEST_PTRN = ./PTRN/test/


all : compile

lib : compile

compile : 
	cd PNLCG  ; $(MAKE) lib
	cd PSTD   ; $(MAKE) lib
	cd PTRN   ; $(MAKE) lib
	cd TRN    ; $(MAKE) lib
	cd PLBFGS ; $(MAKE) lib
	cd LBFGS  ; $(MAKE) lib
	ar rcvs $(LIB)libSEISCOPE_OPTIM.a $(COM)*.o $(KER_PSTD)*.o $(KER_PNLCG)*.o $(KER_LBFGS)*.o $(KER_PLBFGS)*.o $(KER_TRN)*.o $(KER_PTRN)*.o

clean :	
	cd PNLCG  ; $(MAKE) clean
	cd PSTD   ; $(MAKE) clean
	cd PTRN   ; $(MAKE) clean
	cd TRN    ; $(MAKE) clean
	cd PLBFGS ; $(MAKE) clean
	cd LBFGS  ; $(MAKE) clean
	rm $(LIB)*.a

superclean :
	cd PNLCG  ; $(MAKE) clean
	cd PSTD   ; $(MAKE) clean
	cd PTRN   ; $(MAKE) clean
	cd TRN    ; $(MAKE) clean
	cd PLBFGS ; $(MAKE) clean
	cd LBFGS  ; $(MAKE) clean
	rm $(LIB)*.a
	cd $(TEST_PSTD)src ; $(MAKE) clean
	cd $(TEST_PNLCG)src ; $(MAKE) clean
	cd $(TEST_LBFGS)src ; $(MAKE) clean
	cd $(TEST_PLBFGS)src ; $(MAKE) clean
	cd $(TEST_TRN)src ; $(MAKE) clean
	cd $(TEST_PTRN)src ; $(MAKE) clean



