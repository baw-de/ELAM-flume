# Compiler choice (icpc -g, icc, g++, ...)
# flag explanations below / -g
CC = icpc -g
# C++11 (ISO/IEC 14882:2011)
STD = -std=c++11

# Compile-time flags
# Profiling: -profile-functions -profile-loops=all -O1
# Optimization: -O1..3 - O2: no effect / F90: -g
CFLAGS_F90 = -c -fixed -g
CFLAGS_C++ = -c $(STD) -DWM_DP -diag-disable 525 -DWM_LABEL_SIZE=32

# Include directories to search for header files (-I)
# HPC automatix
FOAM_DIR = /automatix/sw/apps/OpenFOAM/OpenFOAM-4.1/OpenFOAM-4.1/src

# Include paths of header files
EXE_INC = \
	-I$(FOAM_DIR)/OpenFOAM/lnInclude \
	-I$(FOAM_DIR)/finiteVolume/lnInclude \
	-I$(FOAM_DIR)/OSspecific/POSIX/lnInclude \
	-I$(FOAM_DIR)/triSurface/lnInclude \
	-I$(FOAM_DIR)/meshTools/lnInclude

# Library path(s) and libraries
EXE_LIBS = \
	-L$(FOAM_LIBBIN) \
	-lOpenFOAM \
	-lfiniteVolume \
	-lfvOptions \
    -lmeshTools \
    -lincompressibleTurbulenceModels

# RULE:
# target: dependencies
#     command(s)

.PHONY: all
all: ELAM-flume

.PHONY: ELAM-flume
ELAM-flume: random.o writeOutput.o vectorRelation.o BehaviorRule.o hydroInterpolation.o sensoryPointCreate.o updateFishLocation.o ELAM-main.o
	$(CC) -o ELAM-flume $(STD) random.o writeOutput.o vectorRelation.o BehaviorRule.o hydroInterpolation.o sensoryPointCreate.o updateFishLocation.o ELAM-main.o $(EXE_INC) $(EXE_LIBS) -lifcore
# -lifcore, -lifport include Intel modules (file handling) for linking - automatically in compiling below

random.o: random.f90
	ifort $(CFLAGS_F90) random.f90

writeOutput.o: writeOutput.f90
	ifort $(CFLAGS_F90) writeOutput.f90

vectorRelation.o: vectorRelation.f90
	ifort $(CFLAGS_F90) vectorRelation.f90

BehaviorRule.o: BehaviorRule.f90
	ifort $(CFLAGS_F90) BehaviorRule.f90

hydroInterpolation.o: hydroInterpolation.cpp hydroInterpolation.h
	$(CC) $(CFLAGS_C++) -DNoRepository hydroInterpolation.cpp $(EXE_INC) $(EXE_LIBS)

sensoryPointCreate.o: sensoryPointCreate.cpp sensoryPointCreate.h
	$(CC) $(CFLAGS_C++) sensoryPointCreate.cpp

updateFishLocation.o: updateFishLocation.cpp
	$(CC) $(CFLAGS_C++) updateFishLocation.cpp $(EXE_INC) $(EXE_LIBS)

ELAM-main.o: ELAM-main.cpp
	$(CC) $(CFLAGS_C++) ELAM-main.cpp $(EXE_INC) $(EXE_LIBS)

# full test for master branch commits - use ASCII instead of binary to enable error searching
test:
	time ./ELAM-flume -case ../eh-rinne/07_lang_pfosten_slot_cut_OF4/ &> output/ELAM-log.txt
	@diff -s ./output/v_TecTrack_ZonesAreTime.dat ./output_ref_401-5_atx/v_TecTrack_ZonesAreTime.dat | head -n 10

# plain run without test
run:
	./ELAM-flume -case ../eh-rinne/07_lang_pfosten_slot_cut_OF4/

clean:
	rm -rf *.o *.so ELAM-flume
