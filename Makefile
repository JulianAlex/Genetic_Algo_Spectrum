# Makefile GeneticFit-Spektrum
#
#CC_FLAGS = -g -ggdb -Wall -m64 -O3     #-pg -Wall -m64 -O3
#
CC_FLAGS = -Wall   

LD_FLAGS =  -lfftw3 -lm

OBJ = fit_main.o spek_correlJw.o spek_dynTheory.o spek_gammaTau.o spek_input.o\
 spek_lineshape.o spek_omegaShift.o spek_output.o spek_spectra.o fit_subworx.o fit_fitnessFkt.o fit_geneticOp.o\
 eispack.o normal.o sorting.o nrutil.o\

SRC = fit_main.c spek_correlJw.c spek_dynTheory.c spek_gammaTau.c spek_input.c\
 spek_lineshape.c spek_omegaShift.c spek_output.c spek_spectra.c fit_subworx.c fit_fitnessFkt.c fit_geneticOp.c\
 eispack.c normal.c sorting.c nrutil.c\

HDR = spek_funk.h spek_head.h fit_head.h
 

CC = gcc     # Compiler

EXECUTABLE = fit



$(EXECUTABLE) : $(SRC) $(HDR)                 #$(OBJ)
#
	$(CC) $(CC_FLAGS) -o fit $(SRC) $(LD_FLAGS)
#
#
##############################################################################
#
# GPROF:
#
# "gcc -g -Wall -o spektrum -pg -lm spek.c -pg"
#
# compile and run program ( => gprof.out )
#
# run gprof:  "gprof ./spektrum > temp.out"
#

