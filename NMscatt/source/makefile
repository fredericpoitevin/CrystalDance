# makefile for NMscatt

#FC = g77
FC = gfortran

FFLAGS = -mcmodel=medium -funroll-all-loops -O3

obj1 = phonon.o
lib  = -lblas -llapack

obj2 = incoh.o

obj3 = coh.o

obj4 = bead.o

all:
	make phonon
	make incoh
	make coh
	make bead

phonon:	$(obj1)
	$(FC) $(FFLAGS) $(obj1) $(lib) -o phonon

incoh:	$(obj2)
	$(FC) $(FFLAGS) $(obj2) -o incoh

coh:	$(obj3)
	$(FC) $(FFLAGS) $(obj3) -o coh

bead:	$(obj4)
	$(FC) $(FFLAGS) $(obj4) -o bead
