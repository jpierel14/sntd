#
# makefile for producing a runable version of MICROLENS 
#              Joachim Wambsganss, January 2008
#
#FC = gfortran -O3 -fdefault-integer-8  -c
FC = gfortran -O3 -c -ffixed-line-length-0
microlens: main.o detko.o input.o massf.o niveau.o \
	numstar.o output.o pos.o poslc.o rands.o raycount.o  \
	setcell.o setmult.o setstar.o  shootcl.o shootfix.o shotayco.o \
	sort.o  cputime.o
	gfortran  -O3 -o microlens detko.o main.o niveau.o pos.o poslc.o \
	rands.o  setcell.o setmult.o setstar.o shootcl.o shootfix.o \
	shotayco.o raycount.o input.o output.o numstar.o sort.o \
	massf.o   cputime.o \
	-ffixed-line-length-132  -std=gnu
#	massf.o   cputime.o  -L.  -lm
#	massf.o -lveclib
#	    
main.o:  main.f
	$(FC) main.f
#	    
detko.o:  detko.f
	$(FC) detko.f
#	    
input.o:  input.f
	$(FC)  input.f
#	    
massf.o:  massf.f
	$(FC) massf.f
#	    
niveau.o:  niveau.f
	$(FC) niveau.f
#	    
numstar.o:  numstar.f
	$(FC) numstar.f
#	    
output.o:  output.f
	$(FC) output.f
#	    
pos.o:  pos.f
	$(FC) pos.f
#	    
poslc.o:  poslc.f
	$(FC) poslc.f
#	    
rands.o:  rands.f
	$(FC) rands.f
#	    
raycount.o:  raycount.f
	$(FC) raycount.f
#	    
setcell.o:  setcell.f
	$(FC) setcell.f
#	    
setmult.o:  setmult.f
	$(FC) setmult.f
#	    
setstar.o:  setstar.f
	$(FC) setstar.f
#	    
shootcl.o:  shootcl.f
	$(FC) shootcl.f
#	    
shootfix.o:  shootfix.f
	$(FC) shootfix.f
#	    
shotayco.o:  shotayco.f
	$(FC)  shotayco.f
#
sort.o:  sort.f
	$(FC) sort.f
#
cputime.o:  cputime.c
	cc  -c cputime.c
