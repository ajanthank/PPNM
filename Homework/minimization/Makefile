CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt Higgs.pdf

Higgs.pdf: out.txt Makefile
	echo '\
	set terminal pdf;\
	set out "$@";\
	set xlabel "E[GeV]";\
	set ylabel "sigma";\
	set tics out;\
	set key out;\
	set title "Higgs peak" ;\
	plot "$<" index 1 using 1:4 with lines title "Fit"\
	,"$<" index 1 using 1:2:3 with yerrorbars title "Data"\
	'|gnuplot

out.txt: main
	./$< 1> $@ 2>log

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~ log*
