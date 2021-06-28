CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt HarmO.pdf

out.txt: main
	./$< 1> $@ 2>log

HarmO.pdf: Plotdat.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "t";\
	set xrange [0:3];\
	set ylabel "y";\
	set title "Harmonic Oscillator";\
	plot \
	"$<" index 0 using 1:2 with lines title "Integration"\
	,"$<" index 0 using 1:3 with lines title "Exact"\
	'|gnuplot


.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~ log*
