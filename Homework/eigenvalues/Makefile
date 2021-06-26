CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt Energy.pdf Eigenfunc.pdf

Energy.pdf: Energy.txt Makefile
	echo '\
	set terminal pdf;\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Particle in a box: Eigenvelues" ;\
	plot \
	"$<" using 1:2 with lines title "Calculated"\
	,"$<" using 1:3 with lines title "Exact"\
	'|gnuplot


Eigenfunc.pdf: Eigenfunc.txt Makefile
	echo '\
	set terminal pdf;\
	set out "$@";\
	set xlabel "Time t (days)";\
	set ylabel "log(activity) arb. u.";\
	set tics out;\
	set key out;\
	set title "Particle in box: Eigenfunctions" ;\
	plot \
	"$<" index 0 using 1:2 with lines title "0."\
	,"$<" index 1 using 1:2 with lines title "1."\
	,"$<" index 2 using 1:2 with lines title "2."\
	,"$<" index 3 using 1:2 with lines title "3."\
	'|gnuplot


out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~
