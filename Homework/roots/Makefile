CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt SE.pdf

SE.pdf: out.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "r[a_0]";\
	set ylabel "Psi(r)";\
	set title "Wavefuntion Solution";\
	plot \
	"$<" index 1 using 1:2 with lines title "Theoretical"\
	,"$<" index 1 using 1:3 with lines title "Calculated"\
	'|gnuplot

out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *~ *.pdf
