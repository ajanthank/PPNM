CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt Error.pdf

Error.pdf: Data.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set title "Error scaling";\
	plot \
	"$<" using 1:2 with lines title "Pseudo-random"\
	,"$<" using 1:3 with lines title "Quasi-random"\
	'|gnuplot

out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~
