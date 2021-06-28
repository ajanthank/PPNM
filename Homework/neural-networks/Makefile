CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt Plot.pdf

Plot.pdf: out.txt Makefile
	echo '\
	set terminal pdf;\
	set out "$@";\
	set xlabel "x";\
	set xrange [-2:2];\
	set ylabel "f(x)";\
	set tics out;\
	set key out;\
	set title "Ann Curvefitting" ;\
	plot\
	"$<" index 0 using 1:2 with points pointtype 6 title "Data"\
	,"$<" index 1 using 1:2 with lines title "ANN exp(-x*x)"\
	,"$<" index 1 using 1:3 with lines title "Exact exp(-x*x)"\
	,"$<" index 1 using 1:4 with lines title "Ann d/dx exp(-x*x)"\
	,"$<" index 1 using 1:5 with lines title "Exact d/dx exp(-x*x)"\
	,"$<" index 1 using 1:6 with lines title "Ann int dx exp(-x*x)"\
	,"$<" index 1 using 1:7 with lines title "Exact int dx exp(-x*x)"\
	'|gnuplot

out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~
