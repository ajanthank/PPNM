CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt Plot.pdf

Plot.pdf: fitdat.txt Makefile
	echo '\
	set terminal pdf;\
	set out "$@";\
	set xlabel "Time t (days)";\
	set ylabel "log(activity) arb. u.";\
	set tics out;\
	set key out;\
	set title "L-S fit to Rutherford/Soddy data" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "Fit"\
	,"$<" index 1 using 1:3 with lines title "+∆c_0"\
	,"$<" index 1 using 1:4 with lines title "-∆c_0"\
	,"$<" index 2 using 1:3 with lines title "+∆c_1"\
	,"$<" index 2 using 1:4 with lines title "-∆c_1"\
	'|gnuplot

out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *~ *.pdf
