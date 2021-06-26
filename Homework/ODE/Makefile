CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt HarmO.pdf S.pdf I.pdf R.pdf

HarmO.pdf: Plotdat.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "t";\
	set ylabel "y";\
	set title "Harmonic Oscillator";\
	plot \
	"$<" index 0 using 1:2 with lines title "Integration"\
	,"$<" index 0 using 1:3 with lines title "Exact"\
	'|gnuplot

S.pdf: Plotdat.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "Time[days]";\
	set ylabel "Number of People";\
	set title "Susceptible Individuals";\
	plot \
	"$<" index 1 using 1:2 with lines title "T_C = 2/3 per day"\
	,"$<" index 2 using 1:2 with lines title "T_C = 1 per day"\
	,"$<" index 3 using 1:2 with lines title "T_C = 2 per day"\
	'|gnuplot

I.pdf: Plotdat.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "Time[days]";\
	set ylabel "Number of People";\
	set title "Infectious Individuals";\
	plot \
	"$<" index 1 using 1:3 with lines title "T_C = 2/3 per day"\
	,"$<" index 2 using 1:3 with lines title "T_C = 1 per day"\
	,"$<" index 3 using 1:3 with lines title "T_C = 2 per day"\
	'|gnuplot

R.pdf: Plotdat.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "Time[days]";\
	set ylabel "Number of People";\
	set title "Removed Individuals";\
	plot \
	"$<" index 1 using 1:4 with lines title "T_C = 2/3 per day"\
	,"$<" index 2 using 1:4 with lines title "(I)T_C = 1 per day"\
	,"$<" index 3 using 1:4 with lines title "(R)T_C = 2 per day"\
	'|gnuplot

out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~
