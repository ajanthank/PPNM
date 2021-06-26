CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt ExerciseC.pdf


out.txt: main
	./$< > $@

ExerciseC.pdf: Time.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N^3";\
	set ylabel "Time[s]";\
	set title "CPU time";\
	plot \
	"$<" using 1:2 with points title "My Time"\
	,"$<" using 1:3 with points title "Gsl Time"\
	'|gnuplot

.PHONEY: clean

clean:
	$(RM) *.o main *.txt ExerciseC.pdf *~
