CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc
.PHONEY: default

defualt: out.txt erfplot.png

out.txt: main
	./$< >$@
erfplot.png:DataB.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "Erf(x)";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf"\
	'|gnuplot

.PHONEY: clean

clean:
	$(RM) main *.txt *png
