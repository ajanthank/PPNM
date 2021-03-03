CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc
.PHONEY: default

defualt: out.txt erfplot.png Besselplot.png

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

Besselplot.png: DataBessel.txt Makefile
	echo '\
	set terminal png;\
        set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "J_0 (x)";\
	set title "Bessel function";\
	plot \
	"$<" using 1:2 with line title "J_0"\
	'|gnuplot

.PHONEY: clean

clean:
	$(RM) main *.txt *png
