CFLAGS = $(shellA gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: plot.pdf plot2.pdf ms.pdf

plot2.pdf : out.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "f(x)";\
	set title "Exponential function";\
	plot \
	"$<" using 1:2 with line title "math.h"\
	,"$<" using 1:3 with line title "Approximation"\
	'|gnuplot

plot.pdf : out.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "f(x)";\
	set title "Exponential function";\
	plot \
	"$<" using 1:2 with line title "math.h"\
	,"$<" using 1:4 with line title "Approximation"\
	'|gnuplot

out.txt: main
	./$< >$@

ms.pdf: ms.tex
	pdflatex $<
	pdflatex $<

.PHONEY: clean

clean:
	$(RM) *.txt *.pdf main *~ *.aux *.out *.log
