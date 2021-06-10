CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt Lin.pdf LinInteg.pdf

Lin.pdf: XY.txt LinData.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "cos(x)";\
	set title "Linear Interpolation";\
	f(x)=cos(x);\
	plot \
	"$<" using 1:2 with points title "Original XY-data"\
	,f(x) with lines title "Exact Cosine"\
	,"LinData.txt" using 1:2 with line title "Interpolated data"\
	'|gnuplot

LinInteg.pdf: LinIntegData.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "sin(x)";\
	set title "Integrated interpoletion";\
	f(x)=sin(x);\
	plot \
	"$<" using 1:2 with lines title "Integral of interpolation"\
	, f(x) with lines title "Exact sin(x)"\
	'|gnuplot

out.txt: main
	./$< > $@

main: main.c interp.c
	$(CC) main.c interp.c -o $@

.PHONEY: clean

clean:
	$(RM) *.txt *.pdf main
