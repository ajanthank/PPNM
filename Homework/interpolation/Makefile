CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt Lin.pdf Quad.pdf

Lin.pdf: XY.txt  Makefile LinData.txt LinIntegData.txt
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set xrange [0:9.5];\
	set ylabel "f(x)";\
	set title "Linear Interpolation";\
	f(x)=cos(x);\
	g(x)=sin(x);\
	plot\
	"$<" index 0 using 1:2 with points title "Original XY-data"\
	,"LinData.txt" using 1:2 with line title "Interpolated data"\
	,"LinIntegData.txt" using 1:2 with lines title "Integral of interpolation"\
	,f(x) with lines title "cos(x)"\
	,g(x) with lines title "sin(x)"\
	'|gnuplot

Quad.pdf: XY.txt out.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set xrange [0:9.5];\
	set ylabel "f(x)";\
	set title "Quadratic interpoletion";\
	plot \
	"$<" index 0 using 1:2 with points title "Original XY-data"\
	,"out.txt" using 1:2 with lines title "Interpolated data"\
	,"out.txt" using 1:3 with lines title "Integral of interpolation"\
	,"out.txt" using 1:4 with lines title "Derivative of interpolation"\
	,"out.txt" using 1:5 with lines title "cos(x)"\
	,"out.txt" using 1:6 with lines title "sin(x)"\
	,"out.txt" using 1:7 with lines title "-sin(x)"\
	'|gnuplot

out.txt: main
	./$< > $@

main: main.c interp.c
	$(CC) main.c interp.c -o $@

.PHONEY: clean

clean:
	$(RM) *.txt *.pdf main *~
