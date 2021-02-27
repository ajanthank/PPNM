CFLAGS = $(shellA gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: erfplot.png gammaplot.png lgammaplot.png cgammaplot.png

erfplot.png : dataset.txt Makefile tabulateddata.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "Erf";\
	set title "Error Function";\
	plot \
	"$<" using 1:2 with line title "Erf from math.h"\
	,"$<" using 1:3 with line title "Gsl-sf-erf"\
	,"$<" using 1:4 with line title "The Homemade Erf"\
	,"tabulateddata.txt" using 1:2 with points title "Tabulated data"\
	'|gnuplot

gammaplot.png : data2.txt Makefile tablegamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "Gamma";\
	set title "Gamma Function";\
	plot \
	"$<" using 1:2 with line title "Gamma function from math.h"\
	,"$<" using 1:3 with line title "gsl_sf_gamma"\
	,"$<" using 1:4 with line title "Approximation"\
	,"tablegamma.txt" with points title "Tabulated data"\
	'|gnuplot

lgammaplot.png : data3.txt Makefile tablelngamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
        set key top left;\
        set tics out;\
        set xlabel "x";\
        set ylabel "ln(Gamma)";\
        set title "ln(Gamma) Function";\
        plot \
        "$<" using 1:2 with line title "log Gamma function from math.h"\
        ,"$<" using 1:3 with line title "gsl_sf_lngamma"\
        ,"$<" using 1:4 with line title "Approximation"\
        ,"tablelngamma.txt" with points title "Tabulated data"\
        '|gnuplot

cgammaplot.png : data4.txt Makefile
	echo '\
	set title "Complex Gamma-function";\
	set out "$@";\
	set terminal png;\
	set xlabel "Re(z)";\
	set ylabel "Im(z)";\
	set zlabel "|Γ(z)|";\
	set zrange [0:5];\
	set tics out;\
	set key off;\
	set pm3d interpolate 2,2;\
	set dgrid3d 50,50 qnorm 5;\
	set hidden3d;\
	splot \
	"$<" using 1:2:3 with pm3d\
	'|gnuplot


dataset.txt :main
	./$< > $@

.PHONEY: data2.txt

main: main.o erf.o gamma.o lngamma.o complex_gamma.o

.PHONEY:clean

clean:
	$(RM) *.o main *.png data*
