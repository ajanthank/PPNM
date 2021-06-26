CC=gcc
CFLAGS= -std=gnu11 -O1 -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs) -lm -Wall

.PHONEY: default

default: out.txt 


out.txt: main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) *.o main *.txt *.pdf *~
