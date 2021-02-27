CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

out.txt: main
	./$< >$@

.PHONEY: clean

clean:
	$(RM) main out.txt
