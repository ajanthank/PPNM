CC = gcc
CFLAGS = -O -Wall
LDLIBS = -lm

default: out.txt out_file.txt assignment3

out.txt error.txt: read-cmdline input.txt
	./read-cmdline $$(cat input.txt) 1> out.txt 2> error.txt

out_file.txt error.txt: input.txt read-stdin Makefile
	./read-stdin < input.txt 1> out_file.txt 2>>error.txt

assignment3: read-file
	./read-file

input.txt: Makefile
	cat /dev/null
	for x in $$(seq 1 0.5 5); do echo "$$x" >> $@; done

.PHONEY : clean

clean :
	$(RM) read-cmdline input.txt out* error.txt read-stdin read-file
