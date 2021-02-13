CC = gcc
RM = rm -f
CFFLAGS = -std-gnu99 -O -Wall
LDLIBs = -lm
cat: out.txt
	@echo -n 'The results of the program is:'
	@cat out.txt

out.txt: math
	@echo 'The run time of the project is:'
	time ./math > out.txt

math: exercise-math.c
	$(CC) exercise-math.c -o math

clean:
	$(RM) math out.txt
