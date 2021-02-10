CC = gcc
RM = rm -f
CFFLAGS = -std-gnu99 -O -Wall
LDLIBs = -lm
cat: out.txt
	@echo -n 'The results of the program is:'
	@cat out.txt

out.txt: combined
	@echo 'The run time of the project is:'
	time ./combined > out.txt

combined: epsilon.c epsilon2.c exercise-math.c
	$(CC) epsilon.c epsilon2.c exercise-math.c -o combined

clean:
	$(RM) combined out.txt
