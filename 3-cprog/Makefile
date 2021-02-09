CC = gcc
RM = rm -f

cat: out.txt
	echo -n 'The results of the program is:'
	cat out.txt

out.txt: epsilon
	echo -n 'The run time of the program is:\n'
	time ./epsilon>out.txt

epsilon: epsilon.c epsilon2.c
	$(CC) epsilon.c epsilon2.c -o epsilon

clean:
	$(RM) epsilon out.txt
