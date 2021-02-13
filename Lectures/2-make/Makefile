CC= gcc
RM= rm -f


default: out.txt
	cat out.txt

out.txt: hello
	./hello > out.txt

hello: hello.c
	$(CC) hello.c -o hello 

#hello.o: hello.c
#	$(CC) -c hello.c

clean:
	$(RM) out.txt hello.o hello

test:
	echo $(CC)
	echo $(RM)
