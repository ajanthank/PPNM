LDLIBS = -lm

output.txt error.txt out.txt : intput.txt stdio Makefile
	./stdio < input.txt 1> output.txt 2> error.txt

input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY: clean
clean:
	$(RM) input.txt out.txt out.file.txt stdio
