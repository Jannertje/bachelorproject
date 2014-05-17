CFLAGS = -Wall $$(gsl-config --cflags) -g -I.
LIBS = $$(gsl-config --libs)
DEPS = main.h error.h tree.h macros.h types.h sorter.h algo.h
OBJ = main.o error.o tree.o sorter.o algo.o
default: $(OBJ)
	gcc -o main $^ $(LIBS) $(CFLAGS)

%.o: %.c $(DEPS)
	gcc -c -o $@ $< $(CFLAGS)

.PHONY: clean valgrind

clean:
	rm -f *.o main
	rm -rf *.dSYM

valgrind: main
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes --dsymutil=yes ./main 
