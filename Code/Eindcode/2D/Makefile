CFLAGS = -Wall $$(gsl-config --cflags) -g3 -I. -fopenmp
LIBS = $$(gsl-config --libs)
DEPS = tri.h pkd.h pair.h edge.h workspace.h tree.h types.h partition.h error.h triangulate.h treegen.h main.h helper.h
OBJ = pkd.o tri.o pair.o edge.o workspace.o tree.o partition.o error.o triangulate.o treegen.o main.o helper.o
CC = gcc-4.9
default: $(OBJ)
	$(CC) -o main $^ $(LIBS) $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean valgrind qcachegrind

clean:
	rm -f *.o main
	rm -rf *.dSYM

valgrind: main
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes --dsymutil=yes ./main

qcachegrind: main
	valgrind --tool=callgrind -v --dump-every-bb=1000000000 ./main 400
