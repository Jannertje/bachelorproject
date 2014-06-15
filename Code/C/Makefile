CFLAGS = -Wall $$(gsl-config --cflags) -g3 -I.
LIBS = $$(gsl-config --libs)
DEPS = tri.h pkd.h pair.h edge.h workspace.h tree.h types.h partition.h
OBJ = pkd.o tri.o pair.o edge.o workspace.o tree.o partition.o
default: $(OBJ)
	gcc -o main $^ $(LIBS) $(CFLAGS)

%.o: %.c $(DEPS)
	gcc -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o main
	rm -rf *.dSYM
