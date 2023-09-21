COMPILER = mpicc
CFLAGS = -Wall -pedantic
COBJS = gaussianLib.o qdbmp.o utilities.o
CEXES =  main

all: ${CEXES}

gaussian: gaussian.c ${COBJS}
	${COMPILER} ${CFLAGS} gaussian.c ${COBJS} -o gaussian -lm

%.o: %.c %.h  makefile
	${COMPILER} ${CFLAGS} -lm $< -c

clean:
	rm -f *.o *~ ${CEXES}
