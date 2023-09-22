COMPILER = mpicc
CFLAGS = -Wall -pedantic
COBJS = gaussianLib.o qdbmp.o utilities.o
CEXES =  main

all: ${CEXES}

run:
	mpiexec -np 4 ./main pencils.bmp out.bmp 10

main: main.c ${COBJS}
	${COMPILER} ${CFLAGS} main.c ${COBJS} -o main -lm

%.o: %.c %.h  makefile
	${COMPILER} ${CFLAGS} -lm $< -c

clean:
	rm -f *.o *~ ${CEXES}
