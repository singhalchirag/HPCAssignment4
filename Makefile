CC=mpicc
FLAGS=-O3
EXECS=mpi sort omp_solved2  omp_solved3 omp_solved4 omp_solved5 omp_solved6 omp_solved7 omp_solved1

all: ${EXECS}

mpi:mpi.c
	${CC} ${FLAGS} -lm -o mpi mpi.c

sort:sort.c
	${CC} ${FLAGS} -o sort sort.c

omp_solved2: mpi_solved2.c
	${CC} ${FLAGS} mpi_solved1.c -o omp_solved2

omp_solved3: mpi_solved3.c
	${CC} ${FLAGS} mpi_solved3.c -o omp_solved3

omp_solved4: mpi_solved4.c
	${CC} ${FLAGS} mpi_solved4.c -o omp_solved4

omp_solved5: mpi_solved5.c
	${CC} ${FLAGS} mpi_solved5.c -o omp_solved5

omp_solved6: mpi_solved6.c
	${CC} ${FLAGS} mpi_solved6.c -o omp_solved6

omp_solved7: mpi_solved7.c
	${CC} ${FLAGS} mpi_solved7.c -o omp_solved7

omp_solved1: mpi_solved1.c
	${CC} ${FLAGS} mpi_solved1.c -o omp_solved1


clean:
	rm -f ${EXECS}
