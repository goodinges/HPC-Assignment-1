all: a b

a:
mpicc int_ring.c -lrt -lm -o int_ring

b:
mpicc jacobi-mpi.c -lrt -lm -o jacobi-mpi
