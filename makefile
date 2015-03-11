a: int_ring.c
	mpicc int_ring.c -lrt -lm -o int_ring
	mpicc jacobi-mpi.c -lrt -lm -o jacobi-mpi
