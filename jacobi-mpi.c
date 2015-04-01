#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main ( int argc, char *argv[] )
{
	int mpisize, rank, tag;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 3) {
		fprintf(stderr, "Wrong number of arguments!\n");
		MPI_ABORT(MPI_COMM_WORLD,1);
	}

	long N = atoi(argv[1]);
	if (N%mpisize != 0)
	{
		fprintf(stderr, "N should be a multiple of p\n");
		MPI_ABORT(MPI_COMM_WORLD,1);
	}
	long maxIterations = atoi(argv[2]);
	int uSize = N/mpisize;

	double *residuals = calloc(N, sizeof(double));
	int i;
	double *residuals_reduced;
	if(rank==0)
	{
		residuals_reduced = calloc(N, sizeof(double));
	}
	int finished = 0;
	tag = 99;

	double h = 1;
	h = h/(N+1);

	MPI_Status status;

	double *u = calloc(uSize, sizeof(double));

	//double** A;
	//A = (double**)malloc((N) * sizeof(double*));

/*	long j;
	for(i=0;i<N;i++){
		A[i] = (double*)malloc((N) * sizeof(double));

		double h2 = h*h;
		double d1 = 2/h2;
		double d2 = -1/h2;
		for(j=0;j<N;j++){
			if(i==j){
				A[i][j] = d1;
			}else if(i==j+1||i==j-1){
				A[i][j] = d2;
			}else{
				A[i][j] = 0;
			}
		}
	}
	*/
	/*
	printf("Let A be:\n");
	for (i = 0; i < N ; i++)
	{
		for(j = 0; j < N; j++)
			printf("%5.2f ", A[i][j]);
		printf("\n");
	}
*/
	double initial_residual;
	double threshold;
	if(rank == 0)
	{
		initial_residual = 0;
		//double h2 = h*h;
		//double d1 = 2/h2;
		//double d2 = -1/h2;
		//double initial_residuals[N];
/*		for(i=0;i<N;i++){
			//double s = 0;
			for(j=0;j<N;j++){
				s += A[i][j]*u0[j];
				if(i==j){
					A[i][j] = d1;
				}else if(i==j+1||i==j-1){
					A[i][j] = d2;
				}else{
					A[i][j] = 0;
				}
			}

			//initial_residuals[i] = s - f[i];
		}
*/
		for(i=0;i<N;i++){
			initial_residual += 1;//initial_residuals[i]*initial_residuals[i];
		}
		initial_residual = sqrt(initial_residual);
		//printf("%6.4f\n",initial_residual);
		threshold = initial_residual;
		threshold /= 1000000;

	}

	//Jacobi
	double *pre_u = calloc(uSize + 2, sizeof(double));

	long k;
	double h2 = h*h;
	double diag = 2/h2;
	double residual;
	for(k=0;k<maxIterations ;k++){
		for(i=0;i<uSize;i++){
			double sum = 1;
			sum -= (-1)*pre_u[i]/h2 + (-1)*pre_u[i+2]/h2;
			u[i] = sum/diag;
		}

		for(i=0;i<uSize+2;i++){
			pre_u[i+1] = u[i];
		}


//		MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(finished == 1)
		{
			break;
		}
		else
		{
//			for (i=0;i<N;i++)
//			{
//				residuals[i] = 0;
//			}
			if(rank!=0)
			{
				MPI_Send(&u[0], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
			}
			if(rank!=mpisize-1)
			{
				MPI_Send(&u[uSize-1], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
				MPI_Recv(&pre_u[uSize+1], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
			}
			if(rank!=0)
			{
				MPI_Recv(&pre_u[0], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
			}
		}

	}

	//Computing individual residuals
	for(i=0;i<uSize;i++)
	{
		int realI = rank*uSize + i;
		if(realI-1>=0)
		{
			residuals[realI-1] -= u[i];
		}
		residuals[realI] += 2*u[i];
		if(realI+1<N)
		{
			residuals[realI+1] -= u[i];
		}
	}
	/*	for(i=0;i<N;i++){
		printf("%d: r:%d: %12.10f\n", rank, i, residuals[i]);
		}*/
	MPI_Reduce(residuals, residuals_reduced, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	/*if(rank==0)
	  for(i=0;i<N;i++){
	  printf("reduced: %d: r:%d: %12.10f\n", rank, i, residuals_reduced[i]);
	  }*/

	if(rank==0)
	{

		for(i=0;i<N;i++){
			residuals_reduced[i] = residuals_reduced[i]/h2 - 1;
			//	printf("%12.10f\n", residuals_reduced[i]);
		}
		residual = 0;
		for(i=0;i<N;i++){
			residual += (residuals_reduced[i])*(residuals_reduced[i]);
		}
		residual = sqrt(residual);
		//printf("%ld %12.10f\n",k,residual);
		if(residual<=threshold){
			finished = 1;
		}
	}
	
	if(rank==0){
		printf("%ld iterations\n",k);
		printf("Residual: %12.10f\n",residual);
	}

	MPI_Finalize();
}
