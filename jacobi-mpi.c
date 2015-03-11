#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main ( int argc, char *argv[] )
{
	int mpisize, rank, tag, origin, destination;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 2) {
		fprintf(stderr, "Please specify how often the messege should be sent around the ring as an argument!\n");
		MPI_ABORT(MPI_COMM_WORLD,1);
	}

	long N = atoi(argv[1]);
	if (N%mpisize != 0)
	{
		fprintf(stderr, "N should be a multiple of p\n");
		MPI_ABORT(MPI_COMM_WORLD,1);
	}
	int uSize = N/mpisize;

	double residuals[N];
	int i;
	for (i=0;i<N;i++)
	{
		residuals[i] = 0;
	}
	double residuals_reduced[N];
	for (i=0;i<N;i++)
	{
		residuals_reduced[i] = 0;
	}
	int finished = 0;
	int message_out;
	int message_in;
	tag = 99;

	double h = 1;
	h = h/(N+1);

	MPI_Status status;

	double u[uSize];
	double upperNeighbor;
	double lowerNeighbor;

	double** A;
	A = (double**)malloc((N) * sizeof(double*));

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
	double pre_u[uSize + 2];
	for(i=0;i<uSize + 2;i++){
		pre_u[i] = 0;
	}

	long k;
	double h2 = h*h;
	double diag = 2/h2;
	for(k=0;;k++){
		for(i=0;i<uSize;i++){
			double sum = 1;
			sum -= (-1)*pre_u[i]/h2 + (-1)*pre_u[i+2]/h2;
			u[i] = sum/diag;
		}

		/*for(i=0;i<uSize;i++){
			printf("%12.10f\n", u[i]);
		}*/

/*		for(i=0;i<N+2;i++){
			printf("%12.10f\n", pre_u[i]);
		}
*/
		for(i=0;i<uSize+2;i++){
			pre_u[i+1] = u[i];
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
				if(rank==1)
				residuals[realI+1] -= u[i];
			}
		}
		MPI_Reduce(residuals, residuals_reduced, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(rank==0)
		{

			for(i=0;i<N;i++){
				residuals[i] = residuals[i]/h2 - 1;
			}
			double residual = 0;
			for(i=0;i<N;i++){
				residual += (residuals[i])*(residuals[i]);
			}
			residual = sqrt(residual);
			//printf("%ld %12.10f\n",k,residual);
			if(residual<=threshold){
				printf("%ld iterations\n",k);
				finished = 1;
			}
			MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
		       	if(finished == 1)
			{
				break;
			}
			else
			{
				for (i=0;i<N;i++)
				{
					residuals[i] = 0;
				}
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
		/*
		for(i=0;i<N;i++){
			double s = 0;
			int j = rank*uSize;
			if(i>=j-1 && i<=j+uSize)
			{
				if(i-1>=j)
				{
					s += (-1);
				}
				s += (-1)*
				for(j=0;j<N;j++){
					s += A[i][j]*u[j];
				}
			}
			residuals[i] = s - f[i];
		}
		double residual = 0;
		for(i=0;i<N;i++){
			residual += residuals[i]*residuals[i];
		}
		residual = sqrt(residual);
		//printf("%i %12.10f\n",k,residual);
		if(residual<=threshold){
			printf("%ld iterations\n",k);
			break;
		}
		*/
	}
}


