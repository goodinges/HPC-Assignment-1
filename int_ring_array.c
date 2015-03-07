/* Ring Communication:
 * Exchange between messages in a ring
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "util.h"
#include <time.h>
#include <stdlib.h>

int main( int argc, char *argv[])
{
  int mpisize, rank, tag, origin, destination;
  long N;
  MPI_Status status;
  timestamp_type time1, time2;
  
  int arraySize = 2*1024*1024/sizeof(int);
  int *a;
  a = (int *) malloc(2*1024*1024);
  srand(time(NULL));
  int j=0;
  for(; j < arraySize ; j++)
  {
	  a[j] = rand();
  }

  char hostname[1024];
  gethostname(hostname, 1024);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc != 2) {
    fprintf(stderr, "Please specify how often the messege should be sent around the ring as an argument!\n");
    MPI_ABORT(MPI_COMM_WORLD,1);
  }

  N = atol(argv[1]);

  int message_out;
  int message_in = -1;
  tag = 99;

  if(rank == mpisize-1)
  {
	  destination = 0 ;
  }
  else
  {
	  destination = rank + 1 ;
  }

  if(rank == 0)
  {
    origin = mpisize - 1 ;

    message_out = 0 ;
  }
  else
  {
    origin = rank - 1;
  }

  int i = 0;
  get_timestamp(&time1);
  for(; i < N ; i++ )
  {
	  if(rank != 0 || i != 0 )
	  {
		  MPI_Recv(a, arraySize, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
	  }

    MPI_Send(a, arraySize, MPI_INT, destination, tag, MPI_COMM_WORLD);
  }

  if(rank == 0 )
  {
	  MPI_Recv(a, arraySize, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
	  get_timestamp(&time2);
	  double elapsed = timestamp_diff_in_seconds(time1,time2);
	  printf("Rank %d hosted on %s received the last message.\n", rank, hostname);

	  printf("Total time spent on communicating: %f seconds.\n", elapsed);
	  printf("Communication per second: %f\n", N*mpisize/elapsed);
	  printf("MB data communication per second: %f MB/s\n", N*mpisize*2/elapsed);
  }

  MPI_Finalize();
  return 0;
}
