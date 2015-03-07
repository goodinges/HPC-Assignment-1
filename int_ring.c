/* Ring Communication:
 * Exchange between messages in a ring
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "util.h"

int main( int argc, char *argv[])
{
  int mpisize, rank, tag, origin, destination;
  long N;
  MPI_Status status;
  timestamp_type time1, time2;

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
		  MPI_Recv(&message_in,  1, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
		  message_out = message_in + rank ;
	  }

    MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
  }

  if(rank == 0 )
  {
	  MPI_Recv(&message_in,  1, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
	  get_timestamp(&time2);
	  double elapsed = timestamp_diff_in_seconds(time1,time2);
	  printf("Rank %d hosted on %s received the last message: %d which is ", rank, hostname, message_in);
	  if( message_in == N * ( (mpisize - 1) * mpisize / 2 ) )
	  {
		  printf("right!\n");
	  }
	  else
	  {
		  printf("wrong!\n");
	  }
	  printf("Total time spent on communicating: %f seconds.\n", elapsed);
	  printf("Communication per second: %f\n", N*mpisize/elapsed);
	  printf("KB data communication per second: %f KB/s\n", N*mpisize*sizeof(int)/1024/elapsed);
  }

  MPI_Finalize();
  return 0;
}
