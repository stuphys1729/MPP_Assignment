/* This is a file with which to perform the MPI communications between processes. */

#include <mpi.h>

#include "precision.h"

#define MAX_DIMS 2
#define TRUE 1
#define FALSE 0

void initialise_MP(int* cart_comm, int* rank, int* dims) {

	int size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, rank);

	// Initialise a cartesian topology
	int* dims[MAX_DIMS];
	int* periods[MAX_DIMS];

	for (int i = 0; i < MAX_DIMS; i++) {
		dims[i] = 0;
		if (i == 0) periods[i] = TRUE;
		else periods[i] = FALSE;
	}

	MPI_Dims_create(size, MAX_DIMS, dims);
	MPI_Cart_create(MPI_COMM_WORLD, MAX_DIMS, dims, periods, TRUE, cart_comm);

}