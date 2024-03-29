#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

double *allocVector(int n)
{
	double *vector = (double *)malloc(n * sizeof(double));
	if(vector) {
		return vector;
	} else {
		exit(1);
	}
}

void init(double *out, int n, int heat, int isStart)
{
	for (int i = 1; i < n; i++) {
		out[i] = 0;
	}
	if (isStart) {
		out[0] = heat;
	}
}

void relax(double* in, double* out, int from, int n, int N) {
	int start = from + (from == 0 ? 1 : 0);
	int end = from + n - (from + n == N ? 1 : 0);
	for(int i = start; i < end; i++) {
		out[i - from] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
	}
}

int isStable(double* in, double* out, int n, double eps) {
	for(int i = 1; i < n-1; i++) {
		if(fabs(out[i] - in[i]) > eps) return 0;
	}
	return 1;
}

int main(int argc, char *argv[])
{
	int N = 100;            // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%d", &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	MPI_Init(&argc, &argv);

	double start = MPI_Wtime();

	int num_ranks, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Make sure that N is divisible by the number of ranks - 1
	N -= N % (num_ranks - 1);
	int blocksize = N / (num_ranks - 1);

	double *a = allocVector(N);
	init(a, N, HEAT, 1);

	int stable;

	if (my_rank == 0 && num_ranks != 1) {
		double *b = allocVector(N);
		init(b, N, HEAT, 1);
		double *tmp;
		do {
			for (int dest = 1; dest < num_ranks; dest++) {
				MPI_Send(a,N, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
			}
			// Receive blocks and combine them
			for (int dest = 1; dest < num_ranks; dest++) {
				MPI_Recv(b + (dest - 1) * blocksize,blocksize, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			stable = isStable(a,b,N,EPS);
			// print(b,N);
			tmp = a;
			a = b;
			b = tmp;
		} while(!stable);
		// Communicate to other processes that the should not wait for more messages
		a[0] = 0;
		for (int dest = 1; dest < num_ranks; dest++) {
			MPI_Send(a,N, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
		}
		a[0] = 1;
		free(b);
	} else {
		int from = blocksize * (my_rank - 1);
		double *b = allocVector(blocksize);
		init(b, blocksize, HEAT, my_rank == 1);
		do {
			MPI_Recv(a, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			stable = !a[0];

			if (!stable) {
				relax(a,b,from,blocksize,N);
				MPI_Send(b,blocksize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
		} while (!stable);
		free(b);
	}

	double end = MPI_Wtime();
	MPI_Finalize();

	if(my_rank == 0) {
		printf("%f\n", end - start);
	}

	free(a);

	return 0;
}
