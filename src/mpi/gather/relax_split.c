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

void init(double *out, int n, int heat)
{
	for (int i = 1; i < n; i++) {
		out[i] = 0;
	}
	out[0] = heat;
}

void relax(double* in, double* out, int n) {
	for(int i = 1; i < n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
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
	double N = 100000.0;    // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%lf", &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	int num_ranks, my_rank;
	MPI_Init(&argc, &argv);
	double start = MPI_Wtime();

	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int L = N / num_ranks;

	printf("TEST%d\n", my_rank);
	printf("TEST2%d\n", num_ranks);

	double *a = allocVector(N);
	double *b = allocVector(N);
  double *sub = allocVector(L);

	init(a, N, HEAT);
	init(b, N, HEAT);
  init(sub, L, 0);
  if(my_rank == 0) sub[0] = HEAT;

	double *tmp;
  int stable;
	do {
		tmp = a;
		a = b;
		b = tmp;
        relax(a+(L*my_rank), sub, L);
		int part_stable = isStable(a+(L*my_rank), sub, L, EPS);
        MPI_Allgather(sub, L, MPI_DOUBLE, b, L, MPI_DOUBLE, MPI_COMM_WORLD);
        stable = 1;
        MPI_Allreduce(&part_stable, &stable, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	} while(!stable);

	double end = MPI_Wtime();
	MPI_Finalize();

	if(my_rank == 0) {
		printf("%f\n", end - start);
	}

	free(a);
	free(b);
	free(sub);

	return 0;
}
