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

int relaxAndStable(double *in, double *out, int n, double eps) {
	int stable = 1;
	for(int i = 1; i < n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
		stable = stable && (fabs(out[i] - in[i]) <= eps);
	}
	return stable;
}

void print(double *out, int n) {
	int i;

	fprintf(stderr, "<");
	for( i=0; i<n; i++) {
		fprintf(stderr, " %f", out[i]);
	}
	fprintf(stderr, ">\n");
}

int main(int argc, char *argv[])
{
	double N = 100.0;    // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%lf", &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	MPI_Init(&argc, &argv);
	double start = MPI_Wtime();

	int num_ranks, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank == 0) {
		double *a = allocVector(N);
		double *b = allocVector(N);

		init(a, N, HEAT);
		init(b, N, HEAT);

		double *tmp;
		do {
			tmp = a;
			a = b;
			b = tmp;
		} while(!relaxAndStable(a, b, N, EPS));

		free(a);
		free(b);
	}

	double end = MPI_Wtime();
	MPI_Finalize();

	if(my_rank == 0) {
		printf("%f\n", end - start);
	}

	return 0;
}