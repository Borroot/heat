#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

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
		if(in[i] == 0.0) return;
	}
}

void print(double *out, int n)
{
	int i;

	fprintf(stderr, "<");
	for( i=0; i<n; i++) {
		fprintf(stderr, " %f", out[i]);
	}
	fprintf(stderr, ">\n");
}

int isStable(double* in, double* out, int n, double eps) {
	for(int i = 1; i < n-1; i++) {
		if(fabs(out[i] - in[i]) > eps) return 0;
	}
	return 1;
}

int main(int argc, char *argv[])
{
	int N = 1000000;    // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%d",  &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	fprintf(stderr, "size   : %f M (%d MB)\n", N/1000000.0, (int)(N*sizeof(double) / (1024*1024)));
	fprintf(stderr, "epsilon: %f\n", EPS);
	fprintf(stderr, "heat   : %f\n", HEAT);

	double start = omp_get_wtime();

	double *a = allocVector(N);
	double *b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	int iterations = 0;
	double *tmp;
	do {
		for (int i = 0; i < 3; i++) {
			tmp = a;
			a = b;
			b = tmp;
			relax(a, b, N);
			iterations++;
		}
	} while(!isStable(a, b, N, EPS));

	double end = omp_get_wtime();

	//print(b, N);

	printf("%f\n", end - start);
	fprintf(stderr, "Iterations: %d\n", iterations);

	return 0;
}
