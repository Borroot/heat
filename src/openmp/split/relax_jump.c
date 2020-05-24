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
	#pragma omp parallel for schedule(static)
	for (int i = 1; i < n; i++) {
		out[i] = 0;
	}
	out[0] = heat;
}

void relax(double *in, double *out, int n)
{
	#pragma omp parallel for schedule(static)
	for(int i = 1; i < n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
	}
}

int isStable(double *out, double *in, int n, double eps, int count)
{
	if(count % 4 == 0) {
		for(int i = 0; i < n; i++) {
			if(fabs(out[i] - in[i]) > eps) {
				return 0;
			}
		}
		return 1;
	}
	return 0;
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

	double start = omp_get_wtime();

	double *a = allocVector(N);
	double *b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	double *tmp;
	int count = 0;
	do {
		tmp = a;
		a = b;
		b = tmp;
		relax(a, b, N);
		count++;
	} while(!isStable(a, b, N, EPS, count));

	double end = omp_get_wtime();
	printf("%f\n", end - start);

	return 0;
}
