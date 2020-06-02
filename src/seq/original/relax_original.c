#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

double *allocVector(int n)
{
	double *v;
	v = (double *)malloc( n*sizeof(double));
	return v;
}

void init(double *out, int n, int h)
{
	int i;

	for( i=1; i<n; i++) {
		out[i] = 0;
	}
	out[0] = h;

}

void relax(double *in, double *out, int n)
{
	int i;
	for( i=1; i<n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
	}
}

bool isStable(double *old, double *new, int n, double eps)
{
	int i;
	bool res=true;

	for( i=1; i<n-1; i++) {
		res = res && ( fabs( old[i] - new[i]) <= eps);
	}
	return res;
}

int main(int argc, char *argv[])
{
	double N = 1000000.0;   // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%lf", &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	double start = omp_get_wtime();
	double *a,*b, *tmp;

	a = allocVector(N);
	b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	do {
		tmp = a;
		a = b;
		b = tmp;
		relax(a, b, N);
	} while(!isStable(a, b, N, EPS));

	double end = omp_get_wtime();
	printf("%f\n", end - start);

	free(a);
	free(b);

	return 0;
}
