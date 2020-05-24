#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

//
// allocate a vector of length "n"
//
double *allocVector(int n)
{
	double *v;
	v = (double *)malloc( n*sizeof(double));
	return v;
}

//
// initialise the values of the given vector "out" of len gth "n"
//
void init(double *out, int n, int h)
{
	int i;

	for( i=1; i<n; i++) {
		out[i] = 0;
	}
	out[0] = h;

}

//
// print the values of a given vector "out" of length "n"
//
void print(double *out, int n)
{
	int i;

	printf( "<");
	for( i=0; i<n; i++) {
		printf( " %f", out[i]);
	}
	printf( ">\n");
}

//
// individual step of the 3-point stencil
// computes values in vector "out" from those in vector "in"
// assuming both are of length "n"
//
void relax(double *in, double *out, int n)
{
	int i;
	for( i=1; i<n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
	}
}

//
// checks the convergence criterion:
// true, iff for all indices i, we have |out[i] - in[i]| <= eps
//
bool isStable(double *old, double *new, int n, double eps)
{
	int i;
	bool res=true;

	for( i=1; i<n-1; i++) {
		res = res && ( fabs( old[i] - new[i]) <= eps);
	}
	return res;
}

// Commandline arguments:
// <N> <EPS> <HEAT>
int main(int argc, char *argv[])
{
	double N = 1000000.0;  // length of the vectors
	double EPS = 0.1;       // convergence criterium
	double HEAT = 100.0;    // heat value on the boundary

	if (argc == 4) {
		sscanf(argv[1], "%lf", &N);
		sscanf(argv[2], "%lf", &EPS);
		sscanf(argv[3], "%lf", &HEAT);
	}

	printf("size   : %f M (%d MB)\n", N/1000000.0, (int)(N*sizeof(double) / (1024*1024)));
	printf("epsilon: %f\n", EPS);
	printf("heat   : %f\n", HEAT);

	double start = omp_get_wtime();

	double *a,*b, *tmp;
	int iterations = 0;

	a = allocVector(N);
	b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	do {
		tmp = a;
		a = b;
		b = tmp;
		relax(a, b, N);
		// print(b, n);
		iterations ++;
	} while(!isStable(a, b, N, EPS));

	double end = omp_get_wtime();
	printf("%f\n", end - start);
	printf("Number of iterations: %d\n", iterations);

	return 0;
}
