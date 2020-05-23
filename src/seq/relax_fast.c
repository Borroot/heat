#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

//
// allocate a vector of length "n"
// malloc nullptr is not being checked
double *allocVector(int n)
{
	double *v;
	v = (double *)malloc( n*sizeof(double));
	return v;
}

//
// initialise the values of the given vector "out" of length "n"
//
void init( double *out, int n, int heat)
{
	int i;

	for( i=1; i<n; i++) {
		out[i] = 0;
	}
	out[0] = heat;

}

//
// print the values of a given vector "out" of length "n"
//
void print( double *out, int n)
{
	int i;

	printf( "<");
	for( i=0; i<n; i++) {
		printf( " %f", out[i]);
	}
	printf( ">\n");
}

bool relaxAndStable(double *in, double *out, int n, double eps) {
	bool notstable = true;
	for(int i = 1; i < n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
		if(notstable) {
			notstable = notstable && (fabs(out[i] - in[i]) <= eps);
		}
	}
	return notstable;
}

// Commandline arguments:
// <N> <EPS> <HEAT>
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

	double *a,*b, *tmp;
	int iterations = 0;

	a = allocVector(N);
	b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	printf("size   : %f M (%d MB)\n", N/1000000, (int)(N*sizeof(double) / (1024*1024)));
	printf("heat   : %f\n", HEAT);
	printf("epsilon: %f\n", EPS);

	do {
		tmp = a;
		a = b;
		b = tmp;
		iterations++;
	} while(!relaxAndStable(a, b, N, EPS));

	printf("Number of iterations: %d\n", iterations);

	return 0;
}
