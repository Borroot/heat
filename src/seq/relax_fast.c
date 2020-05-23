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

//
// individual step of the 3-point stencil
// computes values in vector "out" from those in vector "in"
// assuming both are of length "n"
//
void relax( double *in, double *out, int n)
{
	for(int i=1; i<n-1; i++) {
		out[i] = 0.25*in[i-1] + 0.5*in[i] + 0.25*in[i+1];
	}
}

//
// checks the convergence criterion:
// true, iff for all indices i, we have |out[i] - in[i]| <= eps
//
bool isStable( double *old, double *new, int n, double eps)
{
	for(int i=1; i<n-1; i++) {
		if(!(fabs(old[i] - new[i]) <= eps)) return false;
	}
	return true;
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
    double N = 100000.0;  // length of the vectors
    double EPS = 0.1;       // convergence criterium
    double HEAT = 100.0;    // heat value on the boundary

    if (argc == 4) {
        sscanf(argv[1], "%lf", &N);
        sscanf(argv[2], "%lf", &EPS);
        sscanf(argv[3], "%lf", &HEAT);
    }

	double *a,*b, *tmp;
	int n;
	int iterations = 0;

	a = allocVector(N);
	b = allocVector(N);

	init(a, N, HEAT);
	init(b, N, HEAT);

	n = N;

	printf("size   : %d M (%d MB)\n", n/1000000, (int)(n*sizeof(double) / (1024*1024)));
	printf("heat   : %f\n", HEAT);
	printf("epsilon: %f\n", EPS);

	#if RELAX_AND_STABLE == 1
		do {
			tmp = a;
			a = b;
			b = tmp;
			iterations++;
		} while(!relaxAndStable(a, b, n, EPS));
	#else
		do {
			tmp = a;
			a = b;
			b = tmp;
			relax(a, b, n);
			iterations ++;
		} while(!isStable(a, b, n, EPS));
	#endif

	//print(b, n);
	printf("Number of iterations: %d\n", iterations);

	return 0;
}
