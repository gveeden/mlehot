#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#define PR(format,value) printf("value = %format\t",(value))
#define NL putchar('\n')
#define size_t unsigned
#define THINK_C_CL 0

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define MAX(x, y) ( (x)>(y) ? (x) : (y) )
#define PRINT1(f,x1) PR(f,x1), NL
#define PR_MEAN_VAR(pf,string,var) fprintf(pf,"string %lf ( %e\
)\n",var.sum,var.sumsq)

#define INPUT(for,x) fprintf(stderr," x=?\n"), fscanf(stdin," %for",&(x))

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

struct runsum {
	int count;
	double sum;
	double sumsq;
	} ;
