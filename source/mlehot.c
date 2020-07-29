/* mlehot.c - Modified from programs developed by Hudson (2001, Genetics) 
to estimate pairwise composite likelihoods under a simple hotspot model 
similar to LDhot.  stdin is the output from seqhot.
Usage (e.g.):

ms0 2 30 -t 5. | seqhot 30 hotbeg hotend | mlehot fname gridin

fname = file name of pairwise likelihood file from el1.c (30elouta for n=30)
gridin = input file of rho and hotspot values for likelihood calculations.
  3 1. 2. 4 3 1. 3. 5.
  considers all combinations of background (total) rho equals 1., 2. and 4.
  and hotspot rate 1, 3 and 5 times the background rate

To compile: cc -O -o mlehot mlehot.c lnlikehot.c getprobu.c getprobmat1.c -lm

Output consists of a list of

npairs mle_co mle1 mle2 lik_tot lik_co

for each input data set, where

npairs = number of pairs of sites in the data set
mle_co = MLE under constant rho
mle1 = background MLE of rho under hotspot model
mle2 = MLE of hotspot multiplier
lik_tot = maximum likelihood under hotspot model
lik_co = maximum likelihood under constant rho model

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int npairs;
struct pairst {
  int n;
  int n1;
  int n2;
  int n11;
  double  d12;
  double hot12;
	} *pair ;	

char *fname, *gridin;
double cmax ;
double ****prob, *recrates , *ppoly ;
int nsitesp;
double *co_val, *hot_val;

main(int argc, char *argv[] )
{
  int i,j,k,  n, n1, n2,  nsites,  howmany  ;
  int n00, n01, n10, n11, s1, s2, s3 ;
  double pmle, findmax(int ), c, d12, r, theta;
  double lnlikehot (double,double,int,double *,double ****,double *);
  double minusone = -1.0 , alpha ;
  FILE *pout, *fopen(const char *, const char *), *pfrho;
  unsigned maxpair = 11000;
  double  ****getprobmat(char *, int *, int *, double **, double **);
  double cbest, bestlike, cl ;

  int nloci, covals, gcvals, tractvals, npop, nsam, *config, a, b,d;
  int mle1, mle2, mle3, seed1, seed2, seed3, ntrials, mle_c1, numwin, dummy;
  double **likmat, mle_co, mle_tot, alphag, t0, crate, clen, likact;
  double **getlikmat (char *, int *, int *), nsovert, mig_rate;
  double rtemp, *mletemp, rhotemp;

  int hotvals;
  double hot12;

  if (argc != 3) {
    printf("Usage: mlehot fname gridin\n");
    exit(0);
  }

  fname = argv[1];
  gridin = argv[2];

  prob = getprobmat(fname, &n,  &nsitesp, &recrates , &ppoly) ;
  pair = (struct pairst *) malloc( maxpair*sizeof( struct pairst) ) ; 
  likmat = getlikmat (gridin, &covals, &hotvals);
  
  for ( ; ; ) {
    for (i=0; i<covals; ++i)
      for (j=0; j<hotvals; ++j)
	likmat[i][j] = 0.0;
    if (scanf("%d\n", &npairs) == EOF )
	exit(0);
    if (npairs == 0)
      printf("0.\n");
    else if (npairs > 0) {
      
      if( (npairs+1)  >= maxpair ){
	maxpair = npairs + 10 ;
	pair = (struct pairst *) realloc( pair, maxpair*sizeof( struct pairst) ) ;
      }
      for( i=0; i<npairs ; i++ ) {
	scanf(" %d %d %d %d %lf %lf", &n, &n1, &n2, &n11, &d12, &hot12);
	pair[i].n = n ;
	pair[i].n1 = n1 ;
	pair[i].n2 = n2 ;
	pair[i].n11 = n11 ;
	pair[i].d12 = d12 ;
	pair[i].hot12 = hot12;
      }

      for (i=0; i<covals; ++i)
	for (j=0; j<hotvals; ++j) {
	  likmat[i][j] = lnlikehot (co_val[i], hot_val[j], nsitesp, 
				    recrates, prob, ppoly);
	  printf("");
	}
    }

    mle_tot = mle_co = likmat[0][0];
    mle1 = mle2 = mle3 = mle_c1 = 0;
    for (i=0; i<covals; ++i)
      for (j=0; j<hotvals; ++j)
	if (likmat[i][j] > mle_tot) {
	    mle_tot = likmat[i][j];
	    mle1 = i; mle2 = j;
	}

    for (i=0; i<covals; ++i) 
      if (likmat[i][0] > mle_co) {
	mle_c1 = i;
	mle_co = likmat[i][0];
      }

    printf("%d %lf %lf %lf %lf %lf\n", npairs, co_val[mle_c1], co_val[mle1], 
	   hot_val[mle2], mle_tot, mle_co);
      
  }
}
	
double **getlikmat (char *gridname, int *pco_count, int *phot_count)
{
  FILE *fgrid, *fopen();
  int co, hot, a, b, d;
  double **likmat;

  fgrid = fopen (gridname, "r");

  fscanf (fgrid, "%d ", &co);

  co_val = (double *) malloc (co * sizeof (double));
  for (a=0; a<co; ++a)
    fscanf (fgrid, "%lf ", &co_val[a]);
                                        
  fscanf (fgrid, "%d ", &hot);
  hot_val = (double *) malloc (hot * sizeof (double));
  for (a=0; a<hot; ++a)
    fscanf (fgrid, "%lf ", &hot_val[a]);
  fclose (fgrid);

  *pco_count = co;
  *phot_count = hot;

  likmat = (double **) malloc (co * sizeof (double **));
  for (a=0; a<co; ++a) 
    likmat[a] = (double *) malloc (hot * sizeof (double *));
   
  for (a=0; a<co; ++a)
    for (b=0; b<hot; ++b)
      likmat[a][b] = 0.0;

  return (likmat);
}
