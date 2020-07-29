/*
hotspotliks3.c - Takes as input a list of likelihoods from mlehot (actual 
data) and outputs the probability of observing a difference in log likelihood
(normalized by the number of pairs of sites included in the likelihood 
calculation) greater than or equal to what is observed.  Simulations can match
S and R up to pre-specified tolerances:

wS <= s <= xS and yR <= r <= zR

Usage:

hotspotliks3 grid simlikfile w x y z < actuallikfile

grid = file of background rho and hotspot multiplier rates
simlikfile = file (from mlehot) of simulations under null model
actuallikfile = file (from mlehot) from actual data to be analyzed

For each line of input from stdin, the program outputs

S rho_hat rho_back hotspot den num p

S = number of segregating sites
rho_hat = estimated rho under constant rho model
rho_back = estimated background rho under hotspot model
hotspot = ratio of hotspot rho to background rho
den = number of simulations with rho_hat ($2) and S ($1) roughly equal to 
      value in input file
num = number of simulations in den with ($5-$6) >= actual value
p = probability = max (num/den, 1/den)

If den = 0, output is 0 0 1.

To compile, type cc -O -o hotspotliks3 hotspotliks3.c

The array ***logliks contains likelihood ratio statistics from simlikfile,
with logliks[a][b] corresponding to values with S=a and rho/10 = b.

Direct questions to Jeff Wall <wallj@humgen.ucsf.edu>
3/19/13
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SMAX 1000
#define RMAX 160
#define REPMAX 10000

int main (int argc, char *argv[])
{
  int a, b, num, den, pairs, S, **repcount, R, s, r, gridsize;
  double rho, ***logliks, lik1, lik2, lik3, lik4, testlik, w, x, y, z;
  double *rgrid;
  FILE *simlikfile, *grid;

  if (argc != 7 ) {
    printf("Usage: hotspotliks3 grid simlikfile w x y z < actuallikfile\n");
    exit(0);
  }

  w = atof (argv[3]);
  x = atof (argv[4]);
  y = atof (argv[5]);
  z = atof (argv[6]);

  rgrid = (double *) malloc (RMAX * sizeof (double));
  grid = fopen (argv[1], "r");
  fscanf (grid, "%d ", &gridsize);
  for (b=0; b<gridsize; ++b)
    fscanf (grid, "%lf ", &rgrid[b]);
  fclose (grid);

  logliks = (double ***) malloc ((SMAX+1) * sizeof (double **));
  for (a=1; a<=SMAX; ++a) {
    logliks[a] = (double **) malloc (RMAX * sizeof (double *));
    for (b=1; b<RMAX; ++b)
      logliks[a][b] = (double *) malloc (REPMAX * sizeof (double));
  }
  repcount = (int **) malloc ((SMAX+1) * sizeof (int *));
  for (a=1; a<=SMAX; ++a) {
    repcount[a] = (int *) malloc (RMAX * sizeof (int));
    for (b=0; b<RMAX; ++b)
      repcount[a][b] = 0;
  }

  simlikfile = fopen (argv[2], "r");

  while (fscanf (simlikfile, "%d %lf %*f %*f %lf %lf\n", &pairs, &rho, &lik1,
		 &lik2) != EOF) {

    S = (int) floor (sqrt (2.*pairs));
    for (R=0; R<RMAX; ++R)
      if (rho<(rgrid[R]+1E-08))
	break;
    ++R;

    if (S>SMAX)
      printf("Error. S = %d\n", S);
    else if (R>=RMAX)
      printf("Error. R = %d\n", R);
    else {
      logliks[S][R][repcount[S][R]] = (lik1-lik2) / (double) pairs;
      ++repcount[S][R];
      if (repcount[S][R]>=(REPMAX-2)) {
	printf ("Error: repcount too high\n");
	exit(0);
      }
    }
  }
  fclose (simlikfile);

  while (scanf ("%d %lf %lf %lf %lf %lf\n", &pairs, &rho, &lik1, &lik2, &lik3,
		&lik4) != EOF) {
    S = (int) floor (sqrt (2.*pairs));
    for (R=0; R<RMAX; ++R)
      if (rho<(rgrid[R]+1E-08))
        break;
    ++R;

    num = den = 0;
    if (S>0 && S<=SMAX && R>0 && R<RMAX) {
      for (s=1; s<=SMAX; ++s)
	for (r=1; r<RMAX; ++r)
	  if (s >= w*S && s <= x*S && r >= y*R && r <= z*R) {
	    den += repcount[s][r];
	    testlik = (lik3-lik4) / pairs;
	    for (a=0; a<repcount[s][r]; ++a)
	      if (testlik <= logliks[s][r][a])
		++num;
	  }
      if (den>0 && num>0)
	printf("%d %.1f %.1f %.1f %d %d %f\n", S, rho, lik1, lik2, den,
	       num, (double) num / (double) den);
      else if (den>0 && num==0)
	printf("%d %.1f %.1f %.1f %d %d %f\n", S, rho, lik1, lik2, den,
	       num, 1. / (double) den);
      else
	printf("%d %.1f %.1f %.1f 0 0 1.\n", S, rho, lik1, lik2);
    }
    else 
      printf("%d %.1f %.1f %.1f 0 0 1.\n", S, rho, lik1, lik2);
  }
}
