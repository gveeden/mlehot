/* lnlikehot.c - Modified from lnlikegc.c to allow for a single hotspot
in the center.  To be used with mlehot.c with no gene conversion */

#include <math.h>

struct pairst {
  int n;
  int n1;
  int n2;
  int n11;
  double  d12;
  double hot12;
};

extern struct pairst *pair ;
 
extern npairs;


double lnlikehot (double co, double hot, int nsites, double *recrates, double ****prob, double *ppoly  )
{
  int i, n, n1, n2, n11 ;
  double getprob(int,int,int,int, double, int,double *, double ****, double *);
  double sum, d12, hot12, trec;


  sum = 0.0 ;
  for( i =0 ; i< npairs; i++){
    n = pair[i].n; n1 = pair[i].n1 ; n2 = pair[i].n2 ; n11 = pair[i].n11; 
    d12 = pair[i].d12; hot12 = pair[i].hot12;
    trec = d12 * co + hot12 * co * hot;
    /*    printf("%d %d %d %d %lf %d\n", n, n1, n2, n11, trec, nsites); */

    sum += getprob(n,n1,n2,n11, trec, nsites,recrates, prob, ppoly ) ;
  }
  /*
  printf("sum = %lf\n", sum);
    exit(0);
  */
  return( sum ) ;
}

