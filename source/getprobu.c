#include "my.h"
/*
void main(void)
{
	int  n1, n2, n12, n ;
	double  getprob( int, int, int, int, double) ;
	double r ;
	
	
	while( 1){
		printf (" n n1 n2 n12 r?\n");
		scanf("%d %d %d %d %lf",&n, &n1, &n2, &n12, &r);
		printf(" %e\n", getprob(n, n1, n2, n12, r) ) ;
		}
		


}
*/




double
getprob(int n, int n1, int n2, int n11, double r, int nsites, double *recrates,
	 double ****prob, double *ppoly )
{
	int m1, m2,  k  ;
	double lnval( double), extrap( double, double, double, double, double ) ; 
	double nadprob(int, int, int, int, int, double **** ) ;

 
	m1 = MAX( n1, n2);
	m2 = MIN( n1, n2) ;

	for( k=0; k<nsites; k++){
	  if( recrates[k] == r ) return( lnval( nadprob(n, k,m1,m2,n11, prob)/ppoly[k]  )  ) ;
	  if( recrates[k] >r )  return( extrap(nadprob(n,k-1,m1,m2,n11, prob)/ppoly[k-1] ,
				nadprob(n,k,m1,m2,n11,prob)/ppoly[k] , recrates[k-1],recrates[k],r));
	  }
	return( lnval( nadprob(n,nsites-1,m1,m2,n11,prob)/ppoly[nsites-1] ) ) ;

}

	double
nadprob(int n,  int k, int m1, int m2, int n11, double ****prob )
{
	int n00, n01, n10 ;
	double sum  ;

	n10 = m1 - n11 ;
	n01 = m2 - n11 ;
	n00 = n - n10 - n01 - n11 ;
	if( (n00 == n01) && ( n01 == n10 ) && ( n10 == n11 ) ) return( prob[k][m1][m2][n11] ) ;
	else if( (n00 == n01) && ( n10 == n11) ) 
	  return(  prob[k][m1][m2][n11] + prob[k][m2][n-m1][n01] ) ;
	else if( ( n00 == n10 ) && ( n01 == n11 ) ) 
	  return(  prob[k][m1][m2][n11] + prob[k][n-m2][m1][n10] ) ;
        else if( ( n00 == n11 ) && ( n01 == n10 ) ) 
          return(  prob[k][m1][m2][n11] + prob[k][n-m2][n-m1][n00] ) ;
        else {
	sum = prob[k][m1][m2][n11] ;
	sum += prob[k][ MAX( m2, n-m1) ][ MIN( m2, n-m1)][n01] ;
	sum += prob[k][ MAX( m1, n-m2) ][ MIN( m1, n-m2)][n10] ;
	sum += prob[k][ MAX( n-m1, n-m2) ][ MIN( n-m1, n-m2)][n10] ;
	 
  return(  sum) ;
	}

}

	double
lnval( double p2)
{
	if( p2 > 0.0 ) return( log(p2 ) ) ;
	else return( -9999. ) ;
}


	double 
extrap(  double p1, double p2, double r1, double r2, double r )
{
	double linextrap( double, double, double, double, double), lnval( double) ;

  if( (r1 == 0 ) && ( (p1 == 0) || (p2 == 0.0) ) ) {
	    return( lnval( linextrap(p1, p2, r1, r2, r) ) ) ;
	    }
  else if( r1 == 0.0 )  return( linextrap( log(p1), log(p2), r1, r2, r ) ) ;
  else if( (p1 == 0.0) || ( p2 == 0.0 ) ) return( lnval( linextrap( p1, p2, log(r1),log(r2), log(r))));  
  else return(  linextrap( log(p1), log(p2), log(r1), log(r2), log(r) ) ) ;
}


	double 
linextrap(  double p1, double p2, double r1, double r2, double r )
{
	return(  p1 + (p2-p1)*(r-r1)/(r2-r1) ) ;  

}
