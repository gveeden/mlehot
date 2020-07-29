
#include "my.h"




	double ****getprobmat(char *fname, int *pnsam, int *pnsites, double **precrates, double **pppoly )
{

  int  i ; 
  double     alpha  ;
  FILE *pf, *fopen() ;
  int kmax, k2,   j, k, m  ;
  double ****prob;
  int nsam, nsites, npop, *config, pop;
  double alphag;
  double x , y , *pp ;
  int d1, d2 , ng ;
  
  pf = fopen( fname, "r");
  
  fscanf(pf, " %d", &npop);
  /*  printf("npop = %d\n", npop); */
  config = (int *)malloc((unsigned)((npop+1)*sizeof(int)));
  for( pop=0  ; pop<npop; pop++) {
    fscanf(pf, " %d", config+pop ) ;
    if( config[pop] < 0 ) break;
  }
  fscanf(pf," %d", &nsites) ;
  *pnsites = nsites ;
  *precrates = (double *)malloc( (unsigned)nsites*sizeof(double) ) ;
  pp  = (double *)malloc( (unsigned)nsites*sizeof(double) ) ;
  *pppoly = pp ;
  for( i=0; i<nsites; i++){
    fscanf(pf," %lf", *precrates + i ) ;
  }
  /*
    for (i=0; i<nsites; ++i)
    printf("%lf ", *precrates[i]);
    printf("\n");
  */
  fscanf(pf," %lf", &alphag ) ;
  if( alphag != 0.0 ) fscanf(pf, " %*lf" ) ;
  if( npop > 1 )  fscanf(pf," %*lf" ) ;
  fscanf(pf, " %*d" ) ;
  for(i=0;i<3;i++){
    fscanf(pf," %*d");
  }
  
  while( pop<npop ){ config[pop++]=0 ;  }
  nsam = 0 ;
  for(i=0;i<npop; i++) nsam += config[i] ;
  *pnsam = nsam ;
  /*  printf("nsam = %d\n", nsam); */
  
  prob = (double ****)malloc( (unsigned)nsites*sizeof( double *** ) ) ;
  for( i=0; i<nsites; i++){
    prob[i] = (double ***)malloc( (unsigned)nsam*sizeof( double **) ) ;
    for( j=0; j<nsam; j++) {
      prob[i][j] = (double **)malloc( (unsigned)(j+1)*sizeof( double *) ) ;
      for( k=0; k<=j ; k++) {
	
	prob[i][j][k] = (double *)malloc( (unsigned)( MIN(j,k) - MAX( 0, (j+k-nsam) )+1 )*sizeof( double) ) ;
	if( prob[i][j][k] == 0 ) printf( " whoaa!\n" ) ;
	prob[i][j][k] -= MAX( 0, j+k-nsam ) ;
      }
    }
  }

  for( m=0; m<nsites; m++) pp[m] = 0.0 ;

  for( i=1; i<nsam; i++)
    for( j=1; j<= i ; j++){
      for (m=0; m<nsites; m++) {
	fscanf(pf," freq: %d %d", &d1, &d2) ;
	ng =  fscanf(pf," rec rate: %lf", &y) ;
	for( k= MAX(0,(i+j-nsam))  ; k<=  MIN(i,j)  ; k++) {
	  fscanf(pf," %le", &( prob[m][i][j][k]) ) ;
	  /*	  printf("%lf %d %d %d %d\t", prob[m][i][j][k],m,i,j, k); */
	  if( i != j ) pp[m] += 2.0*prob[m][i][j][k] ;
	  else pp[m] += prob[m][i][j][k] ;
	}
	
      }
    }
  fclose(pf) ;

/*

	for( i=1; i<= nsam/2; i++){
	  for( j=1; j<=i; j++){
		printf(" %d %d\n", i, j ) ;
		kmax = MIN( i, j) - MAX( 0, (i+j-nsam) ) ;
	    for( m=0; m<nsites; m++){
		printf(" %lf ", (*precrates)[m]) ;
		k2 = kmax ;
          for( k= 0 ; k<= ( MIN(i,j) - MAX(0,(i+j-nsam))) ; k++) {
	      if( (i == j)&&(i == nsam/2) ) {
		printf(" %.4e", prob[m][i][j][k]+prob[m][i][j][kmax-k]); } 
	      else if( (i == j)&&(i != nsam/2) ) {
		printf(" %.4e", prob[m][i][j][k] + prob[m][nsam-j][nsam-i][k]+
			2.0*prob[m][nsam-i][j][kmax-k]);}
	     else {
 printf(" %.4e", prob[m][i][j][k] +
			prob[m][nsam-j][nsam-i][k] + prob[m][nsam-j][i][k2] +
			prob[m][ MAX( nsam-i , j)][MIN( nsam-i,j)][k2] ) ;
		   k2--;
		 }
	     }
		printf("\n");
		}}}

*/

	return( prob ) ;
}



