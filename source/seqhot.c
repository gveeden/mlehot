/* seqhot.c - Input is from program ms0 or equivalent.  Output is to be used 
by mlehot for pairwise composite likelihood estimation of rho.  Usage:

seqhot n hotbeg hotend

n = haploid sample size
hotbeg = position of left boundary of putative hotspot (scaled from 0 to 1)
hotend = position of right boundary of putative hotspot (scaled from 0 to 1)

Output: First line is number of pairs.  After that, each output consists of

n n2 n1 n12 x1 x2

n2 = derived allele frequency of second site in pair
n1 = derived allele frequency of first site in pair
n12 = number of haplotypes with both derived alleles
x1 = distance between sites of non-hotspot sequence
x2 = distance between sites of putative hotspot sequence

Compilation: cc -O -o seqhot seqhot.c

Note: File my.h must be present for compilation. 

Modified from a program of R. Hudson.
*/


#include "my.h"
int maxsites = 100 ;

main(int argc,char *argv[])
{
  int nsam, j ,nsites, i, ns, howmany, npop, *config, pop, count ;
  double   r, mig_rate,p, ave, x ;
  char **list, **cmatrix(), allele,na ;
  FILE *pf, *fopen(), *pfin ;
  double told, ttot, theta ,*posit, *tmrca, *tmut, minseg, t  ;
  int   segsites, freq, freq2(), fsite, fsite2, site, site1, f ;
  int n1, n2, n12, temp , allele1, allele2, seed1, seed2, seed3 ;
  double t0, alpha  ;
  
  double crate, clen;
  double hotbeg, hotend, b1, b2, x1, x2;
  
  THINK_C_CL ;
  pfin = stdin ;

  if (argc != 4) {
    printf ("Usage: seqhot n hotbeg hotend\n");
    exit(0);
  }
  nsam = atoi (argv[1]);
  hotbeg = atof (argv[2]);
  hotend = atof (argv[3]);

  list = cmatrix(nsam,maxsites+1);
  posit = (double *)malloc( maxsites*sizeof( double ) ) ;

  while (fscanf(pfin," %d", &segsites ) != EOF) {
    if( segsites >= maxsites){
      maxsites = segsites + 10 ;
      posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
      biggerlist(nsam,maxsites, list) ;
    }
    /*	   readsamples */
    if( segsites > 0) {
      for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
      for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
    }
       
    printf(" %d\n",( segsites*(segsites-1) )/2 ) ;
    allele2 = '1' ;
    for( fsite = 0 ;  fsite<(segsites-1); fsite++ ){
      allele1 = '1' ;
      n1 = frequency(allele1,fsite,nsam, list);
      for( fsite2 = fsite+1; fsite2<segsites; fsite2++){
	n2 = frequency(allele2,fsite2,nsam, list);
	n12 = freq2(allele1,allele2,fsite, fsite2, nsam, list);
	x = posit[fsite2] - posit[fsite] ; 
	
	b1 = (posit[fsite2] < hotend) ? posit[fsite2] : hotend;
	b2 = (posit[fsite] > hotbeg) ? posit[fsite] : hotbeg;
	if (b1 < b2)
	  x2 = 0.;
	else
	  x2 = b1 - b2;
	x1 = x - x2;

	printf("%d %d %d %d %lf %lf\n", nsam, n1, n2, n12, x1, x2 ); 
      }
    }
  }
}

	

	int
frequency(allele,site,nsam, list)
	int site, nsam;
	char **list, allele ;
{
	int i, count=0;
	for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}

	int
freq2(allele1,allele2,site,site2,nsam, list)
	int site, nsam;
	char **list, allele1, allele2 ;
{
	int i, count=0;
	for( i=0; i<nsam; i++) {
		if(  ( list[i][site] == allele1 ) && (  list[i][site2] == allele2)  ) count++ ;
		}
	return( count);
}



/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}

        int
biggerlist(nsam, nmax, list )
        int nsam ;
        unsigned nmax ;
        char ** list ;
{
        int i;

        maxsites = nmax  ;
        for( i=0; i<nsam; i++){
           list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
           if( list[i] == NULL ) perror( "realloc error. bigger");
           }
}                        
