/*Functions for the fitting of Latent Class Analysis models using MCMC
	methods. Two implementations of the model are included in the Bayesian
	formulation: collapsed and not collapsed.
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie
			
	Last modification of this code: Mon 09 May 2016 15:38:03 IST   */
	
#include "BLCA_utils.h"
	
double BLCA_get_max(double *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	double max=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]>max)
				max = x[i];
		}
	
	}
	return(max);
}

int BLCA_get_imax(int *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	int max=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]>max)
				max = x[i];
		}
	
	}
	return(max);
}

double BLCA_get_min(double *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	double min=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]<min)
				min = x[i];
		}
	
	}
	return(min);
}

int BLCA_get_imin(int *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	int min=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]<min)
				min = x[i];
		}
	
	}
	return(min);
}



int BLCA_sample_discrete( double *weights, int len )
{
	/*sample once from a multinomial distribution with weights*/
	
	int i=0;
	double w , u;
	
	u = runif(0.0,1.0) ;
	
	w = weights[0];
	
	while( w < u && i < len )
	{
		i++ ;
		w += weights[i];		
	}
	
	return(i);	
}

int BLCA_random_integer( int n )
{
	//generate a random integer between 0 , ..., n-1
	int i;
	double *w =  calloc(n,sizeof(double));
	for(i=0;i<n;i++) w[i] = 1.0/n ;
	i = BLCA_sample_discrete( w, n);
	free(w);
	return(i);
}


/* a util to permute a vector  */
void BLCA_random_ranshuffle( int *a, int n )
{
	//randomly permute the first n elements of a using Fisher-Yates algorithm
	//  using RNG from R
	
	int i, j, x;
	
	for( i=n-1; i>0; i-- )
	{
		j = (int)( i * runif(0.0,1.0) ) ;
		x = a[i] ;
		a[i] = a[j] ;
		a[j] = x;
	}
	
	return;
}

double BLCA_get_log_sum_exp( double *x, int len )
{
	//note that this function also modifies the x argument

	int i;
	double s=0., max=-DBL_MAX;

	max = BLCA_get_max( x, len );
	
	for( i=0; i<len; i++ )
	{
	  x[i] -= max;
	  x[i] = exp(x[i]);
	  s += x[i];
	}
	
	for( i=0; i<len; i++ ) x[i] /= s;
	
	return( max + log(s) );
}

double BLCA_get_log_sum_exp_all( double *x, int len )
{
  //note that this function also modifies the x argument
  
  int i;
  double s=0., max=-DBL_MAX;
  
  max = BLCA_get_max( x, len );
  
  for( i=0; i<len; i++ )
  {
    x[i] -= max;
    x[i] = exp(x[i]);
    s += x[i];
  }
  
  for( i=0; i<len; i++ ) x[i] /= s;
  
  return( max + log(s) );
}


