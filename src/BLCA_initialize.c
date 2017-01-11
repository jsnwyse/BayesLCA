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


#include "BLCA_initialize.h"

void BLCA_set_prior_on_number_of_components(struct mix_mod *mixmod,int type)
{

	int i;
	
	switch(type)
	{
	
		case RICHARDSON_AND_GREEN:
		
			for(i=0;i<mixmod->maxgroups+1;i++){
				mixmod->log_prior_G[i] = 0.;
			}
		
		break;
		
		case NOBILE_AND_FEARNSIDE:
		
			for(i=1;i<mixmod->maxgroups+1;i++){
				mixmod->log_prior_G[i] = -lgamma((double)i + 1.);
			}
		
		break;
	
	}
	

	return;

}

void BLCA_initialize_data( struct mix_mod *mixmod, int *Y )
{
	int i, j, x, nobs = mixmod->n, nvar = mixmod->d ;
	
	for(i=0;i<nobs;i++)
	{
		for(j=0;j<nvar;j++)
		{
			x = Y[i + j*(nobs)];
			mixmod->Y[j][i] = x;
			mixmod->Yobs[i][j] = x;
		}
	}
}

int BLCA_initialize_simple(struct mix_mod *mixmod,int numgroups)
/*gives a very simple initialization of the model by just dumping the
	first n/G observations in the first group, second n/G in second group...*/
{
	int i,j,jj,k,G=numgroups,d=mixmod->d,n=mixmod->n;

	int m = n/G; /*gives the number of segments*/
	for(k=0;k<G-1;k++){
	
		for(i=0;i<d;i++){
			mixmod->varindicator[i] = 1;
			for(j=0;j<mixmod->ncat[i];j++)
				mixmod->components[k]->N[i][j] = 0;
		}
	
		/*cycle through appropriate data*/
		for(i=k*m;i<(k+1)*m;i++){
			mixmod->z[i] = k;
			mixmod->components[k]->n_g = m;
			for(j=0;j<d;j++){ 
				mixmod->components[k]->N[j][ mixmod->Y[j][i] ] += 1; 
			}
		}
	}
	
	/*special case for last group*/

	for(i=0;i<d;i++){
		for(j=0;j<mixmod->ncat[i];j++)
			mixmod->components[G-1]->N[i][j] = 0;
	}
	
	for(i=(G-1)*m;i<n;i++){
		mixmod->z[i] = G-1;
		mixmod->components[G-1]->n_g = n-(G-1)*m;
		for(j=0;j<d;j++){
			mixmod->components[G-1]->N[j][ mixmod->Y[j][i] ] += 1;
		}
	}

	
	/*compute the log_prob for each of the components*/
	for(k=0;k<G;k++){
		BLCA_recompute_marginal_likelihood_component(mixmod->components[k],mixmod);
	}	
	
	for(i=0;i<n;i++){
		for(j=0;j<d;j++){
			mixmod->undiscriminating->N[j][ mixmod->Y[j][i] ] += 1;
		}
	} 
	
	BLCA_recompute_marginal_likelihood_undiscriminating_variables(mixmod->undiscriminating,mixmod);

	for(k=0;k<G;k++){
		mixmod->whereis[k] = k;
	}
	
	if(FALSE)
	{
		//initialize the weights and vectors of probabilities also
		double s = 0.;
		
		
	
	}

	return(TRUE);
}

void BLCA_initialize_EM( struct mix_mod *mixmod , double *group_weights, double *prob_variables)
{
	int g, j, c, p=0, gap_;
	double s=0.;
	
	//initialize everything randomly
	for( g=0; g<mixmod->G; g++ ) 
	{
		mixmod->weights[g] = rgamma( 1. , 1. ) ;
		s += mixmod->weights[g];
	}
	for( g=0; g<mixmod->G; g++ ) mixmod->weights[g] /= s;
	
	for( g=0; g<mixmod->G; g++ )
	{
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				//initialize the variable probabilities randomly
				s = 0.;
				for( c=0; c<mixmod->ncat[j]; c++ )
				{
					mixmod->components[g]->prob_variables[j][c]  = rgamma( 1., 1.);
					s += mixmod->components[g]->prob_variables[j][c] ;
				}
				for( c=0; c<mixmod->ncat[j]; c++ )
					mixmod->components[g]->prob_variables[j][c] /= s ;
			}
		}
	}
	return;
}


void BLCA_initialize_VB( struct mix_mod *mixmod , double *alpha_ud, double *beta_ud )
{
	int g, j, c, p=0, gap_;
	double s=0.;
	
	//initialize everything randomly
	for( g=0; g<mixmod->G; g++ ) 
		mixmod->alpha_ud[g] = rgamma( 1. , 1. ) ;
	
	
	for( g=0; g<mixmod->G; g++ )
	{
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				for( c=0; c<mixmod->ncat[j]; c++ )
					mixmod->components[g]->beta_ud[j][c]  = rgamma( 1., 1.);
			}
		}
	}
	return;
}



