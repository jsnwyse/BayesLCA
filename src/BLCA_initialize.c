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
			
	Last modification of this code: Tue 08 Aug 2017 21:04:52 IST   */


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


void BLCA_initialize_priors( struct mix_mod *mixmod, double *alpha_prior, double *beta_prior, int type )
{
	int k, j, c, p, G = mixmod->maxgroups, d = mixmod->d ;
	
	// three different types- type 0 is the simple initialization
	if( type == 0 )
	{
		// one parameter for each; these values have already been set above
		//mixmod->alpha = alpha_prior[0];
		//mixmod->beta = beta_prior[0];
	}
	
	if( type == 1 )
	{
		//all vectors
		//all values take the same
		for( k=0; k<G; k++ ) mixmod->alpha_prior[k] = alpha_prior[k];
		p = 0;
		for( j=0;  j<d; j++ )
		{
			for( k=0; k<G; k++ )
			{
				for( c=0; c<mixmod->ncat[j]; c++ ) mixmod->beta_prior[k][j][c] = beta_prior[ p + c ];
			}
			p += mixmod->ncat[j] ;
		}	
	}
	
	if( type == 2 )
	{
		//all take different values
		for( k=0; k<G; k++ ) mixmod->alpha_prior[k] = alpha_prior[k];
		p = 0;
		for( j=0;  j<d; j++ )
		{
			for( k=0; k<G; k++ )
			{
				for( c=0; c<mixmod->ncat[j]; c++ ) mixmod->beta_prior[k][j][c] = beta_prior[ p*G + k*mixmod->ncat[j] + c ];
			}
			p += mixmod->ncat[j] ;
		}	
	}
	
	//trying to fix problem with eject absorb moves when prior_type==1
	mixmod->prior_type = 1;
	
	return;
	
}

void BLCA_initialize_Gibbs_sampler( int init_type, struct mix_mod *mixmod  )
{
	int g, k, i, j, c, p=0, gap_;
	double s=0.;
	
	// go through three possible initialization types making the s 
	// matrix in each case
	
	// option 'single' and membership initialization for init_type == 3
	
	if( mixmod->collapsed || init_type == 0 || init_type == 2 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			//generate a random group 0:G-1
			g = (int)( mixmod->G * runif(0.,1.) ) ;
			mixmod->z[k] = g;
			BLCA_add_to_component( mixmod->components[g], mixmod->Yobs[k] , mixmod, 1 );
			/*mixmod->components[g]->n_g += 1;
			for( j=0; j<mixmod->d; j++ )
			{
				if( mixmod->varindicator[j] )
					mixmod->components[g]->N[ j ][ mixmod->Y[j][k] ] += 1 ;
			}*/
			//to give consistency with previous version
			if( !mixmod->collapsed )
			{
				mixmod->s[k][g] = 1.;
			}	
		}	
	}
	
	// option 'across'
	if( init_type == 1 )
	{
		// generate a portion of a random group
		for( k=0; k<mixmod->n; k++ )
		{
			s = 0.;
			for( g=0; g<mixmod->G; g++ ) 
			{
				mixmod->s[k][g] = runif( 0., 1.);
				s += mixmod->s[k][g];
			}
			for( g=0; g<mixmod->G; g++ ) mixmod->s[k][g] /= s; 
		}			
	}
	
	
	if( init_type == 0 || init_type == 1 )
	{
		for( g=0; g<mixmod->G; g++ )
		{
			for( k=0; k<mixmod->n; k++ )
			{
				mixmod->weights[g] += mixmod->s[k][g] ;
			}
		}
		
		for( g=0; g<mixmod->G; g++ )
		{
			for( j=0; j<mixmod->d; j++ )
			{
				if( mixmod->varindicator[j] )
				{
					for( k=0; k<mixmod->n; k++ )
						mixmod->components[g]->prob_variables[j][ mixmod->Y[j][k]  ] += mixmod->s[k][g] ;
					
					
					for( c=0; c<mixmod->ncat[j]; c++ ) 
						mixmod->components[g]->prob_variables[j][c] /= mixmod->weights[g] ; 
					
				}
			}
		}
		
		for( g=0; g<mixmod->G; g++ ) mixmod->weights[g] /= mixmod->n; 
	
	}
	
	// option 'prior'
	
	if( init_type == 2 )
	{
		//generate the probabilities from the prior
		
		//component weights
		s = 0.;
		for(i=0;i<mixmod->G;i++){
			mixmod->weights[i] = rgamma( mixmod->alpha_prior[i] , 1.  ) ;
			s += mixmod->weights[i];
		}
		for(i=0;i<mixmod->G;i++) mixmod->weights[i] /= s;
		
		//variable probabilities
		for(i=0;i<mixmod->G;i++){
			for(j=0;j<mixmod->d;j++){
			
				if( mixmod->varindicator[j] )
				{
					//only sample if the variable is included in the model
					s = 0.;
					for(k=0;k<mixmod->ncat[j];k++){
						 mixmod->components[i]->prob_variables[j][k] =  rgamma( mixmod->beta_prior[i][j][k] , 1. ) ;
						 s += mixmod->components[i]->prob_variables[j][k];
					}
					for(k=0;k<mixmod->ncat[j];k++) 
						mixmod->components[i]->prob_variables[j][k] /= s;

					
				}
			}
		}
		
	}
	
	if( mixmod->collapsed )
	{
		mixmod->component_compute = 0;
	
		//compute log marginal likelihood for each component
		for( g=0; g<mixmod->G; g++ )
			BLCA_recompute_marginal_likelihood_component( mixmod->components[g], mixmod );
	
		//fill up the part for the undiscriminating variables
		for( k=0; k<mixmod->n; k++ ) BLCA_add_to_component( mixmod->undiscriminating, mixmod->Yobs[k], mixmod, 1 );
		
		/*for( j=0; j<mixmod->d; j++ )
		{
			for( k=0; k<mixmod->n; k++ )
				mixmod->undiscriminating->N[j][ mixmod->Y[j][k] ] += 1;
		} */
	
		BLCA_recompute_marginal_likelihood_undiscriminating_variables( mixmod->undiscriminating, mixmod );
	}

	for( g=0; g<mixmod->G; g++ )
		mixmod->whereis[g] = g;
	
	return;
}


void BLCA_initialize_EM( int init_type, double *init_vals, struct mix_mod *mixmod , double *group_weights, double *prob_variables)
{
	int g, k, j, c, p=0, gap_;
	double s=0.;
	struct component *comp;
	
	//initialize the membership matrix using either across/single or prespecified
	
	if( init_type == 0 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			//generate random group 0:G-1
			g = (int)( mixmod->G * runif( 0.,1.) ) ;
			mixmod->s[k][g] = 1.;
		}
	}
	
	if( init_type == 1 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			s = 0.;
			for( g=0; g<mixmod->G; g++ ) 
			{
				mixmod->s[k][g] = runif( 0., 1.);
				s += mixmod->s[k][g];
			}
			for( g=0; g<mixmod->G; g++ ) mixmod->s[k][g] /= s; 
		}
	}
	
	if( init_type == 2 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			for( g=0; g<mixmod->G; g++ ) mixmod->s[k][g] = init_vals[k + g * mixmod->n ] ; 
		}
	
	}
	
	return;
}


void BLCA_initialize_VB( int init_type, double *init_vals, struct mix_mod *mixmod , double *alpha_ud, double *beta_ud )
{
	int k, g, j, c, p=0, gap_;
	double s=0.;
	struct component *comp;
	
	if( init_type == 0 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			//generate random group 0:G-1
			g = (int)( mixmod->G * runif( 0.,1.) ) ;
			mixmod->s[k][g] = 1.;
		}
	}
	
	if( init_type == 1 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			s = 0.;
			for( g=0; g<mixmod->G; g++ ) 
			{
				mixmod->s[k][g] = runif( 0., 1.);
				s += mixmod->s[k][g];
			}
			for( g=0; g<mixmod->G; g++ ) mixmod->s[k][g] /= s; 
		}
	}
	
	if( init_type == 2 )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			for( g=0; g<mixmod->G; g++ ) mixmod->s[k][g] = init_vals[k + g * mixmod->n ] ; 
		}
	
	}
	
	//initialize everything randomly
	mixmod->di_sum_alpha_ud = 0.;
	for( g=0; g<mixmod->G; g++ ) 
	{
		mixmod->alpha_ud[g] = rgamma( 1. , 1. ) ;
		mixmod->di_sum_alpha_ud += mixmod->alpha_ud[g];
		mixmod->di_alpha_ud[g] = digammaRN( mixmod->alpha_ud[g] );
	}
	mixmod->di_sum_alpha_ud = digammaRN( mixmod->di_sum_alpha_ud );
	
	for( g=0; g<mixmod->G; g++ )
	{
		comp = mixmod->components[g] ;
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				comp->di_sum_beta_ud[j] = 0.;
				for( c=0; c<mixmod->ncat[j]; c++ )
				{
					comp->beta_ud[j][c]  = rgamma( 1., 1.);
					comp->di_sum_beta_ud[j] += comp->beta_ud[j][c] ;
					comp->di_beta_ud[j][c] = digammaRN( comp->beta_ud[j][c] );
				}
				comp->di_sum_beta_ud[j] = digammaRN( comp->di_sum_beta_ud[j] ) ;
			}
		}
	}
	return;
}

void BLCA_initialize_Boot( struct mix_mod *mixmod, double *group_weights, double *prob_variables  )
{

	int i, j, k, c, p, gap_;
	
	struct component *comp;
	
	for( i=0; i<mixmod->G; i++ ) mixmod->weights[i] = group_weights[i]; 
	
	p = 0;
			
	for(j=0;j<mixmod->d;j++){
		gap_ = p * mixmod->G  ;
		for(i=0;i<mixmod->G;i++){
			for(k=0;k<mixmod->ncat[j];k++){
				mixmod->components[i]->prob_variables[j][k] = 
				  prob_variables[ //long index expression
							
					gap_ + i*mixmod->ncat[j] + k 
								
					] ;
			}
		}
		p += mixmod->ncat[j];
	}

	return;
	
}



