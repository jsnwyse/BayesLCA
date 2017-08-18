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


#include "BLCA_component.h"

void BLCA_allocate_component(struct component *component,struct mix_mod *mixmod)
{
	/*allocate memory for a component*/
	int j;
	
	component->N = calloc(mixmod->d,sizeof(int *));
	
	for(j=0;j<mixmod->d;j++){
		component->N[j] = calloc(mixmod->ncat[j],sizeof(int));
	}
	
	if( !mixmod->collapsed )
	{
		component->prob_variables = calloc(mixmod->d,sizeof(double *));
		for(j=0;j<mixmod->d;j++){
			component->prob_variables[j] = calloc(mixmod->ncat[j],sizeof(double));
		}
	}
	
	if( mixmod->VB )
	{
		component->beta_ud = calloc( mixmod->d, sizeof(double *));
		component->di_beta_ud = calloc( mixmod->d, sizeof(double *));
		component->di_sum_beta_ud = calloc( mixmod->d, sizeof(double));
		for(j=0;j<mixmod->d;j++){
			component->beta_ud[j] = calloc(mixmod->ncat[j],sizeof(double));
			component->di_beta_ud[j] = calloc(mixmod->ncat[j],sizeof(double));
		}
	}
	
	return;
}

void BLCA_free_component(struct component *component,struct mix_mod *mixmod)
{
	
	int j;
	
	for(j=0;j<mixmod->d;j++){
		free(component->N[j]);
	}
	free(component->N);
	
	if( !mixmod->collapsed )
	{
		for(j=0;j<mixmod->d;j++){
			free(component->prob_variables[j]);
		}
		free(component->prob_variables);
	}
	
	if( mixmod->VB )
	{
		for(j=0;j<mixmod->d;j++){
			free(component->beta_ud[j]);
			free(component->di_beta_ud[j]);
		}	
		free(component->beta_ud);
		free(component->di_beta_ud);
		free(component->di_sum_beta_ud);
	}
	
	return;
	
}

void BLCA_copy_component(struct component *component_original,struct component *component_target,struct mix_mod *mixmod)
/*copy the contents of the first argument into the second component argument*/
{

	int i,j;
	
	component_target->n_g = component_original->n_g;
	for(i=0;i<mixmod->d;i++){
		for(j=0;j<mixmod->ncat[i];j++){
			component_target->N[i][j] = component_original->N[i][j];
		}
	}
	
	if( !mixmod->collapsed )
	{
		for( i=0; i<mixmod->d; i++ )
		{
			for( j=0; j<mixmod->ncat[i]; j++ )
				component_target->prob_variables[i][j] = component_original->prob_variables[i][j];
		}
	}

	component_target->log_prob = component_original->log_prob;
	
	return;

}


void BLCA_add_to_component( struct component *component, int *x, struct mix_mod *mixmod, int sgn )
{
	int c, j;
	
	component->n_g += sgn ;
	
	for( j=0; j<mixmod->d; j++ ) component->N[j][ x[j] ] += sgn ;
	
	return;
}


void BLCA_recompute_sufficient_statistics_for_components(struct mix_mod *mixmod)
{

	int i,j,k;
	
	/*undiscriminating variables*/
	
	
	
	/*discriminating variables*/
	
	for(k=0;k<mixmod->G;k++){	
		/*reset*/
		for(j=0;j<mixmod->d;j++){
			for(i=0;i<mixmod->ncat[j];i++){
				mixmod->components[k]->N[j][i] = 0;
				mixmod->undiscriminating->N[j][i] = 0;
			}
		}
		mixmod->components[k]->n_g = 0;
	}
	
	int z;
	
	for(i=0;i<mixmod->n;i++){
	
		z = mixmod->z[i];	
		
		mixmod->components[ z ]->n_g += 1;
		
		for(j=0;j<mixmod->d;j++){
			mixmod->components[ z ]->N[j][ mixmod->Y[j][i] ] += 1;
			mixmod->undiscriminating->N[j][ mixmod->Y[j][i] ] += 1;
		}
	
	
	}
	
	return;
}


double BLCA_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{

	int i,j,c,I,k=mixmod->component_compute;
	double log_prob = 0., s, *beta;
	
	if( mixmod->prior_type == 0 )
	{
	
		for(j=0;j<mixmod->d;j++){
	
			if(mixmod->varindicator[j]){
	
				log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma( component->n_g+1. + mixmod->ncat[j]*mixmod->beta );
		
				for(i=0;i<mixmod->ncat[j];i++){
		
					I = (x[j] == i) ?  1 : 0 ;
		
					log_prob += lgamma(component->N[j][i] + I + mixmod->beta);
		
				}
		
			}
		
		}
	
	}else if( mixmod->prior_type == 1 ){
		
		for( j=0; j<mixmod->d; j++ ){
			
			if( mixmod->varindicator[j]){
			
				s =  0.; 
				
				beta = mixmod->beta_prior[k][j] ;
				
				for( c=0; c<mixmod->ncat[j]; c++ ){
					
					I = (x[j] == c) ? 1 : 0 ;
					log_prob += lgamma( component->N[j][c] + I + beta[c] ) - lgamma( beta[c] ) ;
					
					s += beta[c] ;
				}
				
				log_prob += lgamma(s) - lgamma( component->n_g + 1. + s ) ; 
			
			}
		
		}
	
	}

	return(log_prob);
}

double BLCA_compute_log_marginal_likelihood_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{
	/*function used for class prediction*/

	int i,j,I;
	double log_prob = 0.;
	
	log_prob = lgamma(component->n_g+1.+mixmod->alpha);
	
	for(j=0;j<mixmod->d;j++){
	
		if(mixmod->varindicator[j]){
	
			log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma( component->n_g + 1. + mixmod->ncat[j]*mixmod->beta );
		
			for(i=0;i<mixmod->ncat[j];i++){
		
				I = (x[j] == i) ?  1 : 0 ;
		
				log_prob += lgamma(component->N[j][i] + I + mixmod->beta);
		
			}
		
		}
		
	}

	return(log_prob);
}


double BLCA_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod)
{

	int i,j,c,k=mixmod->component_compute;
	double log_prob = 0., s, *beta;
	
	if( mixmod->prior_type ==  0 )
	{
	
		for(j=0;j<mixmod->d;j++){
	
			if(mixmod->varindicator[j]){
		
				log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + mixmod->ncat[j]*mixmod->beta );
		
				for(i=0;i<mixmod->ncat[j];i++){
		
					log_prob += lgamma(component->N[j][i] + mixmod->beta);
		
				}
		
			}
		
		}
	
	}else if( mixmod->prior_type == 1 ){
	
		
		for( j=0; j<mixmod->d; j++ ){
		
			if( mixmod->varindicator[j] ){
			
				s = 0.;
				
				beta = mixmod->beta_prior[k][j];
				
				for( c=0; c<mixmod->ncat[j]; c++ ){
					
					log_prob += lgamma( component->N[j][c] + beta[c] ) - lgamma( beta[c] ) ;
					
					s += beta[c] ;
				}
				
				log_prob += lgamma(s) - lgamma( component->n_g + 1. + s ) ; 
			
			}
		
		}
	
	
	}
	
	/*if(isnan(log_prob)){
		mixmod_warning(1);
		printf("\nFrom compute_log_data_probability_component\n");
	}*/

	return(log_prob);	


}

double BLCA_return_log_marginal_likelihood_component( struct component *component, struct mix_mod *mixmod )
{
	// these quantities will save some time for the Gibbs update of the labels
	int i,j,c, k, wis_k;
	double log_prob, s, r, a, *beta;
	
	if( mixmod->prior_type == 0 )
	{
	
		log_prob = lgamma(component->n_g+mixmod->alpha);
	
		if(component->n_g > 0){
	
			for(j=0;j<mixmod->d;j++){
		
				if(mixmod->varindicator[j]){
			
					log_prob += lgamma(mixmod->ncat[j] *mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + mixmod->ncat[j]*mixmod->beta);
			
					for(i=0;i<mixmod->ncat[j];i++){
			
						log_prob += lgamma(component->N[j][i] + mixmod->beta);
			
					}
			
				}
			
			}
		
		}	
		
	}else if( mixmod->prior_type == 1 ){
	
		k = mixmod->component_compute;
		
		log_prob = lgamma( component->n_g + mixmod->alpha_prior[k] ) /*- lgamma( mixmod->alpha_prior[k] )*/  ;
		
		if(component->n_g > 0){
	
			for(j=0;j<mixmod->d;j++){
		
				if(mixmod->varindicator[j]){
				
					s = 0.; r = 0.;
					beta = mixmod->beta_prior[k][j] ;
					for( c=0; c<mixmod->ncat[j]; c++ )
					{
						a = component->N[j][c] + beta[c] ;
						log_prob += lgamma( a ) - lgamma( beta[c] ) ;
						s += a ;
						r += beta[c] ; 
					}
					log_prob += lgamma(r) - lgamma(s) ;
					
				}	
			}
		
		}
		
	
	}
	
	return( log_prob );
	
}


void BLCA_recompute_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod)
{

	double log_prob = BLCA_return_log_marginal_likelihood_component( component, mixmod );
	
	component->log_prob = log_prob;

	return;
}


double BLCA_return_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod)
{

	int i,j;
	double log_prob, s=0., r=0., *beta;
	
	log_prob = 0.;
	
	if( mixmod->prior_type == 0 )
	{
	
		for(j=0;j<mixmod->d;j++){
	
			if(!mixmod->varindicator[j]){ /*variables not in mixture model*/
		
				log_prob += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(undiscriminating->n_g + mixmod->ncat[j]*mixmod->beta);
			
				for(i=0;i<mixmod->ncat[j];i++){
			
					log_prob += lgamma(undiscriminating->N[j][i] + mixmod->beta);
			
				}
		
			}
		
		}
	
	}else if( mixmod->prior_type == 1 ){
		
		for( j=0; j<mixmod->d; j++ ){
			
			if(!mixmod->varindicator[j]){
				s=0.; r=0.;
				beta = mixmod->beta_prior[0][j];
				
				for( i=0; i<mixmod->ncat[j]; i++ ){
					r += undiscriminating->N[j][i] + beta[i] ;
					log_prob += lgamma( undiscriminating->N[j][i] + beta[i] ) - lgamma( beta[i] );
					s +=  beta[i] ;
				}	
				log_prob += lgamma(s) - lgamma(r) ; 
			}
			
		}
		
	}
	
	return(log_prob);
	
}

void BLCA_recompute_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod)
{
	undiscriminating->log_prob = BLCA_return_marginal_likelihood_undiscriminating_variables( undiscriminating, mixmod );
	return;
}




