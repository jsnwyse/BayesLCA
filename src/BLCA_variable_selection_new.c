#include "BLCA_variable_selection_new.h"

void BLCA_update_model_by_variable_include_exclude( struct mix_mod *mixmod, int *accepted, int *proposed, int k )
{
	
	int i, j, g, S=0, D=0, Sn, Dn, sgn, *wis;
	double logP=0., logPn, cstr;
	double *lprobs, lp_undisc, lr, lp_curr, lp_prop, lr_prior;
	struct component *comp;
	
	for( j=0; j<mixmod->d; j++ )
	{
		if( mixmod->varindicator[j] ) 
		{
			S += mixmod->ncat[j];
			logP += log( mixmod->ncat[j] );
			D += 1;
		}
	}
	
	if( mixmod->varindicator[k] ) sgn = -1; else sgn = 1;
	
	Sn = S + sgn * mixmod->ncat[k] ;
	logPn = logP + sgn * log( mixmod->ncat[k] ) ;
	Dn = D + sgn ; 
	
	cstr = logPn - log( ( Sn - Dn + 1 ) * mixmod->G ) ;
	
	// if constraint not satisfied return
	if( cstr < 0 ) return; 
	
	mixmod->varindicator[k] = 1 - mixmod->varindicator[k] ; 
	
	lr_prior = BLCA_l_prior_variable_include( D, mixmod )
	 				- BLCA_l_prior_variable_include( D + sgn, mixmod ); 
	
	
	*proposed += 1;
	
	lprobs = calloc( mixmod->G, sizeof(double) );
	
	wis = mixmod->whereis;
	
	lp_curr = 0.; lp_prop = 0.;
	
	for( g=0; g<mixmod->G; g++ )
	{
		mixmod->component_compute = g;
		comp = mixmod->components[ wis[g] ] ; 
		lprobs[g] =  BLCA_return_log_marginal_likelihood_component( comp, mixmod );
		lp_prop += lprobs[g];
		lp_curr += comp->log_prob;
	}
	
	mixmod->component_compute = 0;
	lp_undisc = BLCA_return_marginal_likelihood_undiscriminating_variables( mixmod->undiscriminating, mixmod ) ; 
	
	lr = lp_prop - lp_curr + lp_undisc - mixmod->undiscriminating->log_prob + lr_prior ;
	
	Rprintf("\n lr is %lf ", lr );
	
	if( log(runif(0.,1.)) < lr )
	{
		*accepted += 1;
		for( g=0; g<mixmod->G; g++ )
		{
			mixmod->components[ wis[g] ]->log_prob = lprobs[g];
		}
		mixmod->undiscriminating->log_prob = lp_undisc;
	}else{
		mixmod->varindicator[k] = 1 - mixmod->varindicator[k] ; 
	}
	
	free( lprobs );
	
	return;

}
