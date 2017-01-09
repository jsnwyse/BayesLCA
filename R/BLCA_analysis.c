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
			
	Last modification of this code: Thu 04 Aug 2016 21:02:24 IST    */
	
#include "BLCA_analysis.h"

struct results *BLCA_analysis_MCMC_collapsed( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int fixedG, int onlyGibbs, int selectVariables, int *group_memberships, int *variable_inclusion_indicator, int *n_groups, double *log_joint_posterior, double *prior_include, int *var_pattern )
/*fixedG takes the value either TRUE or FALSE as defined in the macros*/
{

	int i, j, k, l, itmod, d_in, ej_case, vs_case, maxgroups = mixmod->maxgroups;
	
	int gap = (num_iteration - num_burnin)/(thin_by);
	
	//writeToFile = FALSE;
	
	//struct results res,*results;
	//results = &res;

	double pr_ej_G,pr_ej_Gm1,pr_ej_Gp1;
	
	struct results *results = (struct results *)malloc(sizeof(struct results));
	BLCA_allocate_results( results, num_iteration, num_burnin, thin_by, mixmod->n, mixmod->d  );
	
	/*initialize with all variables in model*/
	//pass in the pattern from R
	for(j=0;j<mixmod->d;j++){
	   mixmod->varindicator[j] = var_pattern[j];
	}
	
	for(l=0;l<num_iteration;l++){
	
		R_CheckUserInterrupt();
		
		/*check for violation of identifiability constraint*/
		int s = 0;
		double log_p = 0;
		int in = 0;
		
		for(j=0;j<mixmod->d;j++){
		  
		  if(mixmod->varindicator[j]){
		    
		    s += mixmod->ncat[j];
		    log_p += log(mixmod->ncat[j]);
		    in += 1;
		    
		  }
		}
		
		if(!(log_p > log((s - in + 1)*mixmod->G)) && !(l>0)){
		  //Rprintf("\n***Warning***: identifiability constraint not satisfied; log_p = %d, s=%d, in=%d, G=%d",log_p,s,in,mixmod->G);
		}
	
		BLCA_update_allocations_with_gibbs(mixmod);
		
		if(!onlyGibbs){
			BLCA_update_allocations_with_metropolis_move_1(mixmod,&(results->accepted_m1),&(results->proposed_m1));
			BLCA_update_allocations_with_metropolis_move_2(mixmod,&(results->accepted_m2),&(results->proposed_m2));
			BLCA_update_allocations_with_metropolis_move_3(mixmod,&(results->accepted_m3),&(results->proposed_m3));
		}
	
		if(!fixedG){

			if(mixmod->G == 1){
				ej_case = 0;
			}else if(mixmod->G == maxgroups){
				ej_case = 1;
			}else if(mixmod->G == 2){
				ej_case = 2;
			}else if(mixmod->G == maxgroups-1){
				ej_case = 3;
			}else{
				ej_case = 4;
			}
		
		
			switch(ej_case)
			{
				case 0:
					pr_ej_G = 1.;
					pr_ej_Gp1 = .5;
					pr_ej_Gm1 = 0.;
				break;
			
				case 1:
					pr_ej_G = 0.;
					pr_ej_Gp1 = 0.;
					pr_ej_Gm1 = .5;
				break;
			
				case 2:
					pr_ej_G = .5;
					pr_ej_Gp1 = .5;
					pr_ej_Gm1 = 1.;
				break;
			
				case 3:
					pr_ej_G = 0.5;
					pr_ej_Gp1 = 0.;
					pr_ej_Gm1 = 0.5;
				break;
			
				case 4:
					pr_ej_G = 0.5;
					pr_ej_Gp1 = 0.5;
					pr_ej_Gm1 = 0.5;
				break;
			}
		
		
			if( runif(0.0,1.0) < pr_ej_G){		
				BLCA_update_allocations_with_ejection_move(mixmod,&(results->accepted_eject),&(results->proposed_eject),pr_ej_G,pr_ej_Gp1);
			}else{
				BLCA_update_allocations_with_absorb_move(mixmod,&(results->accepted_absorb),&(results->proposed_absorb),pr_ej_G,pr_ej_Gm1);
			}
		
		
		}
		
		if(selectVariables){
			j = BLCA_random_integer(mixmod->d);
			BLCA_update_model_by_variable_include_exclude(mixmod,&(results->accepted_include_exclude_variable),&(results->proposed_include_exclude_variable),j);
			if(mixmod->hprior_model){
				//update the prior probaility variable inclusion using hyperprior
				j = 0;
				for(i=0;i<mixmod->d;i++) j += mixmod->varindicator[i];
				mixmod->prior_prob_variable_include = rbeta( (double)j + mixmod->hprior_model_a0, (double)(mixmod->d-j) + mixmod->hprior_model_b0);
			}
		}
		
		
		if(l > num_burnin-1 && (l+1-num_burnin)%thin_by == 0){
		
			itmod = ((l+1-num_burnin)/thin_by)-1;
		
			if(selectVariables){
			
				for(j=0;j<mixmod->d;j++){
					variable_inclusion_indicator[ itmod + j*gap ] = mixmod->varindicator[j];//results->variable_indicator[i][j];
				}
			
				for(i=0;i<mixmod->d;i++){
				 	if(mixmod->varindicator[i]) results->variable_prob_inclusion[i] += 1.;
				}
			}
			
			results->ngroups[itmod] = mixmod->G;
			
			n_groups[itmod] = mixmod->G;
			
			prior_include[itmod] = mixmod->prior_prob_variable_include;	
				
			for(j=0;j<mixmod->n;j++){
				group_memberships[itmod + j*gap ] = mixmod->z[j];
			}
					
		}
		
		
		/*store the value of the log posterior here- full value (incl prior for no. components)*/
		log_joint_posterior[l] = BLCA_get_full_log_posterior(mixmod);

	}
	
	
	if(selectVariables){
		for(i=0;i<mixmod->d;i++){
			results->variable_prob_inclusion[i] /= ((num_iteration-num_burnin)/thin_by);
		}
	}
	
	return(results);
	
}


void BLCA_analysis_MCMC_Gibbs_sampler( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int *group_memberships, double *group_weights, double *variable_probabilities, double *log_joint_posterior, int verbose, int verbose_update )
{

	//this is the non-collapsed form of the model which uses Gibbs updates for all unknowns...
	// in this model there is no search for the number of groups, and no variable selection moves
	
	int i, j, k, l, c, p, idx, g_new, ind, itmod, *order;
	
	double *w,s,m;
	w = calloc(mixmod->G,sizeof(double));
	
	order = calloc(mixmod->n,sizeof(int));
	
	double **v; 
	v = calloc(mixmod->d,sizeof(double *));
	for(i=0;i<mixmod->d;i++) v[i] = calloc(mixmod->ncat[i],sizeof(double));
	
	int gap = num_iteration/thin_by, gap_;
	
	if( verbose ) Rprintf("\nInitializing sampler... starting burn-in.\n");
	
	for(l=0;l<num_iteration+num_burnin;l++){
	
		R_CheckUserInterrupt();
	
		//sample the weights
		s = 0.;
		for(i=0;i<mixmod->G;i++){
			w[i] = rgamma( (double) mixmod->components[i]->n_g + mixmod->alpha , 1.  ) ;
			s += w[i];
		}
		for(i=0;i<mixmod->G;i++) mixmod->weights[i] = w[i]/s;
		
		//sample the variable probabilities
		for(i=0;i<mixmod->G;i++){
			for(j=0;j<mixmod->d;j++){
			
				if( mixmod->varindicator[j] )
				{
					//only sample if the variable is included in the model
					s = 0.;
					for(k=0;k<mixmod->ncat[j];k++){
						 v[j][k] =  rgamma( (double)mixmod->components[i]->N[j][k] + mixmod->beta , 1. ) ;
						 s += v[j][k];
					}
					for(k=0;k<mixmod->ncat[j];k++) mixmod->components[i]->prob_variables[j][k] = v[j][k]/s;
				}
			}
		}
	
		//sample the memberships in a random order
		
		for(i=0;i<mixmod->n;i++){
			order[i] = i;
		}
		
		BLCA_random_ranshuffle( order, mixmod->n );
		
		for(k=0;k<mixmod->n;k++){
			
			idx = order[k];
			
			for(i=0;i<mixmod->G;i++){
				w[i] = log(mixmod->weights[i]);
				for(j=0;j<mixmod->d;j++){
				
					if( mixmod->varindicator[j] ){
						c = mixmod->Y[j][idx];
						w[i] += log(mixmod->components[i]->prob_variables[j][c]);	
					}	
				}
			}
			
			m = w[0];
			for(i=1;i<mixmod->G;i++){
				if(w[i]>m) m = w[i];
			}
			
			s = 0.;
			for(i=0;i<mixmod->G;i++){
				 w[i] -= m;
				 w[i] = exp(w[i]);
				 s += w[i];
			}
			for(i=0;i<mixmod->G;i++) w[i] /= s;
			
			
			
			g_new = BLCA_sample_discrete( w, mixmod->G );
			
			
			if(g_new != mixmod->z[idx])
			{
				//take out of component counts in old and put into new
				for(j=0;j<mixmod->d;j++)
				{
					if( mixmod->varindicator[j] )
					{
						c = mixmod->Y[j][idx];
						mixmod->components[ mixmod->z[idx] ]->N[j][c] -= 1;
						mixmod->components[ g_new ]->N[j][c] += 1;
					}
				}
				//reduce the count of the old component
				mixmod->components[ mixmod->z[idx] ]->n_g -= 1;
				mixmod->components[ g_new ]->n_g += 1;
				mixmod->z[idx] = g_new;
			}
			
		}
		
		p = 0;
		
		//storage of results
		if(l > num_burnin-1 && (l+1-num_burnin)%thin_by == 0){
		
			if( verbose && l == num_burnin-1 ) 
				Rprintf("\nBurn-in completed...\n");
				
		
			itmod = ((l+1-num_burnin)/thin_by)-1;
			
			
			for(j=0;j<mixmod->n;j++){
				group_memberships[itmod + j*gap ] = mixmod->z[j];
			}		
			
			//there is potentially a bug here-- check this out...
				
			for(j=0;j<mixmod->d;j++){
				gap_ = p*mixmod->G*gap;
				for(i=0;i<mixmod->G;i++){
					for(k=0;k<mixmod->ncat[j];k++){
						variable_probabilities[ //long index expression
							
							gap_ + itmod*mixmod->ncat[j]*mixmod->G + i*mixmod->ncat[j] + k 
								
							] = mixmod->components[i]->prob_variables[j][k];
					}
				}
				p += mixmod->ncat[j];
			}
				
			for(i=0;i<mixmod->G;i++) group_weights[itmod*mixmod->G + i] = mixmod->weights[i];
			
			log_joint_posterior[ itmod ] = BLCA_get_full_log_posterior_x2(mixmod);
			
			if( verbose && (l+1-num_burnin)%verbose_update == 0 )
				Rprintf("\n%d of %d samples completed....\n",l+1-num_burnin,num_iteration);
					
		}	
		
	}

	free(w);
	for(i=0;i<mixmod->d;i++) free(v[i]);
	free(v);
	free(order);
	
	return;

}


void BLCA_analysis_EM( struct mix_mod *mixmod, int max_iterations, int iterations, double *membership_probabilities, double *group_weights, double *variable_probabilities, double *log_likelihood, int MAP, double tol ) 
{
	
	int i, j, g, k, c, p, iter = 0 ;
	double cond = DBL_MAX, llike_new, llike_old = -DBL_MAX ;
	
	while( cond > tol  && iter < max_iterations )
	{
		
		BLCA_E_step( mixmod );
		
		BLCA_M_step( mixmod );
		
		llike_new = BLCA_get_log_likelihood( mixmod ) ;
		
		cond = fabs( mixmod->log_like - llike_old ) ;
		
		llike_old = mixmod->log_like;
		
		log_likelihood[ iter ] = llike_new ;
		
		iter++;
		
	}
	
	//store the results before returning
	
	for( g=0; g<mixmod->G ; g++ ) group_weights[g] = mixmod->weights[g] ;
	
	for( k=0; k<mixmod->n; k++ ) 
	{
		for( g=0; g<mixmod->G; g++ )
			membership_probabilities[ (mixmod->n)*g + k ] = mixmod->s[k][g] ;
	}
	
	p = 0;
	int gap_;
	
	for( j=0; j<mixmod->d; j++ )
	{
		gap_ = p * mixmod->G;
		
		for( g=0; g<mixmod->G; g++ )
		{
			for( c=0; c<mixmod->ncat[j]; c++ )
			{
				variable_probabilities[ gap_ + g * mixmod->ncat[j] + c ] = mixmod->components[g]->prob_variables[j][c] ;
			}
		} 
		
		p += mixmod->ncat[j] ;
	
	}
	
	return;
}



void BLCA_E_step( struct mix_mod *mixmod )
{
	int k, g, j, c;
	
	mixmod->log_like = 0.;
	
	for( k=0; k<mixmod->n; k++ ){
	
		for( g=0; g<mixmod->G; g++ ){
			
			mixmod->s[k][g] = log( mixmod->weights[g] );
			
			for( j=0; j<mixmod->d; j++ ){
				
				if( mixmod->varindicator[j] )
				{
					c = mixmod->Y[j][k] ;
					mixmod->s[k][g] += log( mixmod->components[g]->prob_variables[j][c] ) ;
				}
			
			}	
		}
		
		//BLCA_log_sum_exp also modifies the s[k] arg and turns to weights
		mixmod->log_like += BLCA_get_log_sum_exp( mixmod->s[k] , mixmod->G ) ;
	}

	return;
}


void BLCA_M_step( struct mix_mod *mixmod )
{
	int k, g, j, c;
	
	// get the sum of the probabilities
	for( g=0; g<mixmod->G; g++ )
	{
		mixmod->weights[g] = 0.;
		for( k=0; k<mixmod->n; k++ )
		{
			mixmod->weights[g] += mixmod->s[k][g] ;
		}
	}
	
	// update the variable probabilities
	for( g=0; g<mixmod->G; g++ )
	{
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				for( c=0; c<mixmod->ncat[j]; c++ )
					mixmod->components[g]->prob_variables[j][c] = 0.;
				
				for( k=0; k<mixmod->n; k++ )
				{
					c = mixmod->Y[j][k] ;
					mixmod->components[g]->prob_variables[j][c] += mixmod->s[k][g] ;
				}
				
				for( c=0; c<mixmod->ncat[j]; c++ )
					mixmod->components[g]->prob_variables[j][c] /= mixmod->weights[g] ;
			}
		}
		
		mixmod->weights[g] /= mixmod->n ;
	}
	
	return;
}




