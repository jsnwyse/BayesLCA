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
			
	Last modification of this code: Thu 04 Aug 2016 21:02:24 IST */
	
#include "BLCA_analysis.h"

struct results *BLCA_analysis_MCMC_collapsed( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int fixedG, int onlyGibbs, 
			int selectVariables, int *group_memberships, int *variable_inclusion_indicator, int *n_groups, double *log_joint_posterior, 
			double *prior_include, int *var_pattern, int verbose, int verbose_update )
/*fixedG takes the value either TRUE or FALSE as defined in the macros*/
{

	int i, j, k, l, itmod, d_in, ej_case, vs_case, maxgroups = mixmod->maxgroups;
	
	int gap = (int)( num_iteration / thin_by) ;
	
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
	
	if( verbose ) Rprintf("\nInitializing sampler... starting burn-in...");
	
	for(l=0;l<num_iteration + num_burnin;l++){
	
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
			BLCA_update_model_by_variable_include_exclude_old(mixmod,&(results->accepted_include_exclude_variable),&(results->proposed_include_exclude_variable),j);
			if(mixmod->hprior_model){
				//update the prior probability variable inclusion using hyperprior
				j = 0;
				for(i=0;i<mixmod->d;i++) j += mixmod->varindicator[i];
				mixmod->prior_prob_variable_include = rbeta( (double)j + mixmod->hprior_model_a0, (double)(mixmod->d-j) + mixmod->hprior_model_b0);
			}
		}
		
		if( verbose && l == num_burnin-1 ) 
		{
				Rprintf("\nBurn-in completed...");
		}
		if( l == num_burnin-1 ) BLCA_reset_results( results );
		
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
				group_memberships[ itmod + j*gap ] = mixmod->z[j];
			}
			
			if( verbose && ( (l+1-num_burnin)/thin_by )%verbose_update == 0 )
				Rprintf("\n%d of %d samples completed....",(l+1-num_burnin)/thin_by,gap);
					
		}
		
		
		/*store the value of the log posterior here- full value (incl prior for no. components)*/
		log_joint_posterior[itmod] = BLCA_get_full_log_posterior(mixmod);

	}
	
	if( verbose ) Rprintf("\n");
	
	
	if(selectVariables){
		for(i=0;i<mixmod->d;i++){
			//results->variable_prob_inclusion[i] /= ((num_iteration-num_burnin)/thin_by);
		}
	}
	
	return(results);
	
}


void BLCA_analysis_MCMC_Gibbs_sampler( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int *group_memberships, 
						double *group_weights, double *variable_probabilities, double *log_joint_posterior, int sample_missing_data, int n_missing, 
						int *imputed_missing_values, int *position_missing, int verbose, int verbose_update )
{

	//this is the non-collapsed form of the model which uses Gibbs updates for all unknowns...
	// in this model there is no search for the number of groups, and no variable selection moves
	
	int i, j, k, l, c, p, idx, g_new, y_old, y_new, ind, itmod, *order;
	
	double *w,s,m;
	w = calloc(mixmod->G,sizeof(double));
	
	order = calloc(mixmod->n,sizeof(int));
	
	double **v; 
	v = calloc(mixmod->d,sizeof(double *));
	for(i=0;i<mixmod->d;i++) v[i] = calloc(mixmod->ncat[i],sizeof(double));
	
	int gap = (int)( num_iteration/thin_by ) , gap_;
	
	if( verbose ) Rprintf("\nInitializing sampler... starting burn-in.\n");
	
	for(l=0;l<num_iteration+num_burnin;l++){
	
		R_CheckUserInterrupt();
	  
		
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
				BLCA_add_to_component( mixmod->components[ mixmod->z[idx] ], mixmod->Yobs[idx], mixmod, -1 );
				BLCA_add_to_component( mixmod->components[ g_new ], mixmod->Yobs[idx], mixmod, 1 );
				
				mixmod->z[idx] = g_new;
			}
			
		}		
		
		//sample the weights
		s = 0.;
		for(i=0;i<mixmod->G;i++){
			w[i] = rgamma( (double) mixmod->components[i]->n_g + mixmod->alpha_prior[i] , 1.  ) ;
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
						 v[j][k] =  rgamma( (double)mixmod->components[i]->N[j][k] + mixmod->beta_prior[i][j][k] , 1. ) ;
						 s += v[j][k];
					}
					for(k=0;k<mixmod->ncat[j];k++) mixmod->components[i]->prob_variables[j][k] = v[j][k]/s;
				}
			}
		}
		
		p = 0;
		
		// impute the missing data before carrying out the update step
	  
	  if( sample_missing_data )
	  {
	    for( k=0; k<n_missing; k++ )
	    {
	       i = position_missing[ 2 * k ]; // this is where the problem is!!!
	       j = position_missing[ 2 * k + 1 ];
	      
	       g_new = mixmod->z[i];
	       
	       // take old vector out of component
	       BLCA_add_to_component( mixmod->components[ g_new ], mixmod->Yobs[ i ], mixmod, -1 );
	       
	       y_new = BLCA_sample_discrete( mixmod->components[g_new]->prob_variables[j], mixmod->ncat[j] );
          mixmod->Yobs[i][j] = y_new;
          mixmod->Y[j][i] = y_new;
         
         // add new vector to the component
         BLCA_add_to_component( mixmod->components[ g_new ], mixmod->Yobs[ i ], mixmod, 1 );
	    }
	  }
		
		//storage of results
		if(l > num_burnin-1 && (l+1-num_burnin)%thin_by == 0){
		
			if( verbose && l == num_burnin ) 
				Rprintf("\nBurn-in completed...");
				
		
			itmod = ((l+1-num_burnin)/thin_by)-1;
			
			
			for(j=0;j<mixmod->n;j++){
				group_memberships[itmod + j*gap ] = mixmod->z[j];
			}		
			
				
			for(j=0;j<mixmod->d;j++){
				gap_ = p * mixmod->G * gap;
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

			// if any missing data 
			
			if( sample_missing_data )
			{
			  for( k=0; k<n_missing; k++ )
			  {
			    i = position_missing[ 2 * k ];
			    j = position_missing[ 2 * k + 1 ];
			    
			    imputed_missing_values[ itmod*n_missing + k ] = mixmod->Yobs[i][j];
			  }
			}			
			
			log_joint_posterior[ itmod ] = BLCA_get_full_log_posterior_x2(mixmod);
			
			if( verbose && (l+1-num_burnin)%verbose_update == 0 )
				Rprintf("\n%d of %d samples completed....",l+1-num_burnin,num_iteration);
					
		}	
		
	}
	
	if(verbose) Rprintf("\n");

	free(w);
	for(i=0;i<mixmod->d;i++) free(v[i]);
	free(v);
	free(order);
	
	return;

}


void BLCA_analysis_EM( struct mix_mod *mixmod, int max_iterations, int *iterations, double *membership_probabilities, 
				double *group_weights, double *variable_probabilities, double *log_likelihood, int MAP, double tol, double *eps, int *converged ) 
{
	
	int i, j, g, k, c, p, iter = 0 ;
	double cond = DBL_MAX, llike_new, llike_old = -DBL_MAX ;
	
	//need to make a modification to this for when MAP == TRUE (not done yet)
	
	while( cond > tol  && iter < max_iterations )
	{
	
		R_CheckUserInterrupt();
		
		BLCA_M_step( mixmod );
		
		BLCA_E_step( mixmod );
		
		//llike_new = BLCA_get_log_likelihood( mixmod ) ;
		
		if( mixmod->EM_MAP ) mixmod->log_like += mixmod->log_prior;
		
		cond = fabs( mixmod->log_like - llike_old ) ;
		
		llike_old = mixmod->log_like;
		
		log_likelihood[ iter ] = mixmod->log_like;//llike_new ;
		
		iter++;
		
	}
	
	//Rprintf("\n Log post is %lf ", log_likelihood[iter-1]);
	
	if( cond < tol ) *converged = TRUE; else *converged = FALSE;
	
	//number of iterations to converge (if converged)
	*iterations = iter;
	*eps = cond;
	
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
					if( mixmod->EM_fit ) c = mixmod->Y[j][k] ;
					if( mixmod->BOOT ) c = mixmod->Y[j][ mixmod->boot_idx[k] ] ; 
					mixmod->s[k][g] += log( mixmod->components[g]->prob_variables[j][c] ) ;
				}
			
			}	
		}
		
		//BLCA_log_sum_exp also modifies the s[k] arg and turns to weights
		mixmod->log_like += BLCA_get_log_sum_exp( mixmod->s[k] , mixmod->G ) ;
	}

	//add prior terms if EM_MAP
	/*if( mixmod->EM_MAP ) 
	{
		mixmod->log_like += mixmod->log_prior ; 
	}*/

	return;
}


void BLCA_M_step( struct mix_mod *mixmod )
{
	int k, g, j, c;
	double a, s = 0., r=0., eps=1E-10;
	
	// eps is a small value to prevent numerical instability for small probabilities
	
	struct component *comp;
	
	//separate cases if the mixmod->EM_MAP is TRUE
	
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
	//update the prior term as going
	
	if( mixmod->EM_MAP ) mixmod->log_prior = 0.;
	
	for( g=0; g<mixmod->G; g++ )
	{
	
		comp = mixmod->components[g];
	
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
			
				for( c=0; c<mixmod->ncat[j]; c++ ) comp->prob_variables[j][c] = 0.;
				
				for( k=0; k<mixmod->n; k++ )
				{
					if( mixmod->EM_fit ) c = mixmod->Y[j][k] ; 
					if( mixmod->BOOT ) c = mixmod->Y[j][ mixmod->boot_idx[k] ] ;
					comp->prob_variables[j][c] += mixmod->s[k][g] ;
				}
				
				s = 0.;
				
				for( c=0; c<mixmod->ncat[j]; c++ )
				{
					if( mixmod->EM_MAP )
					{
						comp->prob_variables[j][c] += mixmod->beta_prior[g][j][c] - 1. ;
						s += mixmod->beta_prior[g][j][c] ; 
					}
					else
					{
						comp->prob_variables[j][c] /= mixmod->weights[g] ;
					}
				}
				
				if( mixmod->EM_MAP )
				{
					for( c=0;  c<mixmod->ncat[j]; c++ ) 
					{
						comp->prob_variables[j][c] /= ( mixmod->weights[g] + s - mixmod->ncat[j] );
						a = comp->prob_variables[j][c] < eps ? eps : comp->prob_variables[j][c] ;
						mixmod->log_prior += ( mixmod->beta_prior[g][j][c] - 1. ) * log( a ) ;
						mixmod->log_prior -= lgamma( mixmod->beta_prior[g][j][c] ) ;  
					}
					mixmod->log_prior += lgamma( s ) ;
				}
				
			}
		}
		
		if( mixmod->EM_MAP )
		{
			r += mixmod->alpha_prior[g] ; 
			mixmod->weights[g] +=  mixmod->alpha_prior[g] - 1. ;
		}
		else
		{
			mixmod->weights[g] /= mixmod->n ;
		}
		
	}
	
	if( mixmod->EM_MAP )
	{
		for( g=0; g<mixmod->G; g++ ) 
		{
			mixmod->weights[g] /= ( mixmod->n + r - mixmod->G );
			a = mixmod->weights[g] < eps ? eps : mixmod->weights[g] ; 
			mixmod->log_prior += ( mixmod->alpha_prior[g] - 1.) * log( a ) ; 
			mixmod->log_prior -= lgamma( mixmod->alpha_prior[g] ); 
		}
		
		 mixmod->log_prior += lgamma( r ) ;
	}
	
	
	return;
}

/*-------------------------------------- VB algorithm -----------------------------------*/

void BLCA_analysis_VB( struct mix_mod *mixmod, int max_iterations, int *iterations, double *group_probabilities, double *group_weights, 
				double *sd_group_weights, double *vb_pars_group_weights, double *prob_variables, double *sd_prob_variables, double *vb_pars_prob_variables, double *lower_bound, double tol, int *converged, double *log_post )
{
	int i, j, g, k, c, p, iter = 0 ;
	double cond = DBL_MAX, bound_new, bound_old = -DBL_MAX ;
	
	//need to make a modification to this for when MAP == TRUE (not done yet)
	
	while( cond > tol  && iter < max_iterations )
	{
	
		R_CheckUserInterrupt();
		
		BLCA_VB_phi_step( mixmod );
		
		BLCA_VB_alpha_beta_step( mixmod );
		
		bound_new = BLCA_get_VB_bound( mixmod ) ;
		
		cond = fabs( bound_new - bound_old ) ;
		
		bound_old = bound_new;
		
		lower_bound[ iter ] = bound_new ;
		
		iter++;
		
	}
	
	if( cond < tol ) *converged = TRUE; else *converged = FALSE;
	
	//number of iterations to converge (if converged)
	*iterations = iter;
	
	//store the results before returning
	
	double sum_alpha_ud = 0., var;
	for( g=0; g<mixmod->G ; g++ ) sum_alpha_ud += mixmod->alpha_ud[g] ;
	
	for( g=0; g<mixmod->G ; g++ )  
	{
		group_weights[g] = mixmod->alpha_ud[g] / sum_alpha_ud ;
		vb_pars_group_weights[g] = mixmod->alpha_ud[g] ;
		mixmod->weights[g] = group_weights[g] ;
		var = mixmod->alpha_ud[g] * ( sum_alpha_ud - mixmod->alpha_ud[g] ) / ( sum_alpha_ud * sum_alpha_ud * ( sum_alpha_ud + 1. ) );
		sd_group_weights[g] = sqrt( var );
	}
	
	for( k=0; k<mixmod->n; k++ ) 
	{
		for( g=0; g<mixmod->G; g++ )
			group_probabilities[ (mixmod->n)*g + k ] = mixmod->s[k][g] ;
	}
	
	p = 0;
	int gap_;
	
	//convert the beta_ud into expect values of the theta's
	
	struct component *comp;
	double sum_b = 0., ***var_prob_variables ;
	
	var_prob_variables = calloc( mixmod->G, sizeof(double**) );
	for( g=0; g<mixmod->G; g++ )
	{
		var_prob_variables[g] = calloc( mixmod->d, sizeof(double *) );
		for( j=0; j<mixmod->d; j++ ) 
			var_prob_variables[g][j] = calloc( mixmod->ncat[j], sizeof(double) );
	}
	
	for( g=0; g<mixmod->G; g++ )
	{
		comp = mixmod->components[g] ;
		
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				sum_b = 0.;
				for( c=0; c<mixmod->ncat[j]; c++ ) sum_b += comp->beta_ud[j][c] ;
				for( c=0; c<mixmod->ncat[j]; c++ ) 
				{
					comp->prob_variables[j][c] = comp->beta_ud[j][c] / sum_b;
					var_prob_variables[g][j][c] = comp->beta_ud[j][c] * ( sum_b - comp->beta_ud[j][c] ) / ( sum_b * sum_b * ( sum_b + 1. )  ) ;
				}
			}
		}
	
	}
	
	for( j=0; j<mixmod->d; j++ )
	{
		gap_ = p * mixmod->G;
		for( g=0; g<mixmod->G; g++ )
		{
			for( c=0; c<mixmod->ncat[j]; c++ )
			{
				prob_variables[ gap_ + g * mixmod->ncat[j] + c ] = mixmod->components[g]->prob_variables[j][c] ;
				sd_prob_variables[ gap_ + g * mixmod->ncat[j] + c ] = sqrt( var_prob_variables[g][j][c] ) ;
				//vb pars
				vb_pars_prob_variables[ gap_ + g * mixmod->ncat[j] + c ] = comp->beta_ud[j][c] ;
			}
		} 
		p += mixmod->ncat[j] ;
	}
	
	for( g=0; g<mixmod->G; g++ )
	{
		for( j=0; j<mixmod->d; j++ ) free( var_prob_variables[g][j] );
		free( var_prob_variables[g] );
	}
	free( var_prob_variables );
	
	//compute the log of the posterior... 
	
	*log_post = BLCA_get_log_likelihood( mixmod );
	
	return;
}

void BLCA_VB_phi_step( struct mix_mod *mixmod )
{
	int k, g, j, c;
	double lse;
	struct component *comp;

	for( k=0; k<mixmod->n; k++ ){
	
		for( g=0; g<mixmod->G; g++ ){
			
			comp = mixmod->components[g];
			
			mixmod->s[k][g] = mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ;
			
			for( j=0; j<mixmod->d; j++ ){
				
				if( mixmod->varindicator[j] )
				{
					c = mixmod->Y[j][k] ;
					mixmod->s[k][g] += comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ;
				}
			
			}	
		}
		
		lse = BLCA_get_log_sum_exp( mixmod->s[k], mixmod->G );
	}

	return;	
}

void BLCA_VB_alpha_beta_step( struct mix_mod *mixmod )
{
	
	//update the alphas first
	int k, g, j, c;
	double *colsums = calloc( mixmod->G, sizeof(double) );
	struct component *comp;

	for( g=0; g<mixmod->G; g++ )
	{
		for( k=0; k<mixmod->n; k++ )
		{
			colsums[g] += mixmod->s[k][g] ;
		}
	}
	
	//alpha update 
	mixmod->di_sum_alpha_ud = 0.;
	for( g=0; g<mixmod->G; g++ )
	{
		mixmod->alpha_ud[g] = colsums[g] + mixmod->alpha_prior[g]; 
		mixmod->di_alpha_ud[g] = digammaRN( mixmod->alpha_ud[g] );

		mixmod->di_sum_alpha_ud += mixmod->alpha_ud[g];
		
		//initialize the beta update
		comp = mixmod->components[g];
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				for( c=0; c<mixmod->ncat[j]; c++ )
					comp->beta_ud[j][c] = mixmod->beta_prior[g][j][c] ;
			}
		}
		
	}
	mixmod->di_sum_alpha_ud = digammaRN( mixmod->di_sum_alpha_ud );
	
	//beta update 

	for( g=0; g<mixmod->G; g++ )
	{
		comp = mixmod->components[g];
		for( k=0; k<mixmod->n; k++ )
		{
			for( j=0; j<mixmod->d; j++ )
			{
				if( mixmod->varindicator[j] )
				{
					c = mixmod->Y[j][k] ;
					comp->beta_ud[j][c] += mixmod->s[k][g] ;
				}
			}
		}
	}
	
	for( g=0; g<mixmod->G; g++ )
	{
		comp = mixmod->components[g];	
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				comp->di_sum_beta_ud[j] = 0.;
				for( c=0; c<mixmod->ncat[j]; c++ )
				{
					comp->di_beta_ud[j][c] = digammaRN( comp->beta_ud[j][c] ) ;
					comp->di_sum_beta_ud[j] += comp->beta_ud[j][c] ;
				}
				comp->di_sum_beta_ud[j] = digammaRN( comp->di_sum_beta_ud[j] ) ;
			}
		}
	}

	free( colsums );
	return;
}


void BLCA_analysis_Boot( struct mix_mod *mixmod, int boot_samples, int max_iterations, 
									double *group_weights, double *prob_variables, double *log_likelihood, double tol, int verbose, int verbose_update )
{

	int b, iter, i, j, k, l, c, p, idx, g_new, ind, gap = boot_samples, gap_, *boot_idx;
	double cond = DBL_MAX, llike_new, llike_old = -DBL_MAX ;
	
	//store the probabilities & weights from original EM run
	struct mix_mod *mixmod_0 = BLCA_clone_mixmod( mixmod );
	
	for( b=0; b<boot_samples; b++ )
	{
		// generate indexes for a bootstrap sample
		for( k=0; k<mixmod->n; k++ ) 
		{
			mixmod->boot_idx[k] = (int) ( runif(0.0,1.0) * mixmod->n ) ;
			//Rprintf("\n boot sample %d is %d: %d", b, k, mixmod->boot_idx[k] ) ;	
		}
		
		
		// reset the component weights and item probabilities from original fit
		for( k=0; k<mixmod->G; k++ ) 
		{
			BLCA_copy_component( mixmod_0->components[k], mixmod->components[k], mixmod );
			mixmod->weights[k] = mixmod_0->weights[k];
		}
		
		cond = DBL_MAX; llike_old = -DBL_MAX; iter = 0;
	
		while( cond > tol  && iter < max_iterations )
		{
	
			R_CheckUserInterrupt();
			
			BLCA_E_step( mixmod );
		
			BLCA_M_step( mixmod );
		
			if( mixmod->EM_MAP ) mixmod->log_like += mixmod->log_prior;
		
			cond = fabs( mixmod->log_like - llike_old ) ;
		
			llike_old = mixmod->log_like;
		
			//log_likelihood[ iter ] = mixmod->log_like;
		
			iter++;
		
		}
		
		log_likelihood[b] = mixmod->log_like;
	
		p = 0;
		
		//storage of results	
			
		for(j=0;j<mixmod->d;j++){
			gap_ = p * mixmod->G * gap ;
			for(i=0;i<mixmod->G;i++){
				for(k=0;k<mixmod->ncat[j];k++){
					prob_variables[ //long index expression
							
						gap_ + b * mixmod->ncat[j]*mixmod->G + i*mixmod->ncat[j] + k 
								
						] = mixmod->components[i]->prob_variables[j][k];
				}
			}
				p += mixmod->ncat[j];
		}
				
		for(i=0;i<mixmod->G;i++) group_weights[ b * mixmod->G + i ] = mixmod->weights[i];
			
		//log_joint_posterior[ itmod ] = BLCA_get_full_log_posterior_x2(mixmod);
			
		if( verbose && (b+1)%verbose_update == 0 )
			Rprintf("\n%d of %d bootstrap samples completed....",b+1,boot_samples);
					
	
	}
	
	BLCA_free_mixmod( mixmod_0 );
	
	return;
	
}


