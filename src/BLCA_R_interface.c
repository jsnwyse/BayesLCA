/*Functions written in C for the fitting of Latent Class Analysis models using MCMC
	methods. Two implementations of the model are included in the Bayesian
	formulation: collapsed and not collapsed.
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie

	Last modification of this code: Mon 18 Feb 2019 11:56:42 IST 			
*/

#include "BLCA_mixmod.h"

static void BLCA_VS(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *fixed_groups, int *just_gibbs_updates, int *init_n_groups, int *max_n_groups, int *n_iterations, int *n_burn_in, int *thin_by, int *n_gibbs, int *group_memberships, int *n_groups, double *prior_G, int *variable_select, int *variable_inclusion_indicator, double *prior_prob_include, double *log_joint_posterior, int *hprior_model, double *prior_include, int *var_pattern, int *verbose, int *verbose_update, double *acc_rts );

static void BLCA_VS_COMPUTE_POST_HOC_PARAMETER_ESTIMATES(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *n_groups, int *n_sample, int *memberships, int *variable_inclusion_indicator, int *variable, double *estimate, double *sd_estimate, double *classprob_estimate, double *sd_classprob_estimate ) ; 

static void BLCA_GIBBS_SAMPLER(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *init_type, int *n_groups, int *n_iterations, int *n_burn_in, int *thin_by, int *group_memberships, double *group_weights, double *prob_variables, int *var_in, double *log_joint_posterior, int *sample_missing_data, int *n_missing, int *imputed_missing_values, int *position_missing, int *verbose, int *verbose_update );

static void BLCA_EM_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *n_groups, int *init_type, double *init_vals, int *max_iterations, int *iterations, double *group_probabilities, double *group_weights, double *prob_variables, int *var_in, double *log_likelihood, int *MAP, double *tol, double *eps, int *converged, double *max_log_likelihood );

static void BLCA_VB_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *n_groups, int *init_type, double *init_vals, int *max_iterations, int *iterations, double *group_probabilities, double *group_weights, double *sd_group_weights, double *par_group_weights, double *prob_variables, double *sd_prob_variables, double *par_prob_variables, int *var_in, double *lower_bound, double *tol, int *converged, double *log_post );

static void BLCA_BOOT_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *boot_samples, int *n_groups, double *init_group_weights, double *init_group_probabilities, double *out_group_weights, double *out_prob_variables, double *out_log_likelihood, int *var_in, int *MAP, double *tol, int *max_iterations, int *verbose, int *verbose_update );

static void BLCA_RELABEL( int *n_obs, int *n_sample, int *n_groups, int *labels_in, int *labels_out, int *permutation );

static const R_CMethodDef cMethods[] = {
	{ "BLCA_VS", (DL_FUNC) &BLCA_VS , 29 },
	{ "BLCA_VS_COMPUTE_POST_HOC_PARAMETER_ESTIMATES", (DL_FUNC) &BLCA_VS_COMPUTE_POST_HOC_PARAMETER_ESTIMATES , 14 },
	{ "BLCA_GIBBS_SAMPLER", (DL_FUNC) &BLCA_GIBBS_SAMPLER , 24 },
	{ "BLCA_EM_FIT", (DL_FUNC) &BLCA_EM_FIT , 23 },
	{ "BLCA_VB_FIT", (DL_FUNC) &BLCA_VB_FIT , 25 },
	{ "BLCA_BOOT_FIT", (DL_FUNC) &BLCA_BOOT_FIT , 21 },
	{ "BLCA_RELABEL", (DL_FUNC) &BLCA_RELABEL, 6 },
	NULL
};

void R_init_BayesLCA( DllInfo *info )
{
	R_registerRoutines( info, cMethods, NULL, NULL, NULL );
	R_useDynamicSymbols( info, FALSE );
}



static void BLCA_VS(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *fixed_groups, int *just_gibbs_updates, int *init_n_groups, int *max_n_groups, int *n_iterations, int *n_burn_in, int *thin_by, int *n_gibbs, int *group_memberships, int *n_groups, double *prior_G, int *variable_select, int *variable_inclusion_indicator, double *prior_prob_include, double *log_joint_posterior, int *hprior_model, double *prior_include, int *var_pattern, int *verbose, int *verbose_update, double *acc_rts )
{

	// I have just put the vectors for the priors in here... these have to be passed to the initialize function

	int i,j,n,d,inG,mxG,nit,nburn,fixed,justgibbs;
	struct mix_mod *mixmod;
	struct results *results;
	
	mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *max_n_groups, *init_n_groups, hparam, ncat, TRUE, FALSE, FALSE, FALSE, FALSE );
	
	//change prior on G into a double vector so that it can be passed directly from R
	
	for( i=1; i<*max_n_groups+1; i++ ) mixmod->log_prior_G[i] = log( prior_G[i-1] );
	
	/*if(*prior_G == 0){
		BLCA_set_prior_on_number_of_components(mixmod,RICHARDSON_AND_GREEN);
	}else{
		BLCA_set_prior_on_number_of_components(mixmod,NOBILE_AND_FEARNSIDE);
	}*/
	
	BLCA_initialize_data( mixmod, Y );
	
	BLCA_initialize_priors( mixmod, alpha_prior, beta_prior, *prior_init_type );
	
	mixmod->n_gibbs = (*n_gibbs);
	
	//set prior for variable inclusion-- say 50% variables relevant for explaining clustering
	
	mixmod->prior_prob_variable_include = (*prior_prob_include);
	
	//should there be a hyperprior put on this value?
	
	mixmod->hprior_model = (*hprior_model);
	
	if(mixmod->hprior_model)
	{
		mixmod->hprior_model_a0 = 1.;
		mixmod->hprior_model_b0 = (double)mixmod->d/4.; //use d/4 as a default value- this favours more parsimony in clustering variables
	}
	
	//initialize all variables to be in model
	
	for( j=0; j<mixmod->d; j++ )
		mixmod->varindicator[j] = 1;

	GetRNGstate();

	//set first arg to -1 as collapsed sampler
	BLCA_initialize_Gibbs_sampler( -1, mixmod );
	
	//this call does the MCMC
	
	results = BLCA_analysis_MCMC_collapsed(mixmod,*n_iterations,*n_burn_in,*thin_by,*fixed_groups,*just_gibbs_updates,*variable_select,
												group_memberships, variable_inclusion_indicator, n_groups, log_joint_posterior, prior_include, var_pattern,
												*verbose, *verbose_update );
												
	PutRNGstate();
	
	//record the acceptance rates
	if( !(*just_gibbs_updates) )
	{
		acc_rts[0] = 100.*(double) results->accepted_m1 / results->proposed_m1;
		acc_rts[1] = 100.*(double) results->accepted_m2 / results->proposed_m2;
		acc_rts[2] = 100.*(double) results->accepted_m3 / results->proposed_m3;
	}
	if( !(*fixed_groups) )
	{
		acc_rts[3] = 100.*(double) results->accepted_eject / results->proposed_eject;
		acc_rts[4] = 100.*(double) results->accepted_absorb / results->proposed_absorb;
	}
	if( *variable_select )
	{
		acc_rts[5] = 100.*(double) results->accepted_include_exclude_variable / results->proposed_include_exclude_variable ;
	}

	BLCA_free_results(results,*n_iterations,*n_burn_in,*thin_by);

	BLCA_free_mixmod(mixmod);
	
	return;

}




static void BLCA_VS_COMPUTE_POST_HOC_PARAMETER_ESTIMATES(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *n_groups, int *n_sample, int *memberships, int *variable_inclusion_indicator, int *variable, double *estimate, double *sd_estimate, double *classprob_estimate, double *sd_classprob_estimate )
{

	int i,j,n,d,inG;
	struct mix_mod *mixmod;
	struct results *input;
	double **Estimate,**SE_Estimate;	

	mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *n_groups, *n_groups, hparam, ncat, TRUE, FALSE, FALSE, FALSE, FALSE );
	
	//mixmod->prior_prob_variable_include = (*prior_prob_include);
	
	input = (struct results *)malloc(sizeof(struct results));
	BLCA_allocate_results_x2(input,*n_sample,0,1,*nobs,*nvar,FALSE);

	int x;
	
	for(i=0;i<*nobs;i++){
		for(j=0;j<*nvar;j++){
			x = Y[i + j*(*nobs)];
			mixmod->Y[j][i] = x;
			mixmod->Yobs[i][j] = x;
		}
	}	
		
	for(i=0;i<*n_sample;i++){
		for(j=0;j<*nobs;j++){
			input->memberships[i][j] = memberships[ j*(*n_sample) + i ];
		}
		for(j=0;j<*nvar;j++){
			input->variable_indicator[i][j] = variable_inclusion_indicator[ j*(*n_sample) + i ];
		}
	}
	
	/*allocate the space to store estimated parameters*/
	
	Estimate = calloc(*n_groups,sizeof(double *));
	SE_Estimate = calloc(*n_groups,sizeof(double *));
	for(i=0;i<*n_groups;i++){
		Estimate[i] = calloc(mixmod->ncat[*variable],sizeof(double));
		SE_Estimate[i] = calloc(mixmod->ncat[*variable],sizeof(double));
	}

	BLCA_compute_post_hoc_parameter_estimates_for_variable( mixmod, input, *n_sample, *n_groups, *variable, Estimate, SE_Estimate );

	
	for(i=0;i<*n_groups;i++){
		for(j=0;j<mixmod->ncat[*variable];j++){
			//Rprintf("\nThe value of i = %d, the value of j = %d",i,j);
			estimate[j*(*n_groups) + i] = Estimate[i][j];
			sd_estimate[j*(*n_groups) + i] = SE_Estimate[i][j];
		}
	}
	
	BLCA_compute_post_hoc_parameter_estimates_for_class_probabilities( mixmod, input, *n_sample, *n_groups, 0, classprob_estimate, sd_classprob_estimate );


	for(i=0;i<*n_groups;i++){
		free(Estimate[i]);
		free(SE_Estimate[i]);
	}
	free(Estimate);
	free(SE_Estimate);
	
	BLCA_free_results_x2(input,*n_sample,0,1,FALSE);
	BLCA_free_mixmod(mixmod);

	return;

}

/*------------------------------------ Gibbs sampler (not collapsed)-------------------------------*/

static void BLCA_GIBBS_SAMPLER(int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *init_type, int *n_groups, int *n_iterations, int *n_burn_in, int *thin_by, int *group_memberships, double *group_weights, double *prob_variables, int *var_in, double *log_joint_posterior, int *sample_missing_data, int *n_missing, int *imputed_missing_values, int *position_missing,  int *verbose, int *verbose_update )
{

	int i,j,k,l,n,d,inG,mxG,nit,nburn,fixed,justgibbs,y_new;
  //double *u_prob;
	struct mix_mod *mixmod;
	
	//u_prob = (double *)calloc( 100, sizeof(double) );
	
	mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *n_groups, *n_groups, hparam, ncat, FALSE, TRUE, FALSE, FALSE, FALSE );
	
	if( *verbose ) Rprintf("\n  Allocation complete ... ");

	for(j=0;j<mixmod->d;j++)
		mixmod->varindicator[j] = var_in[j];
	
	BLCA_initialize_data( mixmod, Y );
	
	if( *verbose ) Rprintf("\n  Initialisation complete ... ");
	
	BLCA_initialize_priors( mixmod, alpha_prior, beta_prior, *prior_init_type );
	
	if( *verbose ) Rprintf("\n  Prior initialisation complete ... ");
	
	GetRNGstate();
	
	//initialize the algorithm
	
	BLCA_initialize_Gibbs_sampler( *init_type, mixmod );
	
	if( *verbose ) Rprintf("\n  Gibbs sampler initialisation  complete ... ");
	
	BLCA_analysis_MCMC_Gibbs_sampler(mixmod,*n_iterations,*n_burn_in,*thin_by,
												group_memberships, group_weights, prob_variables, log_joint_posterior, *sample_missing_data, *n_missing, 
												imputed_missing_values, position_missing, *verbose, *verbose_update );
	
	PutRNGstate();
	
	BLCA_free_mixmod(mixmod);
		
	return;
}

/*-----------------------------------------fit using EM algorithm----------------------------------------*/

static void BLCA_EM_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *n_groups, int *init_type, double *init_vals, int *max_iterations, int *iterations, double *group_probabilities, double *group_weights, double *prob_variables, int *var_in, double *log_likelihood, int *MAP, double *tol, double *eps, int *converged, double *max_log_likelihood )
{

		struct mix_mod *mixmod;
		int i, j;
		//if MAP == 1 then a prior regularizer is used in the EM through the hparam vector
		
		mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *n_groups, *n_groups, hparam, ncat, FALSE, TRUE, *MAP, FALSE, FALSE );
		
		for(j=0;j<mixmod->d;j++) 
			mixmod->varindicator[j] = var_in[j];
		
		GetRNGstate();
		
		//initialize by using the values in group_weights and prob_variables
		
		BLCA_initialize_EM( *init_type, init_vals , mixmod , group_weights, prob_variables);
		
		PutRNGstate();
		
		BLCA_initialize_data( mixmod, Y );	
		
		BLCA_initialize_priors( mixmod, alpha_prior, beta_prior, *prior_init_type );

		//important... need a function to initialize the algorithm...a
		//can most likely pass the  initialization in through the arg--
		//initialize_EM(mixmod,*n_groups);	
		
		BLCA_analysis_EM( mixmod, *max_iterations, iterations, group_probabilities, group_weights, prob_variables, log_likelihood, *MAP, *tol, eps, converged ) ;
		
		if( mixmod->EM_MAP ) *max_log_likelihood = mixmod->log_like - mixmod->log_prior ; else *max_log_likelihood = mixmod->log_like ; 
		
		BLCA_free_mixmod(mixmod);
		
		return;
}

/*-----------------------------------------fit using VB approach----------------------------------------*/

static void BLCA_VB_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *n_groups, int *init_type, double *init_vals, int *max_iterations, int *iterations, double *group_probabilities, double *group_weights, double *sd_group_weights, double *par_group_weights, double *prob_variables, double *sd_prob_variables, double *par_prob_variables, int *var_in, double *lower_bound, double *tol, int *converged, double *log_post )
{

		struct mix_mod *mixmod;
		int i, j;
		//if MAP == 1 then a prior regularizer is used in the EM through the hparam vector
		
		mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *n_groups, *n_groups, hparam, ncat, FALSE, TRUE, FALSE, TRUE, FALSE );

		for(j=0;j<mixmod->d;j++) 
			mixmod->varindicator[j] = var_in[j];
		
		GetRNGstate();
		
		//initialize by using the values in group_weights and prob_variables
		
		BLCA_initialize_VB( *init_type, init_vals, mixmod , group_weights, prob_variables );
		
		PutRNGstate();

		BLCA_initialize_data( mixmod, Y );
		
		BLCA_initialize_priors( mixmod, alpha_prior, beta_prior, *prior_init_type );
		
		BLCA_analysis_VB( mixmod, *max_iterations, iterations, group_probabilities, group_weights, sd_group_weights, par_group_weights, prob_variables, sd_prob_variables, par_prob_variables, lower_bound, *tol, converged, log_post );
		
		BLCA_free_mixmod( mixmod );
		
		return;
		
}

/*-----------------------------------------fit using Bootstrap algorithm----------------------------------------*/

static void BLCA_BOOT_FIT( int *Y, int *nobs, int *nvar, int *ncat, double *hparam, int *prior_init_type, double *alpha_prior, double *beta_prior, int *boot_samples, int *n_groups, double *init_group_weights, double *init_prob_variables, double *out_group_weights, double *out_prob_variables, double *out_log_likelihood, int *var_in, int *MAP, double *tol, int *max_iterations, int *verbose, int *verbose_update )
{
		struct mix_mod *mixmod;
		int i, j;
		//if MAP == 1 then a prior regularizer is used in the Bootstrap through the hparam vector
		
		mixmod = BLCA_allocate_mixmod( *nobs, *nvar, *n_groups, *n_groups, hparam, ncat, FALSE, FALSE, *MAP, FALSE, TRUE ); 
		
		for(j=0;j<mixmod->d;j++) 
			mixmod->varindicator[j] = var_in[j];
		
		GetRNGstate();
		
		//initialize by using the values in group_weights and prob_variables
		
		BLCA_initialize_Boot( mixmod, init_group_weights, init_prob_variables );
		
		BLCA_initialize_data( mixmod, Y );	
		
		BLCA_initialize_priors( mixmod, alpha_prior, beta_prior, *prior_init_type );
		
		BLCA_analysis_Boot( mixmod, *boot_samples, *max_iterations, out_group_weights, out_prob_variables, out_log_likelihood, *tol, *verbose, *verbose_update ) ;
		
		PutRNGstate();
		
		BLCA_free_mixmod(mixmod);
		
		return;
}



/*-----------------------------------------Label Switching algorithm----------------------------------------*/

static void BLCA_RELABEL( int *n_obs, int *n_sample, int *n_groups, int *labels_in, int *labels_out, int *permutation )
{

	int n,g,**raw,**relab,**summary,
		i,j,k,t,N,T,**cost,*lab;
	
	//n and no. of groups
	n= *n_obs; //n_obs[0];
	g= *n_groups;//n_groups[0];
	//N = number of samples to undo label switching for
	N= *n_sample; //n_sample[0];

	T=0;

	//allocate memory and initialize
	raw = imatrix(1,N,1,n);
	relab = imatrix(1,N,1,n);
	summary = imatrix(1,g,1,n);
	cost = imatrix(1,g,1,g+1);
	lab = ivector(1,g);


	//copy from labels_in to raw_r
	for(i=0;i<N;i++){
		for(j=0;j<n;j++){
			raw[i+1][j+1] = labels_in[ i + j*N ];
		}	
	}


	for(i=1;i<g+1;i++){
		for(j=1;j<n+1;j++){
			summary[i][j]=0;
		}
	}

	//use the first allocation in raw as the first labelling

	for(i=1;i<n+1;i++){
		summary[raw[1][i]][i]+=1;
	}

	for(i=1;i<n+1;i++){
		relab[1][i] = raw[1][i];
	}
	

	for(i=1;i<g+1;i++){
		permutation[(i-1)*N + 0] = i; //this is the identity mapping as nothing done.
	}

	for(t=2;t<N+1;t++){

	//row t
	//compute the cost matrix
	for(i=1;i<g+1;i++){
	
		for(j=1;j<g+1;j++){
	
			cost[i][j]=0;
		
			for(k=1;k<n+1;k++){
		
				if(raw[t][k]==j){
					cost[i][j]+=summary[i][k];
				}
			
			} 
		
			cost[i][j]=n*(t-1)-cost[i][j];
		
		}
		cost[i][g+1]=0;
	
		}

		T=0;
		
		assct(g,cost,lab,&T);

		//store the permutation
		for(i=1;i<g+1;i++){
			permutation[ (i-1)*N + (t-1) ] = lab[i];
		}

		//relabel based on output from assct
		//update the summary matrix
		for(i=1;i<n+1;i++){
			relab[t][i]=lab[raw[t][i]];
			summary[ lab[raw[t][i]] ][i]+=1;
		}

	}


	for(i=0;i<N;i++){
		for(j=0;j<n;j++){
			labels_out[ i + j*N ] = relab[i+1][j+1];
		}	
	}	
	
	free_imatrix(raw,1,N,1,n);
	free_imatrix(relab,1,N,1,n);
	free_imatrix(summary,1,g,1,n);
	free_imatrix(cost,1,g,1,g+1);
	free_ivector(lab,1,g);

	return;
}




