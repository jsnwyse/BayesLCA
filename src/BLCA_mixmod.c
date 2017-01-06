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


#include "BLCA_mixmod.h"

struct mix_mod *allocate_mixmod(int datasize, int datadimension, int maxgroups, int initgroups,double *prior_hparams,int *ncat, int collapsed)
/*this function allocates and returns a pointer to a mixmod structure...*/
{

	int i,j,k;
	FILE *fp;
	
	struct mix_mod *mixmod = (struct mix_mod *)malloc(sizeof(struct mix_mod));

	mixmod->maxgroups = maxgroups;
	mixmod->G = initgroups;
	mixmod->n = datasize;
	mixmod->d = datadimension;
	mixmod->collapsed = collapsed;
	
	mixmod->Y = calloc(datadimension,sizeof(int *));
	for(i=0;i<datadimension;i++){
		mixmod->Y[i] = calloc(datasize,sizeof(int));
	}
	
	mixmod->Yobs = calloc(datasize,sizeof(int *));
	for(i=0;i<datasize;i++){
		mixmod->Yobs[i] = calloc(datadimension,sizeof(int *));
	}
	
	mixmod->z = calloc(datasize,sizeof(int));
	
	mixmod->ncat = calloc(datadimension,sizeof(int));
	for(i=0;i<datadimension;i++){
		mixmod->ncat[i] = ncat[i];
	}
	
	mixmod->varindicator = calloc(datadimension,sizeof(int));
	
	/*allocate this memory-- only initialize what is needed for initial conditions*/
	mixmod->components = (struct component **)malloc(sizeof(struct component *)*maxgroups);

	for(i=0;i<maxgroups;i++){
		mixmod->components[i] = (struct component *)malloc(sizeof(struct component));
	}
	
	mixmod->undiscriminating = (struct component *)malloc(sizeof(struct component));
	allocate_component(mixmod->undiscriminating,mixmod);
	mixmod->undiscriminating->n_g = mixmod->n; /*there is always n elements in here as these variables
																	do not define a clustering...*/

	
	/*allocate whereis*/
	
	mixmod->whereis = calloc(maxgroups,sizeof(int));
	for(i=0;i<maxgroups;i++)
		mixmod->whereis[i] = -1;
	
	/*allocate for the initial number of components*/
	for(i=0;i<initgroups;i++){
		allocate_component(mixmod->components[i],mixmod);
		mixmod->components[i]->in_use = TRUE;
		mixmod->components[i]->n_g = 0;
	}
	for(i=initgroups;i<maxgroups;i++){
		allocate_component(mixmod->components[i],mixmod);
		mixmod->components[i]->in_use = FALSE;
		mixmod->components[i]->n_g = 0;
	}

	
	/*there will be two hparameters in the default model:
		alpha: dirichlet prior on weights symmetric
		beta: dirichlet prior for within component membership probs
		*/
	
	mixmod->alpha = prior_hparams[0];
	mixmod->beta = prior_hparams[1];

	mixmod->log_prior_G = calloc(maxgroups+1,sizeof(double));
	
	/*assign space for the lookup table for a values*/
	
	//mixmod->table_a = calloc(datasize,sizeof(double));
	
	if(!mixmod->collapsed) mixmod->weights = calloc(mixmod->G,sizeof(double));
	
	return(mixmod);

}

void free_mixmod(struct mix_mod *mixmod)
/*frees the memory used by mixmod object*/
{
	int n = mixmod->n,d = mixmod->d, G = mixmod->G,i,j,k;
	
	/*free up components*/
	for(k=0;k<mixmod->maxgroups;k++){
		free_component(mixmod->components[k],mixmod);
		free(mixmod->components[k]);
	}
	free(mixmod->components);
	
	free_component(mixmod->undiscriminating,mixmod);
	free(mixmod->undiscriminating);
	
	/*free wehreis*/
	free(mixmod->whereis);
	
	/*free data*/
	for(i=0;i<d;i++){
		free(mixmod->Y[i]);
	}
	free(mixmod->Y);
	
	for(i=0;i<n;i++){
		free(mixmod->Yobs[i]);
	}
	free(mixmod->Yobs);
	
	/*free others*/
	free(mixmod->z);
	
	free(mixmod->ncat);
	
	free(mixmod->varindicator);
	
	free(mixmod->log_prior_G);
	
	/*free table*/
	//free(mixmod->table_a);
	
	if(!mixmod->collapsed) free(mixmod->weights);
	
	free(mixmod);
	
	return;
}


double log_normalizing_constant_model(int G,struct mix_mod *mixmod)
/*returns the log of the normalizing constant for a model with G components*/
{
	double z;
	
	z =  lgamma(G*mixmod->alpha) - G*lgamma(mixmod->alpha) - lgamma(mixmod->n+G*mixmod->alpha);
	
	return(z); 	
}

double l_prior_variable_include(int D,struct mix_mod *mixmod)
{
	
	double l;
	
	l = D*log(mixmod->prior_prob_variable_include) + (mixmod->d - D)*log(1.-mixmod->prior_prob_variable_include);
	
	return(l);

}

double get_full_log_posterior(struct mix_mod *mixmod)
{

	double log_full_posterior = 0.;
	int i,d;

	/*model normalizing constant*/
	
	log_full_posterior += log_normalizing_constant_model(mixmod->G,mixmod);
	
	/*components - discriminating*/
	
	for(i=0;i<mixmod->G;i++){
		log_full_posterior += mixmod->components[ mixmod->whereis[i] ]->log_prob;
	}
	
	/*undiscriminating*/
	
	log_full_posterior += mixmod->undiscriminating->log_prob;
	
	/*prior on variable inclusion*/
	d = 0;
	for(i=0;i<mixmod->d;i++) d += mixmod->varindicator[i];
	
	log_full_posterior += l_prior_variable_include(d,mixmod);
	
	log_full_posterior += mixmod->log_prior_G[mixmod->G];
	
	return(log_full_posterior);
	
	
}

double get_full_log_posterior_x2(struct mix_mod *mixmod)
{

	double log_full_posterior;
	int j,i,c;

	/*model normalizing constant*/
	
	log_full_posterior = lgamma(mixmod->G*mixmod->alpha) - mixmod->G*lgamma(mixmod->alpha);
	
	for(i=0;i<mixmod->G;i++){
		log_full_posterior += (mixmod->components[i]->n_g + mixmod->alpha -1.)*log(mixmod->weights[i]);
		for(j=0;j<mixmod->d;j++){
			for(c=0;c<mixmod->ncat[j];c++){
				log_full_posterior += (mixmod->components[i]->N[j][c] + mixmod->beta - 1.)*log(mixmod->components[i]->prob_variables[j][c]);
			}
			log_full_posterior += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta);
		}
	}
	
	return(log_full_posterior);
}


struct results *do_mixmod_analysis( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int fixedG, int onlyGibbs, int selectVariables, int *group_memberships, int *variable_inclusion_indicator, int *n_groups, double *log_joint_posterior, double *prior_include )
/*fixedG takes the value either TRUE or FALSE as defined in the macros*/
{

	int i,j,k,l,itmod,d_in,ej_case,vs_case,maxgroups = mixmod->maxgroups;
	
	int gap = (num_iteration - num_burnin)/(thin_by);
	
	//writeToFile = FALSE;
	
	//struct results res,*results;
	//results = &res;

	double pr_ej_G,pr_ej_Gm1,pr_ej_Gp1;
	
	struct results *results = (struct results *)malloc(sizeof(struct results));
	allocate_results( results, num_iteration, num_burnin, thin_by, mixmod->n, mixmod->d  );
	
	/*initialize with all variables in model*/
	for(j=0;j<mixmod->d;j++){
	   mixmod->varindicator[j] = TRUE;
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
	
		update_allocations_with_gibbs(mixmod);
		
		if(!onlyGibbs){
			update_allocations_with_metropolis_move_1(mixmod,&(results->accepted_m1),&(results->proposed_m1));
			update_allocations_with_metropolis_move_2(mixmod,&(results->accepted_m2),&(results->proposed_m2));
			update_allocations_with_metropolis_move_3(mixmod,&(results->accepted_m3),&(results->proposed_m3));
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
				update_allocations_with_ejection_move(mixmod,&(results->accepted_eject),&(results->proposed_eject),pr_ej_G,pr_ej_Gp1);
			}else{
				update_allocations_with_absorb_move(mixmod,&(results->accepted_absorb),&(results->proposed_absorb),pr_ej_G,pr_ej_Gm1);
			}
		
		
		}
		
		if(selectVariables){
			j = random_integer(mixmod->d);
			update_model_by_variable_include_exclude(mixmod,&(results->accepted_include_exclude_variable),&(results->proposed_include_exclude_variable),j);
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
		log_joint_posterior[l] = get_full_log_posterior(mixmod);

	}
	
	
	if(selectVariables){
		for(i=0;i<mixmod->d;i++){
			results->variable_prob_inclusion[i] /= ((num_iteration-num_burnin)/thin_by);
		}
	}
	
	return(results);
	
}


void do_mixmod_analysis_not_collapsed( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int *group_memberships, double *group_weights, double *variable_probabilities, double *log_joint_posterior )
{

	//this is the non-collapsed form of the model which uses Gibbs updates for all unknowns...
	// in this model there is no search for the number of groups, and no variable selection moves
	
	int i,j,k,l,c,p,idx,g_new,ind,itmod,*order;
	
	double *w,s,m;
	w = calloc(mixmod->G,sizeof(double));
	
	order = calloc(mixmod->n,sizeof(int));
	
	double **v; 
	v = calloc(mixmod->d,sizeof(double *));
	for(i=0;i<mixmod->d;i++) v[i] = calloc(mixmod->ncat[i],sizeof(double));
	
	int gap = (num_iteration - num_burnin)/(thin_by), gap_;
	
	for(l=0;l<num_iteration;l++){
	
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
				s = 0.;
				for(k=0;k<mixmod->ncat[j];k++){
					 v[j][k] =  rgamma( (double)mixmod->components[i]->N[j][k] + mixmod->beta , 1. ) ;
					 s += v[j][k];
				}
				for(k=0;k<mixmod->ncat[j];k++) mixmod->components[i]->prob_variables[j][k] = v[j][k]/s;
			}
		}
	
		//sample the memberships in a random order
		
		for(i=0;i<mixmod->n;i++){
			order[i] = i;
		}
		
		random_ranshuffle( order, mixmod->n );
		
		for(k=0;k<mixmod->n;k++){
			
			idx = order[k];
			
			for(i=0;i<mixmod->G;i++){
				w[i] = log(mixmod->weights[i]);
				for(j=0;j<mixmod->d;j++){
					c = mixmod->Y[j][idx];
					w[i] += log(mixmod->components[i]->prob_variables[j][c]);		
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
			
			
			
			g_new = sample_discrete( w, mixmod->G );
			
			
			if(g_new != mixmod->z[idx])
			{
				//take out of component counts in old and put into new
				for(j=0;j<mixmod->d;j++)
				{
					c = mixmod->Y[j][idx];
					mixmod->components[ mixmod->z[idx] ]->N[j][c] -= 1;
					mixmod->components[ g_new ]->N[j][c] += 1;
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
					
		}	
		

		log_joint_posterior[l] = get_full_log_posterior_x2(mixmod);
		
	}

	free(w);
	for(i=0;i<mixmod->d;i++) free(v[i]);
	free(v);
	
	return;

}

