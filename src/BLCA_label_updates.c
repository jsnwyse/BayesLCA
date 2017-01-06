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
	
#include "BLCA_label_updates.h"

int BLCA_update_allocations_with_gibbs(struct mix_mod *mixmod)
{

	//printf("\nwithin Gibbs...\n");

	int i,ii,ind,j,k,g,g_prime,ig,ig_prime,d=mixmod->d,position,G=mixmod->G,*order,*x;
	double *probs,max,z;
	
	probs = calloc(G,sizeof(double));
	order = calloc(mixmod->n,sizeof(int));
	
	/*cycle through the elements in turn... maybe randomize order?*/
	
	/*randomize the order...*/
	for(i=0;i<mixmod->n;i++){
		order[i] = i;
	}
	
	/*use gsl to shuffle*/
	BLCA_random_ranshuffle(order,mixmod->n);
	
	for(i=0;i<mixmod->n_gibbs;i++){
	
		ind = order[i];
	
		/*current group*/
		g = mixmod->z[ind];
		
		/*point to the entry*/
		x = mixmod->Yobs[ind];
		
		/*cycle through remaining groups*/
		for(k=0;k<G;k++){
			
			if(!(k==g)){
			
				/*get the log probability ratio for this */
				probs[k] = BLCA_get_log_relative_probability_for_gibbs(x,mixmod->components[mixmod->whereis[k]],mixmod->components[mixmod->whereis[g]],mixmod);
				
			//	printf("\nprobs[%d] = %.10f",k,probs[k]);
			
			}else{
			
				/*if k==g this is easy*/
			
				probs[k] = 0.0;
			
			}
		
		}
		
		max = BLCA_get_max(probs,G);
		
		z = 0.;
		
		for(k=0;k<G;k++){
			probs[k] -= max;
			probs[k] = exp(probs[k]);
			z += probs[k];
		}
			
		for(k=0;k<G;k++){
			probs[k] /= z;
		}
					
		/*sample allocation from the vector of weights*/
		g_prime = BLCA_sample_discrete(probs,G);
		
		//printf("\nNumber of groups = %d and g_prime = %d",G,g_prime);
		
		
		/*only need to execute this if the labelling is different*/
		if(!(g_prime == g)){
			
			ig = mixmod->whereis[g];
			ig_prime = mixmod->whereis[g_prime];
			
			//printf("\nNumber of groups = %d and ig = %d and ig_prime = %d",G,ig,ig_prime);
			/*need to recompute component[g_prime]->log_prob and component[g]->log_prob and 
				update the relevant stored quantities*/
			
			mixmod->z[ind] = g_prime;
			
			/*udpate component g*/
			mixmod->components[ig]->n_g -= 1;
			mixmod->components[ig_prime]->n_g += 1;
			
			for(j=0;j<d;j++){
			
				/*update the component counts*/
				
				ii = mixmod->Y[j][ind];
					
				mixmod->components[ig]->N[j][ ii ] -= 1;
				mixmod->components[ig_prime]->N[j][ ii ] += 1;
				
			}
			
			/*recompute the marginal likelihood of both components*/
			
			BLCA_recompute_marginal_likelihood_component(mixmod->components[ig],mixmod);
			BLCA_recompute_marginal_likelihood_component(mixmod->components[ig_prime],mixmod);
				
		}
	
	}
	
	
	free(probs);
	free(order);

	//printf("\nleaving Gibbs...\n");
	
	return(TRUE);
}

double BLCA_get_log_relative_probability_for_gibbs(int *x,struct component *component_k,struct component *component_g,struct mix_mod *mixmod)
/*compute log of the relative probability for assigning entry where x is a pointer to the datapoint*/
{
	int i,j;
	double t1,t2,t3,t4,ll,I;

	
	t1 = 0.;
	t2 = 0.;
	
	for(j=0;j<mixmod->d;j++){
	
		if(mixmod->varindicator[j]){
	
			t1 += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) -lgamma(component_k->n_g + 1. + mixmod->ncat[j]*mixmod->beta);
			t2 += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) -lgamma(component_g->n_g - 1. + mixmod->ncat[j]*mixmod->beta);
		
			for(i=0;i<mixmod->ncat[j];i++){
			
					I = (x[j] == i) ? 1 : 0 ;
				
					t1 += lgamma(component_k->N[j][i] + I + mixmod->beta);
					t2 += lgamma(component_g->N[j][i] - I + mixmod->beta);
			
			}
		
		}
		
	}
	
	/*add all up*/
	
	t3 = lgamma(component_k->n_g+1.+mixmod->alpha)+lgamma(component_g->n_g-1.+mixmod->alpha) 
		+ t1 + t2;
	
	ll = t3 - component_k->log_prob - component_g->log_prob;
	
	return(ll);

}

int BLCA_update_allocations_with_metropolis_move_1(struct mix_mod *mixmod,int *accepted,int *proposed)
	/*this performs the move M1 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	//printf("\nwithin Move 1...\n");
	
	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return(TRUE);
	
	}

	*proposed += 1;

	int i,ii,kk,j,k,g1,g2,ig1,ig2,ntot,*indexes,*proposed_alloc;
	double p,log_acceptance,ag1=1.,ag2=1.;
	struct component *component_g1,*component_g2;
	
	/*the integers g1 and g2 give the components, and ig1 and ig2 their whereis value*/
	/*the doubles ag1 and ag2 give the parameters to the Beta distribution for generating p_1
		take default uniform*/
	
	/*sample the two components*/
	g1 = (int) ( runif(0.0,1.0) * mixmod->G ) ;//gsl_rng_uniform_int(r,mixmod->G);
	g2 = g1;
	while(g2 == g1){
		g2 =  (int) ( runif(0.0,1.0) * mixmod->G ) ;
	}
	
	/*find where in mixmod->components g1 and g2 are*/
	ig1 = mixmod->whereis[g1];
	ig2 = mixmod->whereis[g2];
	
	ntot = mixmod->components[ig1]->n_g + mixmod->components[ig2]->n_g;
	
	/*DO NOT PERFORM THIS MOVE UNLESS ntot > 0*/
	if(ntot == 0){
	
		return(TRUE);
		
	}
	
	/*allocate space for proposal stuff...*/
	
	component_g1 =(struct component *)malloc(sizeof(struct component));
	component_g2 =(struct component *)malloc(sizeof(struct component));
	
	BLCA_allocate_component(component_g1,mixmod);
	BLCA_allocate_component(component_g2,mixmod);
	
	/*allocate a vector to keep track of indexes*/
	indexes = calloc(ntot,sizeof(int));
	proposed_alloc = calloc(ntot,sizeof(int));
	
	k=0;
	for(i=0;i<mixmod->n;i++){
		if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
			indexes[k] = i;
			k+=1;
		}
	}
	
	/*check is ok*/
	if(!(k==ntot)){
		//mixmod_warning(2);
		//printf("\nValue of k=%d and ntot = %d, g1 = %d, g2 = %d",k,ntot,g1,g2);
		//for(i=0;i<mixmod->n;i++) printf("%d,",mixmod->z[i]);
	}
	
	/*generate p and begin reallocation*/
	p = rbeta(ag1,ag2);
	
	/*now propose reallocation*/
	
	component_g1->n_g = 0;
	component_g2->n_g = 0;
	
	for(i=0;i<ntot;i++){
	
		ii = indexes[i];
	
		if( runif(0.0,1.0) < p){
		
			/*reallocate to g1*/
			proposed_alloc[i] = g1;
			component_g1->n_g += 1;
			for(j=0;j<mixmod->d;j++){
				component_g1->N[j][ mixmod->Y[j][ii] ] += 1;
			}
			
		}else{
		
			/*reallocate to g2*/
			proposed_alloc[i] = g2;
			component_g2->n_g += 1;
			for(j=0;j<mixmod->d;j++){
				component_g2->N[j][ mixmod->Y[j][ii] ] += 1;
			}
		
		}		
		
	}
	
	/*evaluate log of acceptance probability*/
	
	BLCA_recompute_marginal_likelihood_component(component_g1,mixmod);
	BLCA_recompute_marginal_likelihood_component(component_g2,mixmod);
	
	log_acceptance = component_g1->log_prob + component_g2->log_prob - mixmod->components[ig1]->log_prob - mixmod->components[ig2]->log_prob
	+lgamma(mixmod->alpha + mixmod->components[ig1]->n_g) + lgamma(mixmod->alpha + mixmod->components[ig2]->n_g) 
			- lgamma(mixmod->alpha + component_g1->n_g) - lgamma(mixmod->alpha + component_g2->n_g);;
	
	//printf("\nThe value of log acceptance is %.10f,",log_acceptance);
	
	if(log(runif(0.0,1.0))<log_acceptance){
	
		/*do the swap!!!*/
		
		*accepted += 1;
		
		for(i=0;i<ntot;i++){
			mixmod->z[indexes[i]] = proposed_alloc[i];
		}
		
		/*copy over the accepted components*/
		
		BLCA_copy_component(component_g1,mixmod->components[ig1],mixmod);
		BLCA_copy_component(component_g2,mixmod->components[ig2],mixmod);
	
	}
	

	BLCA_free_component(component_g1,mixmod);
	BLCA_free_component(component_g2,mixmod);
	
	free(component_g1);
	free(component_g2);
	
	free(indexes);
	free(proposed_alloc);

	//printf("\nleaving Move 1...\n");
	
	return(TRUE);
}

int BLCA_update_allocations_with_metropolis_move_2(struct mix_mod *mixmod,int *accepted,int *proposed)
/*this performs the move M2 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return(TRUE);
	
	}

	int i,ii,j,k,g1,g2,ig1,ig2,curr_n_g1,curr_n_g2,m,ntot,c=0;
	int *indexes,*order;
	double log_acceptance;
	struct component *component_g1,*component_g2;
	
	/*sample the two components*/
	g1 =  (int) ( runif(0.0,1.0) * mixmod->G );
	g2 = g1;
	while(g2 == g1){
		g2 =  (int) ( runif(0.0,1.0) * mixmod->G );
	}	

	/*find where in mixmod->components g1 and g2 are*/
	ig1 = mixmod->whereis[g1];
	ig2 = mixmod->whereis[g2];
	
	if(mixmod->components[ig1]->n_g == 0){
		/*cannot perform move in this case... return*/
		return(TRUE);
	}
	
	/*current component sizes*/
	curr_n_g1 = mixmod->components[ig1]->n_g;
	curr_n_g2 = mixmod->components[ig2]->n_g;
	
	/*allocate the candidate components*/
	
	component_g1 =(struct component *)malloc(sizeof(struct component));
	component_g2 =(struct component *)malloc(sizeof(struct component));
	
	BLCA_allocate_component(component_g1,mixmod);
	BLCA_allocate_component(component_g2,mixmod);
	
	
	*proposed += 1;
	
	order = calloc(curr_n_g1,sizeof(int));
	for(i=0;i<curr_n_g1;i++){
		order[i] = i;
	}
	
	/*shuffle the order*/
	BLCA_random_ranshuffle(order,curr_n_g1);
	
	indexes = calloc(curr_n_g1,sizeof(int));
	for(i=0;i<mixmod->n;i++){
		if(mixmod->z[i] == g1){
			indexes[c] = i;
			c += 1;
		}
	}
	
	m = (int)( runif(0.0,1.0) * curr_n_g1 ) ;  //gsl_rng_uniform_int(r,curr_n_g1);
	
	/*try to reallocate the first i=1,...,m in indexes[order[i]]*/
	
	BLCA_copy_component(mixmod->components[ig1],component_g1,mixmod);
	BLCA_copy_component(mixmod->components[ig2],component_g2,mixmod);
	
	component_g1->n_g -= m;
	component_g2->n_g += m;
	
	for(i=0;i<m;i++){
	
		k = order[i];
		ii = indexes[k];
		
		for(j=0;j<mixmod->d;j++){
		
			component_g1->N[j][ mixmod->Yobs[ii][j] ] -= 1;		
			component_g2->N[j][ mixmod->Yobs[ii][j] ] += 1;	
			
		}
		
	}
	
	/*compute the log acceptance probability*/
	
	BLCA_recompute_marginal_likelihood_component(component_g1,mixmod);
	BLCA_recompute_marginal_likelihood_component(component_g2,mixmod);
	
	log_acceptance = component_g1->log_prob + component_g2->log_prob - mixmod->components[ig1]->log_prob - mixmod->components[ig2]->log_prob
							+ log(component_g1->n_g + m) - log(component_g2->n_g) + lgamma(component_g1->n_g+m+1.) + lgamma(component_g2->n_g-m+1.)
								-lgamma(component_g1->n_g + 1.) - lgamma(component_g2->n_g+1.);
	
	if(log( runif(0.0,1.0) )<log_acceptance){
		
		/*do the swap*/
		
		*accepted += 1;
		
		for(i=0;i<m;i++){
			k = order[i];
			ii = indexes[k];
			mixmod->z[ii] = g2;		
		}
		
		/*copy over the accepted components*/
		BLCA_copy_component(component_g1,mixmod->components[ig1],mixmod);
		BLCA_copy_component(component_g2,mixmod->components[ig2],mixmod);
	
	}
	

	BLCA_free_component(component_g1,mixmod);
	BLCA_free_component(component_g2,mixmod);
	
	free(component_g1);
	free(component_g2);

	free(order);
	free(indexes);
	
	return(TRUE);
}



int BLCA_update_allocations_with_metropolis_move_3(struct mix_mod *mixmod,int *accepted,int *proposed)
/*this performs the move M3 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return(TRUE);
	
	}

	int i,ii,j,k,g1,g2,ig1,ig2,curr_n_g1,curr_n_g2,m,ntot,c=0,d = mixmod->d,identify_g1,identify_g2,id;
	int *indexes,*order,*proposed_allocation;
	double w,log_acceptance,log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.,l1,l2,p1,p2,max;
	struct component *component_g1,*component_g2,*component_backward_g1,*component_backward_g2;
	
	*proposed += 1;
	
	/*sample the two components*/
	g1 =  (int) ( runif(0.0,1.0) * mixmod->G );
	g2 = g1;
	while(g2 == g1){
		g2 =  (int) ( runif(0.0,1.0) * mixmod->G );
	}	

	/*find where in mixmod->components g1 and g2 are*/
	ig1 = mixmod->whereis[g1];
	ig2 = mixmod->whereis[g2];
	
	ntot = mixmod->components[ig1]->n_g + mixmod->components[ig2]->n_g;
	
	if(ntot == 0){
		/*cannot perform if both components empty*/
		return(TRUE);
	}
	
	indexes = calloc(ntot,sizeof(int));
	order = calloc(ntot,sizeof(int));
	proposed_allocation = calloc(ntot,sizeof(int));
	
	/*this move can still be done if either component empty*/
	
	for(i=0;i<mixmod->n;i++){	
		if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
			indexes[c] = i;
			c += 1;
		}
	}
	
	for(i=0;i<ntot;i++){
		order[i] = i;
	}
	
	/*shuffle the order*/
	BLCA_random_ranshuffle(order,ntot);
	
	/*allocate the candidate components*/
	
	component_g1 =(struct component *)malloc(sizeof(struct component));
	component_g2 =(struct component *)malloc(sizeof(struct component));
	
	BLCA_allocate_component(component_g1,mixmod);
	BLCA_allocate_component(component_g2,mixmod);
	
	component_g1->n_g = 0;
	component_g2->n_g = 0;

	
	/*randomizing this part should give a 50% acceptance*/

	k = order[0];
	ii = indexes[k];
	
	if(/*gsl_rng_uniform(r) < 0.5*/mixmod->z[ii] == g1){
	
		/*put it into component_g1*/
		
		component_g1->n_g += 1;
		
		for(j=0;j<d;j++){
			component_g1->N[j][ mixmod->Y[j][ii] ] += 1;
		}
		
		//identify_g1 = mixmod->z[ii];
	
		proposed_allocation[0] = g1;
		
		log_transition_z_to_zprime = log(0.5) - log(1.);
	
	}else{
	
		/*put it into component_g2*/
		
		component_g2->n_g += 1;
		
		for(j=0;j<d;j++){
			component_g2->N[j][ mixmod->Y[j][ii] ] += 1;
		}
	
		proposed_allocation[0] = g2;

		log_transition_z_to_zprime = log(0.5) - log(1.);

	}
	
	/*identify this as being in the same component as original
		...this allows the computing of the backwards probability*/
	
	log_transition_zprime_to_z += log(0.5) - log(1.);
		
	
	for(i=1;i<ntot;i++){
	
		k = order[i];
		ii = indexes[k];
		
		/*do the proposed component values, then sample and make changes*/
		
		/*compute probability generated from g1*/
		
		l1 = BLCA_compute_log_data_probability_with_inclusion_in_component(mixmod->Yobs[ii],component_g1,mixmod)
				+ BLCA_compute_log_data_probability_component(component_g2,mixmod);
		
		l2 = BLCA_compute_log_data_probability_component(component_g1,mixmod)
				+ BLCA_compute_log_data_probability_with_inclusion_in_component(mixmod->Yobs[ii],component_g2,mixmod);
		
		w = ((mixmod->alpha + component_g1->n_g)/(mixmod->alpha + component_g2->n_g))*exp(l1-l2);
		
		p1 = w/(1.+w);
		
		/*make a draw*/
		
		if( runif(0.0,1.0) < p1){
		
			/*put it in g1*/
			
			component_g1->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g1->N[j][ mixmod->Y[j][ii] ] += 1;
			}
			
			if(!(mixmod->z[ii] == g1)){
				log_transition_z_to_zprime += log(p1);
				log_transition_zprime_to_z += log(1.-p1);
			}
			
			proposed_allocation[i] = g1;
			
		}else{
		
		
			/*put it in g2*/
			
			component_g2->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g2->N[j][ mixmod->Y[j][ii] ] += 1;
			}
			
			if(!(mixmod->z[ii] == g2)){
				log_transition_z_to_zprime += log(1.-p1);
				log_transition_zprime_to_z += log(p1);
			}
		
			proposed_allocation[i] = g2;
		}
		
	}
	

	/*compute the acceptance probability*/
	
	BLCA_recompute_marginal_likelihood_component(component_g1,mixmod);
	BLCA_recompute_marginal_likelihood_component(component_g2,mixmod);
	
	log_acceptance = component_g1->log_prob + component_g2->log_prob 
						- mixmod->components[ig1]->log_prob - mixmod->components[ig2]->log_prob
						+ log_transition_zprime_to_z - log_transition_z_to_zprime;
	
	//printf("\n%.4f \t %.4f \t %.4f ",log_acceptance,log_transition_zprime_to_z,log_transition_z_to_zprime );
						
	if(log( runif(0.0,1.0) )<log_acceptance){
	
		*accepted += 1;
	
		/*accept the move and update all quantities*/
		
		/*allocations first*/
	
		BLCA_copy_component(component_g1,mixmod->components[ig1],mixmod);
		BLCA_copy_component(component_g2,mixmod->components[ig2],mixmod);
		
		
		for(i=0;i<ntot;i++){
				
			k = order[i];
			ii = indexes[k];
			mixmod->z[ii] = proposed_allocation[i];
			
		}
	
	}
	
	/*free up all memory*/
	
	BLCA_free_component(component_g1,mixmod);
	BLCA_free_component(component_g2,mixmod);
	
	free(component_g1);
	free(component_g2);

	free(indexes);
	free(order);
	free(proposed_allocation);
	
	return(TRUE);

}


