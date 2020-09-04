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
	
#include "BLCA_results.h"


void BLCA_allocate_results(struct results *results,int iterations,int burn_in,int thin_by,int len,int d)
/*allocates space to store post burn-in iterations*/
{
	
	int i, N = (int)(iterations/thin_by);
	
	//results->ngroups = calloc(N,sizeof(int));
	//results->variable_prob_inclusion = calloc(d,sizeof(double));
	
	BLCA_reset_results( results );
	
	return;
}

void BLCA_reset_results(struct results *results)
{
	results->proposed_m1 = 0; results->accepted_m1 = 0;
	results->proposed_m2 = 0; results->accepted_m2 = 0;
	results->proposed_m3 = 0; results->accepted_m3 = 0;
	results->proposed_eject = 0; results->accepted_eject = 0;
	results->proposed_absorb = 0; results->accepted_absorb = 0;
	results->proposed_add_variable = 0; results->accepted_add_variable = 0;
	results->proposed_remove_variable = 0; results->accepted_remove_variable = 0;
	results->proposed_include_exclude_variable = 0; results->accepted_include_exclude_variable = 0;
	return;
}

void BLCA_free_results(struct results *results,int iterations,int burn_in,int thin_by )
{
	
	//free(results->ngroups);
	//free(results->variable_prob_inclusion);
	
	free(results);
	
	return; 
}

void BLCA_allocate_results_x2(struct results *results,int iterations,int burn_in,int thin_by,int len,int d,int writeToFile)
/*allocates space to store post burn-in iterations*/
{
	
	int i,N = (iterations-burn_in)/thin_by;
	
	if(!writeToFile)
	{
	
		results->memberships = calloc(N,sizeof(int *));
		for(i=0;i<N;i++){
			results->memberships[i] = calloc(len,sizeof(int));
		}
	
		results->MAP_memberships = calloc(len,sizeof(int));
	
		results->variable_indicator = calloc(N,sizeof(int *));
		for(i=0;i<N;i++){
			results->variable_indicator[i] = calloc(d,sizeof(int));
		}
	
		results->log_posterior = calloc(iterations,sizeof(double));
	}
	
	results->ngroups = calloc(N,sizeof(int));
	results->variable_prob_inclusion = calloc(d,sizeof(double));
	
	results->proposed_m1 = 0;
	results->accepted_m1 = 0;
	results->proposed_m2 = 0;
	results->accepted_m2 = 0;
	results->proposed_m3 = 0;
	results->accepted_m3 = 0;
	results->proposed_eject = 0;
	results->accepted_eject = 0;
	results->proposed_absorb = 0;
	results->accepted_absorb = 0;
	results->proposed_add_variable = 0;
	results->accepted_add_variable = 0;
	results->proposed_remove_variable = 0;
	results->accepted_remove_variable = 0;
	results->proposed_include_exclude_variable = 0;
	results->accepted_include_exclude_variable = 0;
	
	return;
}

void BLCA_free_results_x2(struct results *results,int iterations,int burn_in,int thin_by,int writeToFile)
{
	int i,N = (iterations-burn_in)/thin_by;
	
	if(!writeToFile)
	{
	
		for(i=0;i<N;i++){
			free(results->memberships[i]);
			free(results->variable_indicator[i]);
		}
		free(results->memberships);
		free(results->variable_indicator);
	
		free(results->MAP_memberships);
	
		free(results->log_posterior);
	
	}
	free(results->ngroups);
	free(results->variable_prob_inclusion);
	
	free(results);
	
	return; 
}

int write_out_results(struct results *results,int N,int datasize,int datadimension,int onlyGibbs,int fixedG,int selectVariables)
/*write results out to file*/
{

	int i,j,k,Gmax,*freqG;
	
	if(!onlyGibbs){
	
		/*print acceptance rate to screen*/
		Rprintf("\n The acceptance rate for metropolis move 1 was : %.10f \n",(double)results->accepted_m1/(double)results->proposed_m1);
		Rprintf("\n The acceptance rate for metropolis move 2 was: %.10f \n",(double)results->accepted_m2/(double)results->proposed_m2);
		Rprintf("\n The acceptance rate for metropolis move 3 was: %.10f \n",(double)results->accepted_m3/(double)results->proposed_m3);
	
	}
	
	if(!fixedG){
	
		Rprintf("\n The acceptance rate for ejections was: %.10f \n",(double)results->accepted_eject/(double)results->proposed_eject);
		Rprintf("\n The acceptance rate for absorb moves was: %.10f \n",(double)results->accepted_absorb/(double)results->proposed_absorb);
	
	
		/*count up frequencies of numbers of components*/
		k = N;//results->niterations-results->nburnin; /*no. of elements in vector*/
	
		Gmax = BLCA_get_imax(results->ngroups,k);
	
		/*allocate a vector for counts*/
		freqG = calloc(Gmax,sizeof(int));
		for(i=0;i<k;i++){
			freqG[results->ngroups[i]-1] += 1;
		}
	
		/*print out*/
		Rprintf("\n\n\t -- Frequency in each number of components --\n\n\t G \t\t\t Posterior probability\n\n");
		for(i=0;i<Gmax;i++){
		Rprintf("\t %d \t\t\t %.5f \n",i+1,(double)freqG[i]/(double)(k));
		}
	
	
	
		/*for(i=0;i<N;i++){ with
			for(j=0;j<datasize-1;j++){
				fprintf(fp1,"%d,",results->memberships[i][j]);
			}
			fprintf(fp1,"%d\n",results->memberships[i][datasize-1]);
		}
	
		for(i=0;i<N;i++){
			fprintf(fp2,"%d\n",results->ngroups[i]);
		}*/

		free(freqG);
	
	}
	
	if(selectVariables){
	
		Rprintf("\n The acceptance rate for adding/removing variables was: %.10f \n",(double)results->accepted_include_exclude_variable/(double)results->proposed_include_exclude_variable);
		Rprintf("\n\nThe inclusion probability of each of the variables was\n\n");
		for(i=0;i<datadimension;i++){
			Rprintf("\nVariable %d \t\t\t %.5f",i+1,results->variable_prob_inclusion[i]);
		}
		Rprintf("\nNote that this calculation is based only on an empirical calculation\n using the MCMC output of occurances of a variable being in/out");
		

	}

	return(TRUE);
}



