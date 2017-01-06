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
			
	Last modification of this code: Wed 11 May 2016 04:03:05 IST    */

#include "BLCA_post_hoc.h"

void BLCA_compute_post_hoc_parameter_estimates_for_variable( struct mix_mod *mixmod, struct results *input, int n_sample, int n_groups, int variable, double **Estimate, double **SE_Estimate )
{

	int i,j,g,k,z,c,s,*items,**counts,*ng;
	double ***expected,***variance;
	
	ng = calloc(mixmod->G,sizeof(int));
	
	expected = calloc(n_sample,sizeof(double **));
	variance = calloc(n_sample,sizeof(double **));
	for(i=0;i<n_sample;i++){
		expected[i] = calloc(mixmod->G,sizeof(double *));
		variance[i] = calloc(mixmod->G,sizeof(double *));
		for(j=0;j<mixmod->G;j++){
			expected[i][j] = calloc(mixmod->ncat[variable],sizeof(double));
			variance[i][j] = calloc(mixmod->ncat[variable],sizeof(double));
		}
	}
	
	items = calloc(mixmod->n,sizeof(int));
	counts = calloc(mixmod->G,sizeof(int *));
	for(g=0;g<mixmod->G;g++){
		counts[g] = calloc(mixmod->ncat[variable],sizeof(int));
	}
	
	for(s=0;s<n_sample;s++){
	
		for(g=0;g<mixmod->G;g++){
			for(j=0;j<mixmod->ncat[variable];j++){
				counts[g][j] = 0;
			}
			ng[g] = 0;
		}
	
		for(i=0;i<mixmod->n;i++){
		
			z = input->memberships[s][i];
			
			c = mixmod->Yobs[i][variable];
			
			counts[z][c] += 1;
			
			ng[z] += 1;
		
		}
		
		for(g=0;g<mixmod->G;g++){
			for(j=0;j<mixmod->ncat[variable];j++){
				expected[s][g][j] = (counts[g][j] + mixmod->beta)/(ng[g]+(mixmod->ncat[variable]-1.)*mixmod->beta);
				variance[s][g][j] = (counts[g][j]+mixmod->beta)*(ng[g] - counts[g][j] + mixmod->ncat[variable]*mixmod->beta)/(pow(ng[g] + (mixmod->ncat[variable]-1.)*mixmod->beta,2.)*(ng[g]+(mixmod->ncat[variable]-1.)*mixmod->beta + 1.));
				Estimate[g][j] += expected[s][g][j];
				SE_Estimate[g][j] += variance[s][g][j];
			}
		}
		
		
		
	}
	
	
	/*now begin to do the estimates, cycling through samples...*/
	
	for(g=0;g<mixmod->G;g++){
		for(j=0;j<mixmod->ncat[variable];j++){
			Estimate[g][j] /= n_sample;
			SE_Estimate[g][j] /= n_sample;
		}
	}
	
	for(s=0;s<n_sample;s++){
		for(g=0;g<mixmod->G;g++){
			for(j=0;j<mixmod->ncat[variable];j++)
					expected[s][g][j] -= Estimate[g][j];		
		}
	}
	
	double *se;
	se = calloc(mixmod->ncat[variable],sizeof(double));
	
	for(g=0;g<mixmod->G;g++){
	
		for(j=0;j<mixmod->ncat[variable];j++){
			se[j] = 0.;
		}
		
		for(s=0;s<n_sample;s++){
			
			for(j=0;j<mixmod->ncat[variable];j++)
				se[j] += pow(expected[s][g][j],2.);
		
		}
		
		for(j=0;j<mixmod->ncat[variable];j++){
			se[j] /= (n_sample - 1.);
			SE_Estimate[g][j] += se[j];
			SE_Estimate[g][j] = sqrt(SE_Estimate[g][j]);
		}
		
	}	
	
	
	for(i=0;i<n_sample;i++){
		for(j=0;j<mixmod->G;j++){
			free(expected[i][j]);
			free(variance[i][j]);
		}
		free(expected[i]);
		free(variance[i]);
	}

	free(expected);
	free(variance);
	
	for(i=0;i<mixmod->G;i++)
		free(counts[i]);
	
	free(items);
	free(counts);
	free(se);

	return;

}


void BLCA_compute_post_hoc_parameter_estimates_for_class_probabilities( struct mix_mod *mixmod, struct results *input, int n_sample, int n_groups, int variable, double *Estimate, double *SE_Estimate)
{

	int i,j,k,s,g,z,*counts;
	double **expected;
	
	counts = calloc(mixmod->G,sizeof(int));
	expected = calloc(n_sample,sizeof(double *));
	for(i=0;i<n_sample;i++){
		expected[i] = calloc(mixmod->G,sizeof(double));
	}
	

	for(s=0;s<n_sample;s++){
	
		for(g=0;g<mixmod->G;g++){
			counts[g] = 0;
		}
	
		for(i=0;i<mixmod->n;i++){
		
			z = input->memberships[s][i];
			
			counts[z] += 1;
		
		}
		
		for(g=0;g<mixmod->G;g++){
			expected[s][g] = (counts[g]+mixmod->alpha)/(mixmod->n + mixmod->G*mixmod->alpha);
			Estimate[g] += expected[s][g];
			SE_Estimate[g] += (counts[g] + mixmod->alpha)*(mixmod->n - counts[g] + (mixmod->G - 1.)*mixmod->alpha)/(pow(mixmod->n + mixmod->G*mixmod->alpha,2.)*(mixmod->n + mixmod->G*mixmod->alpha + 1.));
		}
		
	}
	
	for(g=0;g<mixmod->G;g++){
		Estimate[g] /= n_sample;
		SE_Estimate[g] /= n_sample;
	}
	
	for(s=0;s<n_sample;s++){
		for(g=0;g<mixmod->G;g++){
			expected[s][g] -= Estimate[g];
			SE_Estimate[g] += pow(expected[s][g],2.)/n_sample;
		}
	}
	
	
	for(g=0;g<mixmod->G;g++){
		SE_Estimate[g] = sqrt(SE_Estimate[g]);
	}


	free(counts);
	free(expected);
	
	return;

}


