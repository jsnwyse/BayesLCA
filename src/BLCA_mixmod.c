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

struct mix_mod *BLCA_allocate_mixmod(int datasize, int datadimension, int maxgroups, int initgroups,double *prior_hparams,int *ncat, int collapsed, int EM_fit, int EM_MAP, int VB )
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
	mixmod->EM_fit = EM_fit;
	mixmod->EM_MAP = EM_MAP;
	mixmod->VB =  VB;
	
	mixmod->Y = calloc(datadimension,sizeof(int *));
	for(i=0;i<datadimension;i++){
		mixmod->Y[i] = calloc(datasize,sizeof(int));
	}
	
	mixmod->Yobs = calloc(datasize,sizeof(int *));
	for(i=0;i<datasize;i++){
		mixmod->Yobs[i] = calloc(datadimension,sizeof(int *));
	}
	
	if( !EM_fit ) mixmod->z = calloc(datasize,sizeof(int));
	
	if( EM_fit || VB ) 
	{
		mixmod->s = calloc(datasize, sizeof(double *));
		for( i=0;i<datasize;i++){
			mixmod->s[i] = calloc(initgroups,sizeof(double));
		}
	}
	
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
	BLCA_allocate_component(mixmod->undiscriminating,mixmod);
	mixmod->undiscriminating->n_g = mixmod->n; /*there is always n elements in here as these variables
																	do not define a clustering...*/

	
	/*allocate whereis*/
	
	mixmod->whereis = calloc(maxgroups,sizeof(int));
	for(i=0;i<maxgroups;i++)
		mixmod->whereis[i] = -1;
	
	/*allocate for the initial number of components*/
	for(i=0;i<initgroups;i++){
		BLCA_allocate_component(mixmod->components[i],mixmod);
		mixmod->components[i]->in_use = TRUE;
		mixmod->components[i]->n_g = 0;
	}
	for(i=initgroups;i<maxgroups;i++){
		BLCA_allocate_component(mixmod->components[i],mixmod);
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
	
	mixmod->log_like = -DBL_MAX;
	
	/*assign space for the lookup table for a values*/
	
	//mixmod->table_a = calloc(datasize,sizeof(double));
	
	if( !mixmod->collapsed ) mixmod->weights = calloc(mixmod->G,sizeof(double));
	
	if( VB ) 
	{
		mixmod->alpha_ud = calloc( mixmod->G, sizeof(double));
		mixmod->di_alpha_ud = calloc( mixmod->G, sizeof(double) );
	}
	
	return(mixmod);

}

void BLCA_free_mixmod(struct mix_mod *mixmod)
/*frees the memory used by mixmod object*/
{
	int n = mixmod->n,d = mixmod->d, G = mixmod->G,i,j,k;
	
	/*free up components*/
	for(k=0;k<mixmod->maxgroups;k++){
		BLCA_free_component(mixmod->components[k],mixmod);
		free(mixmod->components[k]);
	}
	free(mixmod->components);
	
	BLCA_free_component(mixmod->undiscriminating,mixmod);
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
	if( !mixmod->EM_fit ) free(mixmod->z);
	
	if( mixmod->EM_fit || mixmod->VB )
	{
		for(i=0;i<n;i++) free(mixmod->s[i]);
		free(mixmod->s);
	}
	
	free(mixmod->ncat);
	
	free(mixmod->varindicator);
	
	free(mixmod->log_prior_G);
	
	/*free table*/
	//free(mixmod->table_a);
	
	if( !mixmod->collapsed ) free(mixmod->weights);
	if( mixmod->VB ) 
	{
		free(mixmod->alpha_ud);
		free(mixmod->di_alpha_ud);
	}
	
	free(mixmod);
	
	return;
}




