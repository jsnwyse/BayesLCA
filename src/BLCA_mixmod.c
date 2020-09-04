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
 
 Last modification of this code: Wed 29 July 2020  */


#include "BLCA_mixmod.h"

struct mix_mod *BLCA_allocate_mixmod(int n, int d, int G_max, int G, double *prior_hparams,int *ncat, 
                                     int COLLAPSED, int EM_fit, int EM_MAP, int VB, int BOOT )
  /*this function allocates and returns a pointer to a mixmod structure...*/
{
  
  int i,j,k;
  FILE *fp;
  
  struct mix_mod *mm = (struct mix_mod *)malloc(sizeof(struct mix_mod));
  
  mm->maxgroups = G_max;
  mm->G = G;
  mm->n = n;
  mm->d = d;
  mm->collapsed = COLLAPSED;
  mm->EM_fit = EM_fit;
  mm->EM_MAP = EM_MAP;
  mm->VB =  VB;
  mm->BOOT = BOOT;
  
  if( FALSE )
  {
    
    mm->Y = (int **)calloc(d,sizeof(int *));
    for( i=0; i<d; i++ )
      mm->Y[i] = (int *)calloc(n,sizeof(int));
    
    mm->Yobs = (int **)calloc(n,sizeof(int *));
    for( i=0; i<n; i++ )
      mm->Yobs[i] = (int *)calloc(d,sizeof(int));
  }
  
  mm->z = (int *)calloc( n,sizeof(int) );
  
  if( !mm->collapsed ) 
  {
    mm->s = (double **)calloc( n, sizeof(double *) );
    for( i=0; i<n; i++ )
      mm->s[i] = (double *)calloc( G, sizeof(double) );
  }
  
  mm->ncat = (int *)calloc( d, sizeof(int) );
  for( i=0; i<d; i++ ) mm->ncat[i] = ncat[i];
  
  
  mm->varindicator = (int *)calloc( d, sizeof(int) );
  
  /*allocate this memory-- only initialize what is needed for initial conditions*/
  mm->components = (struct component **)malloc( sizeof(struct component *) * G_max );
  
  for( i=0; i<G_max; i++ )
    mm->components[i] = (struct component *)malloc(sizeof(struct component));
  
  
  mm->undiscriminating = (struct component *)malloc(sizeof(struct component));
  BLCA_allocate_component( mm->undiscriminating, mm );
  mm->undiscriminating->n_g = n; /*there is always n elements in here as these variables
   do not define a clustering...*/
  
  //priors
  
  mm->alpha_prior =  (double *)calloc( G_max, sizeof(double) ) ;
  mm->beta_prior = (double ***)calloc( G_max, sizeof(double**) );
  mm->lg_beta_sum = (double **)calloc( G_max, sizeof(double *) );
  mm->lg_sum_beta = (double **)calloc( G_max, sizeof(double *) );
  for( k=0; k<G_max; k++ ) 
  {
    mm->beta_prior[k] = (double **)calloc( d, sizeof(double*) );
    mm->lg_beta_sum[k] = (double *)calloc( d, sizeof(double) );
    mm->lg_sum_beta[k] = (double *)calloc( d, sizeof(double) );
    for( i=0; i<d; i++ )
      mm->beta_prior[k][i] = (double *)calloc( ncat[i], sizeof(double) );
  }
  
  /*allocate whereis*/
  
  mm->whereis = (int *)calloc( G_max, sizeof(int) );
  for( i=0; i<G_max; i++ )
    mm->whereis[i] = -1;
  
  /*allocate for the initial number of components*/
  for( i=0; i<G; i++ )
  {
    BLCA_allocate_component( mm->components[i], mm );
    mm->components[i]->in_use = TRUE;
    mm->components[i]->n_g = 0;
  }
  
  for( i=G; i<G_max; i++ )
  {
    BLCA_allocate_component( mm->components[i], mm );
    mm->components[i]->in_use = FALSE;
    mm->components[i]->n_g = 0;
  }
  
  
  /*there will be two hparameters in the default model:
   alpha: dirichlet prior on weights symmetric
   beta: dirichlet prior for within component membership probs
   */
  
  mm->alpha = prior_hparams[0];
  mm->beta = prior_hparams[1];
  
  mm->log_prior_G = (double *)calloc( G_max + 1 , sizeof(double) );
  
  mm->log_like = -DBL_MAX;
  
  mm->log_prior = -DBL_MAX;
  
  /*assign space for the lookup table for a values*/
  
  //mixmod->table_a = calloc(datasize,sizeof(double));
  
  if( !COLLAPSED ) mm->weights = (double *)calloc( G, sizeof(double) );
  
  if( VB ) 
  {
    mm->alpha_ud = (double *)calloc( G, sizeof(double));
    mm->di_alpha_ud = (double *)calloc( G, sizeof(double) );
  }
  
  if( BOOT )
  {
    mm->boot_idx = (int *)calloc( n, sizeof(int) );
  }
  
  return(mm);
  
}

struct mix_mod *BLCA_clone_mixmod( struct mix_mod *mm )
{
  
  int k;
  double *hparams = (double *)calloc( 2, sizeof(double) );
  hparams[0] = mm->alpha;
  hparams[1] = mm->beta;
  
  struct mix_mod *mm_clone; 
  
  mm_clone = BLCA_allocate_mixmod( mm->n, mm->d, mm->G, mm->G, hparams, mm->ncat, mm->collapsed, mm->EM_fit, mm->EM_MAP, mm->VB, mm->BOOT );
  
  for( k=0; k<mm->G; k++ )
  {
    BLCA_copy_component( mm->components[k], mm_clone->components[k] , mm );
    mm_clone->weights[k] = mm->weights[k]; 
  }
  
  free( hparams );
  
  return( mm_clone );
}


void BLCA_free_mixmod(struct mix_mod *mm)
  /*frees the memory used by mixmod object*/
{
  int n = mm->n,d = mm->d, G = mm->G, i, j, k;
  
  /*free up components*/
  for( k=0; k<mm->maxgroups; k++ )
    {
    BLCA_free_component(mm->components[k], mm);
    free(mm->components[k]);
  }
  free(mm->components);
  
  BLCA_free_component( mm->undiscriminating, mm );
  free(mm->undiscriminating);
  
  // free priors
  free( mm->alpha_prior );
  for( k=0; k<mm->maxgroups; k++ )
  {
    free(mm->lg_beta_sum[k]);
    free(mm->lg_sum_beta[k]);
    for( i=0; i<d; i++ ) free(mm->beta_prior[k][i]);
    free( mm->beta_prior[k] );
  }
  free(mm->beta_prior );
  free(mm->lg_beta_sum);
  free(mm->lg_sum_beta);
  
  /*free wehreis*/
  free(mm->whereis);
  
  /*free data*/
  if( FALSE )
  {
    for( i=0; i<d; i++ )
      free(mm->Y[i]);
    
    free(mm->Y);
    
    for( i=0; i<n; i++ )
      free(mm->Yobs[i]);
    
    free(mm->Yobs);
  }
  /*free others*/
  free(mm->z);
  
  if( !mm->collapsed )
  {
    for( i=0; i<n; i++ ) 
      free(mm->s[i]);
    
    free(mm->s);
  }
  
  free(mm->ncat);
  
  free(mm->varindicator);
  
  free(mm->log_prior_G);
  
  /*free table*/
  //free(mixmod->table_a);
  
  if( !mm->collapsed ) free(mm->weights);
  if( mm->VB ) 
  {
    free(mm->alpha_ud);
    free(mm->di_alpha_ud);
  }
  
  if( mm->BOOT )
    free( mm->boot_idx );
  
  free(mm);
  
  return;
}




