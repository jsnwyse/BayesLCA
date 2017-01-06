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
	
#ifndef __BLCA_UTILS_H__
#define __BLCA_UTILS_H__

#include "BLCA_required_libs.h"

double get_max(double *x,int len) ;

int get_imax(int *x,int len);

double get_min(double *x,int len) ;

int get_imin(int *x,int len) ;

int sample_discrete( double *weights, int len ) ;

int random_integer( int n );

/* a util to permute a vector  */
void random_ranshuffle( int *a, int n ) ;

#endif
