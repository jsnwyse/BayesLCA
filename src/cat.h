/*C headers: solve square assignment problem 
		(translated from Carpaneto & Toth's Fortran implementation)

	Author:	Jason Wyse, 
			School of Computer Science and Statistics,
			Trinity College Dublin,
			Dublin 2, Ireland.
			email: wyseja@tcd.ie
			
Last modified: Sat 16 May 2015 03:14:03 PM IST   */ 


/*assignment_problem_offset.h*/

#ifndef _ass_problem_offset_H_
#define _ass_problem_offset_H_

#include <stdlib.h>
#include "nrutil.h"
#include "R.h"

void assct(int n, int **A,int *c, int *T);

#endif
