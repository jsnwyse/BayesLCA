/*C functions to do solve square assignment problem 
		(translated from Carpaneto & Toth's Fortran implementation)

	Author:	Jason Wyse, 
			School of Computer Science and Statistics,
			Trinity College Dublin,
			Dublin 2, Ireland.
			email: wyseja@tcd.ie
			
Last modified: Sat 16 May 2015 03:14:03 PM IST   */ 



/*
 *  cat.c
 *
 *  
 *  This is a C implementation of the algorithm of 
 *	Carpaneto and Toth for the assignment problem.
 *
 *	This is a direct translation of the Fortran code
 *	given in Carpaneto and Toth (1980).
 *
 *	The call to the function is assct(n,A,c,T)
 *	The meaning of the input parameters is 
 *
 *	n = number of rows and columns of the cost matrix, with
 *		the current dimensions and the maximum value of n is 200
 *	A[i][j] = element in row i and column j of the cost matrix
 *	(at the end of the computation the elements of A are changed)
 *
 *	The meaning of the output parameters is
 *	
 *	c[j] = row assigned to column j (j=1,...,n)
 *	T = cost of the optimal assignment
 *
 *	ALL PARAMETERS ARE INTEGER
 *
 *	The meaning of the local varaibles is
 *	
 *	A[i][j] = element of the cost matrix if A[i][j] is positive,
 *			column of the unassigned zero following in row i (i = 1,...,n)
 *			the unassigned zero of column j (j=1,...,n) if A[i][j] is not
 *			positive
 *	A[i][n+1] = column of the first unassigned zero of row i (i=1,...,n)
 *	ch[i] = column of the next unexplored and unassigned zero of row i (i=1,...,n)
 *	lc[j] = label of column j (j=1,...,n)
 *	lr[i] = label of row i (i = 1,...,n)
 *	lz[i] = column of the last unassigned zero of row i (i=1,...,n)
 *	nz[i] = column of next unassigned zero of row i (i=1,...,n)
 *	rh[i] = unexplored row following the unexplored row i (i=1,...,n)
 *	rh[n+1] = first unexplored row
 *	slc[k] = k-th element contained in the set of the labelled columns
 *	slr[k] = k-th element contained in the set of the lebelled rows
 *	u[i] = unassigned row following the unassigned row i (i=1,...,n)
 *	u[n+1] = first unassigned row
 *
 *	The vectors c, ch, lc, lr, lz, nz, slc, slr must be dimensioned at least at n,
 *	the vectors rh, u at least at n+1,
 *	the matrix A at least at n times n+1
 *
 */

#include "cat.h" 


void assct(int n, int **A,int *c, int *T){

int *ch,*lc,*lr,*lz,*nz,*slc,*slr,*rh,*u;
int h,q,r,s,i,j,l,k,lj,m,lm,kslc,kslr,nl,nm,np1;

np1 = n+1;

/*c = ivector(1,n);*/
ch = ivector(1,n);
lc = ivector(1,n);
lr = ivector(1,n);
lz = ivector(1,n);
nz = ivector(1,n);
slc = ivector(1,n);
slr = ivector(1,n);
rh = ivector(1,np1);
u = ivector(1,np1);

/*INIALIZATION*/
for(i=1;i<np1;i++){
c[i]=0;
lz[i]=0;
nz[i]=0;
u[i]=0;
}

u[np1]=0;
*T=0;

/*REDUCTION OF THE INITIAL COST MATRIX*/
for(j=1;j<np1;j++){

s = A[1][j];

for(l=2;l<np1;l++){
if(A[l][j]<s) s = A[l][j];
}

*T=*T+s;

for(i=1;i<np1;i++) A[i][j] = A[i][j]-s;
}

for(i=1;i<np1;i++){
q = A[i][1];
for(l=2;l<np1;l++){
if(A[i][l]<q) q = A[i][l];
}

*T=*T+q;
l = np1;

for(j=1;j<np1;j++){
A[i][j] = A[i][j]-q;
if(A[i][j]!=0) continue;
A[i][l] = -j;
l=j;
}
}

/*CHOICE OF THE INITIAL SOLUTION*/

k = np1;

for(i=1;i<np1;i++){
lj = np1;
j = -A[i][np1];
goto label8;

label8: 

if(c[j]==0) goto label13;
lj = j;
j = -A[i][j];
if(j!=0) goto label8;
lj = np1;
j = -A[i][np1]; 
goto label9;

label9:

r = c[j];
lm = lz[r];
m = nz[r]; 
goto label10;

label10:

if(m==0) goto label11;
if(c[m]==0) goto label12;
lm = m;
m = -A[r][m];
goto label10;

label11:
lj = j;
j = -A[i][j];
if(j!=0) goto label9;
u[k] = i;
k = i;
continue;

label12:
nz[r] = -A[r][m];
lz[r] = j;
A[r][lm] = -j;
A[r][j] = A[r][m];
A[r][m] = 0;
c[m] = r;
goto label13;

label13:

c[j]=i;
A[i][lj] = A[i][j];
nz[i] = -A[i][j];
lz[i] = lj;
A[i][j] = 0;

}
goto label15;

/*RESEARCH OF A NEW ASSIGNMENT*/
label15:

if(u[np1]==0){ 
free_ivector(ch,1,n);
free_ivector(lc,1,n);
free_ivector(lr,1,n);
free_ivector(lz,1,n);
free_ivector(nz,1,n);
free_ivector(slc,1,n);
free_ivector(slr,1,n);
free_ivector(rh,1,np1);
free_ivector(u,1,np1);
return;
}

for(i=1;i<np1;i++){
ch[i] = 0;
lc[i] = 0;
lr[i] = 0;
rh[i] = 0;
}

rh[np1] = -1;
kslc = 0;
kslr = 1;
r = u[np1];
lr[r] = -1;
slr[1] = r;
if(A[r][np1]==0) goto label22;
goto label17;

label17:
l = -A[r][np1];
if(A[r][l]==0) goto label18;
if(rh[r]!=0) goto label18;
rh[r] = rh[np1];
ch[r] = -A[r][l];
rh[np1] = r;
goto label18;

label18:
if(lc[l]==0) goto label20;
if(rh[r]==0) goto label21;
goto label19;

label19:
l = ch[r];
ch[r] = -A[r][l];
if(A[r][l]!=0) goto label18;
rh[np1] = rh[r];
rh[r] = 0;
goto label18;

label20:
lc[l] = r;
if(c[l]==0) goto label36;
kslc = kslc+1;
slc[kslc] = l;
r = c[l];
lr[r] = l;
kslr = kslr+1;
slr[kslr] = r;
if(A[r][np1]!=0) goto label17;
goto label21;

label21:
if(rh[np1] > 0) goto label35;
goto label22;

/*REDUCTION OF THE CURRENT COST MATRIX*/

label22:

h = 1000000000;
for(j=1;j<np1;j++){
if(lc[j]!=0) continue;
for(k=1;k<kslr+1;k++){
i = slr[k];
if(A[i][j] < h) h = A[i][j];
}
}

*T = *T+h;

for(j=1;j<np1;j++){
if(lc[j]!=0) continue;
for(k=1;k<kslr+1;k++){
i=slr[k];
A[i][j] = A[i][j]-h;
if(A[i][j]!=0) continue;
if(rh[i]!=0) goto label25;
rh[i] = rh[np1];
ch[i] = j;
rh[np1] = i;
goto label25;

label25:

l=n+1;
goto label26;

label26:
nl = -A[i][l];
if(nl==0) goto label27;
l = nl;
goto label26;

label27:
A[i][l] = -j;
continue;
}
}

if(kslc==0) goto label35;

for(i=1;i<np1;i++){
if(lr[i]!=0) continue;
for(k=1;k<kslc+1;k++){
j = slc[k];
if(A[i][j] > 0) goto label32;
l = np1;
goto label30;

label30:
nl = -A[i][l];
if(nl==j) goto label31;
l = nl;
goto label30;

label31:

A[i][l] = A[i][j];
A[i][j] = h;
continue;

label32:
A[i][j] = A[i][j]+h;
continue;
}
}
goto label35;

label35:

r = rh[np1];
goto label19;

/*ASSIGNMENT OF A NEW ROW*/
label36:

c[l] = r;
m = n+1;
goto label37;

label37:
nm = -A[r][m];
if(nm==l) goto label38;
m = nm;
goto label37;

label38:

A[r][m] = A[r][l];
A[r][l] = 0;
if(lr[r] < 0) goto label39;
l = lr[r];
A[r][l] = A[r][np1];
A[r][np1] = -l;
r = lc[l];
goto label36;

label39:
u[n+1] = u[r];
u[r] = 0;
goto label15;


/*free_ivector(c,1,n);*/
free_ivector(ch,1,n);
free_ivector(lc,1,n);
free_ivector(lr,1,n);
free_ivector(lz,1,n);
free_ivector(nz,1,n);
free_ivector(slc,1,n);
free_ivector(slr,1,n);
free_ivector(rh,1,np1);
free_ivector(u,1,np1);

return;
}






















