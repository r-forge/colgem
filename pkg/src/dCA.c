/*
 * @author Erik M Volz
 * @date October 14 2014
 * 
 * Equations for lineage through time and distribution of node heights
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static double *parms; 

//macros
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) < (b)) ? (b) : (a))

#define m (int)parms[0]
#define treeT parms[1]
#define hres parms[2]

#define pend 3
#define F(i, k,l) parms[(int)(pend + hres * (l*m +k ) + i)]
#define fend (int)(pend + hres * pow(m,2))
#define G(i, k,l) parms[fend + (int)(hres * (l*m +k ) + i)]
#define gend (int)(pend + 2 * hres * pow(m,2))
#define Y(i, k)  parms[gend + (int)(hres*k+i)]
#define yend (int)(gend+hres*m)
#define notSampledYet(i,k) parms[yend + (int)(hres*k+i) ]

/* initializers  */
void initfunc_dCA(void (* odeparms)(int *, double *))
{
	int Nparms;
	DL_FUNC get_deSolve_gparms;
	SEXP gparms;
	get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	gparms = get_deSolve_gparms();
	Nparms = LENGTH(gparms);
	parms = REAL(gparms);
}



/* derivatives */
#define A(k) y[k]
#define dA(k) ydot[k]

//~ R version: 
	//~ dA <- function(h, A, parms, ...)
	//~ {
		//~ nsy <- not.sampled.yet(h) 
		//~ with(get.fgy(h), 
		//~ { 
			//~ A_Y 	<- (A-nsy) / .Y
			//~ csFpG 	<- colSums( .F + .G )
			//~ list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
		//~ })
	//~ }
void dCA( int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
	int i =  (int)min( (int)( hres * (*t) / treeT ), hres-1);
	int k,l,z,w;
	
	double a[m]; //normalized nlft 
	double sumA = 0.; 
	for (k = 0; k < m; k++) sumA += A(k);
	for (k = 0; k < m; k++) { 
		dA(k) = 0.;
		if (Y(i,k) > 0) {
			a[k] =  (A(k)-notSampledYet(i,k)) / Y(i,k);
			//~ a[k] =  max(0, min(1, r *  A(k)/Y(i,k)));
			//~ a[k] = max( min(r * A(k)/Y(i,k), 1), 0) ;
		} else{
			a[k] = 1.; //
		} 
	}
	
	//dA
	for (k = 0; k < m; k++){
		for (l = 0; l < m; l++){
			if (k==l){
				dA(k) -= a[l] * (F(i,l,k)) * a[k];
			} else{
				dA(k) += ((1 - a[k]) * F(i,k,l) + G(i,k,l)) * a[l] ;
				dA(k) -= (F(i,l,k) + G(i,l,k)) * a[k];
			}
		}
	}
}
