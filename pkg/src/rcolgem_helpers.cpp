// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma; 

/*
			tree$mstates[extantLines,] <- t( t(Q) %*% t(tree$mstates[extantLines,])  )
			tree$mstates[extantLines,] <- abs(tree$mstates[extantLines,]) / rowSums(abs(tree$mstates[extantLines,]))
			#recalculate A
			A <- colSums(tree$mstates[extantLines,])
*/
//~ TOOD arma version for speed comparison
//~ mat update_mstates(vec extantLines, mat mstates)
//~ {	
//~ } 
//[[Rcpp::export]]
NumericMatrix update_mstates( NumericVector extantLines, NumericMatrix mstates, NumericMatrix Q ){
	//TODO void work ?
	int u, k, l, z, w; 
	int m = mstates.ncol(); 
	NumericVector p_u; 
	for(int iu = 0; iu < extantLines.size(); iu++){
		u = extantLines(iu) ; 
		p_u = NumericVector(mstates(u-1,_)); //clone( mstates(u - 1, _)); 
		mstates(u-1, _) = NumericVector(m, 0.); 
		double sm = 0.; 
		for(k = 0;  k < m; k++){
			for (l = 0; l < m; l++){
				mstates(u -1 , k) += p_u(l) * Q(l,k);
			}
			sm += mstates(u-1, k);
		}
		mstates(u-1, _) = mstates(u-1, _)/ sm; 
	}
	return mstates; 
}

//[[Rcpp::export]]
mat update_mstates_arma(const uvec& extantLines, const mat& Q,  mat mstates)
{
	int u; 
	vec p_u ; 
	double sms; 
	int m = Q.n_rows;
	for (int i = 0; i < extantLines.size(); i++){
		u = extantLines(i); 
		p_u = mstates.col(u-1);
		mstates.col(u-1) = Q * p_u ; 
		//~ for (int k = 0; k < m; k++){
			//~ mstates(u-1, k) = dot( p_u, Q.col(k)); 
		//~ }
		sms = sum(mstates.col(u-1)); 
		mstates.col(u-1) /= sms; 
	}
	return mstates;
}

/*
 { #CPP 
					p_a     <- tree$lstates[alpha,]
					p_i_mat <- tree$mstates[extantLines,]
					A_mat   <- t( matrix( A, nrow=tree$m, ncol=length(extantLines) )  )
					p_a_mat <- t( matrix(p_a, nrow=tree$m, ncol=length(extantLines)) )
					rho_mat <- A_mat /  (A_mat - p_i_mat)
					rterms  <- p_a_mat / (A_mat - p_i_mat)
					lterms  <- rho_mat %*% p_a
					lterms <- matrix(lterms, nrow=length(lterms), ncol=tree$m) #
					new_p_i <- p_i_mat * (lterms - rterms)
					new_p_i <- new_p_i / rowSums(new_p_i)
					corrupted <- is.nan(rowSums(new_p_i))
					new_p_i[corrupted,] <- p_i_mat[corrupted,]
					tree$mstates[extantLines,] <- new_p_i
				}
*/

//[[Rcpp::export]]
mat finite_size_correction(const vec& p_a, const vec& A, const uvec& extantLines, mat mstates)
{
	// NOTE mstates m X n 
	int u; 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	for (int iu = 0; iu < extantLines.size(); iu++){
		u = extantLines(iu); 
		p_u = mstates.col(u-1); 
		rterm = p_a / ( A - p_u); 
		rho = A / (A - p_u ); 
		lterm = dot( rho, p_a); //
		p_u = p_u % (lterm - rterm) ;
		p_u = p_u / sum(p_u ) ; 
		mstates.col(u - 1) = p_u; 
	}
	return mstates; 
}
//~ List finite_size_correction(NumericVector p_a, NumericMatrix mstates, NumericVector A)
