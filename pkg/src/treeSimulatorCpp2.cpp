/* 
 * rcolgem tree simulator 
 * 
 */


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(cpp11)]]



static const double A_EPS = 0.00;
static const double Y_EPS = 1.E-12;


using namespace std;
using namespace Rcpp; 
using namespace arma; 



void incorporateSamples( double h, int samplesAdded, mat lstates, mat mstates, NumericVector sortedSampleHeights, NumericMatrix sortedSampleStates, IntegerVector extant,  NumericVector heights, int m ){
	for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
		if (sortedSampleHeights[i]==h){
			extant[i] = 1; 
			heights[i]  =h;
			for (int k = 0; k < m; k++) {
				lstates(i,k) = sortedSampleStates(i,k);
				mstates(i,k) = sortedSampleStates(i,k);
			}
		}
		if (sortedSampleHeights[i] > h) break;
	}
}

//[[Rcpp::export]]
List simulateTreeCpp2(const NumericVector times, const List Fs, const List Gs, const List Ys
  , const List As
  , NumericVector sortedCoHeights // computed at R level
  , const NumericVector sortedSampleHeights
  , const NumericMatrix sortedSampleStates
  , double maxSampleTime 
  , const int m
  , bool finiteSizeCorrection
)
{
	RNGScope scope;
	
	int n = sortedSampleHeights.size();
	int Nnode = sortedCoHeights.size(); //n-1; 
	// note times decreasing
	double tfin = times[times.size()-1];
	double maxHeight = maxSampleTime -tfin;
	const double deltat  = times[0] - times[1]; 
	
	IntegerVector extant( n + Nnode, 0);
	IntegerVector extantRange(n+Nnode);
	for  (int i = 0; i < n+Nnode; i++) extantRange[i] = i; 
	int numberExtant = 0;
	
	// note root does not have edge
	int ntrees = 1 + ((n -1) - Nnode); 
	NumericMatrix edge(n+Nnode-ntrees, 2); 
	NumericVector edge_length(n+Nnode-ntrees, -1.0);
	NumericVector heights(n+Nnode, -1.0);
	
	mat lstates = zeros(n+Nnode, m);
	mat mstates = zeros(n+Nnode, m);
	
	double t, t0, t1;
	int internalNodeIndex = n; // counter for internal branches
	int internalNodesAdded = 0;
	int samplesAdded = 0; // counter for terminal branches
	int edgesAdded = 0;
	
	vec  A_Y;
	vec Y ;
	vec A;
	mat F = zeros(m,m);
	mat G = zeros(m,m);
	int it, ist, k, l, nco, inco , i, stateOf_a, stateOf_recip;
	int u, v, a, z, w;
	double r, rk, rl, cwk, cwl, twk, twl, h1;
	bool foundv, foundu;
	double rsmsi;
	
	double h;
		
	h = 0.; 
	
	for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
		if (sortedSampleHeights[i]==h){
			extant[i] = 1; 
			heights[i]  =h;
			for (int k = 0; k < m; k++) {
				lstates(i,k) = sortedSampleStates(i,k);
				mstates(i,k) = sortedSampleStates(i,k);
			}
			samplesAdded++;
			numberExtant++; 
		}
		if (sortedSampleHeights[i] > h) break;
	}
	
	const double deltah  = times[0] - times[1]; //note times decreasing
	
	NumericMatrix corates(m,m);
	mat R  = zeros(m, m );
	mat Q  = zeros(m, m );
	vec m_rs_R;
	
	// note times decreasing
	for(it = 0; it < times.size(); it++)
	{
		t =  times[it];
		t1 = t - deltat;
		h = maxSampleTime - t; 
		h1 = maxSampleTime - t1; 

		
		F = as<mat>(Fs[it]);
		G = as<mat>(Gs[it]);
		Y = as<vec>(Ys[it]);
		A = as<vec>(As[it]);
		// transition probabilities 
		R = trans( F + G ) ;
		for ( k = 0; k < Y.size(); k++)
		{
			R.col(k) = R.col(k) / (Y);
		}
		R.diag(0) =   zeros<colvec>(Y.size());
		m_rs_R = -sum( R, 1 );
		// diagonal correction for coalescent events
		A_Y = clamp( A / Y, 0., 1. ); 
		vec diag_xn = zeros<colvec>(m); 
		for (k = 0; k < m; k++){
			diag_xn[k] = ((double)dot(F.row(k) , A_Y))/ Y[k]; 
		}
		R.diag(0) = m_rs_R - diag_xn;
		// transition prob from row to col
		Q = expmat( deltat * R); 
		// renormalise
		for (k = 0; k < m; k++){
			Q.row(k) = Q.row(k) / sum(Q.row(k));
		}
		// </transition probs>
		
		
		//update line states
		for ( i = 0; i < (n + internalNodesAdded); i++)
		{
			if (0!=extant[i]){
				for ( k = 0; k < m; k++)
				{
					mstates(i,k) = std::max(0.0, (double)(dot(Q.col(k), mstates.row(i))) );
				}
				//renormalise
				rsmsi  = sum(mstates.row(i));
				for ( k =0; k < m; k++){
					mstates(i,k) /= rsmsi;
				}
			}
		}
		
		
		// co events
		double corate = 0.;
		for ( k = 0; k < m; k++) {
			//~ corates.row(k) =  A_Y[k] * (F.row(k) % A_Y)
			for ( l = 0; l < m; l++){
				if (k == l){
					corates(k,k) =  A_Y(k) * (std::max(1.0,(A[k]-1))/std::max(1.0,(Y[k]-1)))  * F(k,k);
				} else{
					corates(k,l) =  A_Y(k) * A_Y(l) * F(k,l);
				}
				corate += corates(k,l);
			}
		}
		//~ for (i = internalNodesAdded; i < sortedCoHeights.size(); i++){
		for (int ico = 0; ico < sortedCoHeights.size(); ico++){
			//~ if (sortedCoHeights[i] > h1) break;
			if (sortedCoHeights[ico] > h && sortedCoHeights[ico] <= h1){
				// find k and l
				r = R::runif( 0., corate); 
				double sumcorate = 0.;
				int donordeme = 0;
				int recipdeme = 0;
				bool demes_found = false;
				for (w = 0; w < m; w++){
					for (z = 0; z < m; z++){
						sumcorate += corates(w,z);
						if (sumcorate > r){
							donordeme = w;
							recipdeme = z;
							demes_found = true; 
							break;
						}
					}
					if (sumcorate > r){
						break;
					}
				}

				// sample u & v in proportion to mstates k & l
				foundu = false;
				foundv = false;
				double Ak = 0.;//sum(mstates(_,k));
				double Al = 0.;//sum(mstates(_,l));
				for (i = 0; i < (n+internalNodesAdded); i++){
					if (extant[i]!=0){
						Ak += mstates(i, donordeme);
						Al += mstates(i,recipdeme);
					}
				}
				rk = R::runif(0., Ak); 
				double sumAk =0.;
				double sumAl = 0.;
				for (i = 0; i < (n + internalNodesAdded); i++){
					if (extant[i]!=0){
						sumAk += mstates(i,donordeme);
						if (sumAk > rk){
							u  = i; 
							foundu = true;
							Al -= mstates(i,recipdeme);
							break;
						} 
					}
				}
				rl = R::runif( 0., Al); 
				for (i = 0; i < (n + internalNodesAdded); i++){
					if ((i != u) && (extant[i]!=0)){
						sumAl += mstates(i,recipdeme);
						if (sumAl > rl){
							v = i;
							foundv = true;
							break; 
						} 
					}
				}
				
				// do the co
				if (!foundu || !foundv){
					// just pick two lines at random; 
					IntegerVector uv =  Rcpp::RcppArmadillo::sample( extantRange , 2, false, as<NumericVector>(extant)) ;
					u = uv[0];
					v = uv[1];
				}
				
				a = n + internalNodesAdded;
				internalNodesAdded++;
				heights[a] = sortedCoHeights[ico]; 
				
				edge(edgesAdded, 0) = a;
				edge(edgesAdded, 1) = u;
				edge_length[edgesAdded] = heights[a]-heights[u];
				edgesAdded++;
				edge(edgesAdded, 0) = a;
				edge(edgesAdded, 1) = v;
				edge_length[edgesAdded] = heights[a]-heights[v];
				edgesAdded++;
				
				extant[u] = 0;
				extant[v] = 0;
				extant[a] = 1;
				//~ numberExtant--;
				
				//lstate, mstate of a;
				lstates(a,donordeme) = 1.0;
				mstates(a,donordeme) = 1.0;
				
				//update mstates of lines not in co 
				if (finiteSizeCorrection){
					for (int s = 0; s < internalNodesAdded-1; s++){
						if (extant[s]!=0){
							for (int w = 0; w <m; w++){
								if (w!=k ){
									if (mstates(s,k) > 0){
										mstates(s, w) = mstates(s,w) * (Y[k] / (Y[k] - mstates(s,k)) )  ;
									}
								}
								if (w!=l ){
									if (mstates(s,l) > 0){
										mstates(s, w) = mstates(s,w) * (Y[l] / (Y[l] - mstates(s,l)) )  ;
									}
								}
							}
							mstates(s,k) = mstates(s,k) * std::max(0., ((Y[k]-1.)/Y[k])) ;
							mstates(s,l) = mstates(s,l) * std::max(0., ((Y[l]-1.)/Y[l])) ;
							mstates.row(s) = mstates.row(s) / sum(mstates.row(s));
						}
					}
				}
			}
		}
		
		// find sample times in (t, t1)
		//~ incorporateSamples(   h,  samplesAdded,  lstates,  mstates,  sortedSampleHeights,  sortedSampleStates,  extant,   heights, m);
		// new samples
		for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
			if ((sortedSampleHeights[i]>h) && (sortedSampleHeights[i]<=h1)){
				extant[i] = 1; 
				heights[i]  = sortedSampleHeights[i];
				for (int k = 0; k < m; k++) {
					lstates(i,k) = sortedSampleStates(i,k);
					mstates(i,k) = sortedSampleStates(i,k);
				}
				samplesAdded++;
				numberExtant++;
			}
			if (sortedSampleHeights[i] > h1) break;
		}
		
		
	}
	
	List ret;
	ret["edge"] = edge;
	ret["edge.length"] = edge_length;
	ret["n"] = n;
	ret["Nnode"] = internalNodesAdded;//Nnode
	ret["lstates"] = lstates;
	ret["mstates"] = mstates;
	ret["heights"] = heights;
	return(ret);
}
