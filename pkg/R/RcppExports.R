
dAL <- function(t, y, parms) {
	.Call('sourceCpp_0_dAL',  PACKAGE='rcolgem', t, y, parms)
}

simulateTreeCpp2 <- function(times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , As
		 , sortedCoHeights 
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection) {
	.Call('sourceCpp_2_simulateTreeCpp2',  PACKAGE='rcolgem', times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , As
		 , sortedCoHeights 
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection)
}

updateWCpp <- function( W
  , psi_a 
  , utips
  , vtips
  , utipsW
  , vtipsW )
{
	.Call( 'sourceCpp_0_updateWCpp', PACKAGE='rcolgem'
	  , W
	  , psi_a 
	  , utips
	  , vtips
	  , utipsW
	  , vtipsW
	)
}


#~ sourceCpp_0_update_mstates_arma(SEXP extantLinesSEXP, SEXP QSEXP, SEXP mstatesSEXP) {
update_mstates_arma <- function( extantLines, Q, mstates)
{
	.Call( 'sourceCpp_0_update_mstates_arma', PACKAGE='rcolgem'
	  ,  extantLines, Q, mstates
	)
}
#~ sourceCpp_0_finite_size_correction(SEXP p_aSEXP, SEXP ASEXP, SEXP extantLinesSEXP, SEXP mstatesSEXP) {
finite_size_correction <- function( p_a, A, extantLines, mstates)
{
	.Call( 'sourceCpp_0_finite_size_correction', PACKAGE='rcolgem'
	  ,   p_a, A, extantLines, mstates
	)
}
