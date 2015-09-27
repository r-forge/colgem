
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
