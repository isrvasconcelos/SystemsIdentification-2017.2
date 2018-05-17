suppressMessages(library(control))

###########################################
# Delayer block using 1st order Pade approximation 
DelayerBlock <- function(theta) {

	numDelayer <- c((-theta/2), 1)
	denDelayer <- c((theta/2), 1)

	return(tf(numDelayer, denDelayer))
}

###########################################
# Generic blocks
ApproximationBlock1stOrder <- function(Amax, tau) {

	numPB <- c(0, Amax)
	denPB <- c(tau, 1)

	return(tf(numPB, denPB))
}

ApproximationBlock2ndOrder <- function(nP=1, nPs=0, nPs2=0, dP=1, dPs=0, dPs2=0) {

	numPlant <- c(nPs2, nPs, nP)
	denPlant <- c(dPs2, dPs, dP)

	return(tf(numPlant, denPlant))
}

