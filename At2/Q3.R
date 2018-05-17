suppressMessages(library(control))
suppressMessages(library(Metrics))

###########################################
# Plant Settings
# 3 / (2s + 1)

	nP <- 3
	nPs <- 0
	dP <- 1
	dPs <- 2

	numPlant <- c(nPs, nP)
	denPlant <- c(dPs, dP)

	if (sum(abs(denPlant))==0)
		stop("Invalid parameters")


###########################################
# Delayer
# exp(-theta*s)
	theta = 0.1
	numDelayer <- c((-theta/2), 1)
	denDelayer <- c((theta/2), 1)

###########################################
# State-Space System Model Defaults for Step Reponse
	SetPoint <- 1 # Reference value
	Kc <- 1 # Controller Gain

###########################################
# Response Analysis and Manipulations Defaults
	StepTime <- 35
	Precision <- .01
	TimeInterval <- seq(0, StepTime, Precision)

###########################################
# Assemblying system boxes

	Plant <- series( tf(numPlant, denPlant), tf(numDelayer, denDelayer) )

	System <- feedback(Plant, 1/SetPoint)

	SystemResponse <- step(System, t=TimeInterval) # Running Experiment

	MSR <- SystemResponse$y[1,] # MSR: Magnitude of System Response - (Samples)

	###########################################
	# Response Analysis and Manipulations - Part II

	## First Critical Point : FCP
	FCPIndex <- which( MSR == max(MSR) )[1] # Index to get the First Local Peak
	FCP <- MSR[FCPIndex]
	FCPTimeStamp <- TimeInterval[FCPIndex]

	## Updating sample set to take the next critical point
	NextCritical <- MSR[FCPIndex:length(MSR)]
	NextCriticalTimeStamp <- TimeInterval[FCPIndex:length(MSR)]

	## Second Critical Point : SCP
	SCPIndex <- which(NextCritical == min(NextCritical))[1] # Index to get the First Local Minimum
	SCP <- NextCritical[SCPIndex]
	SCPTimeStamp <- NextCriticalTimeStamp[SCPIndex]

	## Updating sample set to take the next critical point
	NextCritical <- NextCritical[SCPIndex:length(NextCritical)]
	NextCriticalTimeStamp <- NextCriticalTimeStamp[SCPIndex:length(NextCriticalTimeStamp)]

	## Third Critical Point : SCP
	TCPIndex <- which(NextCritical == max(NextCritical))[1] # Index to get the Third Local Peak
	TCP <- NextCritical[TCPIndex]
	TCPTimeStamp <- NextCriticalTimeStamp[TCPIndex]


###########################################
# Estimating
	# Yp1: First Critical Point - First Local Peak
	# Ym1: Second Critical Point - First Local Minimum
	# Yp2: Third Critical Point - Second Local Peak
	
	yInf <- (FCP*TCP - SCP^2)/(FCP + TCP - 2*SCP);
	Dt <- SCPTimeStamp-FCPTimeStamp
	Km <- yInf/(Kc*(SetPoint-yInf))
	Kf <- Kc*Km

	Zeta <- (-log((yInf - SCP)/(FCP - yInf)))/ ( sqrt(pi^2+log((yInf - SCP)/(FCP - yInf))^2) )
	
	Tm <- (Dt/pi)*( Zeta*sqrt(Kf+1) + sqrt(Kf + (Kf+1)*Zeta^2) ) * sqrt( (1-Zeta^2)*(Kf+1) )

	Theta <- (2*Dt*sqrt( (1-Zeta^2)*(Kf+1) )) / ( pi*( Zeta*sqrt(Kf+1) + sqrt(Kf+(Kf+1)*Zeta^2) ))

	TauEst <- sqrt( (Theta*Tm) / (2*(Kf+1)) )
	KEst <- Kf/(Kf+1)


###########################################
# Building Estimated System
	nESP <- KEst
	nESPs <- -(KEst*Theta)/2
	nESPs2 <- 0

	dESP <- 1
	dESPs <- 2*Zeta*TauEst
	dESPs2 <- TauEst

	numES <- c(nESPs2, nESPs, nESP)
	denES <- c(dESPs2, dESPs, dESP)

#	EstimatedSystem <- tf(numES, denES)

#	EstSysResponse <- step(EstimatedSystem, t=TimeInterval)

###########################################
# Plotting and Outputs

	cat(c("First Local Peak: ", FCP))
	cat(c("\nSecond Local Peak: ", TCP))
	cat(c("\nFirst Local Minimum: ", SCP))
#	cat(c("\nMean Square Error: ", mse(SystemResponse$y ,EstSysResponse$y)))

	cat(paste(c("\nyInf: ", yInf, 
	"\nDt: ", Dt , "\nKm: ", Km, "\nZeta: ", Zeta, "\nTm: ", Tm, "\nTheta: ", Theta ,
	"\nTauEst: ", TauEst, "\nKEst: ", KEst, "\n\n")))

	print("Original System:")
	System

setEPS()
	postscript("Q3_Response.eps")
	plot(SystemResponse$t, SystemResponse$y, type='l', col='darkgoldenrod1', lty=1, bty="n")
	legend("bottomright", legend=c("Original System"), lty=1, col=c("darkgoldenrod1"), bty="n")

dev.off()
