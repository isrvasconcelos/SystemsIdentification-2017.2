suppressMessages(library(control))
suppressMessages(library(Metrics))

###########################################
# Plant Settings
# 1 / (s+1)^5 = 1 / (s^5 + 5s^4 + 10s^3 + 10s^2 + 5s + 1)

	nP <- 1
	nPs <- 0
	nPs2 <- 0
	nPs3 <- 0
	dP <- 1
	dPs <- 4
	dPs2 <- 5
	dPs3 <- 2

	numPlant <- c(nPs3, nPs2, nPs, nP)
	denPlant <- c(dPs3, dPs2, dPs, dP)

	if (sum(abs(denPlant))==0)
		stop("Invalid parameters")

###########################################
# Delayer
# exp(-theta*s)
	theta = .5
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

	Plant <- tf(numPlant, denPlant)

	System <- feedback(Plant, 1/SetPoint)

	SystemResponse <- step(System, t=TimeInterval) # Running Experiment

	MSR <- SystemResponse$y[1,] # MSR: Magnitude of System Response - (Samples)

	###########################################
	# Response Analysis and Manipulations - Part II

	## First Critical Point : FCP
	FCPIndex <- min(which( diff(MSR)<0 )) # Index to get the First Local Peak
	FCP <- MSR[FCPIndex]
	FCPTimeStamp <- TimeInterval[FCPIndex]

	## Updating sample set to take the next critical point
	NextCritical <- MSR[FCPIndex:length(MSR)]
	NextCriticalTimeStamp <- TimeInterval[FCPIndex:length(MSR)]

	## Second Critical Point : SCP
	SCPIndex <- min(which( diff(NextCritical)>0 )) # Index to get the First Local Minimum
	SCP <- NextCritical[SCPIndex]
	SCPTimeStamp <- NextCriticalTimeStamp[SCPIndex]

	## Updating sample set to take the next critical point
	NextCritical <- NextCritical[SCPIndex:length(NextCritical)]
	NextCriticalTimeStamp <- NextCriticalTimeStamp[SCPIndex:length(NextCriticalTimeStamp)]

	## Third Critical Point : SCP
	TCPIndex <- min(which(diff( NextCritical)<0 )) # Index to get the Second Local Peak
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

	EstimatedSystem <- tf(numES, denES)

	EstSysResponse <- step(EstimatedSystem, t=TimeInterval)

###########################################
# Plotting and Outputs

	cat(c("First Local Peak: ", FCP))
	cat(c("\nSecond Local Peak: ", TCP))
	cat(c("\nFirst Local Minimum: ", SCP))
	cat(c("\nMean Square Error: ", mse(SystemResponse$y ,EstSysResponse$y)))

	cat(paste(c("\nyInf: ", yInf, 
	"\nDt: ", Dt , "\nKm: ", Km, "\nZeta: ", Zeta, "\nTm: ", Tm, "\nTheta: ", Theta ,
	"\nTauEst: ", TauEst, "\nKEst: ", KEst, "\n\n")))

	print("Original System:")
	System

	print("Estimated System Below:")
	EstimatedSystem

	setEPS()
	postscript("Q2_Response.eps")

	plot(EstSysResponse$t, EstSysResponse$y, 
		xlab='Time (s)', ylab='Step Response', type='l', col='chartreuse4')

	points(SystemResponse$t, MSR, type='l', col='darkgoldenrod1')

	abline(h=0, col="red", lwd=3, lty=2)

	legend("bottomright", legend=c("Estimated System", "Original System"), 
		lty=1, col=c("chartreuse4","darkgoldenrod1"), bty="n")

dev.off()
