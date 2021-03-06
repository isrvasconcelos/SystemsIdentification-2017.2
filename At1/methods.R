source("systemblocks.R")
suppressMessages(library(Metrics))

######################################################################################
## Plot auxiliar

exportMSE_Barplot <- function( data , imgName , title ) {

	ylim <- c( 0 , max(data) )

	setEPS(width=10, height=7)
		postscript(paste(imgName, ".eps", sep='') )

		barplot( data , beside = FALSE , 
			 names.arg=c("Zieg-Nich.","Haggl.","Smith","Sund-Krish.", "Mollen."),
			 col=c("cyan4"), ylab="Mean Square Error", xlab="Method",
			 ylim=ylim,
			 main=title
			)

	dev.off()
}

######################################################################################
## DONE
calibration <- function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	y <- response_filtered
	yraw <- response

	ylim <- c( min(y, yraw) , max(y,yraw) )

	plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)
	points(t, y, type='l', col='red')

	points(x=t, y=tgline, type='l', col='dodgerblue', lty=5)
	abline(v=InflectionPoint, h=y[InflectionIndex], col='lightblue', lty=4) # Reference for inflection point

#	abline(v=t1, col='dodgerblue', lty=3)
#	abline(v=t3, col='dodgerblue', lty=3)

	legend("bottomright",  legend=c("Original System", 
					"Filtered System", 
					"Tangent Line on Inflection"), 
				lty=c(1,1,5), col=c("tan","red", "dodgerblue"), bty="n")

}

######################################################################################
# 1st Order System Identification
## DONE
zieglerNichols <- function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	###########################################
	# Variable fixes
	y <- response_filtered
	yraw <- response

	###########################################
	# Timestamps for evaluations
	t1Index <- which( abs(tgline) < .1 )[1] # Getting the index which tangent line is close to zero

	if(precision == TRUE) # Increase precision for this specific dataset
		t1Index <- which( round( abs(tgline) ,3) < .001 )[1]

	t1 <- t[t1Index]
	if(!length(t1))
		t1=min(t)

	peakIndex <- which( y == max(y) )[1]
	t3Index <- which(  tgline > y[peakIndex] )[1]
	t3 <- t[t3Index]

	t2Index <- which(  tgline > y[peakIndex]*0.632 )[1]
	t2 <- t[t2Index]

	###########################################
	# Parameters for estimation
	Amax <- max(y)

	tau <- t3-t1
	theta <- t1

	###########################################
	# System blocks
	PredictBlock <- series( ApproximationBlock1stOrder(Amax, tau), DelayerBlock(theta) )
	SystemResponse <- step(PredictBlock, t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System","Tangent Line on Inflection","Filtered System","Original System")
	if(!abs(sum(y-yraw))) # Check if response == response_filtered
		legendNames <- c("Estimated System", "Tangent Line on Inflection", "Original System" )

	setEPS()
		imgName <- paste("ziegler-nichols/zn-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)
		points(t, y, type='l', col='red') # Filtered data plot

		points(x=t, y=tgline, type='l', col='dodgerblue', lty=5)
	#	abline(v=InflectionPoint, h=y[InflectionIndex], col='lightblue', lty=4) # Reference for inflection point

		abline(v=t1, col='dodgerblue', lty=3)
		abline(v=t3, col='dodgerblue', lty=3)

		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,5,1,1), col=c("navy","dodgerblue","red","tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))
}


## DONE
hagglund <- function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	###########################################
	# Variable fixes
	y <- response_filtered
	yraw <- response

	###########################################
	# Timestamps for evaluations
	t1Index <- which( abs(tgline) < .1 )[1] # Getting the index which tangent line is close to zero

	if(precision == TRUE) # Increase precision for this specific dataset
		t1Index <- which( round( abs(tgline) ,3) < .001 )[1]

	t1 <- t[t1Index]
	if(!length(t1))
		t1=min(t)

	peakIndex <- which( y == max(y) )[1]

	t2Index <- which(  tgline > y[peakIndex]*0.632 )[1]
	t3 <- t2 <- t[t2Index]

	###########################################
	# Parameters for estimation
	Amax <- max(y)

	tau <- t3-t1
	theta <- t1

	###########################################
	# System blocks
	PredictBlock <- series( ApproximationBlock1stOrder(Amax, tau), DelayerBlock(theta) )
	SystemResponse <- step(PredictBlock, t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System","Tangent Line on Inflection","Filtered System","Original System")
	if(!abs(sum(y-yraw))) # Check if response == response_filtered
		legendNames <- c("Estimated System", "Tangent Line on Inflection", "Original System" )

	setEPS()
		imgName <- paste("hagglund/hag-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)
		points(t, y, type='l', col='red') # Filtered data plot

		points(x=t, y=tgline, type='l', col='dodgerblue', lty=5)
	#	abline(v=InflectionPoint, h=y[InflectionIndex], col='lightblue', lty=4) # Reference for inflection point

		abline(v=t1, col='dodgerblue', lty=3)
		abline(v=t3, col='dodgerblue', lty=3)

		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,5,1,1), col=c("navy","dodgerblue","red","tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))
}

## DONE
sunKris <- function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	###########################################
	# Variable fixes: This procedure evaluates datasets without filtering
	yraw <- response
	y <- yraw 

	###########################################
	# Timestamps for evaluations

	peakIndex <- which( y == max(y) )[1]
	A1 <- y[peakIndex]*0.353
	t1Index <-which(  y > A1 )[1]
	t1 <- t[t1Index]

	A2 <- y[peakIndex]*0.853
	t2Index <-which(  y > A2 )[1]
	t2 <- t[t2Index]

	###########################################
	# Parameters for estimation
	Amax <- max(y)

	tau <- 0.67*(t2-t1)
	theta <- 1.3*t1 - 0.29*t2

	###########################################
	# System blocks
	PredictBlock <- series( ApproximationBlock1stOrder(Amax, tau), DelayerBlock(theta) )
	SystemResponse <- step(PredictBlock, t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System", "Filtered System" ,"Original System")
	if(!abs(sum(y-yraw))) # Check if response == response_filtered
		legendNames <- c("Estimated System", "Original System")

	setEPS()
		imgName <- paste("sun-kris/sk-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)

		abline(v=t1, col='dodgerblue', lty=3)
		abline(v=t2, col='dodgerblue', lty=3)

		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,1), col=c("navy","tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))
}


# DONE
smith1st <-  function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	## DONE
	###########################################
	# Variable fixes: This procedure evaluates datasets without filtering
	yraw <- response
	y <- yraw 

	###########################################
	# Timestamps for evaluations

	peakIndex <- which( y == max(y) )[1]
	A1 <- y[peakIndex]*0.283
	t1Index <-which(  y > A1 )[1]
	t1 <- t[t1Index]

	A2 <- y[peakIndex]*0.632
	t2Index <-which(  y > A2 )[1]
	t2 <- t[t2Index]

	###########################################
	# Parameters for estimation
	Amax <- max(y)

	tau <- 1.5*(t2-t1)
	theta <-t2 - tau

	###########################################
	# System blocks
	PredictBlock <- series( ApproximationBlock1stOrder(Amax, tau), DelayerBlock(theta) )
	SystemResponse <- step(PredictBlock, t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System", "Original System")

	setEPS()
		imgName <- paste("smith/sm1-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)
	#	points(t, y, type='l', col='tan')

	#	points(x=t, y=tgline, type='l', col='dodgerblue', lty=5)
	#	abline(v=InflectionPoint, h=y[InflectionIndex], col='lightblue', lty=4) # Reference for inflection point

		abline(v=t1, col='dodgerblue', lty=3)
		abline(v=t2, col='dodgerblue', lty=3)

		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,1), col=c("navy", "tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))
}

######################################################################################
# TODO: PENDING
# 2nd Order System Identification
smith2nd <-  function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	###########################################
	y <- response_filtered
	yraw <- response

	###########################################
	# Timestamps for evaluations

	peakIndex <- which( y == max(y) )[1]
	A1 <- y[peakIndex]*0.2
	t1Index <-which( y > A1 )[15]
	t1 <- t[t1Index]

	A2 <- y[peakIndex]*0.6
	t2Index <-which( y > A2 )[1]
	t2 <- t[t2Index]

	###########################################
	# Parameters for estimation
	Amax <- max(y)
	zeta <- t1/t2
	tau <- 1 # TODO: Discover how to evaluate tau
	theta <- t2-tau
	SmithBlock <- 0

	if(zeta < 1)
		SmithBlock <- ApproximationBlock2ndOrder( nP=1, dPs2=(tau^2) , dPs=(2*zeta*tau), dP=1 )

	else if(zeta >=1) {
		tau1 <- tau*zeta - tau*((-1 + zeta^2)^1/2)
		tau1 <- tau*zeta + tau*((-1 + zeta^2)^1/2)
		SmithBlock <- series( ApproximationBlock1stOrder(1, tau1), ApproximationBlock1stOrder(1, tau2) )
	}

	###########################################
	# System blocks
	Delayer <- DelayerBlock(theta) 
	Gain <- ApproximationBlock1stOrder(K, 1)

	PredictBlock <- series( Gain , SmithBlock , Delayer )
	SystemResponse <- step( PredictBlock , t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System", "Original System")

	setEPS()
		imgName <- paste("smith/sm2-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)
		points(t, y, type='l', col='red')
		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,1,1), col=c("navy","red", "tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))

}


# DONE
mollenkamp <- function(interval, response, response_filtered, tgline, fileName, precision=FALSE) {

	###########################################
	# Variable fixes: This procedure evaluates datasets without filtering
	yraw <- response
	y <- yraw

	###########################################
	# Timestamps for evaluations

	peakIndex <- which( y == max(y) )[1]
	A1 <- y[peakIndex]*0.15
	t1Index <-which(  y > A1 )[1]
	t1 <- t[t1Index]

	A2 <- y[peakIndex]*0.45
	t2Index <-which(  y > A2 )[1]
	t2 <- t[t2Index]

	A3 <- y[peakIndex]*0.75
	t3Index <-which(  y > A3 )[1]
	t3 <- t[t3Index]

	if(precision==TRUE) {
		t1Index <- which( y/max(y) > 0.15 )[20] # Manual calibration
		t2Index <- which( y/max(y) > 0.45 )[1]
		t3Index <- which( y/max(y) > 0.75 )[1]

		t1 <- t[t1Index]
		t2 <- t[t2Index]
		t3 <- t[t3Index]
	}

	###########################################
	# Parameters evaluation
	MollenkampBlock <- 0
	f2 <- 0
	K <- y[peakIndex]

	x <- (t2-t1) / (t3-t1)
	zeta <-  ( 0.085 - 5.547*((0.475-x)^2) ) / ( x - 0.356 ) 

	if(zeta < 1) { # Underdamped
		f2 <- (0.708)*((2.811)^zeta)

		wn <- f2/(t3-t1)
		f3 <- (0.922)*((1.66)^zeta)
		theta <- t2 - (f3/wn)

		MollenkampBlock <- ApproximationBlock2ndOrder( nP=(wn^2) , dPs2=1 , dPs=(2*zeta*wn) , dP=(wn^2) )
	}

	else if(zeta >=1) { # Overdamped
		f2 <- 2.6*zeta - 0.6

		wn <- f2/(t3-t1)
		f3 <- (0.922)*((1.66)^zeta)
		theta <- t2 - (f3/wn)

		tau1 <- ( zeta + (-1 + zeta^2)^(1/2) ) / wn
		tau2 <- ( zeta - (-1 + zeta^2)^(1/2) ) / wn

		MollenkampBlock <- series( ApproximationBlock1stOrder(1, tau1) , ApproximationBlock1stOrder(1, tau2) )
	}

	###########################################
	# System blocks
	
	Delayer <- DelayerBlock(theta) 
	Gain <- ApproximationBlock1stOrder(K, 1)

	PredictBlock <- series( Gain , MollenkampBlock , Delayer )
	SystemResponse <- step( PredictBlock , t=t) # Running Experiment

	###########################################
	# Plot the data
	ylim <- c( min(y, yraw, SystemResponse$y) , max(y,yraw, SystemResponse$y) )

	legendNames <- c("Estimated System", "Original System")

	setEPS()
		imgName <- paste("mollenkamp/mol-", gsub('.{4}$', '', fileName), ".eps", sep='')
		postscript(imgName)

		plot(t,yraw, type='l', col='tan', xlab="Time (s)", ylab="Amplitude", ylim=ylim)

		abline(v=t1, col='dodgerblue', lty=3)
		abline(v=t2, col='dodgerblue', lty=3)
		abline(v=t3, col='dodgerblue', lty=3)

		points(SystemResponse$t, SystemResponse$y, type='l', col='navy')

		legend("bottomright",  legend=legendNames, lty=c(1,1), col=c("navy","tan"), bty="n")

	dev.off()

	return(mse(SystemResponse$y , y))
}

