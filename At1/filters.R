fftinv <- function( x ) { fft( x, inverse=TRUE ) / length( x ) } # inversa normalizada

getFFTFreqs <- function(data,samplingFreq) {

	Nyq.Freq=samplingFreq/2;
	if ((length(data) %% 2) == 1) { # Nº ímpar de amostras

	FFTFreqs <- c(seq(0, Nyq.Freq, length.out=(length(data)+1)/2), 
		seq(-Nyq.Freq, 0, length.out=(length(data)-1)/2))
	}

	else { # Nº par
	FFTFreqs <- c(seq(0, Nyq.Freq, length.out=length(data)/2), 
		seq(-Nyq.Freq, 0, length.out=length(data)/2))
        }

    return (FFTFreqs)
}

plotFFT <- function(y, samplingFreq, bodePlot = FALSE, title="", ylim=c(0,max(y))) {
	FFTFreqs <- getFFTFreqs(y,samplingFreq)
	FFT <- y

	if(bodePlot==FALSE) {
		if(!is.complex(y[1])) # Aplicar FFT se o dado estiver bruto
			FFT <- fft(y)
		
		FFT <- Mod(FFT)
		FFT[1]=0 # Corrigindo p/ evitar picos no grafico

		if(missing(ylim))
			ylim=c(0,max(FFT))

	}

	FFTdata <- cbind(FFTFreqs, FFT)
	FFTdata2 <- FFTdata[1:nrow(FFTdata)/2,]
	plot(FFTdata2, t="l", pch=20, lwd=2, cex=0.8, main="", xlab="Frequency (Hz)", ylab="Power", ylim=ylim, sub=title)
}


lowpass <- function(fftData=FALSE,samplingFreq=FALSE,threshold=500,plot=FALSE,bodePlot=FALSE,times=1) {

	if((fftData==FALSE && bodePlot==FALSE) || samplingFreq==FALSE) # controle de parametros
		stop("Missing data.")

	if(fftData!=FALSE && !is.complex(fftData[1]))  # Aplicar FFT se o dado estiver bruto
		fftData = fft(fftData)

	if(is.logical(fftData)) # Define valor default p/ plotar diagrama de Bode
		fftData=1:50000

	FFTFreqs <- Mod(getFFTFreqs(fftData,samplingFreq))

	coefs <- 1/(1+(FFTFreqs/threshold)) # vetor com coeficientes multiplicativos para ganho

	if(times>1)
		for(i in 1:times)
			coefs <- coefs*(1/(1+(threshold/FFTFreqs)))

	if(bodePlot==TRUE) 
		plotFFT(coefs, bodePlot=TRUE)

	else if(plot==TRUE) {
		fftData[1]=0 # Corrigindo p/ evitar picos no grafico
		par(mfrow=c(2,1))
		plotFFT(fftData, ylim=c(0,max(Mod(fftData))), title="(No filter)")
		plotFFT(coefs*fftData, ylim=c(0,max(Mod(fftData))), title="(After low-pass filtering)")
	}
	
	else {
		data=coefs*fftData
		data=Re(fftinv(data))
		return(data)
	}
}


lowpass2ndOrder <- function(fftData=FALSE,samplingFreq=FALSE,threshold=500,plot=FALSE,bodePlot=FALSE) {

	if((fftData==FALSE && bodePlot==FALSE) || samplingFreq==FALSE) # controle de parametros
		stop("Missing data.")

	if(fftData!=FALSE && !is.complex(fftData[1]))  # Aplicar FFT se o dado estiver bruto
		fftData = fft(fftData)

	if(is.logical(fftData)) # Define valor default p/ plotar diagrama de Bode
		fftData=1:50000

	s <- FFTFreqs <- Mod(getFFTFreqs(fftData,samplingFreq))

	coefs <- threshold^2/(s^2 + threshold*s/1000 + threshold^2) # vetor com coeficientes multiplicativos para ganho

	if(bodePlot==TRUE) 
		plotFFT(coefs, bodePlot=TRUE)

	else if(plot==TRUE) {
		fftData[1]=0 # Corrigindo p/ evitar picos no grafico
		par(mfrow=c(2,1))
		plotFFT(fftData, ylim=c(0,max(Mod(fftData))), title="(No filter)")
		plotFFT(coefs*fftData, ylim=c(0,max(Mod(fftData))), title="(After high-pass filtering)")
	}
	
	else {
		data=coefs*fftData
		data=Re(fftinv(data))
		return(data)
	}
}


highpass <- function(fftData=FALSE,threshold=500,plot=FALSE,bodePlot=FALSE,times=1) {

	if(fftData==FALSE && bodePlot==FALSE) # controle de parametros
		stop("Missing data.")

	if(fftData!=FALSE && !is.complex(fftData[1]))  # Aplicar FFT se o dado estiver bruto
		fftData = fft(fftData)

	if(is.logical(fftData)) # Define valor default p/ plotar diagrama de Bode
		fftData=1:50000

	FFTFreqs <- Mod(getFFTFreqs(fftData,samplingFreq))

	coefs <- 1/(1+(threshold/FFTFreqs)) # vetor com coeficientes multiplicativos para ganho

	if(times>1) # aplicar o filtro repetidas vezes
		for(i in 1:times)
			coefs <- coefs*(1/(1+(threshold/FFTFreqs)))
	
	if(bodePlot==TRUE) 
		plotFFT(coefs, bodePlot=TRUE)

	else if(plot==TRUE) {
		fftData[1]=0 # Corrigindo p/ evitar picos no grafico
		par(mfrow=c(2,1))
		plotFFT(fftData, ylim=c(0,max(Mod(fftData))), title="(No filter)")
		plotFFT(coefs*fftData, ylim=c(0,max(Mod(fftData))), title="(After high-pass filtering)")
	}
	
	else {
		data=coefs*fftData
		data=Re(fftinv(data))
		return(data)
	}
}

highpass2ndOrder <- function(fftData=FALSE,threshold=500,plot=FALSE,bodePlot=FALSE) {

	if(fftData==FALSE && bodePlot==FALSE) # controle de parametros
		stop("Missing data.")

	if(fftData!=FALSE && !is.complex(fftData[1]))  # Aplicar FFT se o dado estiver bruto
		fftData = fft(fftData)

	if(is.logical(fftData)) # Define valor default p/ plotar diagrama de Bode
		fftData=1:50000

	s <- FFTFreqs <- Mod(getFFTFreqs(fftData,samplingFreq))

	coefs <- s^2/(s^2+threshold*s/1000+threshold^2) # vetor com coeficientes multiplicativos para ganho

	if(bodePlot==TRUE) 
		plotFFT(coefs, bodePlot=TRUE)

	else if(plot==TRUE) {
		fftData[1]=0 # Corrigindo p/ evitar picos no grafico
		par(mfrow=c(2,1))
		plotFFT(fftData, ylim=c(0,max(Mod(fftData))), title="(No filter)")
		plotFFT(coefs*fftData, ylim=c(0,max(Mod(fftData))), title="(After high-pass filtering)")
	}
	
	else {
		data=coefs*fftData
		data=Re(fftinv(data))
		return(data)
	}
}
