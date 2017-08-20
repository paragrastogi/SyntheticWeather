library(forecast)
# Read the de-meaned values. There is no header line
TSin = read.csv(file.path('../GEN/GEN_IWEC_dm.csv'), header = FALSE)
periodicity = 24
ARp <- 4 # AR
MAq <- 4 # MA
SARp <- 1 # SAR
SMAq <- 1 # SMA
n <- dim(TSin)[1] # Length of incoming time series record
s <- dim(TSin)[2] # Number of incoming time series (not currently used)
npaths <- 50 # Number of paths to be simulated
# Pre-allocate out-going time series. In a future implementation, the number of outgoing time series would be equal to "s".
tsOut1 <- array(0, c(n, npaths))
tsOut2 <- array(0, c(n, npaths))

try(mod1 <- Arima(TSin[1], order=c(1, 0, 0), seasonal=list(order=c(0, 0, 0), period=periodicity)))
bic1 <- mod1$bic
try(mod2 <- Arima(TSin[2], order=c(1, 0, 0), seasonal=list(order=c(0, 0, 0), period=periodicity)))
bic2 <- mod2$bic
for (p in 0: ARp) {
	for (q in 0: MAq) {
		for (pp in 0: SARp) {
			for (qq in 0: SMAq) {
				try(mod1 <- Arima(TSin[1], order=c(p, 0, q), seasonal=list(order=c(pp, 0, qq), period=periodicity, method = "CSS-ML")))
				try(mod2 <- Arima(TSin[2], order=c(p, 0, q), seasonal=list(order=c(pp, 0, qq), period=periodicity, method = "CSS-ML")))
				if (!(exists("mod1"))) {
				try(mod1 <- Arima(TSin[1], order=c(p, 0, q), seasonal=list(order=c(pp, 0, qq), period=periodicity, method = "CSS")))
				}
				if (!(exists("mod2"))) {
				try(mod2 <- Arima(TSin[2], order=c(p, 0, q), seasonal=list(order=c(pp, 0, qq), period=periodicity, method = "CSS")))
				}
				if (mod1$bic<bic1) {
				bic1 <- mod1$bic
				 ARmodel1 <- mod1
				}
				if (mod2$bic<bic2) {
				bic2 <- mod2$bic
				 ARmodel2 <- mod2
				}}}}}
for (s in 1:ncol(tsOut1)) { 
	CustomInn1 <- sample(ARmodel1$residuals, n, replace=T)
	CustomInn2 <- sample(ARmodel2$residuals, n, replace=T)
	tsOut1[,s] <- simulate.Arima(ARmodel1, nsim = n, innov=CustomInn1)
	tsOut2[,s] <- simulate.Arima(ARmodel2, nsim = n, innov=CustomInn2)
}
write.table(tsOut1, file="rout1.csv", row.names=FALSE, col.names=FALSE, sep=",")
write(ARmodel1$arma, file="rmodel1.csv", ncolumns=7)
write(ARmodel1$coef, file="rmodel1.csv", ncolumns=sum(ARmodel1$arma), append=TRUE)
write(c(ARmodel1$AIC,ARmodel1$AICc,ARmodel1$BIC,ARmodel1$loglik), file="rmodel1.csv", ncolumns=3, append=TRUE)
write.table(tsOut2, file="rout2.csv", row.names=FALSE, col.names=FALSE, sep=",")
write(ARmodel2$arma, file="rmodel2.csv", ncolumns=7)
write(ARmodel2$coef, file="rmodel2.csv", ncolumns=sum(ARmodel2$arma), append=TRUE)
write(c(ARmodel2$AIC,ARmodel2$AICc,ARmodel2$BIC,ARmodel2$loglik), file="rmodel2.csv", ncolumns=3, append=TRUE)
