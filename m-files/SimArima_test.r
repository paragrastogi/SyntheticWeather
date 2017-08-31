library(forecast)
CustomInn = read.csv(file.path("/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbInnForR_NYC_CPR_TMY-36.csv"))
IncomingX = read.csv(file.path("/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbForR_NYC_CPR_TMY-36.csv"))
periodicity = 0
ARp = 4
MAq = 4
SARp = 0
SMAq = 1
N = 8760
npaths = 50

try( PsiModel <- Arima(IncomingX, order=c(ARp,0,MAq), seasonal=list(order=c(SARp,0,SMAq), period=periodicity), method = 'CSS-ML', n.cond=0) )
try(testsim <- simulate(PsiModel,nsim=N,innov=CustomInn[,1]))
if (!(exists('PsiModel')) | !(exists('testsim'))) {
 try( PsiModel <- Arima(IncomingX, order=c(ARp,0,MAq), seasonal=list(order=c(SARp,0,SMAq), period=periodicity), method = 'CSS', n.cond=0) )
}
if (!(exists('testsim'))) {
  try( PsiModel <- Arima(IncomingX, order=c(ARp,0,MAq), seasonal=list(order=c(0,0,0), period=periodicity), method = 'CSS-ML', n.cond=0) )
}
tsOut <- array(0,c(N,npaths))
for (s in 1:ncol(tsOut)) { 
 tsOut[,s] <- simulate(PsiModel,nsim=N,innov=CustomInn[,s])
}
DataOut = data.frame(tsOut)
write.table(DataOut, file="/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbSimFromR_NYC_CPR_TMY-36.csv", col.names = F, row.names = F, sep = ",")
write.table(PsiModel$residuals, file="/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbResFromR_NYC_CPR_TMY-36.csv", col.names = F, row.names = F, sep = ",")
modtags = c("ar","ma","sar","sma","seasonal","diff","seasonaldiff")
write(modtags,file = "/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbModelFromR_NYC_CPR_TMY-36.csv",sep=",", ncolumns = length(modtags))
write(t(PsiModel$arma), file = "/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbModelFromR_NYC_CPR_TMY-36.csv",sep=",", ncolumns = length(PsiModel$arma), append = TRUE)
write(t(PsiModel$coef), file = "/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbModelFromR_NYC_CPR_TMY-36.csv",sep=",",ncolumns = length(PsiModel$coef), append = TRUE)
PsiModelTags = c("sigma2","loglik","aic","bic")
write(PsiModelTags,file = "/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbModelFromR_NYC_CPR_TMY-36.csv", sep=",",ncolumns = length(PsiModelTags), append=TRUE)
write(c(PsiModel$sigma2,PsiModel$loglik,PsiModel$aic,PsiModel$bic),file = "/home/rasto/Documents/allmycode/SyntheticWeather/m-data/NYC/NYC_36/RinoutN50/tdbModelFromR_NYC_CPR_TMY-36.csv",sep=",", ncolumns = length(PsiModelTags),append=TRUE)
