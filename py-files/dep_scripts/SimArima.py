import subprocess  # This will be used to call R from inside Python
# import numpy as np
import pandas as pd
import sys

__author__ = 'Parag Rastogi'


def callrcmd(pathtsfile, arp, maq, sarp, smaq, npaths):

    # Write an R script using the function SimArima.py
    """
    :type pathtsfile: str
    :param arp: int
    :param maq: int
    :param sarp: int
    :param smaq: int
    :param npaths: int
    """

    # # Size of the incoming de-meaned series, to which the
    # # time series model will be fit.
    # n = demeaned.shape[0]  # Should be 8760
    # # Number of series being brought in.
    # s = demeaned.shape[1]  # Should be more than 1

    # Write the R script
    writetor(pathtsfile, arp, maq, sarp, smaq, npaths)

    # Compose the R command
    rcmd = 'R CMD BATCH SimArima.r '  # + 'ExtraFileR.txt'

    if sys.version_info >= (3, 0):
        subprocess.run(rcmd, shell=True)  # Command for post-3.0
    else:
        subprocess.call(rcmd, shell=True)  # Command for pre-3.0

    # Read back the files that R wrote to disc
    bsout1 = pd.read_csv('rout1.csv', delimiter=',', skiprows=None,
                         header=None, names=None)
    bsout2 = pd.read_csv('rout2.csv', delimiter=',', skiprows=None,
                         header=None, names=None)
    bsout = dict(tdb=bsout1.values, rh=bsout2.values)

    # Read in the ARmodel file
    modlist = dict(tdb=[line.rstrip() for line in open('rmodel1.csv')],
                   rh=[line.rstrip() for line in open('rmodel2.csv')])

    return [bsout, modlist]


# Write an R script to call R
def writetor(pathtsfile, arp, maq, sarp, smaq, npaths):

    """

    :param pathtsfile: basestring
    :param arp: int
    :param maq: int
    :param sarp: int
    :param smaq: int
    :param npaths: int
    """

    f = open('SimArima.r', 'w')

    # f.write('sim_arima <- function() {\n')

    # Load forecast library in R #!/usr/bin/Rscript\n
    f.write('library(forecast)\n')

    # Add path coming from the main script. This file should contain
    # the de-meaned time series.
    f.write('# Read the de-meaned values. There is no header line\n' +
            "TSin = read.csv(file.path(\'{0}\'), header = FALSE)\n".format(
                str.replace(pathtsfile, '\\', '/')))

    # Write all the constants to the script. Initialise the tsOut
    # matrix as well.
    f.write('periodicity = 24\n' +
            'ARp <- ' + str(arp) + ' # AR\n' + 'MAq <- ' + str(maq) +
            ' # MA\n' + 'SARp <- ' + str(sarp) + ' # SAR\n' +
            'SMAq <- ' + str(smaq) + ' # SMA\n')
    f.write('n <- dim(TSin)[1] # Length of incoming time series record\n' +
            's <- dim(TSin)[2] # Number of incoming time series ' +
            '(not currently used)\n' +
            'npaths <- ' + str(npaths) +
            ' # Number of paths to be simulated\n' +
            '# Pre-allocate out-going time series. In a future ' +
            'implementation, the number of outgoing time series ' +
            'would be equal to "s".\n'
            'tsOut1 <- array(0, c(n, npaths))\n' +
            'tsOut2 <- array(0, c(n, npaths))\n\n')

    f.write('try(mod1 <- Arima(TSin[1], order=c(1, 0, 0), ' +
            'seasonal=list(order=c(0, 0, 0), period=periodicity)))\n' +
            'bic1 <- mod1$bic\n' +
            'try(mod2 <- Arima(TSin[2], order=c(1, 0, 0), ' +
            'seasonal=list(order=c(0, 0, 0), period=periodicity)))\n' +
            'bic2 <- mod2$bic\n')

    f.write('for (p in 0: ARp) {\n' +
            '\tfor (q in 0: MAq) {\n' +
            '\t\tfor (pp in 0: SARp) {\n' +
            '\t\t\tfor (qq in 0: SMAq) {\n' +
            '\t\t\t\ttry(mod1 <- Arima(TSin[1], order=c(p, 0, q), ' +
            'seasonal=list(order=c(pp, 0, qq), ' +
            'period=periodicity, method = "CSS-ML")))\n'
            '\t\t\t\ttry(mod2 <- Arima(TSin[2], order=c(p, 0, q), ' +
            'seasonal=list(order=c(pp, 0, qq), ' +
            'period=periodicity, method = "CSS-ML")))\n'
            '\t\t\t\tif (!(exists("mod1"))) {\n'
            '\t\t\t\ttry(mod1 <- Arima(TSin[1], order=c(p, 0, q), ' +
            'seasonal=list(order=c(pp, 0, qq), ' +
            'period=periodicity, method = "CSS")))\n\t\t\t\t}\n'
            '\t\t\t\tif (!(exists("mod2"))) {\n'
            '\t\t\t\ttry(mod2 <- Arima(TSin[2], order=c(p, 0, q), ' +
            'seasonal=list(order=c(pp, 0, qq), ' +
            'period=periodicity, method = "CSS")))\n\t\t\t\t}\n' +
            '\t\t\t\tif (mod1$bic<bic1) {\n' +
            '\t\t\t\tbic1 <- mod1$bic\n' +
            '\t\t\t\t ARmodel1 <- mod1\n\t\t\t\t}\n' +
            '\t\t\t\tif (mod2$bic<bic2) {\n' +
            '\t\t\t\tbic2 <- mod2$bic\n' +
            '\t\t\t\t ARmodel2 <- mod2\n\t\t\t\t}}}}}\n')

    f.write('for (s in 1:ncol(tsOut1)) { \n' +
            '\tCustomInn1 <- sample(ARmodel1$residuals, n, replace=T)\n' +
            '\tCustomInn2 <- sample(ARmodel2$residuals, n, replace=T)\n' +
            '\ttsOut1[,s] <- simulate.Arima(ARmodel1, nsim = n,' +
            'innov=CustomInn1)\n' +
            '\ttsOut2[,s] <- simulate.Arima(ARmodel2, nsim = n,' +
            'innov=CustomInn2)\n' +
            '}\n')
    f.write('write.table(tsOut1, file="rout1.csv", row.names=FALSE,' +
            'col.names=FALSE, sep=",")\n' +
            'write(ARmodel1$arma, file="rmodel1.csv", ncolumns=7)\n' +
            'write(ARmodel1$coef, file="rmodel1.csv",' +
            'ncolumns=sum(ARmodel1$arma), append=TRUE)\n' +
            'write(c(ARmodel1$AIC,ARmodel1$AICc,ARmodel1$BIC,' +
            'ARmodel1$loglik), file="rmodel1.csv", ' +
            'ncolumns=3, append=TRUE)\n')
    f.write('write.table(tsOut2, file="rout2.csv", row.names=FALSE,' +
            ' col.names=FALSE, sep=",")\n' +
            'write(ARmodel2$arma, file="rmodel2.csv", ncolumns=7)\n' +
            'write(ARmodel2$coef, file="rmodel2.csv",' +
            ' ncolumns=sum(ARmodel2$arma), append=TRUE)\n' +
            'write(c(ARmodel2$AIC,ARmodel2$AICc,ARmodel2$BIC,' +
            'ARmodel2$loglik), file="rmodel2.csv", ' +
            'ncolumns=3, append=TRUE)\n')
    f.close()
