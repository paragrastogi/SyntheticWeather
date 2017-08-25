% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function SimArimaR(InnPath, X1path, OutPathSim, ...
    OutPathRes, OutPathModel, periodicity, ARp, ...
    MAq, SARp, SMAq, N, npaths, RsubFolder)

p = inputParser;
p.FunctionName = 'SimArimaR';

addRequired(p,'InnPath',@ischar)
addRequired(p,'X1path',@ischar)
addRequired(p,'OutPathSim',@ischar)
addRequired(p,'OutPathRes',@ischar)
addRequired(p,'OutPathModel',@ischar)
addRequired(p,'periodicity',@isnumeric)
addRequired(p,'ARp',@isnumeric)
addRequired(p,'MAq',@isnumeric)
addRequired(p,'SARp',@isnumeric)
addRequired(p,'SMAq',@isnumeric)
addRequired(p,'N',@isnumeric)
addRequired(p,'npaths',@isnumeric)
addRequired(p,'RsubFolder',@ischar)


parse(p, InnPath, X1path, OutPathSim, OutPathRes, ...
    OutPathModel, periodicity, ARp, MAq, SARp, ...
	SMAq, N, npaths, RsubFolder);

InnPath = p.Results.InnPath;
X1path = p.Results.X1path;
OutPathSim = p.Results.OutPathSim;
OutPathRes = p.Results.OutPathRes;
OutPathModel = p.Results.OutPathModel;
periodicity = p.Results.periodicity;
ARp = p.Results.ARp; MAq = p.Results.MAq;
SARp = p.Results.SARp; SMAq = p.Results. SMAq;
N = p.Results.N; npaths = p.Results.npaths;
RsubFolder = p.Results.RsubFolder;

% This function takes as input: the original time series 
% (in a comma-delimited text file) the custom innovation 
% series (in a comma-delimited text file) path for the 
% output (simulated series) periodicity of the seasonal 
% model AR, MA, SAR, and SMA lags npaths -> number of 
% paths to simulate N -> number of values per path

% The outputs from the function are:
% CSV file with simulated time series. 
% CSV file with residuals from ARIMA fit in R.

RscriptName = 'SimArima.r';
RscriptPath = fullfile(RsubFolder, RscriptName);
fIDwrite = fopen(RscriptPath,'w');
fprintf(fIDwrite,'library(forecast)\r\n');

% R only accepts paths with forward slashes
InnPath = strrep(InnPath,'\','/');
X1path = strrep(X1path,'\','/'); 
OutPathSim = strrep(OutPathSim,'\','/');
OutPathRes = strrep(OutPathRes,'\','/');
OutPathModel = strrep(OutPathModel,'\','/');

% Read in the custom innovations from MATLAB
fprintf(fIDwrite, ['CustomInn = read.csv(file.', ...
	'path("%s"))\r\n'],InnPath);
fprintf(fIDwrite, ['IncomingX = read.csv(file.', ...
	'path("%s"))\r\n'], X1path);

% Initalise the ARIMA model
fprintf(fIDwrite,'periodicity = %d\r\n',periodicity);
fprintf(fIDwrite,'ARp = %d\r\n',ARp);
fprintf(fIDwrite,'MAq = %d\r\n',MAq);
fprintf(fIDwrite,'SARp = %d\r\n',SARp);
fprintf(fIDwrite,'SMAq = %d\r\n',SMAq);
fprintf(fIDwrite,['try( PsiModel <- ', ...
    'Arima(IncomingX, order=c(ARp,0,MAq), ', ...
	'seasonal=list(order=c(SARp,0,SMAq), ', ...
	'period=periodicity), method = ''CSS-ML'', ', ...
	'n.cond=0) )\r\n']);

% For some cases the optimisation algorithm does not
% converge. So, we make two changes to the Arima function:
% method and optim.control. The former is for choosing the
% fitting method and has three options: CSS-ML (default),
% CSS, and ML. CSS is Conditional Sum of Squares, while ML
% is Maximum Likelihood. The option optim.control passes
% parameters to the optimization routine. In this case, we
% change the algorithm
fprintf(fIDwrite,['if (!(exists(''PsiModel''))) ', ...
	'{\r\n try( PsiModel <- Arima(IncomingX, ', ...
	'order=c(ARp,0,MAq), seasonal=list(', ...
    'order=c(SARp,0,SMAq),period=periodicity), ', ...
	'method = ''CSS'', n.cond=0) )\r\n}\r\n']);


% Simulate the model npaths times, for N values in each path
fprintf(fIDwrite,'N = %d\r\n',N);
fprintf(fIDwrite,'npaths = %d\r\n',npaths);
fprintf(fIDwrite,['tsOut <- array(0,c(N,npaths))\r\n', ...
    'for (s in 1:ncol(tsOut)) { \r\n tsOut[,s] <- ', ...
    'simulate(PsiModel,nsim=N,innov=', ...
    'CustomInn[,s])\r\n}\r\n']);

% Put tsOut into a data frame for export
fprintf(fIDwrite,'DataOut = data.frame(tsOut)\r\n');

% Write the Time Series to a file
fprintf(fIDwrite,['write.table(DataOut, file="%s", ', ...
    'col.names = F, row.names = F, sep = ",")\r\n'], ...
	OutPathSim);

% Write the residuals to a file
fprintf(fIDwrite,['write.table(PsiModel$residuals, ', ...
	'file="%s", col.names = F, ', ...
    'row.names = F, sep = ",")\r\n'], ...
	OutPathRes);

% Initalise the matrices that will hold the coefficients
fprintf(fIDwrite,['modtags = c("ar","ma",', ...
    '"sar","sma","seasonal","diff","seasonaldiff")\r\n']);

% Write the model coefficients to a file
fprintf(fIDwrite,['write(modtags,file = ', ...
    '"%s",sep=",", ncolumns = length(modtags))', ...
	'\r\n'],OutPathModel);
fprintf(fIDwrite,['write(t(PsiModel$arma), ', ...
	'file = "%s",sep=",", ncolumns = ', ...
    'length(PsiModel$arma), append = TRUE)\r\n'], ...
	OutPathModel);
fprintf(fIDwrite,['write(t(PsiModel$coef), ', ...
    'file = "%s",sep=",",ncolumns = length(', ...
	'PsiModel$coef), append = TRUE)\r\n'], ...
	OutPathModel);

% Write the model meta-information to the Model file
fprintf(fIDwrite,['PsiModelTags = c("sigma2",', ...
	'"loglik","aic","bic")\r\n']);
fprintf(fIDwrite,['write(PsiModelTags,file = "%s", ',...
	'sep=",",ncolumns = length(PsiModelTags), ', ...
	'append=TRUE)\r\n'],OutPathModel);
fprintf(fIDwrite,['write(c(PsiModel$sigma2,',...
	'PsiModel$loglik,PsiModel$aic,PsiModel$bic),', ...
	'file = "%s",sep=",", ', ...
	'ncolumns = length(PsiModelTags),', ...
	'append=TRUE)\r\n'],OutPathModel);


fclose(fIDwrite);

% Run the R script written above
system(sprintf('R CMD BATCH %s %s', RscriptPath, ...
    fullfile(RsubFolder,'ExtrafileFromRbatch.txt')));

end