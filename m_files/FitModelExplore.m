function varargout = FitModelExplore (ts, varargin)

% This is bit of code to specify a whole bunch of models
% using combinations of the input lag values.
% Generally, one would select the most parsimonious model
% using AIC (Akaike Information Criteria) and BIC (Schwarz
% Bayesian Information Criteria).

p = inputParser;
p.FunctionName = 'FitModelExplore';

% The time series to which the models are fit.
addRequired(p,'ts',@(x) (isvector(x) || isnumeric(x)))

% Non-seasonal lags
addParameter(p,'ar',0,@isnumeric)
addParameter(p,'d',0,@isnumeric)
addParameter(p,'ma',0,@isnumeric)

% Conditional variance model lags
addParameter(p,'archm',0,@isnumeric)
addParameter(p,'archr',0,@isnumeric)

% The variance model
addParameter(p,'varmdl','constant',@ischar)

% Seasonal lags
addParameter(p,'sar',0,@isnumeric)
addParameter(p,'sma',0,@isnumeric)

% seasonality parameter
addParameter(p,'seasonality',0,@isnumeric)

parse(p, ts, varargin{:});

ar = p.Results.ar;
d = p.Results.d;
ma = p.Results.ma;
archm = p.Results.archm;
archr = p.Results.archr;
VarMdlType = p.Results.varmdl;
sar = p.Results.sar;
sma = p.Results.sma;
seasonality = p.Results.seasonality;

% Temporarily turn off rank-deficient warning, since it
% shows up quite a few times during the course of this
warnid = 'MATLAB:rankDeficientMatrix'; 
warnStruct = warning('off',warnid);

% if ~any(sma==0)
% 	sma = [0,sma];
% end
% 
% if ~any(sar==0)
% 	sar = [0,sar];
% end

if strcmpi(VarMdlType,'arch')
    fprintf('Using an ARCH model for variance \r\n')
else
    fprintf('Using constant variance \r\n')
    if exist('archm','var')==1
        fprintf(['Setting archm to Zero, it cannot be', ...
            ' used with constant variance \r\n'])
        archm = 0;
    end
    fprintf('Using constant variance \r\n')
    if exist('archr','var')==1
        fprintf(['Setting archr to Zero, it cannot be', ...
			' used with constant variance \r\n'])
        archr = 0;
    end
end


iter = 0;
for Sea = seasonality
    for P = ar
        for D = d
            for Q = ma
                if strcmpi(VarMdlType,'arch')
                    for M = archm
                        for R = archr
                            iter = iter + 1;
                        end
                    end
                else
                    iter = iter + 1;
                end
            end
        end
    end
end

% Add an iteration for the seasonal AR and MA components, if
% they were passed to the function

if all(sar == 0)
    fprintf('No sar lags specified. \r\n')
else
    iter = iter+1;
    if ~isscalar(sar) && any(sar==0)
        fprintf(['SAR lags are not scalar but one of', ...
			'the elements is zero. ', ...
            'Ignoring that element. \r\n'])
        sar = sar(sar~=0);
    end
end
if all(sma == 0)
    fprintf('No SMA lags specified. \r\n')
else
    iter = iter+1;
    if ~isscalar(sma) && any(sma==0)
        fprintf(['SMA lags are not scalar but one of', ...
			' the elements is zero. ', ...
            'Ignoring that element. \r\n'])
        sma = sma(sma~=0);
    end
end

OptWarnThrown = false(iter,1);
aic = NaN*ones(iter,1);
bic = NaN*ones(iter,1);
LogL = NaN*ones(iter,1);
saveout(iter).MeansMdl = NaN;
saveout(iter).EstMdl = NaN;
saveout(iter).EstParamCov = NaN;
saveout(iter).LogL = NaN;
saveout(iter).info = NaN;
saveout(iter).numParams = NaN;
saveout(iter).params = NaN;
saveout(iter).MdlNo = NaN;
saveout(iter).Const = NaN;
Count = 1;

try
	if seasonality ~= 0
		SeaVector = [0,seasonality];
	else
		SeaVector = seasonality;
	end
	
for Sea = SeaVector
for P = ar
for D = d
for Q = ma
	if P+D+Q == 0
		aic(Count) = NaN; bic(Count) = NaN;
		LogL(Count) = NaN;
		saveout(Count).MeansMdl = NaN;
		saveout(Count).EstMdl = NaN;
		saveout(Count).EstParamCov = NaN;
		saveout(Count).LogL = LogL(Count);
		saveout(Count).info = NaN;
		saveout(Count).numParams = NaN;
		saveout(Count).params = NaN.*ones(1,5);
		saveout(Count).MdlNo = Count;
		saveout(Count).Const = ...
			NaN;
		Count = Count + 1;
		fprintf(['ARIMA model cannot be ',...
			'specified with ALL zero lags \r\n'])
		continue
	end

	if all(sma==0) && all(sar==0)		
		% Intialise the model without Seasonal lags
		fprintf(['ARIMA model specified ', ...
			'with zero lags for sar and sma.', ...
			' Ignoring those two in the ', ...
			'model specification.', ...
			'\r\n'])
		MeansMdl = arima('ARLags',1:P, 'D', D, ...
			'MALags',1:Q, 'seasonality',Sea);
	elseif all(sma==0) && ~all(sar==0)
		% Initialise the model without sma lags
		MeansMdl = arima('ARLags',1:P, 'D', D, ...
			'MALags',1:Q, ...
			'SARLags',sar, 'seasonality',Sea);
	elseif ~all(sma==0) && all(sar==0)
		% Initialise the model without sar lags
		MeansMdl = arima('ARLags',1:P, 'D', D, ...
			'MALags',1:Q, ...
			'SMALags', sma, ...
			'seasonality',Sea);
	else
		% Initialise the model with both sar and sma lags
		MeansMdl = arima('ARLags',1:P, 'D', D, ...
			'MALags',1:Q, 'SARLags',sar, 'SMALags',sma, ...
			'seasonality',Sea);
	end

if strcmpi(VarMdlType,'arch')

	for M = archm
	for R = archr

	try
	varmdl = garch(M,R);

	catch err
	fprintf('%s\n',err.message);
	fprintf('GARCH model cannot be specified ')
	fprintf('with parameters M = %d and ', M)
	fprintf('R = %d \r\n', R)
	aic(Count) = NaN; bic(Count) = NaN;
	LogL(Count) = NaN;
	saveout(Count).MeansMdl = NaN;
	saveout(Count).EstMdl = NaN;
	saveout(Count).EstParamCov = NaN;
	saveout(Count).LogL = LogL(Count);
	saveout(Count).info = NaN;
	saveout(Count).numParams = NaN;
	saveout(Count).params = NaN.*ones(1,5);
	saveout(Count).MdlNo = Count;
	saveout(Count).Const = ...
		NaN;
	Count = Count + 1;
	% following lines: stack
	for e=1:length(err.stack)
		fprintf('%s at %i\n', ...
			err.stack(e).name, ...
			err.stack(e).line);
	end

	continue
	
	end

	MeansMdl.Variance = varmdl;

	try
	[EstMdl_temp, ...
		EstParamCov_temp, ...
		LogL(Count),info] = ...
		estimate(MeansMdl, ...
		ts, 'print', false);
	catch err
	fprintf('%s\n',err.message);
	fprintf(['At parameters ', ...
		'p = %d, d = %d, ', ...
		'q = %d, seasonality ',...
		'= %d \r\n'], ...
		P, D, Q, Sea)
	% following lines: stack
	for e=1:length(err.stack)
		fprintf('%s at %i\n', ...
			err.stack(e).name, ...
			err.stack(e).line);
	end

	aic(Count) = NaN; bic(Count) = NaN;
	LogL(Count) = NaN;
	saveout(Count).MeansMdl = NaN;
	saveout(Count).EstMdl = NaN;
	saveout(Count).EstParamCov = NaN;
	saveout(Count).LogL = LogL(Count);
	saveout(Count).info = NaN;
	saveout(Count).numParams = NaN;
	saveout(Count).params = NaN.*ones(1,5);
	saveout(Count).MdlNo = Count;
	saveout(Count).Const = ...
		NaN;
	Count = Count + 1;
	continue

	end

	[~, warnidLoop] = lastwarn;

	if strcmp(warnidLoop,warnid)
		OptWarnThrown(Count) = true;
	end
	clear warnidLoop

	numParam_temp = sum(any(EstParamCov_temp));
	numObs =  length(ts);

	[aic(Count), bic(Count)] = ...
		aicbic(LogL(Count), ...
		numParam_temp, numObs);

	saveout(Count).MeansMdl = MeansMdl;
	saveout(Count).EstMdl = EstMdl_temp;
	saveout(Count).EstParamCov = EstParamCov_temp;
	saveout(Count).LogL = LogL(Count);
	saveout(Count).info = info;
	saveout(Count).numParams = numParam_temp;
	saveout(Count).params = ...
		([length(EstMdl_temp.AR), ...
		EstMdl_temp.D, ...
		length(EstMdl_temp.MA), ...
		length(EstMdl_temp.Variance.ARCH), ...
		length(EstMdl_temp.Variance.GARCH), ...
		length(EstMdl_temp.SAR), ...
		length(EstMdl_temp.SMA), ...
		Sea]);
	saveout(Count).MdlNo = Count;
	saveout(Count).Const = ...
		EstMdl_temp.Constant;
	Count = Count + 1;
	end
	end

else

	try
	[EstMdl_temp, EstParamCov_temp, ...
	LogL(Count),info] = estimate(...
	MeansMdl, ts, 'print', false);

	catch err
	fprintf('%s\n',err.message);
	fprintf(['At parameters ', ...
	'p = %d, d = %d, ', ...
	'q = %d, seasonality ',...
	'= %d \r\n'], ...
	P, D, Q, Sea)
	% following lines: stack
	for e=1:length(err.stack)
	fprintf('%s at %i\n', ...
		err.stack(e).name, ...
		err.stack(e).line);
	end

	aic(Count) = NaN; bic(Count) = NaN;
	LogL(Count) = NaN;
	saveout(Count).MeansMdl = NaN;
	saveout(Count).EstMdl = NaN;
	saveout(Count).EstParamCov = NaN;
	saveout(Count).LogL = LogL(Count);
	saveout(Count).info = NaN;
	saveout(Count).numParams = NaN;
	saveout(Count).params = NaN.*ones(1,5);
	saveout(Count).MdlNo = Count;
	saveout(Count).Const = ...
	NaN;
	end

	[~, warnidLoop] = lastwarn;
	if strcmp(warnidLoop,warnid)
	OptWarnThrown(Count) = true;
	end
	clear warnidLoop

	numParam_temp = sum(any(EstParamCov_temp));
	numObs =  length(ts);

	[aic(Count), bic(Count)] = aicbic(LogL(...
	Count), numParam_temp, numObs);

	saveout(Count).MeansMdl = MeansMdl;
	saveout(Count).EstMdl = EstMdl_temp;
	saveout(Count).EstParamCov = ...
	EstParamCov_temp;
	saveout(Count).LogL = LogL(Count);
	saveout(Count).info = info;
	saveout(Count).numParams = numParam_temp;
	saveout(Count).params = ...
	([length(EstMdl_temp.AR), ...
	EstMdl_temp.D, ...
	length(EstMdl_temp.MA), ...
	NaN, NaN, ...
	length(EstMdl_temp.SAR), ...
	length(EstMdl_temp.SMA), ...
	Sea]);
	saveout(Count).MdlNo = Count;
	saveout(Count).Const = ...
	EstMdl_temp.Constant;
	Count = Count + 1;
	
end

end
end
end
end

% Output warnings if the diagnostic values have small
% spreads

if nanstd(aic)<=(nanmean(aic)/10)
	fprintf(['AIC standard deviation is %04.2f%%', ...
		' of the Mean .\r\n'], ...
		(nanstd(aic)/(nanmean(aic)/10))*100)
end
if nanstd(bic)<=(nanmean(bic)/10)
	fprintf(['BIC standard deviation is %04.2f%%', ...
		' of the nanmean \r\n'], ...
		(nanstd(bic)/(nanmean(bic)/10))*100)
end
if nanstd(LogL)<=abs(nanmean(LogL)/10)
	fprintf(['Log Likelihood standard deviation ', ...
		'%04.2f%% of the Mean \r\n'], ...
		(nanstd(LogL)/abs(nanmean(LogL)/10))*100)
end

EmptyRows = cellfun(@isempty,{saveout.EstMdl});
NaNRows = isnan(aic);

if isrow(EmptyRows) && iscolumn(NaNRows)
	EmptyRows = EmptyRows';
elseif iscolumn(EmptyRows) && isrow(NaNRows)
	NaNRows = NaNRows';
end


% Truncate to the actual number of iterations
varargout{1} = saveout(~EmptyRows & ~NaNRows);
varargout{2} = aic(~EmptyRows & ~NaNRows);
varargout{3} = bic(~EmptyRows & ~NaNRows);

% Restore disabled warnings
warning('on',warnStruct.identifier);

catch err2
	fprintf('%s\n',err2.message);
	
	% following lines: stack
	for e=1:length(err2.stack)
		fprintf('%s at %i\n', err2.stack(e).name, ...
			err2.stack(e).line);
	end
	% Restore disabled warnings
	warning('on', warnStruct.identifier);
	
end

end