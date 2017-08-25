% Create Synthetic Files

% � All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate 
% Design, 2016 
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function CreateSyntheticFiles(pathEPWfile, namesavefldr, nboot, recdata, ccdata, varargin)

% I haven't uploaded the climate change forecasts or recorded data, so
% the following two variables are set to false:
% ccdata = false; recdata = false;
% The problem is how to pull these data sources when this tool is being
% used online
p = inputParser;
p.FunctionName = 'CreateSyntheticFiles';

addRequired(p, 'pathEPWfile', @ischar)
addRequired(p, 'namesavefldr', @ischar)
addRequired(p, 'nboot', @isnumeric)
addRequired(p, 'recdata', @islogical)
addRequired(p, 'ccdata', @islogical)
addParameter(p, 'nameEPWfolder', @ischar)
addParameter(p, 'ccpath', '', @ischar)
addParameter(p, 'recpath', '', @ischar)

parse(p, pathEPWfile, namesavefldr, nboot, recdata, ccdata, varargin{:})

pathEPWfile = p.Results.pathEPWfile;
namesavefldr = p.Results.namesavefldr;
nboot = p.Results.nboot;
recdata = p.Results.recdata;
ccdata = p.Results.ccdata;
ccpath = p.Results.ccpath;
recpath = p.Results.recpath;

% % Initialise File Paths

[pathEPWfolder, nameEPWfile, ~] = fileparts(pathEPWfile);

nameEPWfolder = pathEPWfolder(end-2:end);

% prompt = ['Please enter a name for the folder in which all source ', ...
%     'and output files are placed. The best is to have a 3-letter ', ...
%     'code representing the weather station. (e.g., GENEVA --> GEN)\n'];
% pathMATsave = fullfile(pathEPWfolder, input(prompt,'s'));

pathMATsave = fullfile(pathEPWfolder, namesavefldr);

if strcmp(pathMATsave,' ')
	pathMATsave = pathEPWfolder;
end

figsavepath = fullfile(pathMATsave, 'SMYfigs');

% prompt = ['Please enter full path to folder containing CC files: \n', ...
%     '(Leave blank if you don''t have them)\n'];
% ccpath = input(prompt,'s');

if strcmp(ccpath,' ')
	ccpath = pathEPWfolder;
end

% if ~contains(figsavepath, nameEPWfolder)
% 	figsavepath = fullfile(figsavepath, nameEPWfolder);
% end

PresentPath = fullfile(figsavepath,'Presentation');

% prompt = 'How many files do you want out? \n';
% nboot = str2double(input(prompt,'s')); 

% Get default colours
DefaultColours
set(0, 'defaultAxesFontName', 'Helvetica')

% % When running special nboot values, change output
% % folders.
% if nboot~=100
% 	figsavepath = [figsavepath, num2str(nboot)];
% 	pathMATsave = [pathMATsave, num2str(nboot)];
% end

% This is the number of climate change files that will be
% generated.
bootlen = nboot;

% Folder for R
RsubFolder = fullfile(pathMATsave, ...
	sprintf('RinoutN%d',nboot));

if exist(figsavepath,'dir')~=7
	mkdir(figsavepath)
end
if exist(pathMATsave,'dir')~=7
	mkdir(pathMATsave)
end
if exist(PresentPath,'dir')~=7
	mkdir(PresentPath)
end
if exist(RsubFolder,'dir')~=7
	mkdir(RsubFolder)
end


% This function is only meant to generate files from
% TMY-style files, so there are only two conditionals here
% for reading the original data.
if ~contains(pathEPWfile,'Meteonorm')
	tmytable = WeatherFileParseEPWMeteonorm(pathEPWfile);
else
	tmytable = WeatherFileParseEPWUSDOE(pathEPWfile);
end


% Initialise time index variables
N = size(tmytable,1); % Should be 8760
t = (1:N)'; % Julian Hour Index
% d = (1:365)'; % h = (1:24)';
nm = 12; % Number of months
mt = datetime(2015,1:nm,1); % Month number
% Number of days in each month
% nm_days = [31;28;31;30;31;30;31;31;30;31;30;31];

% Read the header info to obtain extreme and design
% conditions, in addition to location information.
StationInfo = TMYHeaderReader(pathEPWfile);

% % Dry-bulb Temperature
TDBk = tmytable.TDB + 273.15; % Temperature in Kelvin
% Converting temperature to Kelvin frees us to deal only
% with positive values


% % Global Horizontal Solar Radiation
ghi = tmytable.GHI;
dni = tmytable.DNI;
dhi = tmytable.DHI;
% Limit below which a reading is considered to be ZERO
ghi_limit = 1; % w/m2
% This limit applies to all solar radiation quantities.
% That is to say that if the global horizontal radiation
% was effectively zero, everything else will be zero too.

ghi(ghi <= ghi_limit) = 0;
dhi(ghi <= ghi_limit) = 0;
dni(ghi <= ghi_limit) = 0;

%%
% Use the following convention when putting the three time
% series together: tdb - 1, ghi - 2, w - 3

% Check and save cross-correlations between tdb, ghi, and
% w series
corrOptions.nLags = 7*24; % 7 days
corrOptions.nMA = 0;
corrOptions.nStd = 2;
% Number of standard deviations for confidence bounds

% This is a flag to give the skewness and kurtosis
% functions. The flag indicates that what is being
% calculated is the sample skewness, or kurtosis, which
% means that they are biased by a systematic amount based on
% the sample size. MATLAB will attempt to correct for this
% if flag = 0
shapeflag = 0;

[dStats.tmy.tdb.max, dStats.tmy.tdb.min, ...
	dStats.tmy.tdb.mean, dStats.tmy.tdb.std, ...
	dStats.tmy.tdb.skew, dStats.tmy.tdb.kurt, ...
	dStats.tmy.tdb.Days] = grpstats(tmytable.TDB, ...
	{tmytable.Month, tmytable.Day}, ...
	{@nanmax, @nanmin, @nanmean,@nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});

[dStats.tmy.rh.max, dStats.tmy.rh.min, ...
	dStats.tmy.rh.mean, dStats.tmy.rh.std, ...
	dStats.tmy.rh.skew, dStats.tmy.rh.kurt, ...
	dStats.tmy.rh.Days] = grpstats(tmytable.RH, ...
	{tmytable.Month, tmytable.Day}, ...
	{@nanmax, @nanmin, @nanmean,@nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});

[dStats.tmy.ghi.sum, dStats.tmy.ghi.max, ...
	dStats.tmy.ghi.min, ...
	dStats.tmy.ghi.mean, dStats.tmy.ghi.std, ...
	dStats.tmy.ghi.skew, dStats.tmy.ghi.kurt, ...
	dStats.tmy.ghi.Days] = grpstats(tmytable.GHI, ...
	{tmytable.Month, tmytable.Day}, ...
	{@nansum, @nanmax, @nanmin, @nanmean, @nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});

[dStats.tmy.dni.sum, dStats.tmy.dni.max, ...
	dStats.tmy.dni.min, ...
	dStats.tmy.dni.mean, dStats.tmy.dni.std, ...
	dStats.tmy.dni.skew, dStats.tmy.dni.kurt, ...
	dStats.tmy.dni.Days] = grpstats(tmytable.DNI, ...
	{tmytable.Month, tmytable.Day}, ...
	{@nansum, @nanmax, @nanmin, @nanmean, @nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});

[dStats.tmy.dhi.sum, dStats.tmy.dhi.max, ...
	dStats.tmy.dhi.min, ...
	dStats.tmy.dhi.mean, dStats.tmy.dhi.std, ...
	dStats.tmy.dhi.skew, dStats.tmy.dhi.kurt, ...
	dStats.tmy.dhi.Days] = grpstats(tmytable.DHI, ...
	{tmytable.Month, tmytable.Day}, ...
	{@nansum, @nanmax, @nanmin, @nanmean,@nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});

% Plotting correlations between ghi and tdb statistics
% Calculate the correlation coefficients for TMY data
[Corr.tmy.r,Corr.tmy.pr] = corr(tmytable.TDB, ...
	[tmytable.RH, tmytable.GHI, ...
	tmytable.TDP], ...
	'type', 'Pearson', 'rows','complete');
[Corr.tmy.rho,Corr.tmy.prho] = corr(tmytable.TDB, ...
	[tmytable.RH, tmytable.GHI, ...
	tmytable.TDP], ...
	'type', 'Spearman', 'rows','complete');

[Corr.tmysum.r, Corr.tmysum.pr] = corr( ...
	[dStats.tmy.ghi.sum, ...
	dStats.tmy.tdb.mean], 'type', 'Pearson');
[Corr.tmysum.rho, Corr.tmysum.prho] = corr( ...
	[dStats.tmy.ghi.sum, ...
	dStats.tmy.tdb.mean], 'type', 'Spearman');

% Calculate summary statistics for the TMY data
[mStats.tmy.tdb.max, mStats.tmy.tdb.min, ...
	mStats.tmy.tdb.mean, ...
	mStats.tmy.tdb.std, ...
	mStats.tmy.tdb.skew, ...
	mStats.tmy.tdb.kurt, ...
	mStats.tmy.tdb.months] = ...
	grpstats(tmytable.TDB, tmytable.Month, ...
	{@(x) (quantile(x,0.99)), @(x) (quantile(x,0.01)), ...
	@nanmean,@nanstd, @(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});

[mStats.tmy.rh.max, mStats.tmy.rh.min, ...
	mStats.tmy.rh.mean, ...
	mStats.tmy.rh.std, ...
	mStats.tmy.rh.skew, ...
	mStats.tmy.rh.kurt, ...
	mStats.tmy.rh.months] = ...
	grpstats(tmytable.RH, tmytable.Month, ...
	{@(x) (quantile(x,0.99)), @(x) (quantile(x,0.01)), ...
	@nanmean,@nanstd, @(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});


% Take the first difference of the TMY time series, in
% Celsius
diffs.tmy.tdb.raw = diff(tmytable.TDB);
% Find the min and max of the original and synthetic series
diffs.tmy.tdb.max = max(diffs.tmy.tdb.raw);
diffs.tmy.tdb.min = min(diffs.tmy.tdb.raw);

% Take the first difference of the w time series,
% dimensionless
diffs.tmy.rh.raw = diff(tmytable.RH);
% % Find the min and max of the TMY series
diffs.tmy.rh.max = max(diffs.tmy.rh.raw);
diffs.tmy.rh.min = min(diffs.tmy.rh.raw);

diffs.tmy.ghi.raw = diff(ghi);
diffs.tmy.ghi.max = max(diffs.tmy.ghi.raw);
diffs.tmy.ghi.min = min(diffs.tmy.ghi.raw);

% Calculate summary statistics for the TMY data
[mStats.tmy.ghi.max, mStats.tmy.ghi.min, ...
	mStats.tmy.ghi.mean, ...
	mStats.tmy.ghi.std, ...
	mStats.tmy.ghi.skew, ...
	mStats.tmy.ghi.kurt, ...
	mStats.tmy.ghi.months] = ...
	grpstats(ghi(ghi>0), tmytable.Month(ghi>0), ...
	{@(x) (quantile(x,0.99)), @(x) (quantile(x,0.01)), ...
	@expfit,@iqr, @(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});

close all
%%
% Get the current AMY files, concatenated and de-duplicated,
% in a TABLE

% These are the possible file names and save paths
MATFileNames = {fullfile(pathEPWfolder, ...
	['TruncatedRecordedData', '_', nameEPWfolder,'.mat']);
	fullfile(pathEPWfolder, ...
	['RecordedData', '_', nameEPWfolder, '.mat']);
	fullfile(pathMATsave, ...
	['TruncatedRecordedData', '_', nameEPWfolder,'.mat']);
	fullfile(pathMATsave, ...
	['RecordedData', '_', nameEPWfolder,'.mat'])};
SavedFileExists = cellfun(@(x) (exist(x,'file')==2), ...
	MATFileNames);

% Prefer full file instead of truncated
if SavedFileExists(2)
	load(MATFileNames{2})
else
	if SavedFileExists(4)
		load(MATFileNames{4})
	else
		if SavedFileExists(1)
			load(MATFileNames{1})
		else
			if SavedFileExists(3)
				load(MATFileNames{3})
			else
				recdata = false;
			end
		end
	end
end

% The MAT file loaded above should have a table called
% RecTables
if recdata
	
	% Delete the minutes column - it isn't used anywhere
	RecTables.Minute = [];
	
	% Calculate the daily statistics
	[mprep,dprep,tprep] = prepareSurfaceData(...
		RecTables.Month, RecTables.Day, RecTables.TDB);
	[dStats.rec.tdb.max, dStats.rec.tdb.min, ...
		dStats.rec.tdb.mean, dStats.rec.tdb.std, ...
		dStats.rec.tdb.skew, dStats.rec.tdb.kurt, ...
		dStats.rec.tdb.Days] = grpstats(tprep, ...
		{mprep,dprep}, {@nanmax, @nanmin, ...
		@nanmean,@nanstd, @(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)),'gname'});
	
	[mprep,dprep,wprep] = prepareSurfaceData(...
		RecTables.Month, RecTables.Day, RecTables.RH);
	[dStats.rec.rh.max, dStats.rec.rh.min, ...
		dStats.rec.rh.mean, dStats.rec.rh.std, ...
		dStats.rec.rh.skew, dStats.rec.rh.kurt, ...
		dStats.rec.rh.Days] = grpstats(wprep, ...
		{mprep,dprep}, {@nanmax, @nanmin, ...
		@nanmean,@nanstd, @(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)),'gname'});
	
	% % ghi stats are calculated later
	
	% Actual date time objects
	TimeActual = datenum(RecTables.Year, ...
		RecTables.Month, RecTables.Day, ...
		RecTables.Hour, zeros(size(RecTables,1),1), ...
		zeros(size(RecTables,1),1));
	
	% Step sizes in the recorded data, in hours
	StepSizes = diff(TimeActual);
	
	% Clean the rh values
	RecTables.RH(RecTables.RH>=100 | RecTables.RH<=0) = nan;
	
	% Special cleaning for certain stations
	if ~isempty(strfind(nameEPWfile,'GEN'))
		% GENEVA, which has two very large jumps
		% Between 22:00 on December 31, 1969 and
		% 1:00 on January 1, 1973
		StepSizes(RecTables.Year==1969 & ...
			RecTables.Month==12 & ...
			RecTables.Day==31 & RecTables.Hour==22) = NaN;
		% Between 22:00 on December 31, 1963 and
		% 1:00 on August 11, 1964
		StepSizes(RecTables.Year==1963 & ...
			RecTables.Month==12 & ...
			RecTables.Day==31 & RecTables.Hour==22) = NaN;
		
	elseif ~isempty(strfind(nameEPWfile,'DEL'))
		% While the lowest recorded temperature in DELHI has
		% has been -2.2, the zeros in the record are
		% unfortunately gaps in the record. Also, there are
		% some funky values (e.g., -37�).
		RecTables.TDB(RecTables.TDB<=0) = NaN;
		
		nantdb = isnan(RecTables.TDB);
		
		y = RecTables.TDB; y(nantdb) = NaN;
		RecTables.TDB = inpaint_nans(y);
		
		% The function we have for calculating TDP from
		% rh and TDB is not super-accurate, but then TDP
		% is not a particularly important metric. It is not
		% used in the EPW files. ()
		RecTables.TDP(nantdb) = nan;
		
	end
	
	clear y
	
	% Calculate summary statistics for the recorded data
	[mStats.rec.tdb.max, mStats.rec.tdb.min, ...
		mStats.rec.tdb.mean, ...
		mStats.rec.tdb.std, mStats.rec.tdb.skew, ...
		mStats.rec.tdb.kurt, mStats.rec.tdb.months] = ...
		grpstats(RecTables.TDB, RecTables.Month, ...
		{@(x) (quantile(x,0.99)), ...
		@(x) (quantile(x,0.01)), ...
		@nanmean, @nanstd, ...
		@(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)), 'gname'});
	
	[mStats.rec.rh.max, mStats.rec.rh.min, ...
		mStats.rec.rh.mean, ...
		mStats.rec.rh.std, mStats.rec.rh.skew, ...
		mStats.rec.rh.kurt, mStats.rec.rh.months] = ...
		grpstats(RecTables.RH, RecTables.Month, ...
		{@(x) (quantile(x,0.99)), ...
		@(x) (quantile(x,0.01)), ...
		@nanmean, @nanstd, ...
		@(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)), 'gname'});
	
	StepSizesHours = datevec(StepSizes);
	StepSizesHours = StepSizesHours(:,4);
	% Avoid blowing up the difference
	StepSizesHours(StepSizesHours<1) = NaN;
	% Denominators of less than 1, when passed to the diff
	% function, cause differentials of -+inf.
		
	% Take the first difference of the TMY time series, in
	% Kelvin
	diffs.rec.tdb.raw = diff(RecTables.TDB)./StepSizesHours;
	
	% We found some step sizes to be too large (~15� /
	% hour). These might be result of errors in measurement.
	% Now, we sort out those differences that may have been
	% caused by erroneous measurements.
	
	% Remove NaNs and Infs. Calculate z-score.
	[~,ydata] = prepareCurveData([], diffs.rec.tdb.raw);
	diffs.rec.tdb.zsco = zscore(ydata);
	
	% This parameter can be tuned, depending on how
	% aggressively one is looking for extreme temperatures.
	diffs.rec.tdb.zcut = max([abs(quantile( ...
		diffs.rec.tdb.zsco,0.996)), ...
		abs(quantile(diffs.rec.tdb.zsco,0.004))]);
	
	diffs.rec.tdb.cenabs = abs(diffs.rec.tdb.zsco) >= ...
		diffs.rec.tdb.zcut;
	
	diffs.rec.tdb.raw(diffs.rec.tdb.cenabs) = NaN;
	
	% Find the min and max of the censored differences
	diffs.rec.tdb.max = max(diffs.rec.tdb.raw);
	diffs.rec.tdb.min = min(diffs.rec.tdb.raw);
	% These DO NOT form the limits on the acceptable
	% differences in the synthetic series.
	
	% Repeat for the Humidity Ratio values %
	
	% Take the first difference of the recorded W time
	% series
	diffs.rec.rh.raw = diff(RecTables.RH)./StepSizesHours;
	
	% The Modified Z-score
	diffs.rec.rh.zsco = zscore(diffs.rec.rh.raw);
	
	diffs.rec.rh.zcut = max([abs(quantile( ...
		diffs.rec.rh.zsco,0.99)), ...
		abs(quantile(diffs.rec.rh.zsco,0.01))]);
	
	diffs.rec.rh.cenabs = abs(diffs.rec.rh.zsco) >= ...
		diffs.rec.rh.zcut;
	
	diffs.rec.rh.raw(diffs.rec.rh.cenabs) = NaN;
	
	% % Find the min and max of the censored differences
	diffs.rec.rh.max = max(diffs.rec.rh.raw);
	diffs.rec.rh.min = min(diffs.rec.rh.raw);
	% These DO NOT form the limits on the acceptable
	% differences in the synthetic series since recorded
	% values are not 'available' to the generating function.
	
	% Do the same for the ghi values
	
	% The ghi values have a MUCH SHORTER record than the
	% temperature and humidity values. Therefore, we only
	% consider the 'valid' length.
	
	% First, find the first non-zero ghi value
	GHIfirst = find(RecTables.GHI>0,1,'first');
	% Find the nearest start of the day (i.e. hour 1)
	GHIstartSolar = find(RecTables.Hour(1:GHIfirst)==1, ...
		1, 'last');
	% Set the data before this starting hour to be nan.
	RecTables.GHI(1:GHIstartSolar-1) = NaN;
	% After this, we assume that the ETRH doesn't change
	% with the years, so the ETRH vector from the TMY file
	% can be replicated
	
	% This is the valid record of ghi
	GHIrec = RecTables.GHI(GHIstartSolar:end);
	
	% Assign time stamps with the same year to all valid
	% ghi values (i.e. after the first non-zero value)
	TimeStampsGHI = datenum(2015, ...
		RecTables.Month(GHIstartSolar:end), ...
		RecTables.Day(GHIstartSolar:end), ...
		RecTables.Hour(GHIstartSolar:end), 0,0);
	
	[yprep,mprep,dprep,sprep] = prepareSurfaceData( ...
		RecTables.Year(GHIstartSolar:end), ...
		RecTables.Month(GHIstartSolar:end), ...
		RecTables.Day(GHIstartSolar:end), GHIrec);
	% Calculate GHi stats now that the valid ghi record is
	% known
	[dStats.rec.ghi.max, dStats.rec.ghi.min, ...
		dStats.rec.ghi.mean, dStats.rec.ghi.std, ...
		dStats.rec.ghi.skew, dStats.rec.ghi.kurt, ...
		dStats.rec.ghi.Days] = grpstats(sprep, ...
		{yprep, mprep,dprep}, {@nanmax, @nanmin, ...
		@nanmean,@nanstd, @(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)),'gname'});
	
	clear yprep mprep dprep sprep
	
	% Step sizes in the recorded data, in hours
	StepSizesGHI = diff(TimeStampsGHI);
	StepSizesHoursGHI = datevec(StepSizesGHI);
	StepSizesHoursGHI = StepSizesHoursGHI(:,4);
	StepSizesHoursGHI(StepSizesHoursGHI<1) = NaN;
	% Avoid blowing up the difference
	
	% This is the approximate first derivative
	diffs.rec.ghi.raw = diff(GHIrec)./StepSizesHoursGHI;
	
	% Calculate the median absolute deviation
	% 	diffs.rec.ghi.mad = mad(diffs.rec.ghi.raw,1);
	
	[~,ydata] = prepareCurveData([], diffs.rec.ghi.raw);
	diffs.rec.ghi.zsco = zscore(ydata);
	
	clear ydata
	
	% The handbook quotes Iglewicz and Hoaglin as
	% recommending labelling modified Z-scores with an
	% absolute value of greater than 3.5 as potential
	% outliers.
	
	% We found that to be too conservative in general, since
	% we are looking to induce some variation in the 'mean'
	% data. This parameter can be tuned, depending on how
	% aggressively one is looking for extreme temperatures.
	diffs.rec.ghi.zcut = max([abs(quantile( ...
		diffs.rec.ghi.zsco,0.99)), ...
		abs(quantile(diffs.rec.ghi.zsco,0.01))]);
	
	diffs.rec.ghi.cenabs = abs(diffs.rec.ghi.zsco) >= ...
		diffs.rec.ghi.zcut;
	
	diffs.rec.ghi.raw(diffs.rec.ghi.cenabs) = NaN;
	
	% Find the min and max of the censored differences
	diffs.rec.ghi.max = max(diffs.rec.ghi.raw);
	diffs.rec.ghi.min = min(diffs.rec.ghi.raw);
	
	MonthRec = RecTables.Month(GHIstartSolar:end);
	% Calculate summary statistics for the recorded data
	[mStats.rec.ghi.max, mStats.rec.ghi.min, ...
		mStats.rec.ghi.mean, ...
		mStats.rec.ghi.std, mStats.rec.ghi.skew, ...
		mStats.rec.ghi.kurt, mStats.rec.ghi.months] = ...
		grpstats(GHIrec(GHIrec>0), MonthRec(GHIrec>0), ...
		{@(x) (quantile(x,0.99)), ...
		@(x) (quantile(x,0.01)), ...
		@expfit, @iqr, @(x) (skewness(x,shapeflag)), ...
		@(y) (kurtosis(y,shapeflag)), 'gname'});
	% The ghi statistics are slightly different - the mean
	% is calculated using an appximate exponential
	% distribution fit to the sample. The std is replaced by
	% an IQR. 
	
end

clear yprep mprep dprep sprep tprep wprep
clear StepSizesHours StepSizes
clear StepSizesGHI StepSizesHoursGHI

%%
% Implementing the method of Magnano et al for Dry Bulb
% Temperature (tdb)

% There are two separate fourier series fit to the data.
% FOURIER FIT - LOW FREQUENCY This is to remove the annual
% seasonal component. Instead of removing the calculated
% daily means like Magnano et al, we use the low-fequency
% Fourier fit
% FOURIER FIT - HIGH FREQUENCY This is to remove
% the daily sinusoidal component.

% tdb Fit options are initialised here. The custom function
% is a sum of three fourier series, i.e. one pair of sine
% and cosine terms to account for annual variation, another
% pair for the daily fluctuations, and a third one for 2
% cycles in a year. The third term is unexpected but it is
% used to induce an asymmetry in the Fourier term, which is
% better to explain the unexpectedly high temperature
% generally seen in Autumn. Unexpected if the average
% temperature was perfectly symmetric about peak summer
% in a given year.
FitFunc.tdb.Hourly = fittype(['a0 + ', ...
	'a1*cos(2*pi*x/8760) + b1*sin(2*pi*x/8760) + ' ...
	'a2*cos(2*pi*x/4380) + b2*sin(2*pi*x/4380) + ', ...
	'a3*cos(2*pi*x/24) + b3*sin(2*pi*x/24)'], ...
	'independent', 'x', 'dependent', 'y');

FitFunc.tdb.LowHourly = fittype([ 'a0 + ', ...
	'a1*cos(2*pi*x/8760) + b1*sin(2*pi*x/8760) + ', ...
	'a2*cos(2*pi*x/4380) + b2*sin(2*pi*x/4380)'], ...
	'independent', 'x', 'dependent', 'y');

FitFunc.tdb.HighHourly = fittype( ...
	'a3*cos(2*pi*x/24) + b3*sin(2*pi*x/24)', ...
	'independent', 'x', 'dependent', 'y');

% The humidity ratio series is the only one without a daily
% sinusoidal component
FitFunc.rh.Hourly = fittype(['a0 + ', ...
	'a1*cos(2*pi*x/8760) + b1*sin(2*pi*x/8760) + ', ...
	'a3*cos(2*pi*x/24) + b3*sin(2*pi*x/24)'], ...
	'independent', 'x', 'dependent', 'y');

FitFunc.rh.LowHourly = fittype(['a0 + ', ...
	'a1*cos(2*pi*x/8760) + b1*sin(2*pi*x/8760)'], ...
	'independent', 'x', 'dependent', 'y');

FitFunc.rh.HighHourly = fittype(...
	'a3*cos(2*pi*x/24) + b3*sin(2*pi*x/24)', ...
	'independent', 'x', 'dependent', 'y');

% The method to find the best fit line through the data
FitOpts = fitoptions( 'Method', 'NonlinearLeastSquares' );
FitOpts.Display = 'Off';


% Prepare the data for a curve fit
[xdata.tdb, ydata.tdb] = prepareCurveData([], ...
	TDBk);
[xdata.rh, ydata.rh] = prepareCurveData([], ...
	tmytable.RH);

% Fit model to data.
[fourfits.tdb.fit, fourfits.tdb.gof, ...
	fourfits.tdb.out] = fit(xdata.tdb, ...
	ydata.tdb, FitFunc.tdb.Hourly,  FitOpts);
fourfits.tdb.evalmu = ...
	feval(fourfits.tdb.fit,xdata.tdb);

[fourfits.tdb.fitLow, fourfits.tdb.gofLow, ...
	fourfits.tdb.outLow] = fit(xdata.tdb, ...
	ydata.tdb, FitFunc.tdb.LowHourly,  FitOpts);
fourfits.tdb.evalmuLow = ...
	feval(fourfits.tdb.fitLow,xdata.tdb);

% Before evaluating the high frequency fit, remove the low
% frequency fit from the raw temperature data. In other
% words, fit to the residuals of the earlier fit.

[fourfits.tdb.fitHigh, fourfits.tdb.gofHigh, ...
	fourfits.tdb.outHigh] = fit(xdata.tdb, ...
	fourfits.tdb.outLow.residuals, ...
	FitFunc.tdb.HighHourly,  FitOpts);
fourfits.tdb.evalmuHigh = ...
	feval(fourfits.tdb.fitHigh, xdata.tdb);


[fourfits.rh.fit, fourfits.rh.gof, ...
	fourfits.rh.out] = fit(xdata.rh, ...
	ydata.rh, FitFunc.rh.Hourly,  FitOpts);
fourfits.rh.evalmu = ...
	feval(fourfits.rh.fit, xdata.rh);

[fourfits.rh.fitLow, fourfits.rh.gofLow, ...
	fourfits.rh.outLow] = fit(xdata.rh, ...
	ydata.rh, FitFunc.rh.LowHourly,  FitOpts);
fourfits.rh.evalmuLow = ...
	feval(fourfits.rh.fitLow,xdata.rh);
[fourfits.rh.fitHigh, fourfits.rh.gofHigh, ...
	fourfits.rh.outHigh] = fit(xdata.rh, ...
	fourfits.rh.outLow.residuals, ...
	FitFunc.rh.HighHourly,  FitOpts);
fourfits.rh.evalmuHigh = ...
	feval(fourfits.rh.fitHigh, xdata.rh);


clear xdata ydata

% Save the fourier fits
save(fullfile(pathMATsave, ['FourierFits_', ...
	nameEPWfile,'.mat']), 'fourfits')

close all

% Standardising the residuals does not confer any
% advantage, so we leave that out.
%%


% In the case of New York LaGuardia, models with MA lags are
% creating higher values of autocorrelation in simulated
% residuals. The (P)ACF plots indicate a pure AR model, but
% the log likelihood values don't. 

% % The moral of the story is that choosing a model without human
% % intervention is prone to failure!!

% No seasonal
[hrmodels.tdb, ic.aic.tdb, ic.bic.tdb] = ...
	FitModelExplore(fourfits.tdb.out.residuals, ...
	'ar',0:4, 'ma',0:4, ...
	'sma', 0, 'sar', 0,  'varmdl', 'constant');
[hrmodels.rh, ic.aic.rh, ic.bic.rh] = ...
	FitModelExplore(fourfits.rh.out.residuals, ...
	'ar',0:4, 'ma',0:4, 'sma', 0, ...
	'sar', 0,  'varmdl', 'constant');

% With seasonal
[tdb2, aictdb2, bictdb2] = ...
	FitModelExplore(fourfits.tdb.out.residuals, ...
	'ar',0:4, 'ma',0:4, ...
	'sma', 24, 'sar', 24,  'varmdl', 'constant');
hrmodels.tdb = [hrmodels.tdb, tdb2];
ic.aic.tdb = [ic.aic.tdb;aictdb2];
ic.bic.tdb = [ic.bic.tdb;bictdb2];
% Only SMA
[tdb2, aictdb2, bictdb2] = ...
	FitModelExplore(fourfits.tdb.out.residuals, ...
	'ar',0:4, 'ma',0:4, ...
	'sma', 24, 'sar', 0,  'varmdl', 'constant');
hrmodels.tdb = [hrmodels.tdb, tdb2];
ic.aic.tdb = [ic.aic.tdb;aictdb2];
ic.bic.tdb = [ic.bic.tdb;bictdb2];

% Only SAR
[tdb2, aictdb2, bictdb2] = ...
	FitModelExplore(fourfits.tdb.out.residuals, ...
	'ar',0:4, 'ma',0:4, ...
	'sma', 0, 'sar', 24,  'varmdl', 'constant');
hrmodels.tdb = [hrmodels.tdb, tdb2];
ic.aic.tdb = [ic.aic.tdb;aictdb2];
ic.bic.tdb = [ic.bic.tdb;bictdb2];

[rh2, aicrh2, bicrh2] = ...
	FitModelExplore(fourfits.rh.out.residuals, ...
	'ar',0:4, 'ma',0:4, 'sma', 24, ...
	'sar', 24,  'varmdl', 'constant');
hrmodels.rh = [hrmodels.rh, rh2];
ic.aic.rh = [ic.aic.rh;aicrh2];
ic.bic.rh = [ic.bic.rh;bicrh2];
[rh2, aicrh2, bicrh2] = ...
	FitModelExplore(fourfits.rh.out.residuals, ...
	'ar',0:4, 'ma',0:4, 'sma', 24, ...
	'sar', 0,  'varmdl', 'constant');
hrmodels.rh = [hrmodels.rh, rh2];
ic.aic.rh = [ic.aic.rh;aicrh2];
ic.bic.rh = [ic.bic.rh;bicrh2];
[rh2, aicrh2, bicrh2] = ...
	FitModelExplore(fourfits.rh.out.residuals, ...
	'ar',0:4, 'ma',0:4, 'sma', 0, ...
	'sar', 24,  'varmdl', 'constant');
hrmodels.rh = [hrmodels.rh, rh2];
ic.aic.rh = [ic.aic.rh;aicrh2];
ic.bic.rh = [ic.bic.rh;bicrh2];
% end

save(fullfile(pathMATsave,['HourMdls_', ...
	nameEPWfile,'.mat']), 'hrmodels','ic')

% Hourly models with seasonal autoregressive and moving
% average components work better than those without.


selmdls.tdb = hrmodels.tdb(ic.aic.tdb==min(ic.aic.tdb));
selmdls.rh = hrmodels.rh(ic.aic.rh==min(ic.aic.rh));

% % Infer the residuals for the single selected model
[selmdls.tdb.resid,selmdls.tdb.Var, selmdls.tdb.LogL] = ...
	infer(selmdls.tdb.EstMdl, fourfits.tdb.out.residuals);
[selmdls.rh.resid,selmdls.rh.Var, selmdls.rh.LogL] = ...
	infer(selmdls.rh.EstMdl, fourfits.rh.out.residuals);

% SARMA clearly does better than ARMA

% % Censor outliers in the Residuals
resid.tdb.zsco = zscore(selmdls.tdb.resid);
resid.rh.zsco = zscore(selmdls.rh.resid);
% resid.atm.zsco = zscore(selmdls.atm.resid);

resid.tdb.zcut = max([abs(quantile(resid.tdb.zsco, ...
	0.99)), abs(quantile(resid.tdb.zsco, 0.01))]);
resid.rh.zcut = max([abs(quantile(resid.rh.zsco, ...
	0.99)), abs(quantile(resid.rh.zsco, 0.01))]);


resid.tdb.cenabs = abs(resid.tdb.zsco) >= resid.tdb.zcut;
resid.rh.cenabs = abs(resid.rh.zsco) >= resid.rh.zcut;

y = selmdls.tdb.resid; y(resid.tdb.cenabs) = NaN;
selmdls.tdb.resid = inpaint_nans(y);

y = selmdls.rh.resid; y(resid.rh.cenabs) = NaN;
selmdls.rh.resid = inpaint_nans(y);

clear resid

%%
% Trying block bootstrap. Suggested in Magnano et al and is
% also indicated by parcorr. In parcorr, there are clearly
% significant coefficients at around the 24 and 48 hour
% lags. Trying both a block size of 2*24 and 3*24 since
% Magnano et al suggest the latter.

% A block bootstrap can be done using datasample with
% replacement. Each block of B days is treated as a single
% 'data point'. The resampling weights are reduced for those
% days that are repeated in the reshaping function above.


% Generate a datetime matrix for the synthetic time series
% Try to get nboot years
SynStartYear = 3000;
SynYears = SynStartYear: SynStartYear+(nboot*5);
% Find leap years, by looking to see which of the candidates
% have 29 as the last day of the month for February
E = eomday(SynYears,2);
% Take out leap years
SynYears = SynYears((E ~= 29));
% Only need the first nboot years
SynYears = (SynYears(1:nboot))';


% % Create an year-long datetime vector
yeardate = (datetime(SynYears(1),1,1,1,0,0) : ...
	hours(1) : datetime(SynYears(1),12,31,24,0,0))';
% The increment is one hour.

% % Get the hours, days, months
syntime = [repmat(SynYears(1),N,1), month(yeardate), ...
	day(yeardate), hour(yeardate)];

% Repeat the vector nboot times to signify many years
syntime = repmat(syntime,nboot,1);

% Change the years column
for tt = 2:nboot
	syntime((tt-1)*N:tt*N) = SynYears(1)+(tt-1);
end

% Extract the months
SynMonths = syntime(:,2);

% This value indicates that the block bootstrap is being
% done by splitting the data into 12 subperiods - i.e. a
% month at a time.
SubPeriods = 12;


%%
% Use a 3 day block bootstrap.

b = 3;
% Now resample the residuals using a block size determined
% by the loop index
BlockSize = b*24;
[tdbsyn.resam, tdbResamIdx] = BootstrapCustom( ...
	selmdls.tdb.resid, BlockSize, nboot, SubPeriods);

% The resampled values of the humidity time series are
% resampled with the same indices as the tdb. That is, the
% two series are shuffled in tandem.
rhsyn.resam = selmdls.rh.resid(tdbResamIdx);

clear tdbResamIdx

% Write an R file to be able to simulate the ARIMA model
% with custom innovations

TSname = 'tdb';
% % Write Innovations to a file for R
filenameInn = sprintf('%sInnForR_%s.csv', ...
	TSname, nameEPWfile);
InnPath = fullfile(RsubFolder,filenameInn);
csvwrite(InnPath,tdbsyn.resam,1,0)

% % Original time series
filenameX1 = sprintf('%sForR_%s.csv', TSname, nameEPWfile);
X1path = fullfile(RsubFolder,filenameX1);
csvwrite(X1path,fourfits.tdb.out.residuals,1,0)
% % Output from R
filenameOutSim = sprintf('%sSimFromR_%s.csv', ...
	TSname, nameEPWfile);
OutPathSim = fullfile(RsubFolder,filenameOutSim);
filenameOutRes = sprintf('%sResFromR_%s.csv', ...
	TSname, nameEPWfile);
OutPathRes = fullfile(RsubFolder,filenameOutRes);
filenameOutModel = sprintf('%sModelFromR_%s.csv', ...
	TSname, nameEPWfile);
OutPathModel = fullfile(RsubFolder,filenameOutModel);

% Call a function which writes a custom R file and runs it
% using R. For this to work, R has to be on the system path
SimArimaR(InnPath, X1path, OutPathSim, OutPathRes, ...
	OutPathModel, ...
	length(selmdls.tdb.EstMdl.SAR), ...
	length(selmdls.tdb.EstMdl.AR), ...
	length(selmdls.tdb.EstMdl.MA), ...
	sum(cellfun(@(x) (x~=0), selmdls.tdb.EstMdl.SAR)), ...
	sum(cellfun(@(x) (x~=0), selmdls.tdb.EstMdl.SMA)), ...
	N, nboot, RsubFolder)
% Errors from R will be reported in ExtrafileFromRbatch.txt

% % Read the R output
fspec = repmat('%f',1,bootlen);
fid = fopen(OutPathSim, 'r');
temp = textscan(fid, fspec, ...
	'Delimiter', ',', 'TreatAsEmpty', 'NA');
fclose(fid);
temp2 = cell2mat(temp);
selmdls.tdb.SimResid_R = inpaint_nans(temp2);
clear temp temp2

fspec = repmat('%f',1,bootlen);
fid = fopen(OutPathRes, 'r');
temp = textscan(fid, fspec, ...
	'Delimiter', ',', 'TreatAsEmpty', 'NA');
fclose(fid);
temp2 = cell2mat(temp);
selmdls.tdb.Resid_R = inpaint_nans(temp2);
clear temp temp2


% Now simulate w, with its own bootstrapped noise

TSname = 'rh';
% % Write Innovations to a file for R
filenameInn = sprintf('%sInnForR_%s.csv',TSname, ...
	nameEPWfile);
InnPath = fullfile(RsubFolder,filenameInn);
csvwrite(InnPath,rhsyn.resam,1,0)

% % Original time series
filenameX1 = sprintf('%sForR_%s.csv', TSname, nameEPWfile);
X1path = fullfile(RsubFolder,filenameX1);
csvwrite(X1path,fourfits.rh.out.residuals,1,0)
% % Output from R
filenameOutSim = sprintf('%sSimFromR_%s.csv', TSname, ...
	nameEPWfile);
OutPathSim = fullfile(RsubFolder,filenameOutSim);
filenameOutRes = sprintf('%sResFromR_%s.csv', TSname, ...
	nameEPWfile);
OutPathRes = fullfile(RsubFolder,filenameOutRes);
filenameOutModel = sprintf('%sModelFromR_%s.csv', ...
	TSname, nameEPWfile);
OutPathModel = fullfile(RsubFolder,filenameOutModel);


% Call a function which writes a custom R file and runs it
% using R. For this to work, R has to be on the system path
SimArimaR(InnPath, X1path, OutPathSim, OutPathRes, ...
	OutPathModel, ...
	length(selmdls.rh.EstMdl.SAR), ...
	length(selmdls.rh.EstMdl.AR), ...
	length(selmdls.rh.EstMdl.MA), ...
	sum(cellfun(@(x) (x~=0), selmdls.rh.EstMdl.SAR)), ...
	sum(cellfun(@(x) (x~=0), selmdls.rh.EstMdl.SMA)), ...
	N, nboot, RsubFolder)
% Errors from R will be reported in ExtrafileFromRbatch.txt
% Read the R output
fspec = repmat('%f',1,bootlen);
fid = fopen(OutPathSim, 'r');
temp = textscan(fid, fspec, ...
	'Delimiter', ',', 'TreatAsEmpty', 'NA');
fclose(fid);
temp2 = cell2mat(temp);
selmdls.rh.SimResid_R = inpaint_nans(temp2);
clear temp temp2

zsco = zscore(selmdls.rh.SimResid_R);
zcut = max( ...
	[abs(quantile(zsco,0.99)), ...
	abs(quantile(zsco,0.01))]);
zcens = abs(zsco)>zcut;
y = selmdls.rh.SimResid_R;
y(zcens) = nan;
y = inpaint_nans(y);
selmdls.rh.SimResid_R = y;

fspec = repmat('%f',1,bootlen);
fid = fopen(OutPathRes, 'r');
temp = textscan(fid, fspec, ...
	'Delimiter', ',', 'TreatAsEmpty', 'NA');
fclose(fid);
temp2 = cell2mat(temp);
selmdls.rh.Resid_R = inpaint_nans(temp2);

clear temp temp2 zsco


% ghi is not resampled using a model. Instead, we create
% synthetic series by exploiting the close relationship
% between daily mean tdb and daily sum of ghi.

% Creating synthetic hourly time series
tdbsyn.yearly = bsxfun( @plus, ...
	fourfits.tdb.evalmu, selmdls.tdb.SimResid_R ) - 273.15;
% Reshape for censoring, plotting, and summary statistics
tdbsyn.col = reshape(tdbsyn.yearly,[],1);

rhsyn.yearly = bsxfun( @plus, fourfits.rh.evalmu, ...
	selmdls.rh.SimResid_R );
rhsyn.col = reshape(rhsyn.yearly,[],1);

clear TSname
%%

% Load climate change data

if ccdata

CCdataFilePath = fullfile(ccpath, 'HourlyFutureData.mat');

if exist(CCdataFilePath,'file')==2
	% If concatenated file already exists
	load(CCdataFilePath);
else
	% File does not exist
	CCdataFilePath = CatClimChangeData(ccpath);
	load(CCdataFilePath);
end

% This file should contain two variables:
% FutureWeather - a struct containing the individual
% variables, split by climate change scenario. Within
% each scenario, each variable is organised into a 3D
% matrix of 95x8760x6
% FutureTime - a vector of time for the years of the
% future climate files (2006:2100)

% Creating synthetic hourly time series with climate
% change scenarios
UniqueParams = fieldnames(FutureWeather);

for p = 1:length(UniqueParams)
	FutureWeather.(UniqueParams{p}) = ...
		structfun(@(x) reshape( ...
		permute(x, [2 1 3]), [], 6), ...
		FutureWeather.(UniqueParams{p}), ...
		'UniformOutput',0);
end

FutureYears = unique(FutureTime(:,1));

dVals.tdb.rcp45 = NaN(size(FutureTime,1),1);
dVals.tdb.rcp85 = NaN(size(FutureTime,1),1);
dVals.rh.rcp45 = NaN(size(FutureTime,1),1);
dVals.rh.rcp85 = NaN(size(FutureTime,1),1);
dVals.ghi.rcp45 = NaN(size(FutureTime,1),1);
dVals.ghi.rcp85 = NaN(size(FutureTime,1),1);

% Pick a random model, any one. It doesn't matter which
% model is picked, so long as we are aware that the
% models aer all equivalent. One model is picked for any
% run of the future value generating script. This means
% that each run is different.
PickModel = randi([1 6], 1);

% Monthly selection not used any more since only ONE
% model is used in ONE call to this script.

% Pick the string of 85*8760 values associated with a
% particular model. It is also possible to use a mean
% of all the models.
dVals.tdb.rcp85 = FutureWeather.TDBdmean.rcp85(:, ...
	PickModel);
dVals.tdb.rcp45 = FutureWeather.TDBdmean.rcp45(:, ...
	PickModel);

dVals.rh.rcp85 = WtoRH(FutureWeather.Wdmean.rcp85(:, ...
	PickModel),dVals.tdb.rcp85);
dVals.rh.rcp45 = WtoRH(FutureWeather.Wdmean.rcp45(:, ...
	PickModel),dVals.tdb.rcp45);

dVals.ghi.rcp85 = FutureWeather.GHIdmean.rcp85(:, ...
	PickModel);
dVals.ghi.rcp45 = FutureWeather.GHIdmean.rcp45(:, ...
	PickModel);

clear FutureWeather
clear xdata ydata


% Pick only a certain number of the
% simulated/bootstrapped residuals
PickBoot = randi([1 nboot], bootlen, 1);

tdbsynCC.rcp45.yearly = bsxfun(@plus, ...
	dVals.tdb.rcp45 - 273.15, ...
	repmat(bsxfun(@plus, ...
	fourfits.tdb.evalmuHigh, ...
	selmdls.tdb.SimResid_R(:,PickBoot)), ...
	length(FutureYears) , 1));
tdbsynCC.rcp85.yearly = bsxfun(@plus, ...
	dVals.tdb.rcp85 - 273.15, ...
	repmat(bsxfun(@plus, ...
	fourfits.tdb.evalmuHigh, ...
	selmdls.tdb.SimResid_R(:,PickBoot)), ...
	length(FutureYears) , 1));


rhsynCC.rcp45.yearly = bsxfun(@plus, ...
	dVals.rh.rcp45, ...
	repmat( bsxfun(@minus, ...
	selmdls.rh.SimResid_R(:,PickBoot), ...
	mean(selmdls.rh.SimResid_R(:,PickBoot))), ...
	length(FutureYears) , 1));

rhsynCC.rcp85.yearly = bsxfun(@plus, ...
	dVals.rh.rcp85, ...
	repmat( bsxfun(@minus, ...
	selmdls.rh.SimResid_R(:,PickBoot), ...
	mean(selmdls.rh.SimResid_R(:,PickBoot))), ...
	length(FutureYears) , 1));

clear dVals

% Reshape
tdbsynCC.rcp85.col = reshape( ...
	tdbsynCC.rcp85.yearly,[],1);
tdbsynCC.rcp45.col = reshape( ...
	tdbsynCC.rcp45.yearly,[],1);

% Wait until after the rhsyn series is censored to assign it
% to the future structs.

end

clear PickBoot fourfits

%%

% % CLEANING the tdb series first

% Calculate the first difference of the Synthetic series
diffs.syn.tdb.raw = [diff(tdbsyn.col); 0];
% If the first difference is too large or too small
% (compared to the max and min of the original series), then
% it should be censored.
diffcen.syn.tdb.rtlft = (diffs.syn.tdb.raw > ...
	diffs.tmy.tdb.max) | (diffs.syn.tdb.raw < ...
	diffs.tmy.tdb.min);

diffs.syn.tdb.zsco = zscore(diffs.syn.tdb.raw);

% Chicago seems to need a more aggresive cut for some
% reason.
if ~contains(nameEPWfile,'CHI')
	diffs.syn.tdb.zcut = max( ...
		[abs(quantile(diffs.syn.tdb.zsco,0.99)), ...
		abs(quantile(diffs.syn.tdb.zsco,0.01))]);
else
	% Take the maximum of the 99.9 and 0.01 percentile
	% Z-scores as the cutoff
	diffs.syn.tdb.zcut = max( ...
		[abs(quantile(diffs.syn.tdb.zsco,0.999)), ...
		abs(quantile(diffs.syn.tdb.zsco,0.001))]);
end

% Cut any points with a zscore greater than the cutoff
% (absolute)
diffcen.syn.tdb.topbot = abs(diffs.syn.tdb.zsco) >= ...
	diffs.syn.tdb.zcut;

% The censored series
y = tdbsyn.col;
y(diffcen.syn.tdb.rtlft | diffcen.syn.tdb.topbot) = NaN;
y = reshape(y, N, []);
tdbsyn.yearly = inpaint_nans(y);
tdbsyn.col = reshape(tdbsyn.yearly,[],1);
% Calculate the first difference again - these are the
% censored values
diffs.syn.tdb.cenyear = diff(tdbsyn.yearly,1,1);
diffs.syn.tdb.maxyear = max(diffs.syn.tdb.cenyear);
diffs.syn.tdb.minyear = min(diffs.syn.tdb.cenyear);

% Clean outlandish raw values using the same procedure
rawcen.syn.tdb.zsco = zscore(tdbsyn.col);

if ~contains(nameEPWfile,'CHI')
	rawcen.syn.tdb.zcut = max( ...
		[abs(quantile(diffs.syn.tdb.zsco,0.99)), ...
		abs(quantile(diffs.syn.tdb.zsco,0.01))]);
else
	rawcen.syn.tdb.zcut = max([abs(quantile( ...
		rawcen.syn.tdb.zsco,0.999)), abs(quantile( ...
		rawcen.syn.tdb.zsco,0.001))]);
end

% The outlandish values in the synthetic tdb occur only at
% the 'top' of the distribution, that is to say that some
% values tend to strangely high. But we censor on both
% sides.
rawcen.syn.tdb.topbot = abs(rawcen.syn.tdb.zsco) >= ...
	rawcen.syn.tdb.zcut;

% These are the highest and lowest temperatures ever
% recorded. This is valid for both 'present' weather as
% well as that derived from future climate scenarios.
rawcen.syn.tdb.rechilo = (tdbsyn.col>=56) | ...
	(tdbsyn.col<=-89);

% The censored series
y = tdbsyn.col;
y(rawcen.syn.tdb.topbot | rawcen.syn.tdb.rechilo) = NaN;
y = reshape(y,N,[]);
% Choosing method two (2) since the others seem to produce
% some new outliers, probably due to extrapolation.
tdbsyn.yearly = inpaint_nans(y,2);
% Reshape the vector for some calculations
tdbsyn.col = reshape(tdbsyn.yearly, [], 1);

if ccdata
	
	diffs.rcp45.tdb.raw = [diff(tdbsynCC.rcp45.col); 0];
	diffs.rcp85.tdb.raw = [diff(tdbsynCC.rcp85.col); 0];
	
	diffcen.rcp45.tdb.rtlft = ( ...
		diffs.rcp45.tdb.raw > diffs.tmy.tdb.max) | ...
		(diffs.rcp45.tdb.raw < diffs.tmy.tdb.min);
	diffcen.rcp85.tdb.rtlft = ( ...
		diffs.rcp85.tdb.raw > diffs.tmy.tdb.max) | ...
		(diffs.rcp85.tdb.raw < diffs.tmy.tdb.min);
	
	diffs.rcp45.tdb.zsco = zscore(diffs.rcp45.tdb.raw);
	diffs.rcp85.tdb.zsco = zscore(diffs.rcp45.tdb.raw);
	
	diffs.rcp45.tdb.zcut = max([abs(quantile( ...
		diffs.rcp45.tdb.zsco,0.999)), ...
		abs(quantile(diffs.rcp45.tdb.zsco,0.001))]);
	diffs.rcp85.tdb.zcut = max([abs(quantile( ...
		diffs.rcp85.tdb.zsco,0.999)), ...
		abs(quantile(diffs.rcp85.tdb.zsco,0.001))]);
	
	diffcen.rcp45.tdb.topbot = abs( ...
		diffs.rcp45.tdb.zsco) >= diffs.rcp45.tdb.zcut;
	diffcen.rcp85.tdb.topbot = abs( ...
		diffs.rcp85.tdb.zsco) >= diffs.rcp85.tdb.zcut;
	
	y = tdbsynCC.rcp45.col;
	y(diffcen.rcp45.tdb.rtlft | ...
		diffcen.rcp45.tdb.topbot) = NaN;
	y = reshape(y,N,[]);
	tdbsynCC.rcp45.yearly = inpaint_nans(y,2);
	tdbsynCC.rcp45.col = reshape( ...
		tdbsynCC.rcp45.yearly, [], 1);
	
	
	y = tdbsynCC.rcp85.col;
	y(diffcen.rcp85.tdb.rtlft | ...
		diffcen.rcp85.tdb.topbot) = NaN;
	y = reshape(y,N,[]);
	tdbsynCC.rcp85.yearly = inpaint_nans(y,2);
	tdbsynCC.rcp85.col = reshape( ...
		tdbsynCC.rcp85.yearly, [], 1);
	
	diffs.rcp45.tdb.cenyear = diff( ...
		tdbsynCC.rcp45.yearly,1,1);
	diffs.rcp85.tdb.cenyear = diff( ...
		tdbsynCC.rcp85.yearly,1,1);
	
	diffs.rcp45.tdb.maxyear = max(diffs.rcp45.tdb.cenyear);
	diffs.rcp45.tdb.minyear = min(diffs.rcp45.tdb.cenyear);
	
	diffs.rcp85.tdb.maxyear = max(diffs.rcp85.tdb.cenyear);
	diffs.rcp85.tdb.minyear = min(diffs.rcp85.tdb.cenyear);
	
	rawcen.rcp45.tdb.zsco = zscore(tdbsynCC.rcp45.col);
	rawcen.rcp45.tdb.zcut = max([abs(quantile( ...
		rawcen.rcp45.tdb.zsco,0.999)), ...
		abs(quantile(rawcen.rcp45.tdb.zsco,0.001))]);
	
	rawcen.rcp85.tdb.zsco = zscore(tdbsynCC.rcp85.col);
	rawcen.rcp85.tdb.zcut = max([abs(quantile( ...
		rawcen.rcp85.tdb.zsco,0.999)), ...
		abs(quantile(rawcen.rcp85.tdb.zsco,0.001))]);
	
	rawcen.rcp45.tdb.topbot = abs( ...
		rawcen.rcp45.tdb.zsco) >= rawcen.rcp45.tdb.zcut;
	rawcen.rcp85.tdb.topbot = abs( ...
		rawcen.rcp85.tdb.zsco) >= rawcen.rcp85.tdb.zcut;
	
	rawcen.rcp45.rechilo = ...
		(tdbsynCC.rcp45.col>=56) | ...
		(tdbsynCC.rcp45.col<=-89);
	rawcen.rcp85.rechilo = ...
		(tdbsynCC.rcp85.col>=56) | ...
		(tdbsynCC.rcp85.col<=-89);
	
	y = tdbsynCC.rcp45.col;
	y(rawcen.rcp45.tdb.topbot | ...
		rawcen.rcp45.rechilo) = NaN;
	y = reshape(y,N,[]);
	tdbsynCC.rcp45.yearly = inpaint_nans(y,2);
	tdbsynCC.rcp45.col = ...
		reshape(tdbsynCC.rcp45.yearly, [], 1);
	
	y = tdbsynCC.rcp85.col;
	y(rawcen.rcp85.tdb.topbot | ...
		rawcen.rcp85.rechilo) = NaN;
	y = reshape(y,N,[]);
	tdbsynCC.rcp85.yearly = inpaint_nans(y,2);
	tdbsynCC.rcp85.col = ...
		reshape(tdbsynCC.rcp85.yearly, [], 1);
	
	clear rawcen y
	
end

close all
%%

% % Cleaning the RH series

% We set a limit on the humidity ratio: it must be
% constrained between the minimum and maximum values in the
% TMY file
rawcen.syn.rh.topbot = (rhsyn.col > ...
	max(tmytable.RH)) | ...
	(rhsyn.col < min(tmytable.RH));
% Calculate the first difference of the Synthetic series,
% including nans.
diffs.syn.rh.raw = diff(rhsyn.col);
% If the first difference is too large or too small
% (compared to the max and min of the original series),
% then it should be censored.
diffcen.syn.rh.rtlft = [ ...
	(diffs.syn.rh.raw > diffs.tmy.rh.max) | ...
	(diffs.syn.rh.raw < diffs.tmy.rh.min); true];
% The censored series
y = rhsyn.col;
y(diffcen.syn.rh.rtlft | rawcen.syn.rh.topbot) = NaN;
y = reshape(y, N, []);
y2 = inpaint_nans(y);	
y2(y2>100) = 100; y2(y2<=0) = min(tmytable.RH);
rhsyn.yearly = y2;
rhsyn.col = reshape(rhsyn.yearly, [], 1);

% Calculate the first difference again
diffs.syn.rh.cenyear = diff(rhsyn.yearly,1,1);
diffs.syn.rh.maxyear = max(diffs.syn.rh.cenyear);
diffs.syn.rh.minyear = min(diffs.syn.rh.cenyear);

% if ccdata
% rhsynCC.rcp45.yearly = repmat(rhsyn.yearly,85,1);
% rhsynCC.rcp85.yearly = repmat(rhsyn.yearly,85,1);
% end

if ccdata	
	
	rhsynCC.rcp45.col = reshape(rhsynCC.rcp45.yearly,[],1);
	rhsynCC.rcp85.col = reshape(rhsynCC.rcp85.yearly,[],1);

	rawcen.rcp45.rh.topbot = ...
		(rhsynCC.rcp45.col > max(tmytable.RH)) | ...
		(rhsynCC.rcp45.col < min(tmytable.RH));
	rawcen.rcp85.rh.topbot = ...
		(rhsynCC.rcp85.col > max(tmytable.RH)) | ...
		(rhsynCC.rcp85.col < min(tmytable.RH));
	
	diffs.rcp45.rh.raw = diff(rhsynCC.rcp45.col);
	diffs.rcp85.rh.raw = diff(rhsynCC.rcp85.col);
	
	diffcen.rcp45.rh.rtlft = [ ...
		(diffs.rcp45.rh.raw > diffs.tmy.rh.max) | ...
		(diffs.rcp45.rh.raw < diffs.tmy.rh.min); false];
	diffcen.rcp85.rh.rtlft = [ ...
		(diffs.rcp85.rh.raw > diffs.tmy.rh.max) | ...
		(diffs.rcp85.rh.raw < diffs.tmy.rh.min); false];
	
	y = rhsynCC.rcp45.col;
	y(diffcen.rcp45.rh.rtlft | ...
		rawcen.rcp45.rh.topbot) = NaN;
	y = reshape(y, N, []);
	y2 = inpaint_nans(y);
	y2(y2>100) = 100;
	y2(y2<=0) = min(tmytable.RH);
	rhsynCC.rcp45.yearly = y2;
	rhsynCC.rcp45.col = reshape( ...
		rhsynCC.rcp45.yearly, [], 1);
	
	y = rhsynCC.rcp85.col;
	y(diffcen.rcp85.rh.rtlft | ...
		rawcen.rcp85.rh.topbot) = NaN;
	y = reshape(y, N, []);
	y2 = inpaint_nans(y);	
	y2(y2>100) = 100;
	y2(y2<=0) = min(tmytable.RH);
	rhsynCC.rcp85.yearly = y2;
	rhsynCC.rcp85.col = reshape( ...
		rhsynCC.rcp85.yearly, [], 1);
	
	clear diffcen rawcen y
		
end
%%

% ghi is resampled by exploiting the strong correlation
% between daily sum of ghi and daily mean of temperature.

% First, calculate the daily means of the synthetic
% temperatures
tdbsyn.daily.syn.raw = reshape(tdbsyn.col,24,[]);
tdbsyn.daily.syn.means = mean(tdbsyn.daily.syn.raw);

% Reshape the synthetic and tmy months to separate the means
% by month
SynMonths_res = reshape(syntime(:,2),24,[]);
SynMonths_res = (SynMonths_res(1,:))';
TMYmonths_res = reshape(tmytable.Month,24,[]);
TMYmonths_res = (TMYmonths_res(1,:))';
% It is necessary to do a nearest neighbour bootstrap
% month-by-month since it is quite possible that, for a
% given synthetic day, the closest daily mean in the TMY is
% in a different month.This means that the sunshine hours
% will be considerably different between the synthetic day
% and the neihbour TMY day. So, we get over this problem by
% making an arbitrary separation by month.


% Find the 'k' nearest neighbours in the TMY daily means for
% each synthetic daily mean.
knn = 10;
nIDX = NaN(size(SynMonths_res));

% Cycle through every month
for m = 1:length(mt)
	% Get the synthetic means for this month
	synmeans = tdbsyn.daily.syn.means(SynMonths_res==m)';
	% Get the tmy means for this month
	tmymeans = dStats.tmy.tdb.mean(TMYmonths_res==m);
	% Get the ending day of the current month - this is to
	% make sure that when each month is sampled, the indices
	% stored are for the YEAR, not just the month.
	if m == 1
		% January starts from idx = 1, so IDXadd = 0.
		IDXadd = 0;
	else
		% These are the end-of-month days for all months
		% preceding the current month.
		CurrEOMdays =  eomday(year(mt(1:m-1)),month(mt(1:m-1)));
		% These are the number of days that have passed,
		% i.e. the sum of the end-of-month days
		IDXadd = sum(CurrEOMdays);
	end
	
	% Get the knn nearest neihbours in the TMY record
	% for each daily mean in the synthetic time series.
	[temp1, ~] = knnsearch(tmymeans,synmeans, ...
		'K', knn, 'Distance','mahalanobis');
	
	% Pick one nearest neighbour at random
	nIDX(SynMonths_res==m) = temp1(:,randi([1 10]));
	
	% Add the preceding days of the year
	nIDX(SynMonths_res==m) = nIDX(SynMonths_res==m) + ...
		IDXadd;
	
end

clear m tmymeans synmeans temp1 temp2

% nIDX is the index of a day in the TMY record whose daily
% mean matches some day in the Synthetic series.

% Pick the corresponding daily sums
ghisyn.daily.sums = dStats.tmy.ghi.sum(nIDX);
% Record the daily means as well
ghisyn.daily.means = dStats.tmy.ghi.mean(nIDX);
% Reshape the hourly ghi to be separated by days
ghi_res = reshape(tmytable.GHI,24,[]);
% Pick the days corresponding to each nearest neighbour
ghisyn.daily.raw = ghi_res(:,nIDX);
% Reshape the days to be hourly
ghisyn.col = reshape(ghisyn.daily.raw,[],1);
ghisyn.yearly = reshape(ghisyn.col,N,[]);

% Also copy the other solar radiation values
% Direct Normal Irradiance
% Pick the corresponding daily sums
dnisyn.daily.sums = dStats.tmy.dni.sum(nIDX);
% Reshape the hourly dni to be separated by days
DNI_res = reshape(dni,24,N/24);
% Pick the days corresponding to each nearest neighbour
dnisyn.daily.raw = DNI_res(:,nIDX);
% Reshape the days to be hourly
dnisyn.col = reshape(dnisyn.daily.raw,[],1);
dnisyn.yearly = reshape(dnisyn.col,N,nboot);

% Diffuse Horizontal Irradiance
% Pick the corresponding daily sums
dhisyn.daily.sums = dStats.tmy.dhi.sum(nIDX);
% Reshape the hourly dhi to be separated by days
DHI_res = reshape(dhi,24,N/24);
% Pick the days corresponding to each nearest neighbour
dhisyn.daily.raw = DHI_res(:,nIDX);
% Reshape the days to be hourly
dhisyn.col = reshape(dhisyn.daily.raw,[],1);
dhisyn.yearly = reshape(dhisyn.col,N,nboot);

clear nIDX nD nIDXmonthly
%%

if ccdata
	tdbsynCC.rcp45.daily.raw = reshape( ...
		tdbsynCC.rcp45.col, 24, []);
	tdbsynCC.rcp45.daily.means = (mean( ...
		tdbsynCC.rcp45.daily.raw,1));
	
	tdbsynCC.rcp85.daily.raw = reshape( ...
		tdbsynCC.rcp85.col, 24, []);
	tdbsynCC.rcp85.daily.means = (mean( ...
		tdbsynCC.rcp85.daily.raw,1));
	
	
	% The distribution of future daily means is not
	% different from TMY or recorded data. This means that
	% the series can be treated somewhat interchangeably. We
	% take the future predictions and find the corresponding
	% days from the TMY file. This lets us reconstruct the
	% hourly profile from the TMY and future daily means.
	
	% Reshape the future time (month) vector, and keep only
	% one value per day. Then repeat it a bootlen times,
	% since the Future Daily Means will also be reshaped so
	% that each 'iteration' follows the last. One
	% 'iteration' is 85 years or 31025 values.
	FutMonths_res = reshape(FutureTime(:,2),24,[]);
	FutMonths_res = FutMonths_res(1,:)';
	FutMonths_res = repmat(FutMonths_res, bootlen, 1);
	
	
	% Find the 'knn' nearest neighbours in the TMY daily
	% means for each synthetic daily mean.
	knn = 10;
	% Begin bz pre-allocating the index vectors
	nIDX.rcp45 = NaN(length(FutMonths_res),1);
	nIDX.rcp85 = NaN(length(FutMonths_res),1);
	
	clear CurrEOMdays
	
	% Cycle through every month
	for m = 1:length(mt)
		
		% Get the future values for this month (first
		% index). Only one column from the FutMonths_res
		% variable needs to be used for this comparison
		% since the second index is just a repetition.
		FutureMeans.rcp45 = (tdbsynCC.rcp45.daily.means( ...
			FutMonths_res==m))';
		FutureMeans.rcp85 = (tdbsynCC.rcp85.daily.means( ...
			FutMonths_res==m))';
		
		% Get the tmy means for this month
		tmymeans = dStats.tmy.tdb.mean(TMYmonths_res==m);
		
		% Get the ending day of the current month. When each
		% month is sampled, the indices stored are for the
		% YEAR, not just the month. The way to do this is to
		% add the indices from the previous months.
		if m == 1
			% January starts from idx = 1, so IDXadd = 0.
			IDXadd = 0;
		else
			% These are the end-of-month days for all months
			% preceding the current month.
			CurrEOMdays =  eomday(year(mt(1:m-1)), ...
				month(mt(1:m-1)));
			% These are the number of days that have passed,
			% i.e. the sum of the end-of-month days
			IDXadd = sum(CurrEOMdays);
		end
		
		% Get the knn nearest neihbours in the TMY record of
		% ghi daily means for each daily mean in the future
		% time series. Once again, use only one replicate of
		% the series, since the future values have been
		% replicated bootlen times
		[temp1.rcp45, ~] = knnsearch(tmymeans, ...
			FutureMeans.rcp45, 'K', knn, ...
			'Distance','mahalanobis');
		[temp1.rcp85, ~] = knnsearch(tmymeans, ...
			FutureMeans.rcp85, 'K', knn, ...
			'Distance','mahalanobis');
		
		% Pick one nearest neighbour at random
		RandomNeighbour = randi([1 10]);
		nIDX.rcp45(FutMonths_res==m) = ...
			temp1.rcp45(:, RandomNeighbour);
		nIDX.rcp85(FutMonths_res==m) = ...
			temp1.rcp85(:, RandomNeighbour);
		
		% Add the preceding days of the year
		nIDX.rcp45(FutMonths_res==m) = ...
			nIDX.rcp45(FutMonths_res==m) + IDXadd;
		nIDX.rcp85(FutMonths_res==m) = ...
			nIDX.rcp85(FutMonths_res==m) + IDXadd;
		
	end
	
	clear m tmymeans synmeans temp1 temp2
	clear FutureDaily FutMonths_res TMYmonths_res
	
	% nIDX is the index of a day in the TMY record whose
	% daily mean matches some day in the Synthetic series.
	
	% Pick the daily means corresponding to each nearest
	% neighbour
	ghisynCC.rcp45.daily.means = ...
		dStats.tmy.ghi.mean(nIDX.rcp45);
	ghisynCC.rcp85.daily.means = ...
		dStats.tmy.ghi.mean(nIDX.rcp85);
	% Reshape the hourly ghi to be separated by days
	ghi_res = reshape(ghi,24,[]);
	% Pick the days corresponding to each nearest neighbour
	ghisynCC.rcp45.daily.raw = ghi_res(:,nIDX.rcp45);
	ghisynCC.rcp85.daily.raw = ghi_res(:,nIDX.rcp85);
	% Reshape the days to be hourly
	ghisynCC.rcp45.col = reshape( ...
		ghisynCC.rcp45.daily.raw,[],1);
	ghisynCC.rcp45.yearly = reshape( ...
		ghisynCC.rcp45.col,N,[]);
	ghisynCC.rcp85.col = reshape( ...
		ghisynCC.rcp85.daily.raw,[],1);
	ghisynCC.rcp85.yearly = reshape( ...
		ghisynCC.rcp85.col,N,[]);
		
		
	% Also copy the other solar radiation values
	% Direct Normal Irradiance
	% Reshape the hourly dni to be separated by days
	DNI_res = reshape(dni,24,[]);
	% Pick the days corresponding to each nearest neighbour
	dnisynCC.rcp45.daily.raw = DNI_res(:,nIDX.rcp45);
	dnisynCC.rcp85.daily.raw = DNI_res(:,nIDX.rcp85);
	% Reshape the days to be hourly
	dnisynCC.rcp45.col = reshape( ...
		dnisynCC.rcp45.daily.raw,[],1);
	dnisynCC.rcp45.yearly = reshape( ...
		dnisynCC.rcp45.col,N,[]);
	dnisynCC.rcp85.col = reshape( ...
		dnisynCC.rcp85.daily.raw,[],1);
	dnisynCC.rcp85.yearly = reshape( ...
		dnisynCC.rcp85.col,N,[]);
	
	% Diffuse Horizontal Irradiance
	% Reshape the hourly dhi to be separated by days
	DHI_res = reshape(dhi,24,N/24);
	% Pick the days corresponding to each nearest neighbour
	dhisynCC.rcp45.daily.raw = DHI_res(:,nIDX.rcp45);
	dhisynCC.rcp85.daily.raw = DHI_res(:,nIDX.rcp85);
	% Reshape the days to be hourly
	dhisynCC.rcp45.col = reshape( ...
		dhisynCC.rcp45.daily.raw,[],1);
	dhisynCC.rcp45.yearly = reshape( ...
		dhisynCC.rcp45.col,N,[]);
	dhisynCC.rcp85.col = reshape( ...
		dhisynCC.rcp85.daily.raw,[],1);
	dhisynCC.rcp85.yearly = reshape( ...
		dhisynCC.rcp85.col,N,[]);
	
	
	clear nIDX nIDXmonthly nD
	
end


% % Now calculate some summary statistics

[mStats.syn.tdb.max, mStats.syn.tdb.min, ...
	mStats.syn.tdb.mean, mStats.syn.tdb.std, ...
	mStats.syn.tdb.skew, mStats.syn.tdb.kurt, ...
	mStats.syn.tdb.months] = ...
	grpstats(tdbsyn.col, {syntime(:,1),syntime(:,2)}, ...
	{@(x) (quantile(x,0.99)),@(x) (quantile(x,0.01)), ...
	@nanmean, @nanstd,@(x) (skewness(x, ...
	shapeflag)), @(y) (kurtosis(y,shapeflag)),'gname'});
mStats.syn.tdb = structfun(@(x) reshape(x,12,[]), ...
	mStats.syn.tdb, 'UniformOutput', 0);

[mStats.syn.ghi.max, mStats.syn.ghi.min, ...
	mStats.syn.ghi.mean, mStats.syn.ghi.std, ...
	mStats.syn.ghi.skew, mStats.syn.ghi.kurt, ...
	mStats.syn.ghi.months] = ...
	grpstats(ghisyn.col(ghisyn.col>0), ...
	{syntime((ghisyn.col>0),1), ...
	syntime((ghisyn.col>0),2)}, ...
	{@(x) (quantile(x,0.99)),@(x) (quantile(x,0.01)), ...
	@expfit,@iqr, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});
mStats.syn.ghi = structfun(@(x) reshape(x,12,[]), ...
	mStats.syn.ghi, 'UniformOutput', 0);

[mStats.syn.rh.max, mStats.syn.rh.min, ...
	mStats.syn.rh.mean, ...
	mStats.syn.rh.std, mStats.syn.rh.skew, ...
	mStats.syn.rh.kurt, mStats.syn.rh.months] = ...
	grpstats(rhsyn.col, {syntime(:,1),syntime(:,2)}, ...
	{@(x) (quantile(x,0.99)),@(x) (quantile(x,0.01)), ...
	@nanmean,@nanstd,@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});
mStats.syn.rh = structfun(@(x) reshape(x,12,[]), ...
	mStats.syn.rh, 'UniformOutput', 0);


%%

% Make groups of years for the splitapply functions below
yrgrps = findgroups(syntime(:,1));

% First calculate the extents of the values that will be
% plotted in the probability density functions

if recdata
    
    mines.tdb = round(min([tmytable.TDB; RecTables.TDB; tdbsyn.col]));
    maxes.tdb = round(max([tmytable.TDB; RecTables.TDB; tdbsyn.col]));
    mines.ghi = round(min([tmytable.GHI; RecTables.GHI; ghisyn.col]));
    maxes.ghi = round(max([tmytable.GHI; RecTables.GHI; ghisyn.col]));
    mines.rh = round(min([tmytable.RH; RecTables.RH; rhsyn.col]));
    maxes.rh = round(max([tmytable.RH; RecTables.RH; rhsyn.col]));
    
else
    
    mines.tdb = round(min([tmytable.TDB; tdbsyn.col]));
    maxes.tdb = round(max([tmytable.TDB; tdbsyn.col]));
    mines.ghi = round(min([tmytable.GHI; ghisyn.col]));
    maxes.ghi = round(max([tmytable.GHI; ghisyn.col]));
    mines.rh = round(min([tmytable.RH; rhsyn.col]));
    maxes.rh = round(max([tmytable.RH; rhsyn.col]));
    
end

% Now create a function for each quantity
ecdffun.tdb = @(x) histcounts(x, mines.tdb:2:maxes.tdb, ...
	'Normalization', 'cdf');
ecdffun.ghi = @(x) histcounts(x, mines.ghi:25:maxes.ghi, ...
	'Normalization', 'cdf');
ecdffun.rh = @(x) histcounts(x, mines.rh:5:maxes.rh, ...
	'Normalization', 'cdf');

pdffun.tdb = @(x) histcounts(x, ...
	mines.tdb:2:maxes.tdb, 'Normalization', 'pdf');
pdffun.ghi = @(x) histcounts(x, ...
	mines.ghi:25:maxes.ghi, 'Normalization', 'pdf');
pdffun.rh = @(x) histcounts(x, ...
	mines.rh:5:maxes.rh, 'Normalization', 'pdf');

% Calculate the empiricial CDF of the TMY data
[eCDF.tmy.tdb.e, eCDF.tmy.tdb.x] = ...
	ecdffun.tdb(tmytable.TDB);
[eCDF.tmy.ghi.e, eCDF.tmy.ghi.x] = ...
	ecdffun.ghi(tmytable.GHI(tmytable.GHI>0));
[eCDF.tmy.rh.e, eCDF.tmy.rh.x] = ...
	ecdffun.rh(tmytable.RH);

if recdata
    % Calculate the empiricial CDF of the recorded data
    [eCDF.rec.tdb.e, eCDF.rec.tdb.x] = ...
        ecdffun.tdb(RecTables.TDB);
    [eCDF.rec.ghi.e, eCDF.rec.ghi.x] = ...
        ecdffun.ghi(GHIrec(GHIrec>0));
    [eCDF.rec.rh.e, eCDF.rec.rh.x] = ...
        ecdffun.rh(RecTables.RH);

    % Calculate the PDF of the recorded data
    [ePDF.rec.tdb.y, ePDF.rec.tdb.x] = ...
        pdffun.tdb(RecTables.TDB);
    [ePDF.rec.ghi.y, ePDF.rec.ghi.x] = ...
        pdffun.ghi(GHIrec(GHIrec>0));
    [ePDF.rec.rh.y, ePDF.rec.rh.x] = ...
        pdffun.rh(RecTables.RH);
end

% ECDF of the synthetic hours
[eCDF.syn.tdb.e,eCDF.syn.tdb.x] = ...
	splitapply(@(x) ecdffun.tdb(x),tdbsyn.col,yrgrps);
[eCDF.syn.ghi.e,eCDF.syn.ghi.x] = ...
	splitapply(@(x) ecdffun.ghi(x), ...
	ghisyn.col(ghisyn.col>0),yrgrps(ghisyn.col>0));
[eCDF.syn.rh.e,eCDF.syn.rh.x] = ...
	splitapply(@(x) ecdffun.rh(x),rhsyn.col,yrgrps);

% Calculate the PDF of the TMY data
[ePDF.tmy.tdb.y, ePDF.tmy.tdb.x] = ...
	pdffun.tdb(tmytable.TDB);
[ePDF.tmy.ghi.y, ePDF.tmy.ghi.x] = ...
	pdffun.ghi(tmytable.GHI(tmytable.GHI>0));
[ePDF.tmy.rh.y, ePDF.tmy.rh.x] = ...
	pdffun.rh(tmytable.RH);

% Plotting without zeros
[ePDF.syn.tdb.y, ePDF.syn.tdb.x] = ...
	splitapply(@(x) pdffun.tdb(x),tdbsyn.col,yrgrps);
[ePDF.syn.ghi.y, ePDF.syn.ghi.x] = ...
	splitapply(@(x) pdffun.ghi(x), ...
	ghisyn.col(ghisyn.col>0), yrgrps(ghisyn.col>0));
[ePDF.syn.rh.y, ePDF.syn.rh.x] = ...
	splitapply(@(x) pdffun.rh(x),rhsyn.col,yrgrps);


if ccdata

	% This bit of code 
uniqueyr_t = repmat((ones(8760,1)),1,85);
uniqueyr_t2 = [];
for k = 1:bootlen
	uniqueyr_t2 = [uniqueyr_t2, k*uniqueyr_t];
end
uniqueyr = reshape(uniqueyr_t2, [], 1);
clear uniqueyr_t uniqueyr_t2
FutureTimeRep = [repmat([FutureTime(:,1), ...
	FutureTime(:,2)],bootlen,1), ...
	uniqueyr];

resfunc = @(x) reshape(permute(reshape(x,50,12,[]), ...
	[2,1,3]),12,[]);

[mStats.rcp45.tdb.max, mStats.rcp45.tdb.min, ...
	mStats.rcp45.tdb.mean, mStats.rcp45.tdb.std, ...
	mStats.rcp45.tdb.skew, mStats.rcp45.tdb.kurt, ...
	mStats.rcp45.tdb.months] = ...
	grpstats(tdbsynCC.rcp45.col, ...
	FutureTimeRep, ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), ...
	@nanmean,@nanstd,@(x) (skewness(x, shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});
mStats.rcp45.tdb = structfun(@(x) resfunc(x), ...
	mStats.rcp45.tdb, 'UniformOutput', 0);

[mStats.rcp45.ghi.max, mStats.rcp45.ghi.min, ...
	mStats.rcp45.ghi.mean, mStats.rcp45.ghi.std, ...
	mStats.rcp45.ghi.skew, mStats.rcp45.ghi.kurt, ...
	mStats.rcp45.ghi.months] = grpstats( ...
	ghisynCC.rcp45.col(ghisynCC.rcp45.col>0), ...
	FutureTimeRep(ghisynCC.rcp45.col>0,:), ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), ...
	@expfit, @iqr, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});
mStats.rcp45.ghi = structfun(@(x) resfunc(x), ...
	mStats.rcp45.ghi, 'UniformOutput', 0);

[mStats.rcp45.rh.max, mStats.rcp45.rh.min, ...
	mStats.rcp45.rh.mean, ...
	mStats.rcp45.rh.std, mStats.rcp45.rh.skew, ...
	mStats.rcp45.rh.kurt, mStats.rcp45.rh.months] = ...
	grpstats(rhsynCC.rcp45.col, FutureTimeRep, ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), @nanmean, ...
	@nanstd, @(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});
mStats.rcp45.rh = structfun(@(x) resfunc(x), ...
	mStats.rcp45.rh, 'UniformOutput', 0);

[mStats.rcp85.tdb.max, mStats.rcp85.tdb.min, ...
	mStats.rcp85.tdb.mean, mStats.rcp85.tdb.std, ...
	mStats.rcp85.tdb.skew, mStats.rcp85.tdb.kurt, ...
	mStats.rcp85.tdb.months] = ...
	grpstats(tdbsynCC.rcp85.col, FutureTimeRep, ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), ...
	@nanmean,@nanstd,@(x) (skewness(x, shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});
mStats.rcp85.tdb = structfun(@(x) resfunc(x), ...
	mStats.rcp85.tdb, 'UniformOutput', 0);

[mStats.rcp85.ghi.max, mStats.rcp85.ghi.min, ...
	mStats.rcp85.ghi.mean, mStats.rcp85.ghi.std, ...
	mStats.rcp85.ghi.skew, mStats.rcp85.ghi.kurt, ...
	mStats.rcp85.ghi.months] = grpstats( ...
	ghisynCC.rcp85.col(ghisynCC.rcp85.col>0), ...
	FutureTimeRep(ghisynCC.rcp85.col>0,:), ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), ...
	@expfit, @iqr, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)), 'gname'});
mStats.rcp85.ghi = structfun(@(x) resfunc(x), ...
	mStats.rcp85.ghi, 'UniformOutput', 0);

[mStats.rcp85.rh.max, mStats.rcp85.rh.min, ...
	mStats.rcp85.rh.mean, ...
	mStats.rcp85.rh.std, mStats.rcp85.rh.skew, ...
	mStats.rcp85.rh.kurt, mStats.rcp85.rh.months] = ...
	grpstats(rhsynCC.rcp85.col, FutureTimeRep, ...
	{@(x) (quantile(x,0.99)), ...
	@(x) (quantile(x,0.01)), ...
	@nanmean, @nanstd, ...
	@(x) (skewness(x,shapeflag)), ...
	@(y) (kurtosis(y,shapeflag)),'gname'});
mStats.rcp85.rh = structfun(@(x) resfunc(x), ...
	mStats.rcp85.rh, 'UniformOutput', 0);

% Make groups of years for the splitapply functions
% below
fyrgrps = findgroups(FutureTimeRep(:,1));

% First calculate the extents of the values that will be
% plotted in the probability density functions
mines.tdb = round(min([tmytable.TDB; RecTables.TDB; ...
	tdbsynCC.rcp45.col; tdbsynCC.rcp85.col]));
maxes.tdb = round(max([tmytable.TDB; RecTables.TDB; ...
	tdbsynCC.rcp45.col; tdbsynCC.rcp85.col]));
mines.ghi = round(min([tmytable.GHI; RecTables.GHI; ...
	ghisynCC.rcp45.col; ghisynCC.rcp85.col]));
maxes.ghi = round(max([tmytable.GHI; RecTables.GHI; ...
	ghisynCC.rcp45.col; ghisynCC.rcp85.col]));
mines.rh = round(min([tmytable.RH; RecTables.RH; ...
	rhsynCC.rcp45.col; rhsynCC.rcp85.col]));
maxes.rh = round(max([tmytable.RH; RecTables.RH; ...
	rhsynCC.rcp45.col; rhsynCC.rcp85.col]));

% Now create a function for each quantity
ecdffun.tdb = @(x) histcounts(x, mines.tdb:2:maxes.tdb, ...
	'Normalization', 'cdf');
ecdffun.ghi = @(x) histcounts(x, mines.ghi:25:maxes.ghi, ...
	'Normalization', 'cdf');
ecdffun.rh = @(x) histcounts(x, mines.rh:5:maxes.rh, ...
	'Normalization', 'cdf');

pdffun.tdb = @(x) histcounts(x, ...
	mines.tdb:2:maxes.tdb, 'Normalization', 'pdf');
pdffun.ghi = @(x) histcounts(x, ...
	mines.ghi:25:maxes.ghi, 'Normalization', 'pdf');
pdffun.rh = @(x) histcounts(x, ...
	mines.rh:5:maxes.rh, 'Normalization', 'pdf');


% Calculate the empiricial CDF of the TMY data
[eCDF.tmy.tdb.e, eCDF.tmy.tdb.x] = ...
	ecdffun.tdb(tmytable.TDB);
[eCDF.tmy.ghi.e, eCDF.tmy.ghi.x] = ...
	ecdffun.ghi(tmytable.GHI(tmytable.GHI>0));
[eCDF.tmy.rh.e, eCDF.tmy.rh.x] = ...
	ecdffun.rh(tmytable.RH);

if recdata
    % Calculate the empiricial CDF of the recorded data
    [eCDF.rec.tdb.e, eCDF.rec.tdb.x] = ...
        ecdffun.tdb(RecTables.TDB);
    [eCDF.rec.ghi.e, eCDF.rec.ghi.x] = ...
        ecdffun.ghi(GHIrec(GHIrec>0));
    [eCDF.rec.rh.e, eCDF.rec.rh.x] = ...
        ecdffun.rh(RecTables.RH);
    
    % Calculate the PDF of the recorded data
    [ePDF.rec.tdb.y, ePDF.rec.tdb.x] = ...
        pdffun.tdb(RecTables.TDB);
    [ePDF.rec.ghi.y, ePDF.rec.ghi.x] = ...
        pdffun.ghi(GHIrec(GHIrec>0));
    [ePDF.rec.rh.y, ePDF.rec.rh.x] = ...
        pdffun.rh(RecTables.RH);
end

% Calculate the PDF of the TMY data
[ePDF.tmy.tdb.y, ePDF.tmy.tdb.x] = ...
	pdffun.tdb(tmytable.TDB);
[ePDF.tmy.ghi.y, ePDF.tmy.ghi.x] = ...
	pdffun.ghi(tmytable.GHI(tmytable.GHI>0));
[ePDF.tmy.rh.y, ePDF.tmy.rh.x] = ...
	pdffun.rh(tmytable.RH);


[eCDF.rcp45.tdb.e,eCDF.rcp45.tdb.x] = splitapply( ...
	@(x) ecdffun.tdb(x), tdbsynCC.rcp45.col, fyrgrps);
[eCDF.rcp45.ghi.e,eCDF.rcp45.ghi.x] = splitapply( ...
	@(x) ecdffun.ghi(x), ghisynCC.rcp45.col, fyrgrps);
[eCDF.rcp45.rh.e,eCDF.rcp45.rh.x] = splitapply( ...
	@(x) ecdffun.rh(x), rhsynCC.rcp45.col, fyrgrps);

[eCDF.rcp85.tdb.e, eCDF.rcp85.tdb.x] = splitapply( ...
	@(x) ecdffun.tdb(x), tdbsynCC.rcp85.col, fyrgrps);
[eCDF.rcp85.ghi.e, eCDF.rcp85.ghi.x] = splitapply( ...
	@(x) ecdffun.ghi(x), ghisynCC.rcp85.col, fyrgrps);
[eCDF.rcp85.rh.e, eCDF.rcp85.rh.x] = splitapply( ...
	@(x) ecdffun.rh(x), rhsynCC.rcp85.col, fyrgrps);

[ePDF.rcp45.tdb.y, ePDF.rcp45.tdb.x] = ...
	splitapply(@(x) pdffun.tdb(x),tdbsynCC.rcp45.col,fyrgrps);
[ePDF.rcp45.ghi.y, ePDF.rcp45.ghi.x] = ...
	splitapply(@(x) pdffun.ghi(x),ghisynCC.rcp45.col,fyrgrps);
[ePDF.rcp45.rh.y, ePDF.rcp45.rh.x] = ...
	splitapply(@(x) pdffun.rh(x),rhsynCC.rcp45.col,fyrgrps);

[ePDF.rcp85.tdb.y, ePDF.rcp85.tdb.x] = ...
	splitapply(@(x) pdffun.tdb(x),tdbsynCC.rcp85.col,fyrgrps);
[ePDF.rcp85.ghi.y, ePDF.rcp85.ghi.x] = ...
	splitapply(@(x) pdffun.ghi(x),ghisynCC.rcp85.col,fyrgrps);
[ePDF.rcp85.rh.y, ePDF.rcp85.rh.x] = ...
	splitapply(@(x) pdffun.rh(x),rhsynCC.rcp85.col,fyrgrps);


clear FutureMonthRep

end

% % Plot the summary statistics
clear ePDF eCDF
clear H1 H2 H3 BoxPlotRecMonths
clear nIDX nIDXmonthly nD
clear diffs diffcen

%%
% Calculate the correlation coefficients for synthetic data
[Corr.syn.r,Corr.syn.pr] = corr(tdbsyn.col,...
	[rhsyn.col,ghisyn.col], ...
	'type','Pearson', 'rows','complete');
[Corr.syn.rho,Corr.syn.prho] = corr(tdbsyn.col, ...
	[rhsyn.col,ghisyn.col], ...
	'type','Spearman', 'rows','complete');
Corr.syn.rho = [Corr.syn.rho NaN];
Corr.syn.r = [Corr.syn.r NaN];
Corr.syn.prho = [Corr.syn.prho NaN];
Corr.syn.pr = [Corr.syn.pr NaN];

if recdata
	
	% Real data needs to be cleaned for correlation
	% calculations
	NewDataTable = RecTables(~isnan(RecTables.TDB) & ...
		~isnan(RecTables.GHI) & ~isnan(RecTables.RH), :);
	
	% Calculate the correlation coefficients for
	% measured data
	[Corr.rec.r,Corr.rec.pr] = corr( ...
		NewDataTable.TDB, [NewDataTable.RH, ...
		NewDataTable.GHI, NewDataTable.TDP], ...
		'type','Pearson', 'rows','complete');
	[Corr.rec.rho,Corr.rec.prho] = corr( ...
		NewDataTable.TDB,[NewDataTable.RH, ...
		NewDataTable.GHI, NewDataTable.TDP], ...
		'type','Spearman', 'rows','complete');
	
	clear NewDataTable
else
	Corr.rec.r = zeros(size(Corr.syn.r));
	Corr.rec.pr = zeros(size(Corr.syn.pr));
	Corr.rec.rho = zeros(size(Corr.syn.rho));
	Corr.rec.prho = zeros(size(Corr.syn.prho));
end

if ccdata
	% Climate change data
	[Corr.rcp45.r,Corr.rcp45.pr] = corr( ...
		tdbsynCC.rcp45.col, ...
		[rhsynCC.rcp45.col,ghisynCC.rcp45.col], ...
		'type','Pearson', 'rows','complete');
	[Corr.rcp85.r,Corr.rcp85.pr] = corr( ...
		tdbsynCC.rcp85.col, ...
		[rhsynCC.rcp85.col,ghisynCC.rcp85.col], ...
		'type','Pearson', 'rows','complete');
	
	[Corr.rcp45.rho,Corr.rcp45.prho] = corr( ...
		tdbsynCC.rcp45.col, ...
		[rhsynCC.rcp45.col,ghisynCC.rcp45.col], ...
		'type','Spearman', 'rows','complete');
	[Corr.rcp85.rho,Corr.rcp85.prho] = corr( ...
		tdbsynCC.rcp85.col, ...
		[rhsynCC.rcp85.col,ghisynCC.rcp85.col], ...
		'type','Spearman', 'rows','complete');
else
	Corr.rcp45.r = zeros(size(Corr.syn.r));
	Corr.rcp85.r = zeros(size(Corr.syn.r));
	Corr.rcp45.pr = zeros(size(Corr.syn.pr));
	Corr.rcp85.pr = zeros(size(Corr.syn.pr));
	Corr.rcp45.rho = zeros(size(Corr.syn.rho));
	Corr.rcp85.rho = zeros(size(Corr.syn.rho));
	Corr.rcp45.prho = zeros(size(Corr.syn.prho));
	Corr.rcp85.prho = zeros(size(Corr.syn.prho));
end


FilePathCorr = fullfile(pathMATsave, [nameEPWfile, ...
	'Corr_Perc.txt']);
fIDcorr = fopen(FilePathCorr, 'w');

fprintf(fIDcorr, ['\\toprule\r\n', ...
	'& \\multicolumn{2}{|c|}{W} & ', ...
	'\\multicolumn{2}{|c|}{GHI} \\\\\r\n', ...
	'\\midrule\r\n']);
fprintf(fIDcorr, ['& $r$ & $\\rho$ & $r$ & $\\rho$ ', ...
	'\\\\ \r\n \\midrule\r\n']);

fprintf(fIDcorr, ['Original & %3.2g & %3.2g ', ...
	'& %3.2g & %3.2g \\\\\r\n'], ...
	Corr.rec.r(1), Corr.rec.rho(1), Corr.rec.r(2), ...
	Corr.rec.rho(2));

fprintf(fIDcorr, ['TMY & %3.2g & %3.2g ', ...
	'& %3.2g & %3.2g \\\\\r\n'], ...
	Corr.tmy.r(1), Corr.tmy.rho(1), Corr.tmy.r(2), ...
	Corr.tmy.rho(2));

fprintf(fIDcorr, ['Synthetic & %3.2g & %3.2g ', ...
	'& %3.2g & %3.2g \\\\\r\n'], ...
	nanmean(Corr.syn.r(:,1)), ...
	nanmean(Corr.syn.rho(:,1)), ...
	nanmean(Corr.syn.r(:,2)), ...
	nanmean(Corr.syn.rho(:,2)));

fprintf(fIDcorr, ['RCP4.5 & %3.2g & %3.2g ', ...
	'& %3.2g & %3.2g \\\\\r\n'], ...
	nanmean(Corr.rcp45.r(:,1)), ...
	nanmean(Corr.rcp45.rho(:,1)), ...
	nanmean(Corr.rcp45.r(:,2)), ...
	nanmean(Corr.rcp45.rho(:,2)));

fprintf(fIDcorr, ['RCP8.5 & %3.2g & %3.2g ', ...
	'& %3.2g & %3.2g \\\\\r\n', ...
	'\\bottomrule\r\n'], nanmean(Corr.rcp85.r(:,1)), ...
	nanmean(Corr.rcp85.rho(:,1)), ...
	nanmean(Corr.rcp85.r(:,2)), ...
	nanmean(Corr.rcp85.rho(:,2)));

% Print a few blank lines before the quantiles.
fprintf(fIDcorr, '\r\n\r\n\r\n');

%%

% This part writes the synthetic data in various formats

syndata.Year = syntime(:,1);
syndata.Month = syntime(:,2);
syndata.Day = syntime(:,3);
syndata.Hour = syntime(:,4);
syndata.tdb = tdbsyn.col;
syndata.rh = rhsyn.col;
syndata.ghi = ghisyn.col;
syndata.dni = dnisyn.col;
syndata.dhi = dhisyn.col;

synplottab = structfun(@(x)(reshape(x,N,nboot)), ...
	syndata, 'UniformOutput',0);
% Select a random synthetic series to make plots
synplottab = structfun(@(x)(x(:,randi([1,nboot]))), ...
	synplottab, 'UniformOutput',0);
RawPlots(syndata,figsavepath,nameEPWfile, ...
	'plotser', 'syn', 'visible','off');

% Save to a MAT file in the original EPW folder.
save(fullfile(pathEPWfolder, [nameEPWfile,'_Syn.mat']), ...
	'syndata')

% Save extras to a MAT file
save(fullfile(pathMATsave, ...
	[nameEPWfile,'_SynExtras.mat']), ...
	'tdbsyn', 'rhsyn', 'ghisyn', ...
	'dnisyn', 'dhisyn')

clear syndata


if ccdata
	
	% Write the climate change data as well
	syndata.Year = repmat(FutureTime(:,1),bootlen,1);
	syndata.Month = repmat(FutureTime(:,2),bootlen,1);
	syndata.Day = repmat(FutureTime(:,3),bootlen,1);
	syndata.Hour = repmat(FutureTime(:,4),bootlen,1);
	syndata.tdb = tdbsynCC.rcp45.col;
	syndata.rh = rhsynCC.rcp45.col;
	syndata.ghi = ghisynCC.rcp45.col;
	syndata.dni = dnisynCC.rcp45.col;
	syndata.dhi = dhisynCC.rcp45.col;
	
	% Save to a MAT file
	save(fullfile(pathMATsave, ...
		[nameEPWfile,'_rcp45.mat']), 'syndata')
	
	
	syndata.tdb = tdbsynCC.rcp85.col;
	syndata.rh = rhsynCC.rcp85.col;
	syndata.ghi = ghisynCC.rcp85.col;
	syndata.dni = dnisynCC.rcp85.col;
	syndata.dhi = dhisynCC.rcp85.col;
	
	% Save to a MAT file
	save(fullfile(pathMATsave, ...
		[nameEPWfile,'_rcp85.mat']), 'syndata')
	
	PlotHandle = figure('visible', 'off');
	figname = ['TDBsynAllCC85', '_', nameEPWfile];
	filepath = fullfile(figsavepath, figname);
	h1 = plot(tdbsynCC.rcp85.yearly);
	for k = 1:length(h1)
		h1(k).Marker = '.';
		h1(k).Color = lgrey;
		h1(k).LineStyle = 'none';
	end
	ax = gca; hold(ax, 'on'); ax.Box = 'on';
	h2 = plot(tmytable.TDB);
	h2.LineWidth = 1; h2.Color = red;
	ax.YLabel.Interpreter = 'latex';
	ax.YLabel.String = 'Temperature [\textsuperscript{o}C]';
	ax.FontSize = 22;
	ax.LabelFontSizeMultiplier = 1.25;
	ax.TitleFontSizeMultiplier = 1.25;
	ax.XLim = [0,N];
	ax.XTick = 360:720:8760-360;
	ax.XTickLabel = month(mt,'short');
	ax.XTickLabelRotation = 90;
	leg = legend([h1(1),h2],{'RCP 8.5', 'Typical'}, ...
		'location','south');
	leg.Orientation = 'horizontal';
	leg.FontSize = ax.FontSize;
	SaveThatFig(PlotHandle,filepath,'printpdf',false)
	
	h2.Color = grey;
	
	SaveThatFig(PlotHandle,[filepath,'-BW'], ...
		'printpdf', false, 'printfig',false)
	
	PlotHandle = figure('visible', 'off');
	figname = ['TDBsynAllCC45', '_', nameEPWfile];
	filepath = fullfile(figsavepath, figname);
	h1 = plot(tdbsyn.yearly);
	for k = 1:length(h1)
		h1(k).Marker = '.';
		h1(k).Color = lgrey;
		h1(k).LineStyle = 'none';
	end
	ax = gca; hold(ax, 'on'); ax.Box = 'on';
	h2 = plot(tmytable.TDB);
	h2.LineWidth = 1; h2.Color = orange;
	ax.YLabel.Interpreter = 'latex';
	ax.YLabel.String = 'Temperature [\textsuperscript{o}C]';
	ax.FontSize = 22;
	ax.LabelFontSizeMultiplier = 1.25;
	ax.TitleFontSizeMultiplier = 1.25;
	ax.XLim = [0,N];
	ax.XTick = 360:720:8760-360;
	ax.XTickLabel = month(mt,'short');
	ax.XTickLabelRotation = 90;
	leg = legend([h1(1),h2],{'RCP 4.5', 'Typical'}, ...
		'location','south');
	leg.Orientation = 'horizontal';
	leg.FontSize = ax.FontSize;
	SaveThatFig(PlotHandle,filepath,'printpdf',false)
	
	h2.Color = grey;
	
	SaveThatFig(PlotHandle,[filepath,'-BW'], ...
		'printpdf', false, 'printfig',false)
	
	clear dnisynCC dhisynCC rhsynCC
	clear syndata
end

clear rhsynCC rhsyn

%%


% Magnano et al and Hansen and Driscoll suggest checking for
% episodes (i.e. heat waves or cold snaps). An episode would
% be defined as the number of M consecutive hours occur
% above or below a specified threshold. Magnano et al choose
% 18�, 25�, 34� as their boundaries. This works well for
% their specific application, however we think it might be
% better to work with percentiles for generalisation. For
% example, in the TMY methodology, Wilcox et al use the 67th
% and 33rd percentile for persistence of warm and cold
% spells respectively. Another difference is that while
% Magnano et al look at consecutive hours, Wilcox et al look
% at daily means. Hansen and Driscoll use 90�, 65�, 32�
% (Fahrenheit)

% We take a hybrid approach for the 'spell length' and a
% percentile approach for the temperature boundaries.

% Using ASHRAE DESIGN temperatures is pointless for spell
% length because it never exceeds 2.

% The ashrae design temperatures expressed as quantiles are
quantsASHRAE = 1 - ([0.4, 1.0, 2.0, 50, 98.0, ...
	99.0, 99.6])'/100;

% The data is not stationary, as has been previously shown
% and can be visualised easily. However, mean values are
% often used for meteorological parameters. They serve
% little purpose for design but at least they provide a way
% of classification. In this case, they provide some idea of
% the 'weight' of the distribution. The first and second
% moments are a measure of how much the data is dispersed.
% If we are looking for spells, then we can think of
% capturing spells that are far away from the 'mass' centre
% of the data.

CoverageCheck = @(TS,thresh,func) ...
	(nansum(eval(sprintf('%s(TS,thresh)', ...
	func))) / (numel(TS)));

% Find the values that are higher than the TMY 2 percentile
Coverage.syn.tdb.ge020 = CoverageCheck(tdbsyn.col, ...
	StationInfo.Cool.TDB020, 'ge');
Coverage.tmy.tdb.ge020 = CoverageCheck(tmytable.TDB, ...
	StationInfo.Cool.TDB020, 'ge');

% Find the values that are higher than the TMY 1 percentile
Coverage.syn.tdb.ge010 = CoverageCheck(tdbsyn.col, ...
	StationInfo.Cool.TDB010, 'ge');
Coverage.tmy.tdb.ge010 = CoverageCheck(tmytable.TDB, ...
	StationInfo.Cool.TDB010, 'ge');

% Find the values that are higher than the TMY 0.4
% percentile
Coverage.syn.tdb.ge004 = CoverageCheck(tdbsyn.col, ...
	StationInfo.Cool.TDB004, 'ge');
Coverage.tmy.tdb.ge004 = CoverageCheck(tmytable.TDB, ...
	StationInfo.Cool.TDB004, 'ge');


% Repeat procedure for values that are below the 99.6
% percentile
Coverage.syn.tdb.le996 = 1-CoverageCheck(tdbsyn.col, ...
	StationInfo.Heat.TDB996, 'le');
Coverage.tmy.tdb.le996 = 1-CoverageCheck(tmytable.TDB, ...
	StationInfo.Heat.TDB996, 'le');

% Repeat procedure for values that are below the 99.0
% percentile
Coverage.syn.tdb.le990 = 1-CoverageCheck(tdbsyn.col, ...
	StationInfo.Heat.TDB990, 'le');
Coverage.tmy.tdb.le990 = 1-CoverageCheck(tmytable.TDB, ...
	StationInfo.Heat.TDB990, 'le');

% Find the quantiles, the coverage being assumed on the
% basis of an underlying normal distribution
Quants.syn.tdb = quantile(tdbsyn.col, quantsASHRAE);

if recdata
	Coverage.rec.tdb.ge020 = CoverageCheck( ...
		RecTables.TDB, StationInfo.Cool.TDB020, 'ge');
	Coverage.rec.tdb.ge010 = CoverageCheck( ...
		RecTables.TDB, StationInfo.Cool.TDB010, 'ge');
	Coverage.rec.tdb.ge004 = CoverageCheck( ...
		RecTables.TDB, StationInfo.Cool.TDB004, 'ge');
	Coverage.rec.tdb.le996 = 1-CoverageCheck( ...
		RecTables.TDB, StationInfo.Heat.TDB996, 'le');
	Coverage.rec.tdb.le990 = 1-CoverageCheck( ...
		RecTables.TDB, StationInfo.Heat.TDB990, 'le');
	Quants.rec.tdb = quantile(RecTables.TDB, quantsASHRAE);
end

save(fullfile(pathMATsave,['Coverage_', ...
	nameEPWfile,'.mat']), 'Coverage')

Quants.tmy.tdb = [StationInfo.Cool.TDB004; ...
	StationInfo.Cool.TDB010; ...
	StationInfo.Cool.TDB020; nanmean(tmytable.TDB); ...
	quantile(tmytable.TDB,0.02); ...
	StationInfo.Heat.TDB990; StationInfo.Heat.TDB996];

if ccdata
	Quants.rcp45.tdb = quantile(tdbsynCC.rcp45.col, ...
		quantsASHRAE);
	Quants.rcp85.tdb = quantile(tdbsynCC.rcp85.col, ...
		quantsASHRAE);
end

fprintf(fIDcorr, ['\\toprule\r\n ', ...
	'Percentile & \\multicolumn{5}{c}', ...
	'{Temperatures (\\celsius)} \\\\ \r\n %s & %s ', ...
	'& %s & %s & %s & %s \\midrule \r\n'], '(\\%)', ...
	'Recorded', 'TMY','Synthetic','RCP4.5', 'RCP8.5');

for k = 1:length(Quants.syn.tdb)
	
	if recdata
		QuantRec = Quants.rec.tdb(k);
	else
		QuantRec = 0;
	end
	
	if ccdata
		quant45 = Quants.rcp45.tdb(k);
		quant85 = Quants.rcp85.tdb(k);
	else
		quant45 = 0;
		quant85 = 0;
	end
	
	fprintf(fIDcorr, [ '%3.1f & %3.2f & %3.2f ', ...
		'& %3.2f & %3.2f & %3.2f \\\\\r\n'], ...
		quantsASHRAE(k)*100, Quants.syn.tdb(k), ...
		QuantRec, Quants.tmy.tdb(k), ...
		quant45, quant85);
end
fprintf(fIDcorr, '\\bottomrule\r\n');

fclose(fIDcorr);


if recdata
QuantsSpell = quantile(RecTables.TDB, ...
	[quantsASHRAE(1:3); 0.67; 0.5; ...
	0.33; quantsASHRAE(5:7)]);

for qs = 1:length(QuantsSpell)

	if ismember(qs,[1,2,3,4])
		% These are the 'greater than' quantiles. That
		% is, a heat wave is important if it exceeds
		% these quantiles.
		[dur.syn, durtime.syn] = EpisodeFinder( ...
			tdbsyn.col, ...
			QuantsSpell(qs), 'ge');

		[dur.rec, durtime.rec] = EpisodeFinder( ...
			RecTables.TDB, QuantsSpell(qs), 'ge');

		if ccdata
			[dur.rcp45, durtime.rcp45] = ...
				EpisodeFinder( ...
				tdbsynCC.rcp45.col, ...
				QuantsSpell(qs), 'ge');
			[dur.rcp85, durtime.rcp85] = ...
				EpisodeFinder( ...
				tdbsynCC.rcp85.col, ...
				QuantsSpell(qs), 'ge');
		end

	elseif ismember(qs,[5,6,7,8,9])
		% These are the 'less than' quantiles. That is,
		% a cold snap is important if it goes below
		% these quantiles.
		[dur.syn, durtime.syn] = EpisodeFinder( ...
			tdbsyn.col, ...
			QuantsSpell(qs), 'le');


		[dur.rec, durtime.rec] = EpisodeFinder( ...
			RecTables.TDB, QuantsSpell(qs), 'le');


		if ccdata
			[dur.rcp45, durtime.rcp45] = ...
				EpisodeFinder( ...
				tdbsynCC.rcp45.col, ...
				QuantsSpell(qs), 'le');
			[dur.rcp85, durtime.rcp85] = ...
				EpisodeFinder( ...
				tdbsynCC.rcp85.col, ...
				QuantsSpell(qs), 'le');
		end
	end

	Eps.syn.Dur(qs) = {dur.syn};
	Eps.syn.Times(qs) = {durtime.syn};
	if ~isempty(dur.syn)
		[spellCDF.syn.tdb.e,spellCDF.syn.tdb.x] =  ...
			histcounts(dur.syn, 1000, ...
			'Normalization', 'cdf');
	else
		spellCDF.syn.tdb.e = NaN;
		spellCDF.syn.tdb.x = NaN;
		% 			spellCDF.syn.tdb.flo = NaN;
		% 			spellCDF.syn.tdb.fup = NaN;
	end

	Eps.rec.Dur(qs) = {dur.rec};
	Eps.rec.Times(qs) = {durtime.rec};

	if ~isempty(dur.rec)
		[spellCDF.rec.tdb.e,spellCDF.rec.tdb.x] = ...
			histcounts(dur.rec, 1000, ...
			'Normalization', 'cdf');
	else
		spellCDF.rec.tdb.e = NaN;
		spellCDF.rec.tdb.x = NaN;
	end

	if ismember(qs,[1,2,3,4])
		titleCustom = sprintf( ...
			'Hours above %3.1f ^oC', QuantsSpell(qs));
	elseif ismember(qs,[5,6,7,8,9])
		titleCustom = sprintf( ...
			'Hours below %3.1f ^oC', QuantsSpell(qs));
	end
	figname = sprintf('TDBspellcdf%d',qs);
	filepath = fullfile(figsavepath, ...
		[figname,'_',nameEPWfile]);
	CustomEcdfPlot(spellCDF,'tdb',filepath, ...
		'visible','off', ...
		'xLabelCustom', 'Spell Duration [hrs]', ....
		'plotCI', false, ...
		'ccdata', false, 'titleCustom', titleCustom)

	recempty = cellfun(@isempty, Eps.rec.Dur);
	durlimits.max.rec = cellfun(@max, ...
		Eps.rec.Dur(~recempty));
	durlimits.min.rec = cellfun(@min, ...
		Eps.rec.Dur(~recempty));


	if ccdata
		Eps.rcp45.Dur(qs) = {dur.rcp45};
		Eps.rcp45.Times(qs) = {durtime.rcp45};
		Eps.rcp85.Dur(qs) = {dur.rcp85};
		Eps.rcp85.Times(qs) = {durtime.rcp85};

		[spellCDF.rcp45.tdb.e,spellCDF.rcp45.tdb.x] = ...
			histcounts(dur.rcp45, 1000, ...
			'Normalization', 'cdf');
		[spellCDF.rcp85.tdb.e,spellCDF.rcp85.tdb.x] = ...
			histcounts(dur.rcp85, 1000, ...
			'Normalization', 'cdf');
		figname = sprintf('TDBspellcdf_CC%d',qs);

		filepath = fullfile(figsavepath, ...
			[figname,'_',nameEPWfile]);
		CustomEcdfPlot(spellCDF, 'tdb', ...
			filepath, 'visible', 'off', ...
			'xLabelCustom', 'Spell Duration [hrs]', ...
			'plotCI', false, ...
			'ccdata', true, ...
			'titleCustom', titleCustom)

		durlimits.max.rcp45 = cellfun(@max, ...
			Eps.rcp45.Dur);
		durlimits.max.rcp85 = cellfun(@max, ...
			Eps.rcp85.Dur);
		durlimits.min.rcp45 = cellfun(@min, ...
			Eps.rcp45.Dur);
		durlimits.min.rcp85 = cellfun(@min, ...
			Eps.rcp85.Dur);
	end

	% Don't plot the values corresponding to the median
	if ismember(qs,[4,5,6])
		continue
	end

end

synempty = cellfun(@isempty, Eps.syn.Dur);
durlimits.max.syn = cellfun(@max, ...
	Eps.syn.Dur(~synempty));
durlimits.min.syn = cellfun(@min, ...
	Eps.syn.Dur(~synempty));

end

%%
clear RecTables

close all