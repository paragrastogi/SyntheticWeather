% Function to concatenate climate change prediction files from George
function OutFilePath = CatClimChangeData(CCdataFolderPath, varargin)

p = inputParser;
p.FunctionName = 'CCdataReplace';

addRequired(p, 'CCdataFolderPath', @ischar)
% addParameter(p, 'Plotting', false, @islogical)
addParameter(p, 'FutureYearLen', 95, @isnumeric)
% This is the length of the future record, to be changed based on the
% incoming climate change data.

parse(p, CCdataFolderPath, varargin{:})

CCdataFolderPath = p.Results.CCdataFolderPath;
% This turns plotting of monthly future temperature histograms on or off.
% Plotting = p.Results.Plotting;
FutureYearLen = p.Results.FutureYearLen;


% CCdataFolderPath = 'D:\Weather_Data\ClimateChangeData\GEN';

ModelFolders = dir(CCdataFolderPath);
ModelFolders = ModelFolders ([ModelFolders.isdir]);
ModelFolders = {ModelFolders(cellfun(@isempty, ...
    strfind({ModelFolders.name},'.'))).name};

% Number of models available
NumModels = length(ModelFolders);

% The data are labelled as follows:
% Temperature (tas)
% Minimum daily temperature (tasmin)
% Maximum daily temperature (tasmax)
% Atmospheric pressure (ps)
% Surface wind speed (sfcWind)
% Global horizontal radiation (rsds)
% Specific humidity (huss)

% Each parameter has 3 files dedicated to it. Which means that each folder
% (for each model) will contain 7*3 files.
NumParam = 7;
NumFilesTotal = NumParam*3;

% Make a table to store the climate change data for this folder
CCdata = table();
CCdata.Model = cell(length(ModelFolders)*length(NumFilesTotal),1);
CCdata.SubModel = cell(length(ModelFolders)*length(NumFilesTotal),1);
CCdata.Parameter = cell(length(ModelFolders)*length(NumFilesTotal),1);
CCdata.Data = cell(length(ModelFolders)*length(NumFilesTotal),1);
counter = 1;

for m = 1:NumModels
    
    % Get all the csv files
    FilesList = dir(fullfile(CCdataFolderPath,ModelFolders{m}, '*.csv'));
    FilesList = {FilesList.name};
    
    for f = 1: length(FilesList)
        
        SplitFileName = strsplit(FilesList{f}, {'_', '.'});
        
        % This is the current parameter being read in
        Param = SplitFileName{2};
        if strcmp(Param, 'tas')
            Param = 'TDBdmean';
        elseif strcmp(Param, 'tasmin')
            Param = 'TDBdmin';
        elseif strcmp(Param, 'tasmax')
            Param = 'TDBdmax';
        elseif strcmp(Param, 'ps')
            Param = 'ATMPRdmean';
        elseif strcmp(Param, 'sfcWind')
            Param = 'Wspd_dmean';
        elseif strcmp(Param, 'rsds')
            Param = 'GHIdmean';
        elseif strcmp(Param, 'huss')
            Param = 'Wdmean';
        end
        
        % Period
        SubModel = SplitFileName{3};
        
        ReadData = csvread(fullfile(CCdataFolderPath, ModelFolders{m}, ...
            FilesList{f}));
        
        % Convert specific humidity to humidity ratio
        if strcmp(Param,'Wdmean')
            Wd = -ReadData./(ReadData - 1);
			ReadData = Wd;
        end
        
        if strcmp(SubModel,'historical')
            ReadData = [(1951:2005)', ReadData];
        else
            ReadData = [(2006:2100)', ReadData];
        end
        
        
        CCdata.Model(counter) = ModelFolders(m);
        CCdata.SubModel(counter) = {SubModel};
        CCdata.Parameter(counter) = {Param};
        CCdata.Data(counter) = {ReadData};
        
        counter = counter + 1;
    end
    
end

save(fullfile(CCdataFolderPath,'CollClim.mat'),'CCdata')

UniqueParams = unique(CCdata.Parameter);

temp.rcp85 = NaN(FutureYearLen,365,NumModels);

for p = 1:length(UniqueParams)
	
	% Find the future data for the current parameter
	Tempidx.rcp45 = ( strcmpi(CCdata.Parameter, ...
		UniqueParams{p}) & ...
		strcmpi(CCdata.SubModel,'rcp45') );	
	Tempidx.rcp85 = ( strcmpi(CCdata.Parameter, ...
		UniqueParams{p}) & ...
		strcmpi(CCdata.SubModel,'rcp85') );
	
	% Collect all the future temperature data
	temp.rcp85 = [CCdata.Data{Tempidx.rcp85}];
	temp.rcp45 = [CCdata.Data{Tempidx.rcp45}];
	
	FutureWeather.(UniqueParams{p}).rcp85 = ...
		NaN(FutureYearLen,8760,NumModels); 
	FutureWeather.(UniqueParams{p}).rcp45 = ...
		NaN(FutureYearLen,8760,NumModels);
	
	for r = 1:NumModels
		tempR = temp.rcp85(:,((r-1)*366)+1:(r*366));
		tempR = repmat(tempR(:,2:end),1,1,24);
		tempR2 = reshape(permute(tempR, [1 3 2]), ...
			FutureYearLen,8760); 
		FutureWeather.(UniqueParams{p}).rcp85(:,:,r) = ...
			tempR2;
		
		tempR = temp.rcp45(:,((r-1)*366)+1:(r*366));
		tempR = repmat(tempR(:,2:end),1,1,24); 
		tempR2 = reshape(permute(tempR, [1 3 2]), ...
			FutureYearLen,8760); 
		FutureWeather.(UniqueParams{p}).rcp45(:,:,r) = ...
			tempR2;
	end
		
	% The rows represent years (2006-2100), the columns
	% represent hours (1-8760), and the third dimension
	% represents each model-submodel combination. The daily
	% values are repeated 24 times to create hourly values.
	
	clear temp tempR tempR2
end

FutureYears = 2006:2100;

% % Create an year-long datetime vector
TimeSyn = (datetime(FutureYears(1),1,1,1,0,0) : hours(1) : ...
	datetime(FutureYears(1),12,31,24,0,0))';
% The increment is one hour.

FutureTime = [reshape((repmat(FutureYears, ...
	length(TimeSyn),1)),[],1), ...
	repmat(month(TimeSyn),FutureYearLen,1), ...
	repmat(day(TimeSyn),FutureYearLen,1), ...
	repmat(hour(TimeSyn),FutureYearLen,1)];

% Keep only those years that are the future from now, 2015.
YearCut = FutureYears>2015;

FutureTime = FutureTime(FutureTime(:,1)>2015,:);

for p = 1:length(UniqueParams)
FutureWeather.(UniqueParams{p}) = structfun(@(x) (x(YearCut,:,:)), ...
	FutureWeather.(UniqueParams{p}), 'UniformOutput',0);
end

% Save the results
OutFilePath = fullfile(CCdataFolderPath,'HourlyFutureData.mat');
save(OutFilePath, 'FutureWeather','FutureTime')

end