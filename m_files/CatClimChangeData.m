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
ModelFolders = {ModelFolders(~contains({ModelFolders.name}, '.')).name};

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

scenarios = {'rcp45', 'rcp85'};


for m = 1:NumModels
    
    % Get all the csv files
    FilesList = dir(fullfile(CCdataFolderPath,ModelFolders{m}, '*.csv'));
    FilesList = {FilesList.name};
    
    for f = 1: length(FilesList)
        
        SplitFileName = strsplit(FilesList{f}, {'_', '.'});
        
        % This is the current parameter being read in
        Param = SplitFileName{1};
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
        SubModel = SplitFileName{2};
        
        try
            ReadData = csvread(fullfile(CCdataFolderPath, ModelFolders{m}, ...
                FilesList{f}));
        catch err
            fprintf('Could not read file, probably empty\r\n')
            continue
        end
        
        % Convert specific humidity to humidity ratio
        if strcmp(Param,'Wdmean')
            ReadData(:, end) = -ReadData(:, end)./(ReadData(:, end) - 1);
        end
        
        %         if strcmp(SubModel,'historical')
        %             ReadData = [(1951:2005)', ReadData];
        %         else
        %             ReadData = [(2006:2100)', ReadData];
        %         end
        
        
        CCdata.Model(counter) = ModelFolders(m);
        CCdata.SubModel(counter) = {SubModel};
        CCdata.Parameter(counter) = {Param};
        CCdata.Data(counter) = {ReadData};
        
        counter = counter + 1;
    end
    
end

save(fullfile(CCdataFolderPath,'CollClim.mat'),'CCdata')

UniqueParams = unique(CCdata.Parameter);


for c = 1:length(scenarios)
    
    temp.(scenarios{c}) = NaN(FutureYearLen,365,NumModels);
    
    for p = 1:length(UniqueParams)
        
        % Find the future data for the current parameter
        Tempidx.rcp45 = ( strcmpi(CCdata.Parameter, ...
            UniqueParams{p}) & ...
            strcmpi(CCdata.SubModel,'rcp45') );
        Tempidx.rcp85 = ( strcmpi(CCdata.Parameter, ...
            UniqueParams{p}) & ...
            strcmpi(CCdata.SubModel,'rcp85') );
        
        % Collect all the future temperature data
        temp.(scenarios{c}) = CCdata.Data{Tempidx.rcp85};
        % 	temp.rcp45 = CCdata.Data{Tempidx.rcp45};
        
        FutureWeather.(UniqueParams{p}).(scenarios{c}) = ...
            NaN(FutureYearLen,8760,NumModels);
        %         FutureWeather.(UniqueParams{p}).rcp45 = ...
        %             NaN(FutureYearLen,8760,NumModels);
        
        for r = 1:NumModels
            
            % Find this year of data, for this parameter, submodel,
            % scenario, and model.
            findyear = strcmpi(CCdata.Parameter, UniqueParams{p}) & ...
                strcmpi(CCdata.SubModel, scenarios{c}) & ...
                strcmpi(CCdata.Model, ModelFolders{r});
            temp = CCdata.Data(findyear);
            if isempty(temp)
                % This combination of variables does not exist.
                continue
            end
            temp = temp{1};
            
            % Future years can vary a bit, so if they're 96, the sizes of
            % most variables will increase.
            FutureYears = unique(temp(:,1));
            
            if length(FutureYears)>FutureYearLen
                FutureYearLen = length(FutureYears);
            end
            
            for y = 1:length(FutureYears)
                % For each year.
                
                % Reshape that year's data into the the required shape of
                % y x H x m (year, hours, model).
                temp2 = ...
                    reshape(permute(repmat(temp(temp(:,1)== ...
                    FutureYears(y),:), 1, 1, 24), [1 3 2]), [], 4);
                % 'Standard' hours in this year (all 8760 hours).
                stdhours = datetime(temp2(1,1), 1, 1, 0, 0, 0): ...
                    hours(1):datetime(temp2(1,1), 12, 31, 23, 0, 0);
                stdhours = stdhours(:);
                
                % Sort the vector by month and then day.
                temp2 = sortrows(temp2, [2,3]);
                
                % Actual number of hours that exist for this year.
                acthours = datetime(temp2(:,1), temp2(:,2), ...
                    temp2(:,3), reshape(repmat(0:23, 1, length(unique(temp2(:,2:3), 'rows'))), [], 1), ...
                    zeros(size(temp2,1), 1), zeros(size(temp2,1), 1)) ;
                
                % Eliminate leap years - sorry!
                if any(month(stdhours)==2 & day(stdhours)==29)
                    stdhours(month(stdhours)==2 & day(stdhours)==29) = [];
                end
                
                % Logical index of actual hours that overlap with the
                % expected 'standard' hours.
                hourcommon = ismember(stdhours, acthours);
                
                FutureWeather.(UniqueParams{p}).(scenarios{c})(y,hourcommon,r) = (temp2(:,end))';
            end
        end
        
        % The rows represent years (2006-2100), the columns
        % represent hours (1-8760), and the third dimension
        % represents each model-submodel combination. The daily
        % values are repeated 24 times to create hourly values.
        
        clear temp tempR tempR2
    end
    
    
    
    % % Create an year-long datetime vector
    TimeSyn = (datetime(FutureYears(1),1,1,0,0,0) : hours(1) : ...
        datetime(FutureYears(1),12,31,23,0,0))';
    % The increment is one hour.
    
    FutureTime = [reshape((repmat(FutureYears, ...
        length(TimeSyn),1)),[],1), ...
        repmat(month(TimeSyn),FutureYearLen,1), ...
        repmat(day(TimeSyn),FutureYearLen,1), ...
        repmat(hour(TimeSyn),FutureYearLen,1)];
    
    % % Keep only those years that are the future from now, 2015.
    % YearCut = FutureYears>2015;
    
    % FutureTime = FutureTime(FutureTime(:,1)>2015,:);
    
    % for p = 1:length(UniqueParams)
    % FutureWeather.(UniqueParams{p}) = structfun(@(x) (x(YearCut,:,:)), ...
    % 	FutureWeather.(UniqueParams{p}), 'UniformOutput',0);
    % end
    
    % Save the results
    OutFilePath = fullfile(CCdataFolderPath,'HourlyFutureData.mat');
    save(OutFilePath, 'FutureWeather','FutureTime')
    
end

end
