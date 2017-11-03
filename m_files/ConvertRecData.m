% % This function takes the recorded data for a given location and saves it
% as a MAT file.
function ConvertRecData(EPWfilePath, MATsavePath, StationID, varargin)

p = inputParser;
p.FunctionName = 'ConvertRecData';

% addRequired(p,'EPWfolder',@ischar)
addRequired(p,'EPWfilePath',@ischar)
addRequired(p,'MATsavePath',@ischar)
addRequired(p,'StationID',@ischar)
addParameter(p,'destroyer',false,@islogical)
addParameter(p,'truncate',false,@islogical)
addParameter(p,'amyfmt','',@ischar)

parse(p, EPWfilePath, MATsavePath, StationID, varargin{:})

EPWfilePath = p.Results.EPWfilePath;
MATsavePath = p.Results.MATsavePath;
StationID = p.Results.StationID;
destroyer = p.Results.destroyer;
truncate = p.Results.truncate;
amyfmt = p.Results.amyfmt;

[EPWfolderPath,EPWfileName,~] = fileparts(EPWfilePath);

EPWfolder = EPWfolderPath(end-2:end);

if truncate
    MATfilePath = fullfile(MATsavePath,['TruncatedRecordedData_', ...
        StationID, '.mat']);
else
    MATfilePath = fullfile(MATsavePath,['RecordedData_', ...
        StationID, '.mat']);
end

if ~destroyer
    if exist(MATfilePath, 'file')==2
        fprintf(['Skipping making a new MAT file for incoming ', ...
            'EPW file %s .\r\n'], EPWfileName)
        return
    end
end

FileParserNames = FileParserKeywordLister;


if isempty(amyfmt)
    % No format specified for files, try all.
AMYfiles = [dir(fullfile(EPWfolderPath,'*.txt')); ...
    dir(fullfile(EPWfolderPath,'*.csv')); ...
    dir(fullfile(EPWfolderPath,'*.wy2'))];
else
    AMYfiles = dir(fullfile(EPWfolderPath,['*.', amyfmt]));
end

AMYfiles = {AMYfiles.name};

% Remove the synthetic time series saved in the same folder.
AMYfiles = AMYfiles(cellfun(@isempty, strfind(AMYfiles, 'syn')) & ...
    cellfun(@isempty, strfind(AMYfiles, 'Syn')));

% Exclude any files explicitly labelled as TMY.
AMYfiles = AMYfiles(cellfun(@isempty, strfind(AMYfiles, 'TMY')));


% Also EXCLUDE SODA FILES FOR NOW
AMYfiles = AMYfiles(cellfun(@isempty, strfind(AMYfiles, 'SODA')));

% Keep only the files from the CURRENT STATION
AMYfiles = AMYfiles(~cellfun(@isempty, ...
    cellfun(@(x) (strfind(x,StationID)), AMYfiles,'UniformOutput',0)));

if isempty(AMYfiles)
    fprintf(['No Recorded Data found for current station. ', ...
        'Quitting... (Note that SODA files have been ', ...
        'omitted so far.)\r\n'])
	return
end

% fIDlog = fopen([EPWfileName,'log.txt'],'w');

for a = 1:length(AMYfiles)
    
    RecDataFile = AMYfiles{a};
    RecDataFilePath = fullfile(EPWfolderPath,RecDataFile);
    
    % Split the current file name into its consituents
    SplitName = strsplit(RecDataFile(1:end-4),'_');
    % Only known sources are in AMYKeyword field
    Source = intersect(FileParserNames.AMYKeyword, ...
        SplitName, 'stable');
    if isempty(Source)
        % If a valid source is not found
        fprintf(['Unknown file source in file %s.', ...
            'Trying to read as EPW...\n'], RecDataFile);
        try 
            TempData = WeatherFileParseEPWUSDOE(RecDataFilePath);
        catch
            fprintf('Couldn''t read as EPW either, skipping file...\n');
            continue
        end
    else
        % Make the cell into a string.
        Source = Source{1};
        
        % If a valid source is found
        fprintf('Calling %s function on current file %s\n', ...
            ['WeatherFileParse',Source], RecDataFile );
        [ScreenOutput, TempData] = evalc(['feval(matlab',...
            '.lang.makeValidName([''WeatherFileParse'',', ...
            'Source]), RecDataFilePath)']);
        
        fprintf(ScreenOutput);
                
    end
    
    if a == 1
        % The data table needs to be created first
        RecTables = TempData;
        
    else
        % Joining new data with existing table
        RecTablesOld = RecTables;
        
        RecTables = outerjoin(RecTablesOld, ...
            TempData,'MergeKeys',true);
        
        clear RecTablesOld
        % The mergekeys option ensures that the date stamp is merged.
        % Duplicates are resolved below.
        
    end
end

clear TempData

DateColumnNames = {'Year','Month','Day','Hour','Minute'};


if truncate
    % % Find the first and last valid points in the time series of
    % % GHI and TDB
    
    FirstValidPoints.GHI = find(~isnan(RecTables.GHI),1);
    FirstValidPoints.TDB = find(~isnan(RecTables.TDB),1);
    %         FirstValidPoints(3) = find(~isnan(DataTablesInd.RH),1);
    LastValidPoints.GHI = find(~isnan(RecTables.GHI),1,'last');
    LastValidPoints.TDB = find(~isnan(RecTables.TDB),1,'last');
    %         LastValidPoints(3) = find(~isnan(DataTablesInd.RH),1,'last');
    
    % Take the last FIRST valid point and FIRST last valid point.
    PointsToCut = [structfun(@max,FirstValidPoints), ...
        structfun(@min,LastValidPoints)];
    
    % If the function call asks for truncating between the first and last
    % valid points of the GHI and TDB series, keep data table only between
    % these points.
    RecTables = RecTables(PointsToCut(1):PointsToCut(2),:);
end


% Following EPlus convention, the wind direction is 0� only when
% the data says so. If there is no wind, i.e. WSPD == 0, then wind
% direction is 180�. Wind direction 360� is due north, reset that
% to 0� to avoid unnecessary fractions in direction.
if ismember('WSPD',RecTables.Properties.VariableNames)
    RecTables.WDR(RecTables.WSPD==0) = 180;
end

if ismember('WDR',RecTables.Properties.VariableNames)
    RecTables.WDR(RecTables.WDR==360) = 0;
    RecTables.WDR(isnan(RecTables.WSPD)) = NaN;
end

if ismember('GHI',RecTables.Properties.VariableNames)
    RecTables.GHI(RecTables.GHI<0) = 0;
end

checkcells = (varfun(@iscellstr, RecTables));
if any(checkcells.Variables)
    RecTables(:,checkcells.Variables) = [];
end

DateColumns = RecTables{:, ismember( ...
    RecTables.Properties.VariableNames, DateColumnNames(1:4))};
% Keeping only Year, Month, Day, Hour

try    
    DateColNames = unique(DateColumns,'rows');
    RecTablesTemp = NaN(size(DateColNames,1), size(RecTables,2));
    for c = 6:size(RecTables,2)
        if c == 6
            [RecTablesTemp(:,c), DateColNames] = grpstats( ...
                RecTables{:,c}, DateColumns, {@nanmean, 'gname'});
        else
            RecTablesTemp(:,c) = grpstats( ...
                RecTables{:,c}, DateColumns, @nanmean);
        end
    end    
    
    % Convert date columns to numbers
    DateColNames = cellfun(@str2double, DateColNames);

catch err
    
    fprintf('%s \r\n',err.message);
    for e = 1:length(err.stack)
        fprintf('%s at %i\n', err.stack(e).name, ...
            err.stack(e).line);
    end
    
    fprintf(['Something happened to grpstats. ',...
        'Trying the same thing manually.\r\n']);
    % Take each unique hour
    tic;
    DateColNames = unique(DateColumns,'rows');
    RecTablesTemp = NaN(size(DateColNames,1), size(RecTables,2));
    
    for d = 1:size(DateColNames,1)
        % Find the matching dates in the overall data table
        DateColFind = (DateColNames(d,1)==RecTables{:,1} & ...
            DateColNames(d,2)==RecTables{:,2} & ...
            DateColNames(d,3)==RecTables{:,3} & ...
            DateColNames(d,4)==RecTables{:,4});
        % Assign the unique hour
        RecTablesTemp(d,1:4) = DateColNames(d,:);
        % The first five columns are ignored since they are the date
        % columns.
        if sum(DateColFind)>1
            % More than one matching hour found - take a mean.
            RecTablesTemp(d,6:end) = nanmean(RecTables{DateColFind, 5:end});
        elseif sum(DateColFind)==1
            % Only one matching hour found
            RecTablesTemp(d,6:end) = RecTables{DateColFind, 5:end};
        else
            % No data found for current hour.
            RecTablesTemp(d,6:end) = NaN(1,size(RecTables,2-4));
        end
    end
    toc
end

% If nanmean encounters only NaNs, it will return a NaN. That is,
% no records exist for a particular hour from any source.

% Reassign the values from the temporary array to the overall table
% and delete it

RecTables2 = [DateColNames, zeros(size(RecTablesTemp,1),1), ...
    RecTablesTemp(:,6:end)];

RecTables2 = array2table(RecTables2);
RecTables2.Properties.VariableNames = RecTables.Properties.VariableNames;
RecTables2.Properties.VariableUnits = RecTables.Properties.VariableUnits;
RecTables2.Properties.Description = RecTables.Properties.Description;

RecTables = RecTables2;

clear RecTables2 RecTablesTemp

% % Convert all columns to double
% RecTables{:,~ismember( ...
%     RecTables.Properties.VariableNames, DateColumnNames)} ...
%     = double(RecTables{:,~ismember(...
%     RecTables.Properties.VariableNames,DateColumnNames)});

% There are some instrumentation or documentation errors in all the data
% files. One common error is values of 'ZERO'

save(MATfilePath,'RecTables')

% fclose(fIDlog);

end