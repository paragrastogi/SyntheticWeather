function SMY_Create(stcode, path_syn_fldr, ...
    path_save_fldr, varargin)


% This script will write EPW files from synthetic weather
% files. The synthetic weather files should be in MAT
% format, unless otherwise specified in the optional inputs.

p = inputParser;
p.FunctionName = 'SMY_Create';

% The first input is a cell array of strings specifying
% which folders to run
addRequired(p, 'stcode', @ischar)

% This is where the EPW files should be
addRequired(p, 'path_syn_fldr', @ischar)

% This is where the new files will be written, in a subfolder
% based on the 'stcode' input.
addRequired(p, 'path_save_fldr', @ischar)

% Do you want existing files deleted? They will be left
% unharmed by default, and the corresponding new file (i.e.,
% one with the same name), will not be written.
addParameter(p, 'destroyer', false, @islogical)
addParameter(p, 'nameLOGfile', fullfile('LogFiles', ...
    sprintf('LogFile_SMY_Create_%s.txt',date)),@ischar)
addParameter(p, 'synformat', 'mat', @ischar)
addParameter(p, 'scenario', 'syn', @ischar)
addParameter(p, 'sublabel', 'x', @isnumeric)
addParameter(p, 'srcEPWfile', '', @ischar)

parse(p,stcode, path_syn_fldr, ...
    path_save_fldr, varargin{:})

stcode = p.Results.stcode;
destroyer = p.Results.destroyer;
path_syn_fldr = p.Results.path_syn_fldr;
path_save_fldr = p.Results.path_save_fldr;
nameLOGfile = p.Results.nameLOGfile;
synformat = p.Results.synformat;
scenario = p.Results.scenario;
sublabel = p.Results.sublabel;
srcEPWfile = p.Results.srcEPWfile;

if exist(path_save_fldr, 'dir') ~= 7
    mkdir(path_save_fldr)
end

if exist('LogFiles', 'dir') ~= 7
    mkdir('LogFiles')
end

% % Open a log file to record all output message
fIDlog = fopen(nameLOGfile,'w');
fprintf(fIDlog, ['Commencing SMY_Create script at ', ...
    '%s \r\n'], char(datetime));

if destroyer
    fprintf(fIDlog, ['Destroyer is set to true.\r\n', ...
        'All existing SMY files will be deleted.\r\n', ...
        'I AM BECOME DEATH, THE DESTROYER OF ', ...
        'EXISTING FILES... \r\n']);
else
    fprintf(fIDlog, ['Destroyer is false, so no files will', ...
        ' be deleted. This isn''t so much fun... though', ...
        ' it is probably faster.\r\n']);
end

fprintf(fIDlog, ['WARNING\r\n', ...
    'All input files are treated as text files.\r\n', ...
    'This could cause problems with XLS and XLSX files.\r\n']);


% Skip folders that have not been requested
if exist(path_syn_fldr, 'dir') ~= 7
    fprintf('I could not find incoming files folder %s. \r\n', ...
        path_syn_fldr)
    fclose(fIDlog);
    return
end

if exist(path_save_fldr, 'dir') ~=7
    mkdir(path_save_fldr)
end


fprintf(fIDlog,'Processing folder %s\n', path_syn_fldr);

[~, FileNameMaster, ~] = fileparts(srcEPWfile);

if contains(srcEPWfile, 'Meteonorm')
    % Read the file using Meteonorm file parser.
    mastertable = WeatherFileParseEPWMeteonorm(srcEPWfile);
else
    % Read the file using USDOE file parser.
    mastertable = WeatherFileParseEPWMeteonorm(srcEPWfile);
end

% Add source name to the description
mastertable.Properties.Description = ...
    [mastertable.Properties.Description, '_', FileNameMaster];

fprintf(fIDlog,'Master Data is from %s \n', ...
    FileNameMaster);


% Find the synthetic data file in the current
% folder, corresponding to this station.
if strcmpi(synformat,'mat')
    if strcmpi(scenario, 'syn')
        MATfileSyn = dir(fullfile(path_syn_fldr, ...
            sprintf('%s*Syn.mat', ...
            FileNameMaster(1:end-4))));
    elseif strcmpi(scenario, 'rcp45')
        MATfileSyn = dir(fullfile(path_syn_fldr, ...
            sprintf('%s*rcp45.mat', ...
            FileNameMaster(1:end-4))));
    elseif strcmpi(scenario, 'rcp85')
        MATfileSyn = dir(fullfile(path_syn_fldr, ...
            sprintf('%s*rcp85.mat', ...
            FileNameMaster(1:end-4))));
    end
elseif strcmpi(synformat,'csv')
    MATfileSyn = dir(fullfile(path_syn_fldr, ...
        sprintf('%s*Syn.csv', ...
        FileNameMaster(1:end-4))));
end

% Remove MAT files with the 'extras' suffix from
% the list since they contain, well, extra data.
MATfileSyn = MATfileSyn(cellfun(@isempty,...
    regexp({MATfileSyn.name},'extras')));
% Keep only the names from the list
MATfileSyn = char({MATfileSyn.name});

% If no synthetic data files are found, skip this city
if isempty(MATfileSyn)
    fprintf(fIDlog,['No synthetic data files found for', ...
        ' base file %s\n'], FileNameMaster);
    fclose(fIDlog);
    return
elseif all(size(MATfileSyn)>1)
    fprintf(fIDlog,['Multiple synthetic data files found for', ...
        ' base file %s. Continuing with data file name that matches best.\n'], FileNameMaster);
    for f = 1:size(MATfileSyn, 1)
        if strcmp(MATfileSyn(f, 1:end-8), FileNameMaster)
            temp = MATfileSyn(f,:);
            break
        else
            temp = 'na';
        end
    end
    MATfileSyn = temp;
end


if strcmpi(synformat,'mat')
    % This loads a struct called 'syndata'
    load(fullfile(path_syn_fldr, MATfileSyn), 'syndata');
    
elseif strcmpi(synformat,'csv')
    readsynCSV = csvread(path_syn_fldr, MATfileSyn);
    syndata.tdb = readsynCSV(:,1);
    syndata.rh = readsynCSV(:,2);
    syndata.ghi = readsynCSV(:,3);
    syndata.dni = readsynCSV(:,4);
    syndata.dhi = readsynCSV(:,5);
end

if exist('syndata', 'var')==1
    fprintf(fIDlog, ['Synthetic data struct ', ...
        'successfully loaded.\r\n']);
end

% The length of a year is 8760, and the number
% of synthetic years read in is nboot. (number
% of bootstraps)
N = 8760;

UniqueSynYears= unique(syndata.Year);
PickBootLen = length(syndata.Year) / length(UniqueSynYears) / N;
SynYears = repmat(UniqueSynYears, PickBootLen, 1);

% New batch has a set of sub labels shifted by the value of sublabel.
temp1 = (0:1:(PickBootLen-1)) + sublabel;
temp2 = repmat(temp1, length(UniqueSynYears),1);
temp2 = temp2(:);

SubLabelArray = cellfun(@num2str, num2cell(temp2), 'UniformOutput', 0);

SynYearsNames = strcat(cellfun(@num2str, ...
    num2cell(SynYears), 'UniformOutput', 0), '_', SubLabelArray);

% Path to the new file, i.e. output file.
FilePathNew = fullfile(path_save_fldr, cellfun(@(x) ...
    sprintf('%s_%s_%s.epw', stcode, ...
    scenario, x), SynYearsNames, 'UniformOutput', 0));

% Check if file already exists. If destroyer was true, the existing
% file should have been deleted already.
FileExister = cellfun(@(x) exist(x,'file')==2, FilePathNew);

% Logical to record succesful creation of files
copysuccess = false(1,length(FilePathNew));


% Format string for EPW files obtained from
% USDOE website
% 1990,1,1,1,60,?9?9?9?9E0?9?9?9?9*9?9?9?9?9?9?9?9*
% 9?9?9?9?9,18.7,14.4,76,101000,0,0,340,0,0,0,0,0,0,
% 0,0,0.9,0,0,1.3,77777,9,999999999,0,0.0000,0,0,
% 999.000,999.0,99.0
strbegin = '%4d,%d,%d,%d,%02d,%7d';
strrepeat = ',%.1f';
LineEnder = '\r\n';

FormatEPWData = [strbegin, repmat(strrepeat, ...
    1, size(mastertable,2)-6), LineEnder];
formatEPWhead = '%s\r\n';

% Read in the header from the EPW master file as a series of lines.

% Number of header lines in a regular EPW file
numEPWhead = 8;
% Save the header lines
HeaderSave = cell(numEPWhead,1);

% Open file for low-level read/write operations
fileIDMaster = fopen(srcEPWfile,'r');

% Copy Header lines almost verbatim (except
% COMMENTS 1 field)

new_src_add = upper(['; The following ', ...
    'data columns replaced with synthetic ', ...
    'values: ''GHI'', ''DNI'', ''DHI'', ', ...
    '''RH'', ''TDB''.']);
new_src_add_syn = [new_src_add, upper(['The years are dummy values, ', ...
    'they do not imply an actual ', ...
    'prediction for a specific year. ', ...
    'The hourly values of each parameter are not ', ...
    'predictions either, so do not think of them as', ...
    'the exact values at some specific future hour.'])];
new_src_add_rcp = [new_src_add, upper(['The years represent the ', ...
    'prediction for a specific year in the ', ...
    'future. However, a prediction should ', ...
    'not be taken too literally, since it ', ...
    'comes from an approximate model based', ...
    ' on the state of the art in 2005.', ...
    'All data was downloaded from the CORDEX', ...
    ' project web site with the help of ', ...
    'Georgios Mavormatidis, EMPA Duebendorf ', ...
    '(Zurich) in October 2015.'])];

for el = 1:numEPWhead
    
    % Get current line of master file
    FileLine = fgetl(fileIDMaster);
    
    if strcmpi(FileLine(1:10),'COMMENTS 1')
        % Add name of new source to header line
        % COMMENTS 1
        if strcmpi(scenario,'syn')
            new_src_add = new_src_add_syn;
        else
            new_src_add = new_src_add_rcp;
        end
        FileLineNew = strcat(FileLine, new_src_add);
        FileLine = FileLineNew;
    end
    
    HeaderSave(el) = {FileLine};
    
end

fclose(fileIDMaster);


for f = 1:length(FilePathNew)
    
    if ~destroyer
        % If the existing files should not be
        % destroyed, then skip the current
        % iterate if file already exists.
        if FileExister(f)
            fprintf(fIDlog, ['File %s exists, ', ...
                'not overwritten\r\n']);
        end
    end
    
    FilePathNewCurr = FilePathNew{f};
    
    fileIDNew = fopen(FilePathNewCurr,'w');
    
    % Write the header lines
    cellfun(@(x) fprintf(fileIDNew, formatEPWhead, x), ...
        HeaderSave);
    
    % Get the numeric part of the table and
    % convert it to an array.
    if f == 1
        CurrIdx = 1:N;
    else
        CurrIdx = ((f-1)*N)+1:(f*N);
    end
    
    % Break loop when you are out of the syndata index.
    if any(CurrIdx>numel(syndata.Year))
        break
    end
    
    
    % Copy the entire master table first.
    mtablecurr = mastertable;
    
    % First the dates
    mtablecurr.Year = syndata.Year(CurrIdx);
    mtablecurr.Month = mastertable.Month;
    mtablecurr.Day = mastertable.Day;
    mtablecurr.Hour = mastertable.Hour;
    
    % Convert the quality flags to a meaningless
    % number. This lets the table be converted
    % to
    mtablecurr.QualFlags = repmat(9999999, ...
        size(mastertable,1),1);
    
    % TDB
    mtablecurr.TDB = syndata.tdb(CurrIdx);
    
    % GHI
    mtablecurr.GHI = syndata.ghi(CurrIdx);
    
    % DNI
    mtablecurr.DNI = syndata.dni(CurrIdx);
    
    % DHI
    mtablecurr.DHI = syndata.dhi(CurrIdx);
    
    % RH
    mtablecurr.RH = syndata.rh(CurrIdx);
    
    % Convert to array and rotate for writing.
    mtablecurr = (table2array(mtablecurr))';
    
    % Write the numerical data to the file
    copysuccess(f) = fprintf(fileIDNew, ...
        FormatEPWData, mtablecurr);
    
    fclose(fileIDNew);
    
    if copysuccess(f) == 0
        fprintf(fIDlog,['There was an ', ...
            'error copying line number ', ...
            '%d .\r\n'], LineCounter);
    end
    
    fprintf('Written file number %d\r\n', f);
    
end

fclose(fIDlog);

% Close any remaining files.
fclose all;

end