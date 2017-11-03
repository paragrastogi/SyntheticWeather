% Create Synthetic Files

% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function [SuccessfulCopy, varargout] = FileCopier(FilePathMaster, ...
    FilePathNew, MasterTable, NewTable, varargin)


p = inputParser;
p.FunctionName = 'FileCopier';

addRequired(p,'FilePathMaster', @ischar)
addRequired(p,'FilePathNew',@ischar)
addRequired(p,'MasterTable',@istable)
addRequired(p,'NewTable',@istable)
% File ID to which the function should write the output. If not given, then
% the function writes to screen, i.e., fID = 1
addParameter(p,'fID',1,@isnumeric)
addParameter(p,'OneYearCheck',true,@islogical)

parse(p,FilePathMaster, FilePathNew, MasterTable, NewTable, varargin{:})

FilePathMaster = p.Results.FilePathMaster;
FilePathNew = p.Results.FilePathNew;
MasterTable = p.Results.MasterTable;
NewTable = p.Results.NewTable;
fID = p.Results.fID;
OneYearCheck = p.Results.OneYearCheck;


fprintf(fID,'\r\nBEGINNING of comments from FileCopier function\r\n');
% % This function copies the data contained in NewTable to a new %
% file. The file path of the new file is FilePathNew. The format used is %
% obtained from MasterTable, which is the data contained in %
% FilePathMaster.

% Length of the year, usually 8760
N = size(MasterTable,1);
% Copy year values
ActualYear = unique(NewTable.Year);
% MasterYear = unique(MasterTable.Year);

% This copier only works on equal-sized tables of 8760 rows each
if size(NewTable,1)~=8760
    fprintf(fID,['Incoming Data Table size is not 8760,', ...
        'please call regulariser module first .\r\n']);
    SuccessfulCopy = false;
    return
else
    if size(MasterTable,1)~=size(NewTable,1)
	fprintf(fID, ['The sizes of the incoming master and new ', ...
	'source tables do not match. Returning control to main ', ...
	'script, please use file annualizer or regulariser module.\r\n']);
	SuccessfulCopy = false;
	return
	else
	if OneYearCheck 
	if numel(ActualYear) > 2
		fprintf(fID,['More than two unique years found in ', ...
		'input file. Returning control to main script,', ...
		'please use file annualizer module.\r\n']);
		SuccessfulCopy = false;
		return

	elseif numel(ActualYear) == 2
		Occur1 = sum(NewTable.Year==ActualYear(1));
		Occur2 = sum(NewTable.Year==ActualYear(2));
		if (Occur1 > 1) && (Occur2 > 1)
			fprintf(fID,['Two unique years found in input file ', ...
			'and each has more than one hour in the file.', ...
			'Returning control to main script,', ...
			'please use file regularizer module.\r\n']);
			SuccessfulCopy = false;
			return
		elseif (Occur1 > 1) && (Occur2 == 1)
			ActualYear = ActualYear(1);
		elseif (Occur1 == 1) && (Occur2 > 1)
			ActualYear = ActualYear(2);
		else
			fprintf(fID, ['Two unique years found in input file ', ...
			'and neither has more than one hour in the file.', ...
			'Returning control to main script,', ...
			'please use file regularizer module.\r\n']);
		end
	else
		fprintf(fID,['Only one unique year found in file, ', ...
			'continuing...\r\n']);
	end
	else
		fprintf(fID,['You have disabled the unique year check for this file, ', ...
			'continuing with one or more unique years in file...\r\n']);
		SaveYears = NewTable.Year;
		ActualYear = 2005;
		NewTable.Year = repmat(ActualYear,size(NewTable.Year));
	end        
    end
end


% Check for Leap Year
if OneYearCheck
	LeapYearNess = DetermineLeapYear(ActualYear);
	OriginalYear = ActualYear;
	
	if LeapYearNess
		ActualYear = ActualYear - 1;
	end
else
	OriginalYear = SaveYears;
end

% 'Minute' is not needed in sorting hourly resolution files, so change to
% EPW default value (0 or 60). This will create a 'minute' field in the new
% table if one doesn't already exist. For now, set both to zero. Then
% change to Master value zero.
MasterMinute = unique(MasterTable.Minute);
NewTable.Minute = zeros(height(NewTable),1);
MasterTable.Minute = zeros(height(NewTable),1);

% The datevec function makes midnight hour 0. But the EPW convention makes
% midnight hour 24, so the values should have been changed already.

if any(NewTable.Hour==0)
    OrigDescr = NewTable.Properties.Description;
    NewTable.Properties.Description = 'MATLAB';
    NewTable = DateFormatter(NewTable);
    NewTable.Properties.Description = OrigDescr;
end

% Get unique numbers for each timestamp with the year from the new file
% Real TMY years are not used since the months are from arbitrary years.
% The values from new source file have to replace corresponding values in
% the USDOE EPW file (TMY) without regard for actual year.
DummySeconds = zeros(N,1);
DateNumMaster = datenum(repmat(ActualYear,N,1), ...
    MasterTable.Month, MasterTable.Day, MasterTable.Hour,...
    MasterTable.Minute, DummySeconds);
DateNumNew = datenum(ActualYear, NewTable.Month,...
    NewTable.Day, NewTable.Hour, NewTable.Minute,...
    DummySeconds);

% Find intersection of dates, with corresponding indices
[~, MatchIdxMaster, MatchIdxNew] = intersect(DateNumMaster, ...
    DateNumNew,'stable');
[DateMisMatches, ~, MisMatchIdxN] = ...
    setxor(DateNumMaster, DateNumNew,'stable');

if ~isempty(DateMisMatches)
    
    fprintf(fID, ['Some Date/Time stamps do not match for the files', ...
    ' (using dummy year and second vectors). Non-matching time-', ...
    'stamps between master and new source files number %d . \r\n'], ...
        numel(DateMisMatches));
    
    DateNumNew(MisMatchIdxN) = nan;
    DateNumNew(MatchIdxNew) = DateNumMaster(MatchIdxMaster);
    
    clear DateNumNewTemp
    
else
    % No mismatching dates
    fprintf(fID,'All dates match\n');
end

% If the incoming data table contains a string/non-numeric input, then
% remove it. Here we keep only numeric inputs.
StrInputs = cellfun(@isnumeric, table2cell(NewTable(1,:)));
NewTable = NewTable(:,StrInputs);

% Find the variable names common to both data tables. Names are returned in
% the same order as the Master Data Table.
VarNamesToBeCopied = intersect(...
    MasterTable.Properties.VariableNames,...
    NewTable.Properties.VariableNames,'stable');

if isempty(VarNamesToBeCopied) % None of the names match
    fprintf(fID, ['Found no common variable names to copy, i.e. none', ...
    ' of the variable names match between master and new ', ...
    'file. Returning the master file as the new file. \r\n\r\n']);
    SuccessfulCopy = false;
    return
else
    DateColHeaders = {'Year', 'Month', 'Day', 'Hour', 'Minute'};
    DateIntersect = intersect(VarNamesToBeCopied,DateColHeaders,'stable');
    
    if numel(DateIntersect) == numel(VarNamesToBeCopied)
        SuccessfulCopy = false;
        fprintf(fID, ['Only the date columns match (Y,M,D,H,', ...
        'minutes). No new file created. \r\n']);
        return
    else
        % Throw out date columns from variable names to be copied
        VarNamesToBeCopied = VarNamesToBeCopied(6:end);
    end
    
end

% Fill the output table completely with zeros DataTableOut =
% [array2table(zeros(size(MasterTable,1),5)),...
% MasterTable(:,6),array2table(zeros(size(MasterTable,1),...
% size(MasterTable,2)-6))];

% Copy timestamps from New File to output table
[y,m,d,h,mi,~] = datevec(DateNumMaster);

% The DateNumMaster vector is used since it has been generated using the
% same year as the incoming file and it has not NaNs.
DataArrayOut(:,1:5) = [y,m,d,h,mi];

% The datevec function makes midnight hour 0. But the EPW convention makes
% midnight hour 24, so the values should have been changed already.
if any(DataArrayOut(:,4)==0)
    DataArrayOut = DateFormatter(DataArrayOut);
end

if ~OneYearCheck
	DataArrayOut(:,1) = SaveYears;
end

% Column 6 will be copied directly from master EPW file, since it is the
% column of Quality Flags, which is not checked in this script. Initialise
% with zeros here and add the text column later since an array cannot
% handle text input.
DataArrayOut(:,6) = zeros(size(MasterTable,1),1);

DataArrayOut(:,7:(size(MasterTable,2))) = ...
    zeros(size(MasterTable,1),size(MasterTable,2)-6);

% Declare InsaneColumns and UnitProblems vectors
InsaneColumns = false(1,size(MasterTable,2));
UnitProblems = false(1,size(MasterTable,2));
FinalVarsToKeep = true(size(VarNamesToBeCopied,2),1);

for cidx = 1:size(VarNamesToBeCopied,2)
    
    
    % Find variable index in each file
    VarIdxMaster = find(strcmp(VarNamesToBeCopied{cidx}, ...
        MasterTable.Properties.VariableNames));
    VarIdxNewSrc = find(strcmp(VarNamesToBeCopied{cidx}, ...
        NewTable.Properties.VariableNames));
    
    % Find variable in master file
    MasterCopy = MasterTable.(VarIdxMaster);
    % Find variable with same name in New Source file
    NewCopy = NewTable.(VarIdxNewSrc);
    
    if all(isnan(NewCopy))
        NewCopy = MasterCopy;
        DataArrayOut(:,VarIdxMaster) = NewCopy;
        continue
    end
    
    % If a date doesn't exist in the incoming file, then the corresponding
    % data must be marked as NaN.
    NewCopy(isnan(DateNumNew)) = nan;
    
    % Are they both numerical? Only numerical tables are copied from new
    % files.
    CheckLogicalMaster = isa(MasterTable.(VarIdxMaster),'numeric');
    CheckLogicalNewSrc = isa(NewTable.(VarIdxNewSrc),'numeric');
    
    if ~CheckLogicalMaster
        fprintf(fID, ['The current master column to be overwritten ', ...
        'is not numeric.\r\n Sanity Check is false.\r\n']);
        continue
    elseif ~CheckLogicalNewSrc
        fprintf(fID, ['The current new source column to be ', ...
        'overwritten is not numeric.\r\n Sanity Check is false.\r\n']);
        continue
    else
        % Check to see if units match
        CheckUnits = strcmpi(MasterTable.Properties.VariableNames...
            {VarIdxMaster},NewTable.Properties.VariableNames...
            {VarIdxNewSrc});
        
        % If units don't match, send to unit correcter
        if ~CheckUnits
            [NewCopy, CheckUnits] = UnitConv(...
                NewTable.Properties.VariableUnits{VarIdxNewSrc}, ...
                NewTable.(VarIdxNewSrc), VarIdxNewSrc, ...
                NewTable.Properties.Description,...
                MasterTable.Properties.Description);
            if ~CheckUnits
                % If units still don't match, record name of variable in
                % new table
                UnitProblems(cidx) = true;
                NewCopy = MasterCopy;
            end
        end
        
        if CheckUnits
            % If unit conversion works, or was not required, check order of
            % magnitude
            CheckSanity = SanityCheck(MasterTable.(VarIdxMaster), ...
                NewTable.(VarIdxNewSrc));
        else
            % if unit conversion was false, sanity check is false anyway.
            CheckSanity = false;
            NewCopy = MasterCopy;
            fprintf(fID,['Discrepancy in sanity check. If class ', ...
                'of both input and output were numeric, then check ', ...
                'old and new mean of Master variable "%s" and \r\n', ...
                'New variable "%s" for file \r\n %s. Master copy ', ...
                'unchanged \r\n\r\n'], ...
                MasterTable.Properties.VariableNames{VarIdxMaster}, ...
                NewTable.Properties.VariableNames{VarIdxNewSrc}, ...
                FilePathNew);
            % Record the name of the column
            InsaneColumns(cidx) = true;
            % NewTable.Properties.VariableNames{VarIdxNewSrc}
        end
        
    end
    
    if any(isnan(NewCopy))
        % Take the NaN values and interpolate using a Shape-preserving
        % piecewise cubic interpolation.
        try
            NewCopy(isnan(NewCopy)) = interp1(find(~isnan(NewCopy)), ...
                NewCopy(~isnan(NewCopy)), ...
                find(isnan(NewCopy)),'pchip',nan);
            % Interp1 will NOT extrapolate around the edges due to the last
            % input. Instead it will set values outside domain to NaN.
        catch err_interp
            % Incoming NewCopy does not have a sufficient number of points
            % to be able to interpolate. Remove this particular variable
            % from the list of variable names to be copied.
            %             if
            %             strcmp(NewTable.Properties.VariableNames...
            %                     {VarIdxNewSrc}, 'TDP')
            %             if cidx>1 VarNamesToBeCopied =
            %             [VarNamesToBeCopied(1:cidx-1), ...
            %                 VarNamesToBeCopied(cidx+1:end)];
            %             else
            %                 VarNamesToBeCopied =
            %                 VarNamesToBeCopied(2:end);
            %             end
            FinalVarsToKeep(cidx) = false;
            
            fprintf(fID,'%s \r\n',err_interp.message);
            for e=1:length(err_interp.stack)
                fprintf(fID,'%s at %i\n', err_interp.stack(e).name, ...
                    err_interp.stack(e).line);
            end
            fprintf(fID,'Error in file %s .\r\n', FilePathNew);
            %             end
        end
    end
    
    % Take the remaining NaN values from Newcopy and replace them with
    % corresponding values from Mastercopy. This may create some jumps but
    % it will at least put in realistic values from the Master file.
    if any(isnan(NewCopy))
        NewCopy(isnan(NewCopy)) = MasterCopy(isnan(NewCopy));
    end
    
    % Finally, copy column to new data table, using indices from the master
    % table
    if CheckLogicalMaster && CheckLogicalNewSrc && CheckUnits ...
            && CheckSanity
        DataArrayOut(:,VarIdxMaster) = NewCopy;
    else
        DataArrayOut(:,VarIdxMaster) = MasterCopy;
    end
    
end

CopiedVarsData = cell2mat(cellfun(@(x) ...
    (NewTable.(matlab.lang.makeValidName(x))), ...
    VarNamesToBeCopied,'UniformOutput',0));
ZeroRows = (sum(CopiedVarsData,2)<1);

for cidx = 1:size(VarNamesToBeCopied,2)
    
    % Find variable index in each file
    VarIdxMaster = find(strcmp(VarNamesToBeCopied{cidx}, ...
        MasterTable.Properties.VariableNames));
    %     VarIdxNewSrc = find(strcmp(VarNamesToBeCopied{cidx}, ...
    %         NewTable.Properties.VariableNames));
    
    % Find variable in master file
    MasterCopy = MasterTable.(VarIdxMaster);
    % Find variable with same name in New Source file
    %     NewCopy = NewTable.(VarIdxNewSrc);
    
    DataArrayOut(ZeroRows, VarIdxMaster) = MasterCopy(ZeroRows);
    
end


% Now copy from the master file those variables that were not available in
% the new data file
for ncidx = 7:size(MasterTable,2)
    
    if ismember(MasterTable.Properties.VariableNames{ncidx}, ...
            VarNamesToBeCopied)
        continue
    else
        DataArrayOut(:,ncidx) = MasterTable.(ncidx);
    end
    
end


% % Special consideration of the relationship between TDB, RH, and TDP

% Find indices of RH, TDB, and TDP. These are taken from master table since
% the copier loop also uses the master index.
VarIdxRH = find(strcmp('RH', ...
    MasterTable.Properties.VariableNames));
VarIdxTDB = find(strcmp('TDB', ...
    MasterTable.Properties.VariableNames));
VarIdxTDP = find(strcmp('TDP', ...
    MasterTable.Properties.VariableNames));

if ismember('TDB', VarNamesToBeCopied)
    % If TDB was copied
    if ismember('TDP', VarNamesToBeCopied) && ...
            ismember('RH', VarNamesToBeCopied)
        fprintf(fID, ['TDB, TDP, and RH were present in the ', ...
        'actual data. No changes made to any of the three ', ...
        'variables.\r\n']);
    elseif ismember('TDP', VarNamesToBeCopied) && ...
            ~ismember('RH', VarNamesToBeCopied)
        fprintf(fID, ['TDB and TDP were present in the actual data ', ...
        'but not RH. RH calculated from TDB and TDP.\r\n', ...
        'Using Lawrence''s (02/2005) approximation\n']);
        
        % Replace RH column in output data array
        DataArrayOut(:,VarIdxRH) = 100 - ...
            5.*(DataArrayOut(:,VarIdxTDB)-DataArrayOut(:,VarIdxTDP));
        
    elseif ~ismember('TDP', VarNamesToBeCopied) && ...
            ismember('RH', VarNamesToBeCopied)
        fprintf(fID,['TDB and RH were present in the actual data ', ...
        'but not TDP. TDP calculated from TDB and RH.\r\n', ...
        'Using Lawrence''s (02/2005) approximation\n']);
        
        % Replace TDP column in output data array
        DataArrayOut(:,VarIdxTDP) = DataArrayOut(:,VarIdxTDB) - ...
            ((100-DataArrayOut(:,VarIdxRH))./5);
        
    elseif ~ismember('TDP', VarNamesToBeCopied) && ...
            ~ismember('RH', VarNamesToBeCopied)
        fprintf(fID, ['TDB was present in the actual data but not.', ...
        ' TDP or RH. RH picked from corresponding TMY file and ', ...
        'TDP calculated from TDB and (TMY) RH.\r\n']);
        
        % Get RH column from master file
        DataArrayOut(:,VarIdxRH) = MasterTable.RH;
        
        % Replace TDP column in output data array
        DataArrayOut(:,VarIdxTDP) = DataArrayOut(:,VarIdxTDB) - ...
            ((100-DataArrayOut(:,VarIdxRH))./5);
        
    end
    
end


% Reset minute back to the minute values from the master file
DataArrayOut(:,5) = repmat(MasterMinute,N,1);

% Reset year back to original year.
if OneYearCheck
	DataArrayOut(:,1) = repmat(OriginalYear,N,1);
else
	DataArrayOut(:,1) = OriginalYear;
end
% If original year was not a leap year, this statement doesn't make a
% difference anyway.


% Format string for EPW files obtained from USDOE website
% 1990,1,1,1,60,?9?9?9?9E0?9?9?9?9*9?9?9?9?9?9?9?9*9?9?9?9?9,18.7,14.4,
% 76,101000,0,0,340,0,0,0,0,0,0,0,0,0.9,0,0,1.3,77777,9,999999999,0,
% 0.0000,0,0,999.000,999.0,99.0
StringBegin = '%d,%d,%d,%d,%d,%44s'; StringRepeated = ',%.1f';
if ispc
    LineEnder = '\r\n';
elseif isunix
    LineEnder = '\r\n';
end

FormatEPWData = [StringBegin, repmat(StringRepeated, 1, ...
    size(MasterTable,2)-6), LineEnder];
FormatEPWHeader = '%s\r\n';


% Make a copy of master file
SuccessfulCopy = copyfile(FilePathMaster,FilePathNew);
if ~SuccessfulCopy
    fprintf(fID, ['Unable to copy file for some reason. ', ...
    'Last read file is %s\n'],FilePathOriginal);
end

% Open both files for low-level read/write operations
fileIDMaster = fopen(FilePathMaster,'r');
fileIDNew = fopen(FilePathNew,'w');

LineCounter = 1; % Line Number counter
EPWHeaderLinesCount = 8; % Number of header lines in a regular EPW file

while ~feof(fileIDMaster) % End of file is not reached
    
    FileLine = fgetl(fileIDMaster); % Get current line of master file
    
    % Copy Header lines almost verbatim (except COMMENTS 1 field)
    
    if LineCounter >= 1 && LineCounter <= EPWHeaderLinesCount
        
        if strcmpi(FileLine(1:10),'COMMENTS 1')
            % Add name of new source to header line COMMENTS 1
            NewSourceAddendum = sprintf(['; Some data columns ', ...
                'replaced with values from %s'], ...
                NewTable.Properties.Description);
            FileLineNew = strcat(FileLine, NewSourceAddendum);
            FileLine = FileLineNew;
        end
        
        % Write FileLine string to new file
        SuccessfulPrint = fprintf(fileIDNew, FormatEPWHeader, FileLine);
        % One line lifted from header of master file
        LineCounter = LineCounter + 1;
        
    elseif LineCounter > EPWHeaderLinesCount && ...
            LineCounter <= (size(MasterTable,1) + EPWHeaderLinesCount)
        % Now we write the data
        
        TableLineCounter = LineCounter - EPWHeaderLinesCount;
        % Only write the data one line at a time, until the end of file and
        % last line of table are reached
        SuccessfulPrint = fprintf(fileIDNew, FormatEPWData, ...
            DataArrayOut(TableLineCounter,1:5), ...
            MasterTable.QualFlags{TableLineCounter}, ...
            DataArrayOut(TableLineCounter,7:end));
        clear TempCell
        
        LineCounter = LineCounter + 1;
        
    else
        SuccessfulPrint = 0;
        fprintf(fID,'Line number is greater than 8768.\r\n');
    end
    
    if SuccessfulPrint == 0
        fprintf(fID,'There was an error copying line number %d .\r\n', ...
            LineCounter);
    end
    
end

% Close both files
fclose(fileIDMaster);fclose(fileIDNew);

if exist(FilePathNew,'file')
    SuccessfulCopy = true;
    fprintf(fID,'Successfully created %s \r\n',FilePathNew);
else
    SuccessfulCopy = false;
    fprintf(fID,'Unable to create %s \r\n',FilePathNew);
end

fprintf(fID,'End of comments from FileCopier function \r\n\r\n');

CheckLogicalOverall = CheckLogicalNewSrc && CheckLogicalMaster;

varargout{1} = CheckLogicalOverall; varargout{2} = CheckUnits;
varargout{3} = CheckSanity;

end