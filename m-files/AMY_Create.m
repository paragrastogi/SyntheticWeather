% Create Synthetic Files

% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE
% LAUSANNE, Switzerland, Interdisciplinary Laboratory of
% Performance-Integrate Design, 2016 Parag Rastogi 
% See the LICENSE.TXT file for more details.

% This script will take in files of certain formats and from
% certain sources to create EPW files for simulation using
% EnergyPlus.

% The basis for creating new files has to always be an EPW
% file. This script will replace individual columns as it
% finds them. All replacements are based on timestamps. In
% case of duplicate time stamps, the script will average the
% quantities.

% I'm sorry it currently relies on a very specific folder
% structure. But I think this can be dealt with in due
% course!

function AMY_Create(folderstorun,varargin)

% The first input is a cell array of strings specifying
% which folders to run. The reasoning is that the script
% will have a list of ALL folders and will only run the ones
% you want.

p = inputParser;
p.FunctionName = 'AMY_Create';

% Check inputs

addRequired(p,'folderstorun', @iscellstr)
% destroyer tells the script whether to overwrite existing
% files or not
addParameter(p,'destroyer', false, @islogical)
% Path to the parent or top-level folder, called 'top
% folder'
addParameter(p, 'TopFolderPath', fullfile('..', '..', ...
	'Weather_Data','WeatherFilesForAnalysis'), @ischar)

% If the incoming weather record is smaller than 80% of a
% full year, it will be rejected
cutoff_default = 0.8;
addParameter(p,'cutoff',cutoff_default,@isnumeric)
% Log the screen outputs
addParameter(p,'LogFileName', ...
	sprintf('LogFile_AMY_Create_%s.txt',date), ...
	@ischar)

parse(p,folderstorun,varargin{:})

folderstorun = p.Results.folderstorun;
destroyer = p.Results.destroyer;
TopFolderPath = p.Results.TopFolderPath;
LogFileName = p.Results.LogFileName;
cutoff = p.Results.cutoff;

% % Open a log file to record all output message
fIDlog = fopen(LogFileName,'w');
fprintf(fIDlog,['Commencing AMY_Create script ', ...
	'at %s \r\n'], char(datetime));

if destroyer
	fprintf(fIDlog,['I am become death, the ', ...
		'destroyer of files...\r\n', ...
		'since destroyer is set to true.\r\n', ...
		'All existing AMY files will be deleted.\r\n']);
else
	fprintf(fIDlog,['destroyer is false, no files will be ', ...
		'deleted... Kind of boring if you ask me.\r\n']);
end


% The path to the parent or top-level folder, called 'top
% folder', is an optional input.

% Make a list of files in the top folder
TopFolderDirList = dir(TopFolderPath);
% Dir function picks up everything, including '.', '..', and
% files. In this case, we are only interested in the city
% folders. Remove folder names that start with a dot
TopFolderDirList(strncmp({TopFolderDirList.name}, '.', 1)) = [];
% Remove files
TopFolderDirList = TopFolderDirList([TopFolderDirList.isdir]);
% Keep only the names from the list
TopFolderDirList = {TopFolderDirList.name};

% Make a struct with the keywords for TMY and AMY - usually
% based on the source of the data
FileParserNames = FileParserKeywordLister;

fprintf(fIDlog,['WARNING\r\n', ...
	'All input files are treated as text files.\r\n', ...
	'This could cause problems with XLS and XLSX files.\r\n']);

for CityNumber = 1:length(TopFolderDirList)
	% Cycle through the entire top-level folder
	
	%     try
	
	% Skip folders that have not been requested
	if ~ismember(TopFolderDirList{CityNumber}, folderstorun)
		fprintf(fIDlog, ['Skipping folder %s, because you ', ...
			'didn''t ask for it... \r\n'], ...
			TopFolderDirList{CityNumber});
		continue
	end
	
	if ~ismember(TopFolderDirList{CityNumber}, TopFolderDirList)
		fprintf(fIDlog,'Unknown folder name %s. Also skipped... \r\n', ...
			TopFolderDirList{CityNumber});
		continue
	end
	
	fprintf(fIDlog,'Processing folder %s\n', TopFolderDirList{CityNumber});
	
	% Search directory of current city for all files. If a
	% file other than a weather file of known format exists,
	% it will be handled by the last 'else' condition in the
	% file handling loop below.
	
	% First set the path to city subfolder
	CityDirPath = fullfile(TopFolderPath, TopFolderDirList{CityNumber});
	
	
	% Find only EPW files
	WthrFilesListEPW = dir(fullfile(CityDirPath,'*.epw'));
	% Keep only the names from the list
	WthrFilesListEPW = {WthrFilesListEPW.name};
	
	% Some folders have multiple stations. Files from
	% different stations must be treated separately.
	
	% Split the station names coming from the EPW files
	SplitStations = cellfun(@(x) (strsplit(x,{'_','.'})), ...
		WthrFilesListEPW,'UniformOutput',0);
	% Look to see if the second word in the split names
	% matches any of the TMY keywords
	SecondWord = cell2mat(cellfun(@(x) (strcmp(x{2}, ...
		FileParserNames.TMYKeyword)), SplitStations,'UniformOutput',0));
	if any(SecondWord)
		SplitStations = unique(cellfun(@(x) (x{1}), ...
			SplitStations,'UniformOutput',0));
		NumStations = 1;
	else
		SplitStations = unique(cellfun(@(x) ([x{1},'_',x{2}]), ...
			SplitStations,'UniformOutput',0));
		NumStations = numel(SplitStations);
	end
	
	for StIdx = 1:NumStations
		
		if NumStations>1
			% There is more than one station, which means
			% that the file source is at the third location,
			% i.e. the station name is two words long.
			StNames = 2;
		else
			% Else there is only one word in the station
			% name.
			StNames = 1;
		end
		
		StationID = SplitStations{StIdx};
		
		fprintf(fIDlog,'Processing station %s\n', StationID);
		
		% Find only EPW files
		WthrFilesListEPW = dir(fullfile(CityDirPath,'*.epw'));
		% Keep only the names from the list
		WthrFilesListEPW = {WthrFilesListEPW.name};
		
		% Keep only the files from the CURRENT STATION
		WthrFilesListEPW = WthrFilesListEPW(~cellfun(@isempty, ...
			cellfun(@(x) (strfind(x,StationID)), ...
			WthrFilesListEPW,'UniformOutput',0)));
		
		% Find the existing AMY files. This needs to be
		% modified based on your naming convention. In my
		% case, the word 'and' only exists in an AMY file.
		WthrFilesAMYIdx = ~cellfun(@isempty,(cellfun(@(x) ...
			(regexp(x,'and')), WthrFilesListEPW,'UniformOutput',0)));
				
		if destroyer
			% If the destroyer is true, delete AMY files
			cellfun(@delete,fullfile(CityDirPath, ...
				WthrFilesListEPW(WthrFilesAMYIdx)))
			fprintf(fIDlog,'%d files deleted in folder.\r\n', ...
				sum(WthrFilesAMYIdx));
		end
		
		% List of EPW format files that were not AMY files
		WthrFilesListEPW = WthrFilesListEPW(~WthrFilesAMYIdx);
		
		% Split the names again.
		SplitNames = cellfun(@(x) (strsplit(x,{'_','.'})), ...
			WthrFilesListEPW,'UniformOutput',0);
		
		% The sources of the files.
		EPWFileSources = unique(cellfun(@(x) (x{StNames+1}), ...
			SplitNames, 'UniformOutput',0));
		
		
		if length(EPWFileSources)>1
			% There is more than one TMY file available for
			% the given station. Then exclude METEONORM.
			MasterIdx = find(cell2mat(cellfun(@isempty, ...
				strfind(EPWFileSources,'Meteonorm'), ...
				'UniformOutput',0)),1,'first');
			% Source of master file
			FileSrcMaster = char(SplitNames{MasterIdx}(StNames+1));
			% Master file name
			FileNameMaster = WthrFilesListEPW{MasterIdx};
			% Path to master file
			FilePathMaster = fullfile(CityDirPath,FileNameMaster);
			% Read the file using USDOE parser
			MasterTable = WeatherFileParseEPWUSDOE(FilePathMaster);
			% Add source name to the description
			MasterTable.Properties.Description = ...
				[MasterTable.Properties.Description, ...
				'_', FileSrcMaster];
		else
			% There is only one TMY file available for the
			% given station This file is most likely the
			% Meteonorm file.
			MasterIdx = 1;
			% Source of master file
			FileSrcMaster = char(SplitNames{MasterIdx}(StNames+1));
			% Master file name
			FileNameMaster = WthrFilesListEPW{MasterIdx};
			% Path to master file
			FilePathMaster = fullfile(CityDirPath,FileNameMaster);
			% Read the file using Meteonorm parser
			MasterTable = WeatherFileParseEPWMeteonorm(FilePathMaster);
			% Add source name to the description
			MasterTable.Properties.Description = ...
				[MasterTable.Properties.Description, ...
				'_', FileSrcMaster];
		end
		
		
		fprintf(fIDlog,'Master Data is from %s \n', FileNameMaster);
		
		
		% Find all files in the current folder
		WthrFilesListOthers = dir(CityDirPath);
		% Remove folders
		WthrFilesListOthers = WthrFilesListOthers(~[...
			WthrFilesListOthers.isdir]);
		% Remove EPW files from the list
		WthrFilesListOthers = WthrFilesListOthers(cellfun(@isempty,...
			regexp({WthrFilesListOthers.name},'.epw')));
		% Remove TMY3 files from the list since they are not
		% actual data
		WthrFilesListOthers = WthrFilesListOthers(cellfun(@isempty,...
			regexp({WthrFilesListOthers.name},'TMY3')));
		% Remove MAT files from the list since they are not
		% actual data
		WthrFilesListOthers = WthrFilesListOthers(cellfun(@isempty,...
			regexp({WthrFilesListOthers.name},'.mat')));
		% Remove any files/folders whose name begins with a
		% dot
		WthrFilesListOthers = WthrFilesListOthers(~strncmp({...
			WthrFilesListOthers.name},'.',1));
		% Keep only the names from the list
		WthrFilesListOthers = {WthrFilesListOthers.name};
		% Keep only the files from the CURRENT STATION
		WthrFilesListOthers = WthrFilesListOthers(~cellfun(@isempty, ...
			cellfun(@(x) (strfind(x,StationID)), ...
			WthrFilesListOthers,'UniformOutput',0)));
		
		% If no actual weather data files are found, skip
		% this city
		if isempty(WthrFilesListOthers)
			fprintf(fIDlog, ['No actual data files found for ', ...
				'station %s\n'], FileNameMaster);
			continue
		end
		
		% This will store the names of the sources
		SourceCollect = repmat({''},1,length(WthrFilesListOthers));
		
		% Now cycle through each actual data file to collect
		% the data and write it into an AMY file
		
		for FileNum = 1:length(WthrFilesListOthers)
			
			% Get the name and extension of the current file
			[~,CurrentFileName,CurrentFileExt] = ...
				fileparts(WthrFilesListOthers{FileNum});
			% Ignore the warning here that CurrentFilePath
			% is not used. It IS used - inside an eval
			% function.
			CurrentFilePath = fullfile(CityDirPath, ...
				WthrFilesListOthers{FileNum});
			
			% Split the current file name into its
			% consituents
			SplitName = strsplit(CurrentFileName,'_');
			
			
			if strcmp(CurrentFileExt,'.epw')
				% Skip EPW files
				fprintf(fIDlog,'Skipping EPW file...\n');
				continue
				
			elseif strcmpi(CurrentFileExt,'.csv') || ...
					strcmpi(CurrentFileExt,'.txt') || ...
					strcmpi(CurrentFileExt,'.wy2')
				
				% Only known sources are in AMYKeyword field
				Source = intersect(FileParserNames.AMYKeyword, ...
					SplitName, 'stable');
				% If the encountered source is not in the
				% AMY keywords, that means we don't have a
				% parser for that format. So, skip it for
				% now.
				if isempty(Source)
					% If a valid source is not found
					fprintf(fIDlog,['Unknown file source in file %s.', ...
						'Skipping...\r\n'], CurrentFileName);
					continue
				else
					% Make the cell into a string.
					Source = Source{1};
					
					% If a valid source is found
					fprintf(fIDlog,['Calling %s function on current', ...
						' file %s\n'], ['WeatherFileParse',Source], ...
						CurrentFileName);
					[ScreenOutput, TempData] = evalc(['feval(matlab',...
						'.lang.makeValidName([''WeatherFileParse'',', ...
						'Source]),CurrentFilePath)']);
					fprintf(fIDlog,ScreenOutput);
					
					% Record the name of the source
					SourceCollect{FileNum} = Source;
				end
				
			else
				fprintf(fIDlog,['Unknown file extension %s for file ', ...
					'named %s. Skipping...\n'], ...
					CurrentFileExt, CurrentFileName);
				TempData = nan;
			end
			
			
			% Special cleaning for NCDC files
			if strcmp(Source,'NCDC')
				TempData = NCDCcleaner(TempData);
			end
			
			% Change the hour to match EPW convention. This
			% depends on the incoming weather source.
			if any(TempData.Hour==0)
				TempData = DateFormatter(TempData);
			end
			
			
			if FileNum == 1
				% The data table needs to be created first
				RecTables = TempData;
				
			else
				% Joining new data with existing table
				RecTablesOld = RecTables;
				
				RecTables = outerjoin(RecTablesOld, ...
					TempData, 'MergeKeys', true, 'Type', 'full');
				
			end
			
			
		end
		
		%         % Following EPlus convention, the wind
		%         direction is 0° only when % the data says
		%         so. If there is no wind, i.e. WSPD == 0,
		%         then wind % direction is 180°. Wind
		%         direction 360° is due north, reset that %
		%         to 0° to avoid unnecessary fractions in
		%         direction. if
		%         ismember('WSPD',RecTables.Properties.VariableNames)
		%             RecTables.WDR(RecTables.WSPD==0) =
		%             180;
		%         end
		%
		%         if
		%         ismember('WDR',RecTables.Properties.VariableNames)
		%             RecTables.WDR(RecTables.WDR==360) = 0;
		%             RecTables.WDR(isnan(RecTables.WSPD)) =
		%             NaN;
		%         end
		%
		%         if
		%         ismember('GHI',RecTables.Properties.VariableNames)
		%             RecTables.GHI(RecTables.GHI<10) = 0;
		%         end
		
		% Send all weather data for sanity checks
		RecTables = WeatherCheck(RecTables);
		MasterTable = WeatherCheck(MasterTable);
		
		DateColumnNames = {'Year','Month','Day','Hour','Minute'};
		
		% Sorting the table using dates is straightforward
		% since the resolution of interest is one hour.
		[RecTables, ~] = sortrows(RecTables, DateColumnNames,'ascend');
		% This sorting is useful later since it greatly
		% simplifies the detection of duplicate timestamps.
		
		
		% Outerjoin creates a lot of NaN (not a number)
		% values. These can be replaced by the previous
		% valid record or interpolation. We use
		% nearest-neighbour interpolation here.
		
		% Find the first and last valid points in the time
		% series of GHI and TDB
		KeyVarsPresent = intersect(...
			RecTables.Properties.VariableNames,{'GHI','TDB'});
		
		if numel(KeyVarsPresent)<2
			fprintf(fIDlog,['Either GHI or TDB is missing for ', ...
				'station %s. \r\n NO AMY FILES WRITTEN. \r\n'], StationID);
			SuccessfulCopy = false;
			continue
		end
		
		FirstValidPoints(1) = find(~isnan(RecTables.GHI),1);
		FirstValidPoints(2) = find(~isnan(RecTables.TDB),1);
		
		LastValidPoints(1) = find(~isnan(RecTables.GHI),1,'last');
		LastValidPoints(2) = find(~isnan(RecTables.TDB),1,'last');
		
		% Take the LAST first valid point and FIRST last
		% valid point.
		PointsToCut = [max(FirstValidPoints), min(LastValidPoints)];
		
		% Keep data table only between these points
		RecTables = RecTables(PointsToCut(1):PointsToCut(2),:);
		
		DateColumns = RecTables{:, ismember( ...
			RecTables.Properties.VariableNames, DateColumnNames(1:4))};
		% Keeping only Year, Month, Day, Hour
		
		% UniqueDates = unique(DateColumns, 'rows');
		RecTablesTemp = NaN(size(DateColumns,1), size(RecTables,2));
		
		% Deduplicate by taking a nanmean based on hour
		% Record the date column names only the first time.
		[RecTablesTemp(:,6), DateColNames] = grpstats( ...
			RecTables{:,6}, DateColumns, {@nanmean, 'gname'});
		% Then cycle through the rest of the columns.
		for k = 7:size(RecTables,2)
			RecTablesTemp(:,k) = grpstats(RecTables{:,k}, ...
				DateColumns, @nanmean);
		end
		
		% Convert cell array of strings to matrix of doubles
		DateColNames = cellfun(@str2double, DateColNames);
		
		% If nanmean encounters only NaNs, it will return a
		% NaN. That is, no records exist for a particular
		% hour from any source.
		
		% Reassign the values from the temporary array to
		% the overall table and delete it.
		RecTables2 = [DateColNames, zeros(size(RecTablesTemp,1),1) , ...
			RecTablesTemp(:,6:end)];
		RecTables2 = array2table(RecTables2);
		RecTables2.Properties.VariableNames = ...
			RecTables.Properties.VariableNames;
		RecTables2.Properties.VariableUnits = ...
			RecTables.Properties.VariableUnits;
		RecTables2.Properties.Description = ...
			RecTables.Properties.Description;
		
		RecTables = RecTables2;
		
		clear RecTablesTemp RecTables2
		
		% Convert all non-date columns to double
		RecTables{:, ~ismember( ...
			RecTables.Properties.VariableNames, DateColumnNames)} ...
			= double(RecTables{:, ~ismember(...
			RecTables.Properties.VariableNames,DateColumnNames)});
		
		% Arbitrarily cutting GHI values lower than ~10
		% W/m2. These are usually errors.
		if ismember('GHI', RecTables.Properties.VariableNames)
			RecTables.GHI(RecTables.GHI<=10) = 0;
		end
		
		
		% Special data cleaning for Eplus Input % WIND When
		% the wind speed is zero, reset the wind direction
		% to the EnergyPlus default value (E+ Input-Output
		% Reference, table 44, pg 1741)
		if ismember('WSPD',RecTables.Properties.VariableNames)
			RecTables.WDR(RecTables.WSPD==0) = 180;
		end
		if ismember('WDR',RecTables.Properties.VariableNames)
			RecTables.WDR = floor(RecTables.WDR);
		end
		
		% % ATMOSPHERIC PRESSURE When the atm pressure is
		% zero, reset it to the EnergyPlus default value (E+
		% Input-Output Reference,tab. 44,pg 1741) = STP US
		% Standard Temperature and Pressure taken from
		% ASHRAE Handbook Fundamentals 2009, pg 1.1
		if ismember('ATMPR',RecTables.Properties.VariableNames)
			STP = 101325; % Pa
			% EnergyPlus only accepts pressure values
			% between 31000 Pa and 120000 Pa, so anything
			% below and above that is set to STP.
			RecTables.ATMPR(RecTables.ATMPR <= 31000 | ...
				RecTables.ATMPR >= 120000) = STP;
		end
		
		% Remove blanks from SourceCollect
		SourceCollect = SourceCollect(~cellfun(@isempty,SourceCollect));
		% Add the source name in the overall table
		% description
		if iscellstr(SourceCollect)
			Descriptor = strjoin(unique(SourceCollect),'and');
		else
			Descriptor = char(SourceCollect);
		end
		
		% Assign descriptor to overall table
		RecTables.Properties.Description = Descriptor;
		
		% Send overall table to annualiser
		AnnualisedTable = FileAnnualiser(RecTables);
		
		% Annualiser should return N separate tables nested
		% in a struct, where N is the number of unique years
		% that were present in the input, i.e. RecTables.
		
		% Cycle through all the individual tables in the
		% struct
		for k = 1:numel(AnnualisedTable)
			
			% Assign Variable names and Units
			VarNames = ...
				AnnualisedTable{k}.Properties.VariableNames;
			VarUnits = ...
				AnnualisedTable{k}.Properties.VariableUnits;
			Description = ...
				AnnualisedTable{k}.Properties.Description;
			
			AnnualisedTable{k}.Properties.VariableNames = ...
				VarNames;
			AnnualisedTable{k}.Properties.VariableUnits = ...
				VarUnits;
			AnnualisedTable{k}.Properties.Description = ...
				Description;
			
			clear DataAnnualOut
		end
		
		
		% The master name of the current city is either the
		% same as the station name (if station name is only
		% one word) or the first word in the station name.
		CityMasterName = strjoin(SplitName(1:StNames),'_');
		
		% Logical to record succesful creation of files
		SuccessfulCopy = false(1,length(AnnualisedTable));
		
		% To see if certain years have to be discarded or
		% not
		DiscardYear = false(1,length(AnnualisedTable));
		
		for annidx = 1:length(AnnualisedTable)
			
			% Unique year in the TempData
			UniqueYearTemp = unique(AnnualisedTable{annidx}.Year);
			
			% If the year is a leap year, remove Feb 29
			if pvl_leapyear(UniqueYearTemp)
				AnnualisedTable{annidx} = ...
					LeapRemover(AnnualisedTable{annidx});
			end
			
			% First check the size of the whole table
			SizeOfAMY = size(AnnualisedTable{annidx},1);
			
			if SizeOfAMY < cutoff*size(MasterTable,1)
			% Do not create new file because the table
			% is too small
			fprintf(fIDlog,['AMY file not created for year %d ', ...
				'because the records are less than %d %% of ', ...
				'8760. \r\n'], UniqueYearTemp, cutoff*100);
			SuccessfulCopy(annidx) = false;
			continue

			else
			% Check the quality of only the following
			% key variables
			KeyVariables = {'TDB','GHI','RH'};

			% Cycle through all variables to check the
			% quality of the key variables
			for subidx = 1:size(AnnualisedTable{annidx},2)

			% Take the table in question
			CheckTable = AnnualisedTable{annidx};

			% If the current variable is a key
			% variable
			if ismember(CheckTable.Properties.VariableNames, ...
					KeyVariables)
				% Check for NaN values
				CheckNaNs = isnan(CheckTable{:,subidx});

				if sum(CheckNaNs) > cutoff*length(CheckNaNs)
				% There are more than the
				% tolerable number of NaN
				% values.
				fprintf(fIDlog, ['More than %d %% ', ...
					'NaN values found in year %d for ', ...
					'station %s. \r\n'], ...
					cutoff*100, UniqueYearTemp, StationID);
				fprintf(fIDlog, '%% %% Variable name is %s. %% %%', ...
					CheckTable.Properties.VariableNames{subidx});
				DiscardYear(annidx) = true;
				end
			else
				continue
			end

			end

			if DiscardYear(annidx)
				fprintf(fIDlog,['AMY file not created for year %d', ...
					' of station %s due to an insufficient ', ...
					'number of records. \r\n'], ...
					UniqueYearTemp, StationID);
				SuccessfulCopy(annidx) = false;
				continue
			end


			% If both quality checks were passed, then
			% regularise the table before writing to
			% file

			% If the year is shorter than standard
			% (8760), regularise it,
			if  SizeOfAMY < size(MasterTable,1) && ...
					SizeOfAMY > cutoff*size(MasterTable,1)

				fprintf(fIDlog, ['Calling FILE REGULARISER ', ...
					'function. \r\n']);

				[FileRegulariserOutput, AnnualisedTable{annidx}] = ...
					evalc(['FileRegulariser(AnnualisedTable{', ...
					'annidx},MasterTable)']);
				% The evalc function captures all the
				% output from

				FileRegulariserOutput = strrep( ...
					FileRegulariserOutput, '\','\\');
				FileRegulariserOutput = strrep( ...
					FileRegulariserOutput, '..\\..\\','\\');

				fprintf(fIDlog,FileRegulariserOutput);

			end
				
			end
			
			% Path to the new file, i.e. output file
			FilePathNew = strrep(FilePathMaster,FileNameMaster,...
				[CityMasterName,'_',Descriptor,'_', num2str(unique(...
				AnnualisedTable{annidx}.Year),'%4d'),'.epw']);
			
			% Check if file already exists. If destroyer was
			% true, the existing file should have been
			% deleted already.
			FileExister = exist(FilePathNew,'file')==2;
			
			if FileExister
				% Do not create new file
				fprintf(fIDlog,['File named %s already exists in the', ...
					' target folder, skipping the write function.\n'], ...
					FilePathNew);
				SuccessfulCopy(annidx) = false;
				
			else
				% Create a new file
				
				try
					% The evalc function captures all the
					% output from the evaluated function.
					[FileCopierOutput, SuccessfulCopy(annidx)] = ...
						evalc(['FileCopier(FilePathMaster, FilePathNew,', ...
						'MasterTable, AnnualisedTable{annidx})']);
					
					% Change the backslashes in the output
					% to be double for writing out.
					FileCopierOutput = strrep(FileCopierOutput,'\','\\');
					FileCopierOutput = strrep(FileCopierOutput, ...
						'..\\..\\','\\');
					
				catch err_copy
					SuccessfulCopy(annidx) = false;
					fprintf('%s \r\n',err_copy.message);
					for e=1:length(err_copy.stack)
						fprintf('%s at %i\n', err_copy.stack(e).name, ...
							err_copy.stack(e).line);
					end
					fprintf('Error in folder %s .\r\n', ...
						TopFolderDirList{CityNumber})
				end
				
				if SuccessfulCopy(annidx)
					fprintf(fIDlog,FileCopierOutput);
				else
					fprintf(fIDlog,['There was some error in ', ...
						'copying file %s using master file %s .\n'], ...
						FilePathNew, FilePathMaster);
				end
			end
			
		end
		
	end
	
end

% if any(~SuccessfulCopy)
%     fprintf(fIDlog,'One or more files could not be
%     successfully copied. \r\n');
% end

fclose all;

clear DataTableDuplicates* dupidx Idx2 k nanidx
clear TempData RecTablesOld