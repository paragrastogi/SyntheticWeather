function SMY_Create(folderstorun, pathEPWtopfolder, ...
	pathOUTtopfolder, varargin)


% This script will write EPW files from synthetic weather
% files. The synthetic weather files should be in MAT
% format, unless otherwise specified in the optional inputs.

p = inputParser;
p.FunctionName = 'SMY_Create';

% The first input is a cell array of strings specifying
% which folders to run
addRequired(p, 'folderstorun', @iscellstr)

% This is where the EPW files should be
addRequired(p, 'pathEPWtopfolder', @ischar)

% This is where the new files will be written, in a subfolder
% based on the 'folderstorun' input.
addRequired(p, 'pathOUTtopfolder', @ischar)

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

parse(p,folderstorun, pathEPWtopfolder, ...
	pathOUTtopfolder, varargin{:})

folderstorun = p.Results.folderstorun;
destroyer = p.Results.destroyer;
pathEPWtopfolder = p.Results.pathEPWtopfolder;
pathOUTtopfolder = p.Results.pathOUTtopfolder;
nameLOGfile = p.Results.nameLOGfile;
synformat = p.Results.synformat;
scenario = p.Results.scenario;
sublabel = p.Results.sublabel;
srcEPWfile = p.Results.srcEPWfile;

if exist(pathOUTtopfolder, 'dir') ~= 7
	mkdir(pathOUTtopfolder)
end

if exist('LogFiles', 'dir') ~= 7
	mkdir('LogFiles')
end

% % Open a log file to record all output message
fIDlog = fopen(nameLOGfile,'w');
fprintf(fIDlog,['Commencing SMY_Create script at ', ...
	'%s \r\n'], char(datetime));

if destroyer
	fprintf(fIDlog, ['Destroyer is set to true.\r\n', ...
		'All existing SMY files will be deleted.\r\n', ...
		'I AM BECOME DEATH, THE DESTROYER OF ', ...
		'EXISTING FILES... \r\n']);
else
	fprintf(fIDlog, ['destroyer is false, no files will', ...
		' be deleted.\r\n This isn''t so much fun... though', ...
		' it is probably faster.\r\n']);
end

% Path to the parent or top-level folder, called
% EPWTopFolder Make a list of files in the top folder
TopFolderDirList = dir(pathEPWtopfolder);
% Dir function picks up everything, including '.', '..',
% and files. In this case, we are only interested in the
% city folders. Remove folder names that start with a dot
TopFolderDirList(strncmp({TopFolderDirList.name}, ...
	'.', 1)) = [];
% Remove files
TopFolderDirList = TopFolderDirList([ ...
	TopFolderDirList.isdir]);
% Keep only the names from the list
TopFolderDirList = {TopFolderDirList.name};

% Make a struct with the keywords for TMY and AMY - usually
% based on the source of the data
FileParserNames = FileParserKeywordLister;

fprintf(fIDlog, ['WARNING\r\n', ...
	'All input files are treated as text files.\r\n', ...
	'This could cause problems with XLS and XLSX files.\r\n']);

for CityNumber = 1:length(TopFolderDirList)
	% Cycle through the entire top-level folder
	
	% First set the paths to city subfolders
	CityDirPathIn = fullfile(pathEPWtopfolder, ...
		TopFolderDirList{CityNumber});
	CityDirPathOut = fullfile(pathOUTtopfolder, ...
		TopFolderDirList{CityNumber});
	
	
	% Skip folders that have not been requested
	if ~ismember(TopFolderDirList{CityNumber},folderstorun)
		continue
	end
	
	if ~ismember(TopFolderDirList{CityNumber}, ...
			TopFolderDirList)
		fprintf(fIDlog,'Unknown folder name %s .\r\n', ...
			TopFolderDirList{CityNumber});
		continue
	end
	
	if exist(CityDirPathOut, 'dir') ~=7
		mkdir(CityDirPathOut)
	end
	
% 	try
		
		fprintf(fIDlog,'Processing folder %s\n', ...
			TopFolderDirList{CityNumber});
		
		% Search directory of current city for all files. If
		% a file other than a weather file exists, it will
		% be handled by the last 'else' condition in the
		% file handling loop below.
		
		
		% Find only EPW files
		WthrFilesListEPW = dir(fullfile(CityDirPathIn,'*.epw'));
		% Keep only the names from the list
		WthrFilesListEPW = {WthrFilesListEPW.name};
		
		% Split the station names coming from the EPW files
		SplitStations = cellfun(@(x) (strsplit(x,{'_','.'})), ...
			WthrFilesListEPW,'UniformOutput',0);
		
		% Look to see if the second word in the split names
		% matches any of the TMY keywords
		SecondWord = cell2mat(cellfun(@(x) (strcmp(x{2}, ...
			FileParserNames.TMYKeyword)), SplitStations, ...
			'UniformOutput',0));
		
		if any(SecondWord)
			SplitStations = unique(cellfun(@(x) (x{1}), ...
				SplitStations,'UniformOutput',0));
			NoOfStations = 1;
		else
			SplitStations = unique(cellfun(@(x) ([x{1},'_',x{2}]), ...
				SplitStations,'UniformOutput',0));
			NoOfStations = numel(SplitStations);
		end
		
		for StIdx = 1:NoOfStations
			
			if NoOfStations>1
				% There is more than one station, which
				% means that the file source is at the third
				% location, i.e. there are two station names
				StNames = 2;
			else
				% Else there is only one word in the station
				% name
				StNames = 1;
			end
			
			% This is the ID of the current station.
			StationID = SplitStations{StIdx};
			
			% This is the 'master name' which is the same as
			% the epw folder name.
			CityMasterName = TopFolderDirList{CityNumber};
			
			fprintf(fIDlog,'Processing station %s\r\n', StationID);
			
			% Find only EPW files
			WthrFilesListEPW = dir(fullfile(CityDirPathIn,'*.epw'));
			
			% Remove folders, just in case they were caught
			% by the dir function
			WthrFilesListEPW = WthrFilesListEPW(~[ ...
				WthrFilesListEPW.isdir]);
			
			% Remove any files whose name begins with a dot
			WthrFilesListEPW = WthrFilesListEPW(~strncmp({...
				WthrFilesListEPW.name},'.',1));
			
			% Keep only the names from the list
			WthrFilesListEPW = {WthrFilesListEPW.name};
			
			% Keep only the files from the CURRENT STATION
			WthrFilesListEPW = WthrFilesListEPW(~cellfun(@isempty, ...
				cellfun(@(x) (strfind(x,StationID)), ...
				WthrFilesListEPW,'UniformOutput',0)));
			
			clear TMYfind AMYfind	
            		
			% Search for all TMY keywords
			TMYfind = ~cellfun(@isempty,(strfind(WthrFilesListEPW, ...
				FileParserNames.TMYKeyword{1})));
			for k = 2:length(FileParserNames.TMYKeyword)
				TMYfind = TMYfind | ~cellfun(@isempty,(strfind( ...
					WthrFilesListEPW, FileParserNames.TMYKeyword{k})));
			end
			
			% Special check to remove the NRELIndia files
			% from the list of TMY files. This is because
			% they contain a TMY keyword (NREL), which gets
			% caught by the TMYfind check.
			NRELIndiafind = ~cellfun(@isempty,(strfind( ...
				WthrFilesListEPW, 'NRELIndia')));
			
			TMYfind = TMYfind & ~NRELIndiafind;
            
            if any(TMYfind)
                % Keep only the TMY files from this list
                WthrFilesListEPW = WthrFilesListEPW(TMYfind);
            else
                % The file probably didn't contain any of the expected TMY
                % keywords - just pick one randomly.
                WthrFilesListEPW = WthrFilesListEPW(randi(length(WthrFilesListEPW)));
            end
			
			
			SplitNames = cellfun(@(x) (strsplit(x,{'_','.'})), ...
				WthrFilesListEPW, 'UniformOutput',0);
			
			EPWFileSources = cellfun(@(x) (x{StNames+1}), ...
				SplitNames, 'UniformOutput',0);
			
			
			% If there is more than one TMY file available
			% for a given station, then favour the ones
			% coming from NREL (EnergyPlus website)
			
            if isempty(srcEPWfile)
                if length(EPWFileSources)>1
                    % There is more than one TMY file available
                    % for the given station
                    MasterIdx = find(cell2mat(cellfun(@isempty,...
                        strfind(EPWFileSources,'Meteonorm'), ...
                        'UniformOutput',0)) ,1,'first');
                    % Source of master file
                    FileSrcMaster = char(SplitNames{MasterIdx}(StNames+1));
                    % Master file name
                    FileNameMaster = WthrFilesListEPW{MasterIdx};
                    % Path to master file
                    FilePathMaster = fullfile(CityDirPathIn,FileNameMaster);
                    % Read the file using USDOE parser
                    mastertable = WeatherFileParseEPWUSDOE...
                        (FilePathMaster);
                    % Add source name to the description
                    mastertable.Properties.Description = ...
                        [mastertable.Properties.Description, ...
                        '_', FileSrcMaster];
                else
                    % There is only one TMY file available for
                    % the given station This file is most likely
                    % the Meteonorm file
                    MasterIdx = 1;
                    % Source of master file
                    FileSrcMaster = char(SplitNames{MasterIdx}(StNames+1));
                    % Master file name
                    FileNameMaster = WthrFilesListEPW{MasterIdx};
                    % Path to master file
                    FilePathMaster = fullfile(CityDirPathIn, FileNameMaster);
                    % Read the file using Meteonorm parser
                    mastertable = WeatherFileParseEPWMeteonorm...
                        (FilePathMaster);
                    % Add source name to the description
                    mastertable.Properties.Description = ...
                        [mastertable.Properties.Description, ...
                        '_', FileSrcMaster];
                end
                
            else
                [~, FileNameMaster, ~] = fileparts(srcEPWfile);
                % Path to master file
                FilePathMaster = srcEPWfile;
                % Read the file using Meteonorm parser
                mastertable = WeatherFileParseEPWUSDOE...
                    (FilePathMaster);
                % Add source name to the description
                mastertable.Properties.Description = ...
                    [mastertable.Properties.Description, ...
                    '_', FileNameMaster];
            end
                
			
			fprintf(fIDlog,'Master Data is from %s \n', ...
				FileNameMaster);
			
			
			% Find the synthetic data file in the current
			% folder, corresponding to this station.
			if strcmpi(synformat,'mat')
				if strcmpi(scenario, 'syn')
					MATfileSyn = dir(fullfile(CityDirPathIn, ...
						sprintf('%s*Syn.mat', ...
						FileNameMaster(1:end-4))));
				elseif strcmpi(scenario, 'rcp45')
					MATfileSyn = dir(fullfile(CityDirPathIn, ...
						sprintf('%s*rcp45.mat', ...
						FileNameMaster(1:end-4))));
				elseif strcmpi(scenario, 'rcp85')
					MATfileSyn = dir(fullfile(CityDirPathIn, ...
						sprintf('%s*rcp85.mat', ...
						FileNameMaster(1:end-4))));
				end
			elseif strcmpi(synformat,'csv')
				MATfileSyn = dir(fullfile(CityDirPathIn, ...
					sprintf('%s*Syn.csv', ...
					FileNameMaster(1:end-4))));
			end
			% Remove MAT files with the 'extras' suffix from
			% the list since they contain, well, extra data.
			MATfileSyn = MATfileSyn(cellfun(@isempty,...
				regexp({MATfileSyn.name},'extras')));
			% Keep only the names from the list
			MATfileSyn = char({MATfileSyn.name});
			
			% If no actual weather data files are found,
			% skip this city
            if isempty(MATfileSyn)
                fprintf(fIDlog,['No synthetic data files found for', ...
                    ' base file %s\n'], FileNameMaster);
                continue
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
				load(fullfile(CityDirPathIn, MATfileSyn));
				
			elseif strcmpi(synformat,'csv')
				readsynCSV = csvread(CityDirPathIn, MATfileSyn);
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
			PickBootLen = length(syndata.Year) / ...
				length(UniqueSynYears) / N;
			SynYears = repmat(UniqueSynYears, PickBootLen, 1);
			
			% New batch has a set of sub labels shifted by the value of
			% sublabel.
            temp1 = (0:1:(PickBootLen-1)) + sublabel;
            temp2 = repmat(temp1,length(UniqueSynYears),1);
            temp2 = temp2(:);
            % SubLabelArray = cell(length(UniqueSynYears)*PickBootLen,1);
            SubLabelArray = cellfun(@num2str, num2cell(temp2), 'UniformOutput', 0); % {};
            
            %             [temp1, temp2] = ndgrid(cellfun(@num2str, ...
            %                 num2cell(SynYears), 'UniformOutput', 0), SubLabelArray);
            % 			SynYearsNames = strcat(temp1(:), '-', temp2(:));
            SynYearsNames = strcat(cellfun(@num2str, ...
                num2cell(SynYears), 'UniformOutput', 0), '-', SubLabelArray);
			
			OutFilePath = fullfile(pathOUTtopfolder, CityMasterName);
			% Path to the new file, i.e. output file
			FilePathNew = fullfile(OutFilePath, cellfun(@(x) ...
				sprintf('%s_%s_%08s.epw', StationID, ...
				scenario, x), SynYearsNames, 'UniformOutput', 0));
			
			% Check if file already exists. If destroyer was
			% true, the existing file should have been
			% deleted already.
			FileExister = cellfun(@(x) exist(x,'file')==2, ...
				FilePathNew);
			
			
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
			
			% Read in the header from the EPW master file as
			% a series of lines.
			
			% Number of header lines in a regular EPW file
			numEPWhead = 8;
			% Save the header lines
			HeaderSave = cell(numEPWhead,1);
			
			% Open file for low-level read/write operations
			fileIDMaster = fopen(FilePathMaster,'r');
			
			% Copy Header lines almost verbatim (except
			% COMMENTS 1 field)
			
			for el = 1:numEPWhead
				
				% Get current line of master file
				FileLine = fgetl(fileIDMaster);
				
				if strcmpi(FileLine(1:10),'COMMENTS 1')
					% Add name of new source to header line
					% COMMENTS 1
					newsrcadd = upper(['; The following ', ...
						'data columns replaced with synthetic ', ...
						'values: ''GHI'', ''DNI'', ''DHI'', ', ...
						'''RH'', ''TDB''.']);
					if strcmpi(scenario,'syn')
						newsrcadd = [newsrcadd, ...
							upper(['The years are dummy values, ', ...
							'they do not imply an actual ', ...
							'prediction for a specific year. ', ...
							'The synthetic values are not ', ...
							'predictions either, so they do not ', ...
							'represent the exact value at some ', ...
							'specific future hour.'])];
					else
						newsrcadd = [newsrcadd, ...
							upper(['The years represent the ', ...
							'prediction for a specific year in the ', ...
							'future. However, a prediction should ', ...
							'not be taken too literally, since it ', ...
							'comes from an approximate model based', ...
							' on the state of the art in 2005.', ...
							'All data was downloaded from the CORDEX', ...
							' project web site with the help of ', ...
							'Georgios Mavormatidis, EMPA Duebendorf ', ...
							'(Zurich) in October 2015.'])];
					end
					FileLineNew = strcat(FileLine, newsrcadd);
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
			
		end
		
		
% 	catch err
% 		
% 		fprintf('%s \r\n', err.message);
% 		for e=1:length(err.stack)
% 			fprintf('%s at %i\n', err.stack(e).name, ...
% 				err.stack(e).line);
% 		end
% 		fprintf('Error in folder %s .\r\n', ...
% 			TopFolderDirList{CityNumber})
% 	end
end

fclose(fIDlog);
% Close any remaining files, especially the log file.
fclose all;

end