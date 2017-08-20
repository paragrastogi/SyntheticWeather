% Create Synthetic Files

% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function DataTableOut = FileRegulariser(DataTableIn,DataTableMaster)

% Check the size of the incoming table
if size(DataTableIn,1) == 8760
    fprintf('Incoming table has a height of 8760 rows, are you sure ')
    fprintf('it needed to be sent to the File Regulariser function?\n')
end

% Size of output table
N = 8760;
% Record the actual year of the incoming table
ActualYear = DataTableIn.Year(1);
% This function here determines whether a given year is a leap year or not
LeapYearNess = DetermineLeapYear(ActualYear);
if LeapYearNess
    % Modify the value of year to prevent calculation errors due to Feb 29
    ModYear = ActualYear-1;
    DataTableIn = LeapRemover(DataTableIn);
else
    % Else, the new 'fake' year is the same as the incoming ones
    ModYear = ActualYear;
end


DateIdx(1) = find(strcmp(DataTableIn.Properties.VariableNames,'Year'));
DateIdx(2) = find(strcmp(DataTableIn.Properties.VariableNames,'Month'));
DateIdx(3) = find(strcmp(DataTableIn.Properties.VariableNames,'Day'));
DateIdx(4) = find(strcmp(DataTableIn.Properties.VariableNames,'Hour'));
DateIdx(5) = find(strcmp(DataTableIn.Properties.VariableNames,'Minute'));

DateCols = DataTableIn(:,DateIdx);
[A,~,~]= setxor(DateIdx,1:size(DataTableIn,2));
DataCols = DataTableIn(:,A);
DataLabels = DataTableIn.Properties.VariableNames(A);
DataTableIn(:,1:5) = DateCols;
DataTableIn(:,6:size(DataTableIn,2)) = DataCols;
DateIdx = 1:5;
DataTableIn.Properties.VariableNames = [{'Year','Month','Day',...
    'Hour','Minute'},DataLabels] ;

% Preallocate output table of zeros, and assing it all the metadata from
% the incoming table.
DataTableOut = array2table(zeros(N,size(DataTableIn,2)));
DataTableOut.Properties.VariableNames = DataTableIn.Properties.VariableNames;
DataTableOut.Properties.VariableUnits = DataTableIn.Properties.VariableUnits;
DataTableOut.Properties.VariableDescriptions = ...
    DataTableIn.Properties.VariableDescriptions;
DataTableOut.Properties.Description = DataTableIn.Properties.Description;

if ~any(strcmp(DataTableIn.Properties.VariableNames,'DateNumStamps'))
    DataTableIn.DateNumStamps = datenum(DataTableIn.Year, ...
        DataTableIn.Month, DataTableIn.Hour, DataTableIn.Minute, ...
        zeros(size(DataTableIn,1),1), zeros(size(DataTableIn,1),1));
end

% Beginning date of the year
DateBeg = datenum(ModYear,1,1,0,0,0);
% Create the datenum stamps for the outgoing table
DataTableOut.DateNumStamps = (arrayfun(@(x) (addtodate(DateBeg,x, ...
    'hour')),1:N))';

DataTableOut.(1) = repmat(ModYear,N,1);
DataTableOut.(2) = month(DataTableOut.DateNumStamps);
DataTableOut.(3) = day(DataTableOut.DateNumStamps);
DataTableOut.(4) = hour(DataTableOut.DateNumStamps);
DataTableOut.(5) = zeros(N,1);

DataTableOut.Properties.VariableNames(1:5) = {'Year','Month','Day',...
    'Hour','Minute'} ;

if any(DataTableOut.Hour==0)
    OrigDescr = DataTableOut.Properties.Description;
    DataTableOut.Properties.Description = 'MATLAB';
    DataTableOut = DateFormatter(DataTableOut);
    DataTableOut.Properties.Description = OrigDescr;
end


% Do a setxor of the date stamps of the two tables, without the year.
% SETXOR will find all elements that are not in the intersection of the
% two vectors. In this case, we use three vectors together, so we use the
% 'rows' option. This ensures that each row (of 3 elements) is taken as a
% single entity.
[xOrInOut,idxOutxOr,~] = setxor(DataTableOut{:,2:4}, ...
    DataTableIn{:,DateIdx(2:4)},'rows');

% [~,idxInIs,idxOutIs] = intersect(DataTableIn{:,2:4},DataTableOut{:,2:4},'rows');
% TempInDateNum = arrayfun(@(x) (addtodate(x,ModYear-ActualYear,'year')),...
%     DataTableIn.DateNumStamps);

% Dates that exist in the incoming and outgoing table
[~ ,idxInIs2,idxOutIs2] = intersect(DataTableIn{:,DateIdx(2:4)}, ...
    DataTableOut{:,2:4}, 'rows');

% Variable names that don't match with the master list will be removed
[~,IdxUn] = setdiff(DataTableIn.Properties.VariableNames, ...
    DataTableMaster.Properties.VariableNames);

if isempty(xOrInOut) %&& EndMissing
    %DataTableOut(1:end-1,:) = DataTableIn(1:end,:);
    fprintf('The incoming table seems to have no missing dates.')
    fprintf('Returning the same data table.\n')
    DataTableOut = DataTableIn;
    return
    
else
    if length(idxOutxOr)==1 && (idxOutxOr == N)
        % If only the end of the file was missing, then copy the rest of the table
        % and exit.
        DataTableOut(idxOutxOr,:) = DataTableIn(idxOutxOr-1,:);
        return
        
    else
        
        % Either there is more than one missing index or that index is not
        % equal to N, the last datestamp.
        for misidx = 1:size(DataTableIn,2)
            
            if strcmp(DataTableIn.Properties.VariableNames{misidx}, ...
                    'DateNumStamps')
                continue
            end
            
            if ismember(misidx,DateIdx) || ismember(misidx,IdxUn)
                continue
            end
            
            % Place nan values in the extracted vector where there
            % are holes in the incoming table
            % Fill the rest of the extraced vector with the
            % corresponding values from the incoming table.
            ExtractedVector = nan*ones(N,1);
            ExtractedVector(idxOutIs2) = ...
                DataTableIn{idxInIs2,misidx};
            
            % Index of current variable in incoming table
            CurrentVar = DataTableIn.Properties.VariableNames{misidx};
            % Index of current variable in master table
            VarIdxMaster = (strcmp(CurrentVar, ...
                DataTableMaster.Properties.VariableNames));
            
            % Interpolate over the holes. The inputs for the interpolate
            % function are as follows:
            % 1. The indices which are not holes, ~isnan
            % 2. The values in the incoming table at those (non-hole)
            % indices, values(~isnan)
            % 3. The indices of the holes, idxOut or isnan
            % The output are the values at the holes
            
            NoNaNs = ~isnan(ExtractedVector);
            LastVal = find(NoNaNs,1,'last');
            
            if LastVal<N
                % The last genuine value is not at the end of the year
                % Find the index of the last genuine value
                ExtractedVector(LastVal:N) = ...
                    DataTableMaster{LastVal:N,VarIdxMaster};
            end
            
            % Find NaN values again
            NoNaNs = ~isnan(ExtractedVector);
            
            % First NaN value
            % Find the index of the first genuine value
            FirstVal = find(NoNaNs,1,'first');
            
            if FirstVal>1
                % The first genuine value is not at the beginning of the
                % year
                ExtractedVector(1:FirstVal) = ...
                    DataTableMaster{1:FirstVal,VarIdxMaster};
            end
            
            NoNaNs = ~isnan(ExtractedVector);
            ExNaNs = isnan(ExtractedVector);
            
            if (sum(NoNaNs))<(24*30)
                % The incoming vector has less than a month of genuine
                % values
                DataTableOut{:,misidx} = ...
                    DataTableMaster{:,VarIdxMaster};
            else
                % At least some values are not NaNs
                try
                    DataTableOut{ExNaNs,misidx} = ...
                        interp1(find(NoNaNs), ....
                        ExtractedVector(NoNaNs), ...
                        find(ExNaNs),'pchip');
                catch err
                    disp(err)
                    fprintf('Error thrown by FileRegulariser function ')
                    fprintf('for Master Description %s using PCHIP INTERP.\n', ...
                        DataTableMaster.Properties.Description)
                    try
                        DataTableOut{(isnan(ExtractedVector)),misidx}...
                            = interp1(find(NoNaNs), ....
                            ExtractedVector(NoNaNs), ...
                            find(isnan(ExtractedVector)),'nearest');
                    catch err
                        fprintf('Error thrown by FileRegulariser function ')
                        fprintf('for Master Description %s using NEAREST INTERP.\n', ...
                            DataTableMaster.Properties.Description)
                        disp(err)
                    end
                end
            end
            
            DataTableOut{NoNaNs,misidx} = ...
                ExtractedVector(NoNaNs);
            
            if any(idxOutxOr==N)
                
                % If there were more than one missing index, and the last of
                % those was the end index, then that needs to be handled with
                % this special statement below.
                DataTableOut(N,misidx) = DataTableOut(N-1,misidx);
                % Remove the last index from the list of missing indices
                % idxOut = idxOut(idxOut~=N);
                % xOrInOut = xOrInOut(idxOut~=N);
            end
            
        end
        
        
    end
    
end

% If the incoming year was a leap year, then the actual year was modified
% by subtracting one. This is now re-added.
if LeapYearNess
    DataTableOut.Year = DataTableOut.Year+1;
end

DataTableOut(:,IdxUn) = [];

end