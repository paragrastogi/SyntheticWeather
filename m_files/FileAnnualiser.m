% Create Synthetic Files

% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function DataTableOut = FileAnnualiser(DataTableIn)

ActualYears = unique(DataTableIn.Year);

DataTableOut{1,length(ActualYears)} = array2table(zeros(8760,size(DataTableIn,2)));

for k = 1:length(ActualYears)
    
    TempTable = DataTableIn(DataTableIn.Year==ActualYears(k),:);
    DataTableOut{1,k} = TempTable;

end

end