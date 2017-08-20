% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate 
% Design, 2016 
% Parag Rastogi
% See the LICENSE.TXT file for more details.

% Compose years for an assessment of solar production
% extremes. These years are "solar quantile" years or
% "temperature quantile years". They SHOULD NOT be used for
% energy studies but will work for Solar Energy or Daylight
% applications. Basically, these time series have not
% actually occured - because the hottest july of the last
% 30-odd years is not automatically followed by the hottest
% august. It doesn't make sense to simulate a time-dependent
% / dynamic process with an "impossible" time series)

clear; clc;

% Load recorded data
load('D:\Weather_Data\WeatherFilesForAnalysis\GEN\RecordedData_GEN.mat')

% Delete data points after jan 31, 2016
RecTables(RecTables.Year==2016 & ...
	ismember(RecTables.Month, 2:12),:) = [];

% Collect the statistics, year-by-year for each month
[ghistats.mean, ghistats.sum, ghistats.min, ...
	ghistats.max, ghistats.names] = grpstats(RecTables.GHI, ...
	{RecTables.Year,RecTables.Month}, ...
	{@nanmean, @nansum, @nanmin, @nanmax, 'gname'});

[tdbstats.mean,tdbstats.median, tdbstats.min, ...
	tdbstats.max, tdbstats.names] = grpstats(RecTables.TDB, ...
	{RecTables.Year,RecTables.Month}, ...
	{@nanmean, @nanmedian, @nanmin, @nanmax, 'gname'});
% Convert the names to numbers so it's easier to use them.
ghistats.names = cell2mat(cellfun(@str2double, ...
	ghistats.names, 'UniformOutput', 0));

colltable = [ghistats.names,ghistats.sum,tdbstats.median];
uniyears = unique(ghistats.names(:,1), 'stable');

for k = 1:12
	[temp1, temp2] = max(colltable(colltable(:,2)==k,3));
	maxes.ghisum(k,:) = temp1;
	maxes.ghiidx(k,:) = temp2;	
	tempyear = colltable(colltable(:,2)==k,1);
	maxes.ghiyears(k,:) = tempyear(temp2);
	
	[temp1, temp2] = min(colltable(colltable(:,2)==k,3));
	mines.ghisum(k,:) = temp1;
	mines.ghiidx(k,:) = temp2;	
	tempyear = colltable(colltable(:,2)==k,1);
	mines.ghiyears(k,:) = tempyear(temp2);
	
	[temp1, temp2] = sort(colltable(colltable(:,2)==k,3));
	midpointer = ceil(length(temp1)/2);
	middles.ghisum(k) = temp1(midpointer);
	middles.ghiidx(k) = temp2(midpointer);	
	tempyear = colltable(colltable(:,2)==k,1);
	middles.ghiyears(k) = tempyear(temp2(midpointer));
	
	[temp1, temp2] = max(colltable(colltable(:,2)==k,4));
	maxes.tdbmed(k) = temp1;
	maxes.tdbidx(k) = temp2;	
	tempyear = colltable(colltable(:,2)==k,1);
	maxes.tdbyears(k) = tempyear(temp2);
	
	[temp1, temp2] = min(colltable(colltable(:,2)==k,4));
	mines.tdbmed(k) = temp1;
	mines.tdbidx(k) = temp2;	
	tempyear = colltable(colltable(:,2)==k,1);
	mines.tdbyears(k) = tempyear(temp2);
	
	[temp1, temp2] = sort(colltable(colltable(:,2)==k,4));
	midpointer = ceil(length(temp1)/2);
	middles.tdbmed(k) = temp1(midpointer);
	middles.tdbidx(k) = temp2(midpointer);	
	tempyear = colltable(colltable(:,2)==k,1);
	middles.tdbyears(k) = tempyear(temp2(midpointer));
end

plothand(1) = figure; ax(1) = gca; ax(1).Box = 'on'; hold(ax(1),'on');
g1 = plot(maxes.ghisum, '^:');
g2 = plot(middles.ghisum,'s--');
g3 = plot(mines.ghisum,'o-');
hold(ax(1), 'off')
leg = legend(ax(1),[g1, g2, g3], 'max','mid','min');
ax(1).Title.String = 'GHI';

plothand2 = figure; ax(2) = gca; ax(2).Box = 'on'; hold(ax(2),'on');
t1 = plot(maxes.tdbmed, '^:');
t2 = plot(middles.tdbmed,'s--');
t3 = plot(mines.tdbmed,'o-');
hold(ax(2), 'off')
leg2 = legend(ax(2),[g1, g2, g3], 'max','mid','min');
ax(2).Title.String = 'TDB';

%%

hidx.ghi = RecTables.Year==maxes.ghiyears(1) & ...
		RecTables.Month==1;
midx.ghi = RecTables.Year==middles.ghiyears(1) & ...
		RecTables.Month==1;
lidx.ghi = RecTables.Year==mines.ghiyears(1) & ...
		RecTables.Month==1;
	
hidx.tdb = RecTables.Year==maxes.tdbyears(1) & ...
		RecTables.Month==1;
midx.tdb = RecTables.Year==middles.tdbyears(1) & ...
		RecTables.Month==1;
lidx.tdb = RecTables.Year==mines.tdbyears(1) & ...
		RecTables.Month==1;
	
for k = 2:12
	hidx2 = RecTables.Year==maxes.ghiyears(k) & ...
		RecTables.Month==k;
	hidx.ghi = hidx.ghi | hidx2;
	midx2 = RecTables.Year==middles.ghiyears(k) & ...
		RecTables.Month==k;
	midx.ghi = midx.ghi | midx2;
	lidx2 = RecTables.Year==mines.ghiyears(k) & ...
		RecTables.Month==k;
	lidx.ghi = lidx.ghi | lidx2;
	
	hidx3 = RecTables.Year==maxes.tdbyears(k) & ...
		RecTables.Month==k;
	hidx.tdb = hidx.tdb | hidx3;
	midx3 = RecTables.Year==middles.tdbyears(k) & ...
		RecTables.Month==k;
	midx.tdb = midx.tdb | midx3;
	lidx3 = RecTables.Year==mines.tdbyears(k) & ...
		RecTables.Month==k;
	lidx.tdb = lidx.tdb | lidx3;
end

clear hidx2 midx2 lidx2

%%

HighYear = RecTables(hidx.ghi,:);
MidYear = RecTables(midx.ghi,:);
LowYear = RecTables(lidx.ghi,:);
HighYear = sortrows(HighYear,{'Month','Day'}, ...
	{'ascend','ascend'});
MidYear = sortrows(MidYear,{'Month','Day'}, ...
	{'ascend','ascend'});
LowYear = sortrows(LowYear,{'Month','Day'}, ...
	{'ascend','ascend'});

pathEPWfolder = fullfile('D:','Weather_Data', ...
	'WeatherFilesForAnalysis','GEN');
pathMasterFile = fullfile(pathEPWfolder, ...
	'GEN_Meteonorm.epw');
MasterTable = WeatherFileParseEPWMeteonorm(pathMasterFile);

pathNewFile = fullfile(pathEPWfolder, 'GEN_highrad.epw');
if size(HighYear,1)~=8760
	HighYear(HighYear.Month==2 & HighYear.Day==29,:) = [];
end

Copied = FileCopier(pathMasterFile, pathNewFile, ...
    MasterTable, HighYear, 'OneYearCheck', false);

pathNewFile = fullfile(pathEPWfolder, 'GEN_midrad.epw');
if size(MidYear,1)~=8760
	MidYear(MidYear.Month==2 & ...
		MidYear.Day==29,:) = [];
end

Copied = FileCopier(pathMasterFile, ...
    pathNewFile, MasterTable, MidYear, ...
	'OneYearCheck', false);

pathNewFile = fullfile(pathEPWfolder, 'GEN_lowrad.epw');
if size(LowYear,1)~=8760
	LowYear(LowYear.Month==2 & LowYear.Day==29,:) = [];
end
Copied = FileCopier(pathMasterFile, ...
    pathNewFile, MasterTable, LowYear, ...
	'OneYearCheck', false);

%%

HighYear = RecTables(hidx.tdb,:);
MidYear = RecTables(midx.tdb,:);
LowYear = RecTables(lidx.tdb,:);
HighYear = sortrows(HighYear,{'Month','Day'}, ...
	{'ascend','ascend'});
MidYear = sortrows(MidYear,{'Month','Day'}, ...
	{'ascend','ascend'});
LowYear = sortrows(LowYear,{'Month','Day'}, ...
	{'ascend','ascend'});

pathNewFile = fullfile(pathEPWfolder, 'GEN_hightemp.epw');
if size(HighYear,1)~=8760
HighYear(HighYear.Month==2 & HighYear.Day==29,:) = [];
end
Copied = FileCopier(pathMasterFile, pathNewFile, ...
    MasterTable, HighYear, 'OneYearCheck', false);

pathNewFile = fullfile(pathEPWfolder, 'GEN_midtemp.epw');
if size(MidYear,1)~=8760
MidYear(MidYear.Month==2 & MidYear.Day==29,:) = [];
end
Copied = FileCopier(pathMasterFile, pathNewFile, ...
    MasterTable, MidYear, 'OneYearCheck', false);

pathNewFile = fullfile(pathEPWfolder, 'GEN_lowtemp.epw');
if size(LowYear,1)~=8760
LowYear(LowYear.Month==2 & LowYear.Day==29,:) = [];
end
Copied = FileCopier(pathMasterFile, pathNewFile, ...
    MasterTable, LowYear, 'OneYearCheck', false);