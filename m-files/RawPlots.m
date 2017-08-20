% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function RawPlots(datain,SavePath,wthrfilename, varargin)

p = inputParser;
p.FunctionName = 'RawPlots';

addRequired(p,'datain',@(x) (istable(x) || isstruct(x)))
addRequired(p,'SavePath',@(x) (ischar(x) || iscell(x)))
addRequired(p,'wthrfilename',@(x) (ischar(x) || iscell(x)))
addParameter(p,'plotser','Unknown',@ischar)
addParameter(p,'visible','on',@ischar)

parse(p, datain, SavePath, wthrfilename, varargin{:})

if istable(datain)
	datain = p.Results.datain;
elseif isstruct(datain)
	datain = struct2table(p.Results.datain);
end

if iscell(p.Results.SavePath)
    SavePath = p.Results.SavePath{1};
else
    SavePath = p.Results.SavePath;
end

if iscell(p.Results.wthrfilename)
    wthrfilename = p.Results.wthrfilename{1};
else
    wthrfilename = p.Results.wthrfilename;
end

plotser = p.Results.plotser;
visible = p.Results.visible;

set(0, 'defaultAxesFontName', 'Helvetica')
set(0, 'defaultAxesUnits', 'normalized')

% Convert all variable names to lower case to avoid errors
% later.
datain.Properties.VariableNames = ...
	lower(datain.Properties.VariableNames);

% If the incoming weather file name contains the file
% extension, then remove it.
if strcmpi(wthrfilename(end-4:end),'.epw')
	wthrfilename = wthrfilename(1:end-4);
end

% Should be 8760 or some multiple thereof
N = size(datain,1); 
t = (1:N)'; % Julian Hour Index
nm = 12; % Number of months
mt = datetime(2015,1:nm,1); % MONTH number
DefaultColours

% DO not print EPS files if the incoming time series is too long.
if size(datain.tdb,1)>10000
    EPSprint = false;
else
    EPSprint = true;
end


% Dry-bulb Temperature in Celsius
dstats.means.tdb = grpstats(datain.tdb, ...
	{datain.month, datain.day}, @nanmean);


% % Global Horizontal Solar Radiation
dstats.sums.ghi = grpstats(datain.ghi, ...
	{datain.month, datain.day}, @expfit);

% Humidity Ratio (converted from Relative Humidity in the
% TMY file)
if isfield(datain, 'rh')
	
    dstats.means.rh = grpstats(datain.rh, ...
		{datain.month, datain.day}, @nanmean);
elseif isfield(datain, 'w')

    dstats.means.w = grpstats(datain.w, ...
		{datain.month, datain.day}, @nanmean);
end


gsums_rep = reshape(repmat(dstats.sums.ghi,1,24)',[],1);

t = 1:8760;
if strcmpi(plotser, 'syn')
	ghires = reshape(datain.ghi, 8760, []);
else
	ghires = datain.ghi;
end

PlotHandle = figure('visible',visible);
[ax,line1,L2] = plotyy(t,ghires, t,gsums_rep);
if strcmpi(plotser, 'syn')
	for k = 1:length(line1)
		line1(k).LineStyle = 'none'; line1(k).Marker = '.';
		line1(k).Color = [lgrey, 0.1];
	end
	L2.LineWidth = 1.5; L2.Color = blue;
else
	line1.LineStyle = 'none'; line1.Marker = '.';
	line1.Color = lgrey;
	L2.LineWidth = 1.5; L2.Color = blue;
	
end
ax(1).XLim = [0, t(end)]; ax(2).XLim = [0, t(end)];
ax(1).XTick = 24*15:(24*30):24*345;
ax(1).XTickLabel = month(mt,'shortname');
ax(1).YLabel.String = 'Hourly Power (GHI) [W/m^2]';
ax(2).YLabel.String = 'Daily sum of GHI (Energy) [Wh/m^2]';
set(ax(2).YLabel, 'Units', 'Normalized', ...
	'Position', [1.15, 0.5, 0]);
% Make sure that the axes are black
ax(1).YColor = blackest; 
ax(1).YLabel.Color = ax(1).YColor;
ax(2).YColor = blackest; 
ax(2).YLabel.Color = ax(2).YColor;
ax(2).YLabel.Rotation = -90;
ax(2).YLabel.Position = ax(2).YLabel.Position - ...
	[0.025, 0, 0];

for k = 1:2
ax(k).FontSize = 20;
ax(k).LabelFontSizeMultiplier = 1.5;
ax(k).TitleFontSizeMultiplier = 1.5;
ax(k).YMinorTick = 'on'; 
ax(k).Box = 'on'; 
end

FigName = ['RawSolar', plotser];
FilePath = fullfile(SavePath,[FigName,'_', ...
	wthrfilename]);
SaveThatFig(PlotHandle, FilePath)

%%%%%%%%%%%%%%%

dmeans_rep = reshape(repmat(dstats.means.tdb,1,24)',[],1);

t = 1:8760;
if strcmpi(plotser, 'syn')
	tdbplot = reshape(datain.tdb, 8760, []);
else
	tdbplot = datain.tdb;
end

PlotHandle = figure('visible',visible);
line1 = plot(t,tdbplot);
ax = gca; hold(ax, 'on');
L2 = plot(t,dmeans_rep);

if strcmpi(plotser, 'syn')
	for k = 1:length(line1)
		% line1(k).LineStyle = 'none'; line1(k).Marker = '.';
		line1(k).Color = [lgrey, 0.1];
		line1(k).LineWidth = 0.5; 
		line1(k).LineStyle = '-';
		line1(k).Marker = 'none';
	end
	
	L2.LineWidth = 1.5; L2.Color = blue;
	L2.LineStyle = '-'; L2.Marker = 'none';
	
	title('Synthetic Time Series')
	
	legend([line1(1) L2], 'Hourly', 'Daily Means')
	
else
	line1.LineStyle = 'none'; line1.Marker = '.';
	line1.Color = lgrey;
	L2.LineWidth = 1.5; L2.Color = blue;
	
	title('Typical Time Series')	
	legend([line1 L2], 'Hourly', 'Daily Means')
end

ax.XLim = [0, t(end)];
ax.XTick = 24*15:(24*30):24*345;
ax.XTickLabel = month(mt,'shortname');
ax.YLabel.String = 'Dry bulb temperature (TDB) [^oC]';
ax.YColor = blackest; ax.YLabel.Color = ax.YColor;

for k = 1:1
ax(k).FontSize = 20;
ax(k).LabelFontSizeMultiplier = 1.5;
ax(k).TitleFontSizeMultiplier = 1.5;
ax(k).YMinorTick = 'on'; 
ax(k).Box = 'on'; 
end

FigName = ['RawTemperature', plotser];
FilePath = fullfile(SavePath,[FigName,'_', ...
	wthrfilename]);
SaveThatFig(PlotHandle, FilePath)

%%%%%%%%%%%%%%%%%%

if isfield(datain, 'W')
    
    Wdmeans_rep = reshape(repmat(dstats.means.w,24,1),[],1);
    
    PlotHandle = figure('visible',visible);
    line1 = plot(t,W, t,Wdmeans_rep);
    ax = gca;
    line1(1).LineStyle = 'none'; line1(1).Marker = '.';
    line1(1).Color = lgrey;
    line1(2).LineWidth = 1.5; line1(2).Color = blue;
    ax.XLim = [-48, max(t)+48];
    ax.XTick = 24*15:(24*30):24*345;
    ax.XTickLabel = month(mt,'shortname');
    ax.YTickLabel = cellstr(num2str((min(W):0.0001:max(W))','%3.2g'));
    ax.YLabel.String = 'Humidity Ratio (W) [unitless]';
    ax.YColor = blackest; ax.YLabel.Color = ax.YColor;
	
	for k = 1:1
		ax(k).FontSize = 20;
		ax(k).LabelFontSizeMultiplier = 1.5;
		ax(k).TitleFontSizeMultiplier = 1.5;
		ax(k).YMinorTick = 'on';
		ax(k).Box = 'on';
	end
	
    FigName = 'RawHumidRatio';
    FilePath = fullfile(SavePath, ...
		[FigName,'_',wthrfilename]);
	SaveThatFig(PlotHandle, FilePath)

    
elseif isfield(datain, 'RH')
    
    rhmeans_rep = reshape(repmat( ...
		dstats.means.rh,24,1),[],1);
    
    PlotHandle = figure('visible',visible);
    line1 = plot(t, RH, t,rhmeans_rep);
    ax = gca;
    line1(1).LineStyle = 'none'; line1(1).Marker = '.';
    line1(1).Color = lgrey;
    line1(2).LineWidth = 1.5; 
	line1(2).Color = blue;
    ax.XLim = [-48, max(t)+48];
    ax.XTick = 24*15:(24*30):24*345;
    ax.XTickLabel = month(mt,'shortname');
    ax.YTick = (10:10:100);
    ax.YTickLabel = cellstr(num2str((0:10:100)','%3.2g'));
    ax.YLabel.String = 'Relative Humidity (RH) [%]';
    ax.YColor = blackest; ax.YLabel.Color = ax.YColor;
    ax.Color = blackest;
	
	for k = 1:1
		ax(k).FontSize = 20;
		ax(k).LabelFontSizeMultiplier = 1.5;
		ax(k).TitleFontSizeMultiplier = 1.5;
		ax(k).YMinorTick = 'on';
		ax(k).Box = 'on';
	end
    
	FigName = ['RawRelHumid', plotser];
    FilePath = fullfile(SavePath, ...
		[FigName,'_',wthrfilename]);
	SaveThatFig(PlotHandle, FilePath)
    
end

end