% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function CustomPdfPlot(plotser,plotf,FilePath,varargin)

% % This function brings together the code needed to plot the descriptive
% % plots used in the CreateSyntheticFiles script

p = inputParser;
p.FunctionName = 'CustomPdfPlot';

addRequired(p,'plotser',@isstruct)
addRequired(p,'plotf',@ischar)
addRequired(p,'FilePath',@(x) (ischar(x) | iscell(x)))
addParameter(p,'xlabelCustom','',@(x) (ischar(x) | iscell(x)))
addParameter(p,'titleCustom','',@(x) (ischar(x) | iscell(x)))
addParameter(p,'visible','on',@ischar)
addParameter(p,'plotCI',false,@islogical)
addParameter(p,'CCdata',false,@islogical)

parse(p, plotser, plotf, FilePath, varargin{:})

plotser = p.Results.plotser;

if iscell(p.Results.FilePath)
	FilePath = p.Results.FilePath{1};
else
	FilePath = p.Results.FilePath;
end

if iscell(p.Results.xlabelCustom)
	xlabelCustom = p.Results.xlabelCustom{1};
else
	xlabelCustom = p.Results.xlabelCustom;
end

if iscell(p.Results.titleCustom)
	titleCustom = p.Results.titleCustom{1};
else
	titleCustom = p.Results.titleCustom;
end

visible = p.Results.visible;
CCdata = p.Results.CCdata;
% Plot Confidence Intervals or NOT
plotCI = p.Results.plotCI;
% Convert the incoming plotf string to uppercase
plotf = p.Results.plotf;


DefaultColours
% set(0, 'defaultAxesFontName', 'Helvetica')
set(0, 'defaultAxesUnits', 'normalized')

PlotHandle = figure('visible',visible);
% Get current axes, and turn on HOLD
ax = gca; hold(ax, 'on');

if CCdata
	% Climate change data
	
	for k = 1:size(plotser.rcp45.(plotf).x,1)
	H2(k) = plot(plotser.rcp45.(plotf).x(k,2:end), ...
		plotser.rcp45.(plotf).y(k,:));
	H3(k) = plot(plotser.rcp85.(plotf).x(k,2:end), ...
		plotser.rcp85.(plotf).y(k,:));
		
	H2(k).LineStyle = '-';
	H2(k).LineWidth = 1;
	H2(k).Color = orange;
	H3(k).LineStyle = '-';
	H3(k).LineWidth = 1;
	H3(k).Color = red;
	end

else
	% Synthetic data			
	for k = 1:size(plotser.syn.(plotf).x,1)
	H2(k) = plot(ax,plotser.syn.(plotf).x(k,2:end), ...
		plotser.syn.(plotf).y(k,:));
	H2(k).LineStyle = '-';
	H2(k).LineWidth = 1;
	H2(k).Color = orange;
	end
end


% Recorded data
if isfield(plotser,'rec')
	
H1 = plot(plotser.rec.(plotf).x(2:end), ...
	plotser.rec.(plotf).y);

H1.LineStyle = '-';
H1.LineWidth = 3;
H1.Color = blue;
end

if isfield(plotser,'tmy')
	% TMY data
	H4 = plot(plotser.tmy.(plotf).x(2:end), ...
		plotser.tmy.(plotf).y);
	H4.Color = grey;
	H4.LineStyle = ':';
	H4.LineWidth = 3;
end

if plotCI
	H5(1) = plot(plotser.rec.(plotf).x(2:end), ...
		plotser.rec.(plotf).flo);
	H5(2) = plot(plotser.rec.(plotf).x(2:end), ...
		plotser.rec.(plotf).fup);
end

hold(ax, 'off');

ax.FontSize = 20;
ax.LabelFontSizeMultiplier = 1.125;
ax.TitleFontSizeMultiplier = 1.125;

if isfield(plotser,'rec')
	LegendCustom = {'Recorded'};
	LegendVec = H1;
else
	LegendCustom = {};
	LegendVec = [];
end

if ~CCdata
	
if isfield(plotser,'tmy')
	LegendCustom = [LegendCustom, {'Synthetic','TMY'}];
	LegendVec = [LegendVec, H2(1), H4];
else
	LegendCustom = [LegendCustom, {'Synthetic'}];
	LegendVec = [LegendVec, H2(1)];
end

else
	
if isfield(plotser,'tmy')
	LegendCustom = [LegendCustom, {'RCP45','RCP85','TMY'}];
	LegendVec = [LegendVec, H2(1), H3(1), H4];
else
	LegendCustom = [LegendCustom, {'RCP45','RCP85'}];
	LegendVec = [LegendVec, H2(1), H3(1)];
end

end

leg = legend(ax,LegendVec,LegendCustom, 'Location', 'northeast');
leg.FontSize = ax.FontSize;


if plotCI
	H5(1).Color = lgrey;
	H5(2).Color = lgrey;
	H5(1).LineWidth = 4;
	H5(2).LineWidth = 4;
end

ax.Title.String = titleCustom;

ax.YLabel.String = 'PDF';

if isempty(xlabelCustom)
	if strcmpi(plotf,'TDB')
		ax.XLabel.String = 'Temperature [°C]';
	elseif strcmpi(plotf,'GHI')
		ax.XLabel.String = 'Solar Radiation [W/m^2]';
	elseif strcmpi(plotf,'W')
		ax.XLabel.String = 'Humidity Ratio [unitless]';
	elseif strcmpi(plotf,'RH')
		ax.XLabel.String = 'Relative Humidity [%]';
	end
else
	ax.XLabel.String = xlabelCustom;
end

SaveThatFig(PlotHandle, FilePath)


% Re-colour in Grey and Save again
if isfield(plotser,'rec')
	H1.Color = grey;
end

for k = 1:length(H2)
	H2(k).Color = grey;
end

if CCdata
	% RCP8.5, if it exists. RCP 4.5 takes the properties of
	% the SYN series.
	for k = 1:length(H3)
		H3(k).Color = grey;
		H3(k).LineStyle = '-.';
	end
end

if isfield(plotser,'tmy')
	H4.Color = grey;
end

if plotCI
	H5(1).Color = grey;
	H5(2).Color = grey;
end

SaveThatFig(PlotHandle, [FilePath,'-BW'], ...
	'printfig', false)

end