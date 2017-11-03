% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function CustomFourierPlot(FitIn, DataTable, PlotVar, FilePath, varargin)

p = inputParser;
p.FunctionName = 'CustomFourierPlot';

addRequired(p,'FitIn',@isstruct)
addRequired(p,'DataTable',@isnumeric)
% Optional input specifying the predictor names. This is not needed if the
% TableIn input is a table with the correct predictor names.
addRequired(p,'PlotVar',@ischar)

% Path to save figures
addRequired(p,'FilePath',@ischar)
% Figures visible or not
addParameter(p,'visible','off',@ischar)
addParameter(p,'Silent',true,@islogical)

parse(p, FitIn, DataTable, PlotVar, FilePath, varargin{:})

FitIn = p.Results.FitIn;
DataTable = p.Results.DataTable;
PlotVar = p.Results.PlotVar;
FilePath = p.Results.FilePath;
visible = p.Results.visible;
Silent = p.Results.Silent;

% Call the default colours
DefaultColours
set(0, 'defaultAxesFontName', 'Helvetica')
set(0, 'defaultAxesUnits', 'normalized')

N = size(DataTable,1); % Should be 8760 or some multiple thereof
t = (1:N)'; % Julian Hour Index
nm = 12; % Number of months
mt = datetime(2015,1:nm,1); % Month number


PlotHandle = figure ('visible', visible);
ax = gca; hold(ax, 'on');

if strcmpi(PlotVar, 'TDB')
	% Fits are in Kelvin but plots should prefereably be in Celsius.
	p1 = plot(FourierFits.(PlotVar).evalmu - 273.15);
else
	p1 = plot(FourierFits.(PlotVar).evalmu);
end
p1.Colour = orange;

p2 = plot(DataTable.(PlotVar), '.');
p2.Colour = blue;

p2.MarkerSize = 4;

ax.XLim = [0, t(end)];
ax.XTick = 24*15:(24*30):24*345;
ax.XTickLabel = month(mt, 'shortname');

FindVarIdx = strcmp( DataTable.Properties.VariableNames, PlotVar);

if strcmpi(PlotVar, 'TDB')
	ax.YLabel.String = sprintf('Hourly dry bulb temperature (%s) [%s]', ...
		DataTable.Properties.VariableNames{FindVarIdx}, ...
		DataTable.Properties.VariableUnits{FindVarIdx});
elseif strcmpi(PlotVar, 'W')
	ax.YLabel.String = sprintf('Hourly humidity ratio (%s) [%s]', ...
		DataTable.Properties.VariableNames{FindVarIdx}, ...
		DataTable.Properties.VariableUnits{FindVarIdx});
end

SaveThatFig(PlotHandle, FilePath)
% print(PlotHandle,'-dpng','-cmyk', '-painters', FilePath)
% print(PlotHandle,'-depsc','-cmyk', '-painters', FilePath)
% EPStoPDF(FilePath); delete([FilePath,'.eps']);

p1.Color = rgb2gray(p1.Color);
p2.Color = rgb2gray(p2.Color);

% Print the B&W version
SaveThatFig(PlotHandle, [FilePath, '-BW'], 'printfig', false);
% print(PlotHandle,'-depsc','-cmyk', '-painters', [FilePath,'-BW'])
% EPStoPDF([FilePath,'-BW','.eps']); delete([FilePath,'-BW','.eps']);

end