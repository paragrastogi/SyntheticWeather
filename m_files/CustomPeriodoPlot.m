% © All rights reserved. 
% ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Interdisciplinary Laboratory of Performance-Integrate Design, 2016
% Parag Rastogi
% See the LICENSE.TXT file for more details.

function CustomPeriodoPlot(plotser,filepath,varargin)

p = inputParser;
p.FunctionName = 'CustomPeriodoPlot';


addRequired(p,'plotser',@(x) (isnumeric(x) & isreal(x) & ...
    (isrow(x) | iscolumn(x))))
addRequired(p,'filepath',@(x) (ischar(x) | iscell(x)))

% addParameter(p,'ylabelCustom','',@(x) (ischar(x) | iscell(x)))
addParameter(p,'titleCustom','',@(x) (ischar(x) | iscell(x)))
addParameter(p,'visible','on',@ischar)
addParameter(p,'logplot',true,@islogical)
addParameter(p,'numwaves',1:1:876,@isnumeric)

parse(p, plotser, filepath, varargin{:})

plotser = p.Results.plotser;

if iscell(p.Results.filepath)
    filepath = p.Results.filepath{1};
else
    filepath = p.Results.filepath;
end

logplot = p.Results.logplot;

if iscell(p.Results.titleCustom)
    titleCustom = p.Results.titleCustom{1};
else
    titleCustom = p.Results.titleCustom;
end

visible = p.Results.visible;
numwaves = p.Results.numwaves;

if isrow(numwaves)
	numwaves = (numwaves)';
end
	

DefaultColours
set(0, 'defaultAxesFontName', 'Helvetica')
set(0, 'defaultAxesUnits', 'normalized')

% Initialising the Fourier transform

% Custom frequency vector, indicating the specific
% frequencies we are interested in. The denominator is the
% number of samples, i.e. 8760. At a sampling frequency of 1
% sample/hour, the numerator indicates the number of waves
% (over the sample period) that you are looking at. For
% example, 365 indicates we looking at 365 waves per year,
% approx. wavelength of 1 day.
% The default number of waves plotted are from 1 to 876,
% i.e. periods of 1 year to 10 hours
% numwaves = 1:1:876;
freq = numwaves./8760;

[~,yData] = prepareCurveData([],plotser);

[pxx, ~, ~] = periodogram(yData, [], freq , [], ...
    'ConfidenceLevel', 0.95, 'psd');
% if logplot
% 	pxxcum = cumsum(log10(pxx));
% else
	pxxcum = cumsum(pxx);
% end

PlotHandle = figure('visible',visible);
ax(1) = subplot(2,1,1);

p1 = plot(ax(1),freq,pxx,'-');

if logplot
	ax(1).YAxis.Scale = 'log';
    ax(1).YLabel.String = 'log PSD'; 
else
    ax(1).YLabel.String = 'PSD [dB/Hz]'; 
end

hold(ax(1),'on'); ax(1).Box = 'on';
fmax = freq(pxx==max(pxx));
line1 = plot(ax(1), fmax, 0);
line1.Marker = 'o'; 
hold(ax(1),'off'); ax(1).XMinorTick = 'on';

ax(1).XLim = [-0.001, round(max(freq),3)];
if max(numwaves)==400
	ax(1).XTick = [1/8760, 365/8760];
else
	ax(1).XTick = [1/720, 1./(48:-12:12)];
end
ax(1).XTickLabel = cellstr(num2str(round(1./ax(1).XTick)'));
ax(1).XGrid = 'on'; 
ax(1).XLabel.String = 'Wavelength [hrs]';

line1.Color = lgrey; 
line1.MarkerFaceColor = line1.Color;
line1.MarkerEdgeColor = line1.Color;
% line1.LineWidth = 2;

p1.Color = blue;
p1.LineWidth = 1.5;

%%%%%%%

ax(2) = subplot(2,1,2);
p2 = plot(ax(2), freq, pxxcum); 
hold(ax(2),'on'); ax(2).Box = 'on';
ax(2).XMinorTick = 'on';

ax(2).YLabel.String = 'Cum. PSD';
ax(2).XGrid = 'on'; 
ax(2).XLabel.String = 'Wavelength [hrs]';
ax(2).XLim = [-0.001, round(max(freq),3)];
if max(numwaves)==400
	ax(2).XTick = [1/8760, 365/8760];
else
	ax(2).XTick = [1/720, 1./(48:-12:12)];
end
ax(2).XTickLabel = cellstr(num2str(round(1./ax(2).XTick)'));

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1], ...
    'Box','off','Visible','off', ...
	'Units','normalized', 'clipping' , 'off');
text(0.5, 1,titleCustom, 'HorizontalAlignment','center', ...
    'VerticalAlignment', 'top', 'FontSize', 20)

p2.Color = blue;
p2.LineWidth = 2;

for k = 1:2
	
	ax(k).YColor = blackest;
	ax(k).XColor = blackest;
	ax(k).FontSize = 20;
	ax(k).LabelFontSizeMultiplier = 1.25;
	ax(k).TitleFontSizeMultiplier = 1.25;
	
end

SaveThatFig(PlotHandle, filepath)

 
% Re-colour in Grey and Save again
p1.Color = grey;
p2.Color = grey;

% % Upper Plot
% p1.Color = rgb2gray(p1.Color);
% line1.Color = rgb2gray(line1.Color);
% line1.MarkerFaceColor = rgb2gray(line1.MarkerFaceColor);
% line1.MarkerEdgeColor = rgb2gray(line1.MarkerEdgeColor);
% 
% % Lower Plot
% p2.Color = rgb2gray(p2.Color);

SaveThatFig(PlotHandle, [filepath, '-BW'], ...
	'printfig', false)

end