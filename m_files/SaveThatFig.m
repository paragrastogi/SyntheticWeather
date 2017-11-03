function SaveThatFig(plothand, pathFIGfile, varargin)

p = inputParser;
p.FunctionName = 'SaveThatFig';

addRequired(p,'plothand')
addRequired(p,'pathFIGfile',@ischar)
addParameter(p,'printfig',true, @islogical)
addParameter(p,'printpdf', true, @islogical)
addParameter(p,'printpng', true, @islogical)
addParameter(p,'orient', 'landscape', @ischar)
addParameter(p,'changecolours', true, @islogical)

parse(p, plothand, pathFIGfile, varargin{:})

plothand = p.Results.plothand;
pathFIGfile = p.Results.pathFIGfile;
printfig = p.Results.printfig;
printpdf = p.Results.printpdf;
printpng = p.Results.printpng;
orient = p.Results.orient;
changecolours = p.Results.changecolours;

DefaultColours

% Set some default properties.
if changecolours
plothand.Color = whitest;
set(gca, 'Color', whitest);
set(gca, 'XColor', blackest);
set(gca, 'YColor', blackest);
end

set(gca, 'box', 'on');

set(gca,'FontName', 'Helvetica')

% Save figure file first
if printfig
if ~isempty(strfind(pathFIGfile, 'bw')) || ...
		~isempty(strfind(pathFIGfile, 'BW'))
	fprintf(['Not saving fig file for B/W version ', ...
		'of the figure. \r'])
else
	savefig(plothand, [pathFIGfile,'.fig'], 'compact')
end
end

plothand.Units = 'centimeters';

% Print PNG before changing any settings, force it to be
% portrait.
% plothand.PaperOrientation = 'portrait';
if printpng
% Delete the existing png first, since sometimes MATLAB
% refuses to overwrite existing PNG files.
if exist([pathFIGfile,'.png'],'file')==2
	delete([pathFIGfile,'.png'])
end

if strcmp(orient,'landscape')
plothand.PaperOrientation = 'landscape';
plothand.Position = [0, 0, 29.7, 21];
plothand.PaperPosition = [0, 0, 29.7, 21];
plothand.PaperSize = [29.7, 21];
plothand.PaperPositionMode = 'auto';
elseif strcmp(orient,'portrait')
plothand.PaperOrientation = 'portrait';
plothand.Position = [0, 0, 21, 29.7];
plothand.PaperPosition = [0, 0, 21, 29.7];
plothand.PaperSize = [21, 29.7];
plothand.PaperPositionMode = 'auto';
elseif strcmp(orient,'none')
	fprintf('Not changing paper size or orientation.\r\n')
end
saveas(plothand, pathFIGfile, 'png')
end

% Store the originals so the figure can be reverted back to
% them.
pos = plothand.Position;
paperpos = plothand.PaperPosition;


if printpdf
% Change settings for the PDF.

if strcmp(orient,'landscape')
% This is the default orientation as well
plothand.PaperOrientation = 'landscape';
plothand.Position = [0, 0, 29.7, 21];
plothand.PaperPosition = [0, 0, 29.7, 21];
plothand.PaperPositionMode = 'auto';

elseif strcmp(orient,'portrait')
plothand.PaperOrientation = 'portrait';
plothand.Position = [0, 0, 21, 29.7];
plothand.PaperPosition = [0, 0, 21, 29.7];
plothand.PaperPositionMode = 'auto';

else strcmp(orient,'none')
	fprintf('Not changing paper size or orientation.\r\n')
end

print(plothand,'-dpdf', '-painters', '-fillpage', pathFIGfile)

end


if ismember(orient,{'landscape','portrait'})
% Revert back to the original settings
plothand.Position = pos;
plothand.PaperPosition = paperpos;
end

end