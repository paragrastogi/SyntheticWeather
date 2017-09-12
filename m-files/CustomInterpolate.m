function yOut = CustomInterpolate (yIn, method)

p = inputParser;
p.FunctionName = 'CustomInterpolate';

addRequired(p,'yIn',@isnumeric)
addRequired(p,'method',@ischar)

parse(p, yIn, method);

yIn = p.Results.yIn;
method = p.Results.method;

% Sample points with valid values
x = find(~isnan(yIn));
% The actual values at these points
v = yIn(~isnan(yIn));
% The values that were rejected
xq = find(isnan(yIn));
% Interpolate using the nearest neighbour method
yq = interp1(x, v, xq, method);
% Assign the interpolated values back to the original vector
yIn(isnan(yIn)) = yq;
% Put the interpolated values back into the structure
yOut = yIn;

end