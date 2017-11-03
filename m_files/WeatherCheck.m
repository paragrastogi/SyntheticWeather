function TableOut = WeatherCheck(tabin, varargin)

% This function performs some basic checks on an incoming table of
% weather parameters. If it finds phyiscally unfeasible values, it corrects
% them and then interpolates using the specified interpolation method.

p = inputParser;
p.FunctionName = 'WeatherCheck';

% Incoming table
addRequired(p,'tabin',@istable)
addParameter(p,'Silent',true,@islogical)
addParameter(p,'interpmeth','nearest',@ischar)

parse(p,tabin,varargin{:});

tabin = p.Results.tabin;
Silent = p.Results.Silent;
interpmeth = p.Results.interpmeth;

InfVals = table2array(varfun(@isinf,tabin));
for k = 1:size(InfVals,2)
	tabin{InfVals(:,k),k} = NaN;
end

tabin.Properties.VariableNames = ...
	upper(tabin.Properties.VariableNames);

% Temperature - TDB (dry bulb) and TDP (dew point) in KELVIN

if ismember('TDB', tabin.Properties.VariableNames)
	
	% The temperature might be in Celsius. In this case, the temperature
	% will definitely be always less than 200. In case the temperature is
	% indeed in Kelvin, it is impossible that the temperature will ever be
	% below 200 K, i.e. -73°C.
	
	if any(tabin.TDB<=200)
		tabin.TDB = tabin.TDB + 273.15;
		ConvertBack = true;
		if ~Silent
			fprintf(['Looks like the input temperature vector was ', ...
				'in Celsius. Converting to Kelvin for more checks ', ...
				'and censoring. \r\n'])
		end
	else
		ConvertBack = false;
	end
	
	
	% If any temperature value is below 200 or above 333.15, that would
	% imply a Celsius temperature of -73 or 60 degrees! These values need
	% to be censored.
	tabin.TDB(tabin.TDB<=200 | tabin.TDB>=333.15) = NaN;
	
	% Replace the NaNs using nearest neighbour interpolation.
	NaNfinder = (isnan(tabin.TDB));
	tabin.TDB(NaNfinder) = interp1(find(~NaNfinder), ...
		tabin.TDB(~NaNfinder), find(NaNfinder), interpmeth);
	
	% Convert back to Celsius
	if ConvertBack
		tabin.TDB = tabin.TDB - 273.15;
	end
end


if ismember('TDP', tabin.Properties.VariableNames)
	
	if any(tabin.TDP<=200)
		tabin.TDP = tabin.TDP + 273.15;
		ConvertBack = true;
		if ~Silent
			fprintf(['Looks like the input temperature vector was ', ...
				'in Celsius. Converting to Kelvin for more checks ', ...
				'and censoring. \r\n'])
		end
	else
		ConvertBack = false;
	end
	
	tabin.TDP(tabin.TDP<=200 | tabin.TDP>=333.15) = NaN;
	
	% Replace the NaNs using nearest neighbour interpolation.
	NaNfinder = (isnan(tabin.TDP));
	tabin.TDP(NaNfinder) = interp1(find(~NaNfinder), ...
		tabin.TDP(~NaNfinder), find(NaNfinder),interpmeth);
	
	% Convert back to Celsius
	if ConvertBack
		tabin.TDP = tabin.TDP - 273.15;
	end
end

% Relative Humidity (%) and Humidity Ratio (dimensionless fraction)

if ismember('RH', tabin.Properties.VariableNames)
	
	% If any temperature value is below 200 or above 333.15, that would
	% imply a Celsius temperature of -73 or 60 degrees! These values need
	% to be censored.
	tabin.RH(tabin.RH<=0 | tabin.RH>=100) = NaN;
	
	% Replace the NaNs using nearest neighbour interpolation.
	NaNfinder = (isnan(tabin.RH));
	tabin.RH(NaNfinder) = interp1(find(~NaNfinder), ...
		tabin.RH(~NaNfinder), find(NaNfinder),interpmeth);
	
end

if ismember('W', tabin.Properties.VariableNames)
	
	%if all(isnan(tabin.W))
	
	if sum(~isnan(tabin.W))>(size(tabin,1)/10)
		% If any temperature value is below 200 or above 333.15, that would
		% imply a Celsius temperature of -73 or 60 degrees! These values need
		% to be censored.
		tabin.W(tabin.W<=0 | tabin.W>=1) = NaN;
		
		% Replace the NaNs using nearest neighbour interpolation.
		NaNfinder = (isnan(tabin.W));
		tabin.W(NaNfinder) = interp1(find(~NaNfinder), ...
			tabin.W(~NaNfinder), find(NaNfinder),interpmeth);
	else
		fprintf(['All humidity ratio values in Weather Check ', ...
			'are NaN. Ignoring this.\r\n'])
	end
end

% Following EPlus convention, the wind direction is 0° only when
% the data says so. If there is no wind, i.e. WSPD == 0, then wind
% direction is 180°. Wind direction 360° is due north, reset that
% to 0° to avoid unnecessary fractions in direction.
if ismember('WSPD',tabin.Properties.VariableNames) && ...
		ismember('WDR',tabin.Properties.VariableNames)
	tabin.WDR(tabin.WSPD==0) = 180;
	tabin.WDR(tabin.WDR==360) = 0;
	tabin.WDR(isnan(tabin.WSPD)) = NaN;
end

if ismember('GHI',tabin.Properties.VariableNames)
	tabin.GHI(tabin.GHI<=1) = 0;
end

TableOut = tabin;

end