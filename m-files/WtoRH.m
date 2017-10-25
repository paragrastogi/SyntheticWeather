function RHo = WtoRH(Win,TDBin,varargin)

% % This function was written by Parag Rastogi, 09 July 2015

% % Convert Humidity Ratio to Relative Humidity

% % This function reverses the calculation procedure for converting
% % Relative Humidity to Humidity Ratio described in the ASHRAE
% % fundamentals handbook and implemented in the complementary function
% % RHtoW.

% % All equation and page numbers given below refer to ASHRAE Standard
% % 90.1:2010 Fundamentals (SI Ed.)

p = inputParser;
p.FunctionName = 'WtoRH';

% Relative Humidity
addRequired(p,'Win',@isnumeric)
% Temperature (Absolute preferred)
addRequired(p,'TDBin',@isnumeric)

% US Standard Atmosphere pressure is 101.325 kPa (at sea level)
% ASHRAE 2009 Fundamentals, page 1.1
pres_def = 101.325*1000; % [Pa]
addParameter(p,'Pressure',pres_def,@isnumeric)

addParameter(p,'Silent',true,@islogical)

parse(p,Win,TDBin,varargin{:});

Win = p.Results.Win;
TDBin = p.Results.TDBin;
Pressure = p.Results.Pressure;
Silent = p.Results.Silent;

% In case there is an Inf value in the incoming vector, converting it to
% NaN makes error correction simpler.
TDBin(isinf(TDBin)) = NaN;
Win(isinf(Win)) = NaN;

% The temperature might be in Celsius. In this case, the temperature will
% definitely be always less than 200. In case the temperature is indeed in
% Kelvin, it is unlikely that the temperature will ever be below 200 K,
% i.e. -73°C.
if any(TDBin<=200)
    
    TDBin = TDBin + 273.15;
    
    if ~Silent
        fprintf(['Looks like the input temperature vector was in ', ...
            'Celsius. Converting to Kelvin. \r\n'])
    end
    
end

if length(Win)~=length(TDBin)
    error(['Check lengths of incoming Humidity and Temperature', ...
        ' vectors in call to WtoRH function. \r\n'])
end

% Put the variables in a table
WeatherTable = table(TDBin,Win);
WeatherTable.Properties.VariableNames = {'TDB','W'};
WeatherTable.Properties.VariableUnits = {'K','dimless'};

% Send the variables for cleaning
WeatherTable = WeatherCheck(WeatherTable);

TDBin = WeatherTable.TDB;
Win = WeatherTable.W;

clear WeatherTable

N = length(Win);
if length(Pressure)==1
    Pressure = repmat(Pressure,N,1);
end

%Humidity ratio W, [unitless fraction] % Equation (22), pg 1.8
p_w = ((Win./0.621945).*Pressure) ./ (1 + (Win./0.621945));

% Equations 5 and 6 are for calculating the saturation pressure of water
% vapour.

% Constants for Eq. 5 
% Temperature -200°C to 0°C
C1 = -5.6745359*10^+03; C2 = 6.3925247; 
C3 = -9.6778430*10^-03; C4 = 6.2215701*10^-07;
C5 = 2.0747825*10^-09; C6 = -9.4840240*10^-13;
C7 = 4.1635019;

% Constants for Eq. 6
% Temperature 0°C to 200°C
C_8 =  -5.8002206*10^+03; C_9 =  1.3914993;
C_10 = -4.8640239*10^-02; C_11 = 4.1764768*10^-05;
C_12 = -1.4452093*10^-08; C_13 = 6.5459673;

% This is to distinguish between the two versions of equation 5.
Ice = TDBin<=273.15; 

lnp_ws = NaN(size(TDBin));
% Eq. 5, pg 1.2
lnp_ws(Ice) = C1./TDBin(Ice) + C2 + C3.*TDBin(Ice) + ...
    C4.*TDBin(Ice).^2 + C5.*TDBin(Ice).^3 + C6*TDBin(Ice).^4 + ...
    C7*log(TDBin(Ice));
% Eq. 6, pg 1.2
lnp_ws(~Ice) = C_8./TDBin(~Ice) + C_9 + C_10.*TDBin(~Ice) + ...
    C_11.*TDBin(~Ice).^2 + C_12.*TDBin(~Ice).^3 + ...
    C_13.*log(TDBin(~Ice));  
% Temperature in the above formula must be absolute, i.e. in Kelvin

% Continuing from eqs. 5 and 6
p_ws = exp(lnp_ws); % [Pa]

phi = p_w./p_ws; % [Pa] Formula(24), pg 1.8

RHo = phi.*100; % Relative Humidity from fraction to percentage

% Correct instances of RH exceeding 100% or below 0%
RHo(RHo<=0 | RHo>=100) = nan;
RHo = CustomInterpolate(RHo,'nearest');

end
