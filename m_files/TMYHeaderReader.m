function OutStruct = TMYHeaderReader(FilePath, varargin)

p = inputParser;
p.FunctionName = 'TMYHeaderReader';

addRequired(p,'filename',@(x) (ischar(x) | (iscell(x) && numel(x)==1)))
addParameter(p,'Silent',true,@islogical)

parse(p,FilePath, varargin{:})

% If incoming filename is a cell array of one element, then convert it to
% a character array containing that element
if ischar(FilePath)
    FilePath = p.Results.filename;
elseif ~ischar(FilePath)
    FilePath = p.Results.filename{1};
end

Silent = p.Results.Silent;

% % This function reads the header lines of a EPW format file and returns
% % two types of information:
% % 1. Station metadata
% % 2. Station design data

NumberHeaderEPW = 8; % Number of header lines in EPW files
fileID = fopen(FilePath);
if fileID <3
    OutStruct = struct();
    fprintf('File could not be opened: %s \n',FilePath)
    return
end
HeaderLines = cell(NumberHeaderEPW,1);
for hdridx = 1:NumberHeaderEPW
    HeaderLines{hdridx,1} = fgetl(fileID);
end
fclose(fileID);

% % A complete description of the various fields is available on the
% % EnergyPlus website. The assignments in this script were done based on
% % the instructions in the 'Weather Data Reference'. Specifically, pages
% % 09-13.

% % The description of the header lines in the EnergyPlus 'Weather Data
% % Reference' is outdated and does not correspond exactly to the ASHRAE
% % design conditions actually found in the weather file.

% % The following script is based on some detective work using the 2009
% % ASHRAE fundamentals handbook.

% % For details on the data read below, see the ASHRAE Fundamentals
% % Handbook. Information below taken from 2009 edition, pages 14.1-14.3,
% % Table 1.



% The first header line has the following elements:
% Location Latitude = vector or scalar latitude in decimal degrees
% (positive is northern hemisphere)
% Location Longitude = vector or scalar longitude in decimal degrees
% (positive is east of prime meridian)
% Location Altitude = vector or scalar site elevation (or altitude)
% in meters.
% Split header line using comma as delimiter
LocationInfo = strsplit(HeaderLines{1},',');
% Assign to relevant row in Overall List struct
OutStruct.WMO = (LocationInfo{6});
OutStruct.Latitude = str2double(LocationInfo{7});
OutStruct.Longitude = str2double(LocationInfo{8});
OutStruct.TZ = str2double(LocationInfo{9}); % Time zone, difference from UTC
OutStruct.Altitude = str2double(LocationInfo{10}); % Site altitude



% Split header line using comma as delimiter
DesignInfo = strsplit(HeaderLines{2},',');

% Find the different blocks
for k = 1:length(DesignInfo)    
    if strcmpi(DesignInfo{k},'Heating')
        HeatBreaker = k;
    elseif strcmpi(DesignInfo{k},'Cooling')
        CoolBreaker = k;        
    elseif strcmpi(DesignInfo{k},'Extremes')
        ExtremeBreaker = k;        
    end    
end

% If blocks are not found, output error message.
% Error messages are suppressed if Silent is true
if (exist('HeatBreaker','var')~=1)
    HeatInfo = zeros(1,15);
    if ~Silent
        fprintf('Heating design information not found in file.\n')
    end
else
    HeatInfo = str2double({DesignInfo{HeatBreaker+1:CoolBreaker-1}});
end

if (exist('CoolBreaker','var')~=1)
    CoolInfo = zeros(1,32);
    if ~Silent
        fprintf('Cooling design information not found in file.\n')
    end
else
    CoolInfo = str2double({DesignInfo{CoolBreaker+1:ExtremeBreaker-1}});
end

if (exist('ExtremeBreaker','var')~=1)
    ExtremeInfo = zeros(1,16);
    if ~Silent
        fprintf('Extreme design information not found in file.\n')
    end
else
    ExtremeInfo = str2double({DesignInfo{ExtremeBreaker+1:end}});
end


% % % Annual Heating and Humidification Design Conditions
% MCDB = Mean Coincident Dry Bulb (Temperature)
% MCWB = Mean Coincident Wet Bulb (Temperature)

% Coldest month, (i.e., month with lowest average dry-bulb temperature;
% 1 = January, 12 = December).
OutStruct.Heat.Month = int32(HeatInfo(1));

% Heating Dry Bulb Temperature (99.6%, 99%)
OutStruct.Heat.TDB996 = HeatInfo(2);
OutStruct.Heat.TDB990 = HeatInfo(3);

% Humidification
% Dew-point temperature - 99.6%
OutStruct.Heat.TDP996 = HeatInfo(4);
% Humidity Ratio
OutStruct.Heat.W996 = HeatInfo(5);
% Mean Coincident Dry Bulb Temperature (MCDB)
OutStruct.Heat.TDP_MCDB996 = HeatInfo(6);

% Dew-point temperature - 99%
OutStruct.Heat.TDP990 = HeatInfo(7);
% Humidity Ratio
OutStruct.Heat.W990 = HeatInfo(8);
% Mean Coincident Dry Bulb Temperature (MCDB)
OutStruct.Heat.TDP_MCDB990 = HeatInfo(9);

% Coldest Month WS/MCDB
% Coincident Wind Speed - 4%
OutStruct.Heat.WS004 = HeatInfo(10);
OutStruct.Heat.WS_MCDB004 = HeatInfo(11);
% Coincident Wind Speed - 1%
OutStruct.Heat.WS001 = HeatInfo(12);
OutStruct.Heat.WS_MCDB001 = HeatInfo(13);

% Mean Coincident Wind Speed AND Prevailing Coincident Wind Direction
% corresponding to 99.6% Heating Dry Bulb Temperature 
OutStruct.Heat.TDB_MCWS996 = HeatInfo(14);
OutStruct.Heat.TDB_PCWD996 = HeatInfo(15);


% % % EXTREME Annual Design Conditions

% Annual Extreme Daily Mean Min and Max
OutStruct.Extreme.AnExTDBMeanMin = ExtremeInfo(5);
OutStruct.Extreme.AnExTDBMeanMax = ExtremeInfo(6);

% Annual Extreme Daily Standard Deviation Min and Max
OutStruct.Extreme.AnExTDBStdMin = ExtremeInfo(7);
OutStruct.Extreme.AnExTDBStdMax = ExtremeInfo(8);

% n-year Extreme Daily Mean Min and Max
% 5 years
OutStruct.Extreme.AnExTDB5Min = ExtremeInfo(9);
OutStruct.Extreme.AnExTDB5Max = ExtremeInfo(10);
% 10 years
OutStruct.Extreme.AnExTDB10Min = ExtremeInfo(11);
OutStruct.Extreme.AnExTDB10Max = ExtremeInfo(12);
% 20 years
OutStruct.Extreme.AnExTDB20Min = ExtremeInfo(13);
OutStruct.Extreme.AnExTDB20Max = ExtremeInfo(14);
% 50 years
OutStruct.Extreme.AnExTDB50Min = ExtremeInfo(15);
OutStruct.Extreme.AnExTDB50Max = ExtremeInfo(16);


% % % Annual Cooling, Dehumidification, and Enthalpy Design Conditions

% Hottest month (i.e., month with highest average dry-bulb temperature;
% 1 = January, 12 = December).
OutStruct.Cool.Month = int32(CoolInfo(1));

% Hottest month DB range
OutStruct.Cool.HotMonthRng = CoolInfo(2);

% Cooling Dry Bulb and Mean Coincident Wet Bulb Temperature (0.4%)
OutStruct.Cool.TDB004 = CoolInfo(3);
OutStruct.Cool.TDB_MCWB004 = CoolInfo(4);

% Cooling Dry Bulb and Mean Coincident Wet Bulb Temperature (1.0%)
OutStruct.Cool.TDB010 = CoolInfo(5);
OutStruct.Cool.TDB_MCWB010 = CoolInfo(6);

% Cooling Dry Bulb and Mean Coincident Wet Bulb Temperature (2.0%)
OutStruct.Cool.TDB020 = CoolInfo(7);
OutStruct.Cool.TDB_MCWB020 = CoolInfo(8);

% Cooling Evaporation Wet Bulb and Mean Coincident
% Wet Bulb Temperature (0.4%)
OutStruct.Cool.EWB004 = CoolInfo(9);
OutStruct.Cool.EWB_MCDB004 = CoolInfo(10);

% Cooling Evaporation Wet Bulb and Mean Coincident
% Wet Bulb Temperature (1.0%)
OutStruct.Cool.EWB010 = CoolInfo(11);
OutStruct.Cool.EWB_MCDB010 = CoolInfo(12);

% Cooling Evaporation Wet Bulb and Mean Coincident
% Wet Bulb Temperature (2.0%)
OutStruct.Cool.EWB020 = CoolInfo(13);
OutStruct.Cool.EWB_MCDB020 = CoolInfo(14);

% Mean Coincident Wind Speed (to 0.4% WB)
OutStruct.Cool.MCWS004 = CoolInfo(15);
% Prevailing Coincident Wind Direction (to 0.4% WB) - 0°N , 90°E
OutStruct.Cool.PCWD004 = CoolInfo(16);


% Cooling Dew Point Temperature (0.4%, 1.0%, 2.0%)
OutStruct.Cool.TDP004 = CoolInfo(17);
OutStruct.Cool.TDP010 = CoolInfo(20);
OutStruct.Cool.TDP020 = CoolInfo(23);

% (Cooling Dew Point) Mean Coincident Dry Bulb Temperature (each value
% below corresponds to one of the cooling design temperature values above)
OutStruct.Cool.TDP_MCDB004 = CoolInfo(18);
OutStruct.Cool.TDP_MCDB010 = CoolInfo(21);
OutStruct.Cool.TDP_MCDB020 = CoolInfo(24);

% Coincident Humidity Ratio (0.4%, 1.0%, 2.0%),
OutStruct.Cool.W004 = CoolInfo(19);
OutStruct.Cool.W010 = CoolInfo(22);
OutStruct.Cool.W020 = CoolInfo(25);

% Enthalpy (0.4%, 1.0%, 2.0%)
OutStruct.Cool.H004 = CoolInfo(26);
OutStruct.Cool.H010 = CoolInfo(28);
OutStruct.Cool.H020 = CoolInfo(30);

% (Enthalpy) Coincident Dry Bulb Temp (0.4%, 1.0%, 2.0%)
OutStruct.Cool.H_MCDB004 = CoolInfo(27);
OutStruct.Cool.H_MCDB010 = CoolInfo(29);
OutStruct.Cool.H_MCDB020 = CoolInfo(31);

% Number of hours between 8 AM and 4 PM (inclusive) with dry-bulb
% temperature between 12.8 and 20.6°C.
OutStruct.Cool.HoursASHRAERng = CoolInfo(32);

end

% % There seems to be a typo on page 10 of the Eplus weather data reference
% % document. The field number N23 is repeated - which means that the field
% % numbers of the last four values (in the DESIGN CONDITIONS line) are
% % wrong. I am using the proper numbers here - giving me 27 fields, not
% % 26.

% % Heating Degree Days Base Temperature,
% OutStruct.HDDBT = HeatInfo(end-1);
% % Heating Degree Days,
% OutStruct.HDD = HeatInfo(end);
%
% % Cooling Degree Days Base Temperature,
% OutStruct.CDDBT = CoolInfo(end-1);
% % Cooling Degree Days
% OutStruct.CDD = CoolInfo(end);

% Location
% Design Conditions
% Typical/Extreme Periods
% Ground Temperatures
% Leap Year Indicators
% Comments
