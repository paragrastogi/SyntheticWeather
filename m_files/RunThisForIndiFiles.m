% outfolder = fullfile('..', 'syn_data','NYC');
outfolder = 'F:\WeatherData\SyntheticData\NYC\NYC_36';


% synfile = fullfile('..', 'syn_data','NYC', 'NYC_IWEC_Syn.mat');

% srcEPWfile = fullfile('..', 'm_data', 'NYC', 'NYC_IWEC.epw');
srcEPWfile = 'E:\WeatherDataNew\HistoricalData\NYC\NYC_CPR_TMY-36.epw';

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

% SMY_Create('NYC', synfile, outfolder, ...
%     'sublabel', 50, 'scenario', 'syn', 'srcEPWfile', srcEPWfile)

% synfile = 'F:\WeatherData\SyntheticData\NYC\NYC_36\NYC_CPR_TMY-36_rcp85.mat';
% SMY_Create('NYC', synfile, outfolder, ...
%     'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', srcEPWfile)

synfile = 'F:\WeatherData\SyntheticData\NYC\NYC_36\NYC_CPR_TMY-36_rcp45.mat';
SMY_Create('NYC', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', srcEPWfile)