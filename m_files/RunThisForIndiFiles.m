% outfolder = fullfile('..', 'syn_data','NYC');
outfolder = 'F:\WeatherData\SyntheticData\was';


% synfile = fullfile('..', 'syn_data','NYC', 'NYC_IWEC_Syn.mat');

% srcEPWfile = fullfile('..', 'm_data', 'NYC', 'NYC_IWEC.epw');

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

% SMY_Create('NYC', synfile, outfolder, ...
%     'sublabel', 50, 'scenario', 'syn', 'srcEPWfile', srcEPWfile)

srcEPWfile = 'E:\WeatherData\HistoricalData\was\NIST_TMY3_2015-01_to_2015-12_TMY3.epw';
synfile = 'E:\WeatherData\SyntheticData\was\NIST_TMY3_2015-01_to_2015-12_TMY3_rcp45.mat';
SMY_Create('NIST_2015', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', srcEPWfile)

srcEPWfile = 'E:\WeatherData\HistoricalData\was\NIST_TMY3_2015-01_to_2015-12_TMY3.epw';
synfile = 'E:\WeatherData\SyntheticData\was\NIST_TMY3_2015-01_to_2015-12_TMY3_rcp85.mat';
SMY_Create('NIST_2015', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', srcEPWfile)