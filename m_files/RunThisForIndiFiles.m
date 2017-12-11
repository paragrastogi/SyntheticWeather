% outfolder = fullfile('..', 'syn_data','NYC');
outfolder = '/media/rasto/LargeKaali/WeatherData/SyntheticData/was';


% synfile = fullfile('..', 'syn_data','NYC', 'NYC_IWEC_Syn.mat');

% srcEPWfile = fullfile('..', 'm_data', 'NYC', 'NYC_IWEC.epw');

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

% SMY_Create('NYC', synfile, outfolder, ...
%     'sublabel', 50, 'scenario', 'syn', 'srcEPWfile', srcEPWfile)

srcEPWfile = '/media/rasto/LargeKaali/WeatherData/HistoricalData/was/NIST_TMY3_2016-01_to_2016-12_TMY3.epw';
synfile = '/media/rasto/LargeKaali/WeatherData/SyntheticData/was/NIST_TMY3_2016-01_to_2016-12_TMY3_rcp45.mat';
SMY_Create('NIST_2016', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', srcEPWfile)

srcEPWfile = '/media/rasto/LargeKaali/WeatherData/HistoricalData/was/NIST_TMY3_2016-01_to_2016-12_TMY3.epw';
synfile = '/media/rasto/LargeKaali/WeatherData/SyntheticData/was/NIST_TMY3_2016-01_to_2016-12_TMY3_rcp85.mat';
SMY_Create('NIST_2016', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', srcEPWfile)