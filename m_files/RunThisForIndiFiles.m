% outfolder = fullfile('..', 'syn_data','NYC');
if strcmp(getenv('computername'), 'SATAH')
    prefixpath = 'f:\WeatherData';
elseif isunix && strcmp(getenv('USER'), 'rasto')
    prefixpath = '/media/rasto/LargeKaali/WeatherData';
end

outfolder = fullfile(prefixpath, 'SyntheticData', 'was');

epw_filename = 'USA_VA_Dulles-Washington.Dulles.Intl.AP.724030_TMYx-15';

% synfile = fullfile('..', 'syn_data','NYC', 'NYC_IWEC_Syn.mat');

% pathEPWfile = fullfile('..', 'm_data', 'NYC', 'NYC_IWEC.epw');

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

pathEPWfile = fullfile(prefixpath, 'HistoricalData', 'was', ...
    [epw_filename, '.epw']);

% SMY_Create('NYC', synfile, outfolder, ...
%     'sublabel', 50, 'scenario', 'syn', 'srcEPWfile', pathEPWfile)

synfile = fullfile(prefixpath, 'SyntheticData', 'was', ...
    [epw_filename, '_rcp45.mat']);
SMY_Create(epw_filename, synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', pathEPWfile)

% pathEPWfile = fullfile(prefixpath, 'HistoricalData', 'was', ...
%     'NIST_TMY3_2015-01_to_2015-12_TMY3.epw');
synfile = fullfile(prefixpath, 'SyntheticData', 'was', ...
    [epw_filename, '_rcp85.mat']);
SMY_Create(epw_filename, synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', pathEPWfile)