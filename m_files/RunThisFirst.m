% An example run of the CreateSyntheticFiles function.

% Path to source file.
% pathEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_IWEC.epw');
if strcmp(getenv('computername'), 'SATAH')
    prefixpath = 'd:\WeatherData';
elseif isunix && strcmp(getenv('USER'), 'rasto')
    prefixpath = '/media/rasto/LargeKaali/WeatherData';
end

epw_filename = 'IND_GJ_Ahmedabad-Patel.Intl.AP.426470_ISHRAE2014';

pathEPWfile = fullfile(prefixpath, 'HistoricalData', ...
    'amd', [epw_filename, '.epw']);

% Path to folder where you want all the new stuff saved.
% path_save_fldr = fullfile('..', 'syn_data', 'GEN');
path_save_fldr = fullfile(prefixpath, 'SyntheticData', 'amd');
% if strcmp(getenv('computername'), 'SATAH')
%     path_save_fldr = 'f:\WeatherData\SyntheticData\amd';
% elseif isunix && strcmp(getenv('USER'), 'rasto')
%     path_save_fldr = '/media/rasto/LargeKaali/WeatherData/SyntheticData/amd';
% end

% Number of samples. Change this to get more/less samples.
nboot = 25;

% Is recorded data present? Default is no. If this is present, then the
% script will calculate a bunch of summary statistics. However, you must
% provide a path to the files containing recorded data.
recdata = false;

% Same for climate change 'data' - i.e., outputs from atmospheric models
% of climate change.
% Setting this to true will generate two variants: rcp45 and rcp85,
% corresponding to the RCP 4.5 and RCP 8.5 simulations from the IPCC
% respectively.
ccdata = true;

ccpath = fullfile(prefixpath, 'CCdata', 'amd');
% if strcmp(getenv('computername'), 'SATAH')
%     path_save_fldr = 'f:\WeatherData\CCData\amd';
% elseif isunix && strcmp(getenv('USER'), 'rasto')
%     ccpath = '/media/rasto/LargeKaali/WeatherData/CCdata/amd';
% end

% Random seed. If you use the same random seed for any run with the same
% source data, the random samples produced will be exactly the same.
randseed = 42;

% If this is not the first time you are running these files, the hourly
% Fourier and SARIMA models might be stored. You can specify the path to
% these here.
% hrmdlfile = fullfile(prefixpath, 'SyntheticData', 'amd', ['FourierFits_', epw_filename, '.mat']);
% fourierfile = fullfile(prefixpath, 'SyntheticData', 'amd', ['HourMdls_', epw_filename, '.mat']);
% , ...
%     'hrmdlfile', hrmdlfile, 'fourierfile', fourierfile

CreateSyntheticFiles(pathEPWfile, path_save_fldr, nboot, recdata, ...
    ccdata, 'randseed', randseed, 'ccpath', ccpath)
