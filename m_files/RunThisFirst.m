% An example run of the CreateSyntheticFiles function.

% Path to source file.
% pathEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_IWEC.epw');
pathEPWfile = 'E:\WeatherDataNew\HistoricalData\NYC\NYC_CPR_TMY-36.epw';

% Path to folder where you want all the new stuff saved.
% path_save_fldr = fullfile('..', 'syn_data', 'GEN');
path_save_fldr = 'F:\WeatherData\SyntheticData\NYC\NYC_36';

% Number of samples. Change this to get more/less samples.
nboot = 50;

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

ccpath = 'F:\WeatherData\CCdata\NYC';

% Random seed. If you use the same random seed for any run with the same
% source data, the random samples produced will be exactly the same.
randseed = 1;

% If this is not the first time you are running these files, the hourly
% Fourier and SARIMA models might be stored. You can specify the path to
% these here.
hrmdlfile = 'F:\WeatherData\SyntheticData\NYC\NYC_36\FourierFits_NYC_CPR_TMY-36.mat';
fourierfile = 'F:\WeatherData\SyntheticData\NYC\NYC_36\HourMdls_NYC_CPR_TMY-36.mat';

CreateSyntheticFiles(pathEPWfile, path_save_fldr, nboot, recdata, ...
    ccdata, 'randseed', randseed, 'hrmdlfile', hrmdlfile, ...
    'fourierfile', fourierfile, 'ccpath', ccpath)