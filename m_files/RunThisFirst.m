% An example run of the CreateSyntheticFiles function.

% Path to source file.
pathEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_Meteonorm.epw');

% Path to folder where you want all the new stuff saved.
namesavefldr = fullfile('..', 'syn_data', 'GEN');

% Number of samples.
nboot = 50;

% Is recorded data present? Default is no. If this is present, then the
% script will calculate a bunch of summary statistics. However, you must
% provide a path to the files containing recorded data.
recdata = false;

% Same for climate change 'data' - i.e., outputs from atmospheric models
% of climate change.
ccdata = false;

% Random seed. If you use the same random seed for any run with the same
% source data, the random samples produced will be exactly the same.
randseed = 42;

% If this is not the first time you are running these files, the hourly
% Fourier and SARIMA models might be stored. You can specify the path to
% these here.
% hrmdlfile = ''
% fourierfile = ''

CreateSyntheticFiles(pathEPWfile, namesavefldr, nboot, recdata, ...
    ccdata, 'randseed', randseed)