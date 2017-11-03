% An example run of the CreateSyntheticFiles function.

pathEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_Meteonorm');
namesavefldr = fullfile('..', 'syn_data', 'GEN', 'GEN_Meteonorm');

CreateSyntheticFiles(pathEPWfile, namesavefldr, nboot, recdata, ccdata, varargin)