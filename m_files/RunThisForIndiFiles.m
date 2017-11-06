outfolder = fullfile('..', 'syn_data','GEN');
synfile = fullfile('..', 'syn_data','GEN', 'GEN_IWEC_Syn.mat');
srcEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_IWEC.epw');

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

SMY_Create('GEN', synfile, outfolder, ...
    'sublabel', 50, 'scenario', 'syn', 'srcEPWfile', srcEPWfile)

% SMY_Create('GEN', outfolder, synfile, ...
%     'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', srcEPWfile)
% 
% SMY_Create('GEN', outfolder, synfile, ...
%     'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', srcEPWfile)