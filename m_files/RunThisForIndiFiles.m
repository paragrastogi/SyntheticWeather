outfolder = fullfile('..', 'syn_data','GEN');
synfolder = fullfile('..', 'syn_data','GEN');
srcEPWfile = fullfile('..', 'm_data', 'GEN', 'GEN_IWEC.epw');

% Delete the existing files that match the names of the new synthetic
% files. Use this power wisely.
destroyer = false;

SMY_Create('GEN', outfolder, synfolder, ...
    'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', srcEPWfile)

SMY_Create('GEN', outfolder, synfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', srcEPWfile)