outfolder = 'C:\Users\prastogi\Documents\allmycode\SyntheticWeather\m-data';
synfolder = 'C:\Users\prastogi\Documents\WeatherData\SyntheticData';
pathEPWfile = 'C:\Users\prastogi\Documents\allmycode\SyntheticWeather\m-data\GEN\GEN_IWEC.epw';

SMY_Create({'GEN'}, outfolder, synfolder, ...
    'sublabel', 50, 'scenario', 'rcp85', 'srcEPWfile', pathEPWfile)

SMY_Create({'GEN'}, outfolder, synfolder, ...
    'sublabel', 50, 'scenario', 'rcp45', 'srcEPWfile', pathEPWfile)