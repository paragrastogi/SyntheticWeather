function FileParserKeywordList = FileParserKeywordLister()

% This is just a list of the TMY and AMY keywords used in all of the
% scripts

FileParserKeywordList.TMYKeyword = {'NREL', 'IWEC', 'ISHRAE', 'CWEC', ...
    'IGDG','TMY3','Meteonorm'};
FileParserKeywordList.AMYKeyword = {'NSRDB', 'NRELIndiaSolar', ...
    'NCDC','MeteoSuisse','WY2','NASASaudi'};

end