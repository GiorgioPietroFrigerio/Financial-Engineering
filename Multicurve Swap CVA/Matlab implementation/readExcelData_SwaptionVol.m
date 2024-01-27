function data = readExcelData_SwaptionVol(filename,formatData,N_sheet)
% reads swaption volatilities from market values
%
% INPUT
% filename:    excel file name where data are stored
% formatData:  data format in Excel
% N_sheet:     1 -> 24 June 2022, 2 -> 31 January 2023
%
% OUTPUT
% data:        updated struct containing swaption volatilities, with the relative expiry and tenor

% read implied volatilities
data.vol = xlsread(filename, N_sheet, 'Q17:AE37')*1e-4;

% read expiries
[~, expiries] = xlsread(filename, N_sheet, 'P17:P37');
data.expiry = ones(length(expiries),1);

% settlement date
[~, settlement] = xlsread(filename, N_sheet, 'B2');

% date conversion
settle = datenum(settlement, formatData);

% convert expiries from years and months to number
for i = 1:length(expiries)
    idx = find(isletter(expiries{i})==1);
    if contains(expiries{i}, 'M')
        data.expiry(i) = change_weekends(addtodate(settle,str2num(expiries{i}(1:idx-1)), 'month'));
    else
        data.expiry(i) = change_weekends(addtodate(settle,str2num(expiries{i}(1:idx-1)), 'year'));
    end
end

% read tenors
[~, tenors] = xlsread(filename, N_sheet, 'Q16:AE16');
data.tenor = ones(length(tenors),1);

% convert tenors from years to number
for i = 1:length(tenors)
    idx = find(isletter(tenors{i})==1);
    data.tenor(i) = str2num(tenors{i}(1:idx-1));
end

end