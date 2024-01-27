function data = readExcelData_SwapAnnex(filename,formatData)
% reads swap annex for the amortized contract
%
% INPUT
% filename:    excel file name where data are stored
% formatData:  data format in Excel
%
% OUTPUT
% data:        struct updated with start and end dates of each payment period, with the corresponding notional

% save payment dates
[~, pay_dates] = xlsread(filename, 1, 'A2:A61');
data.PayDates = datenum(pay_dates, formatData);

% save start date of each payment period
[~, start_dates] = xlsread(filename, 1, 'B2:B61');
data.StartDates = datenum(start_dates, formatData);

% save end date of each payment period
[~, end_dates] = xlsread(filename, 1, 'C2:C61');
data.EndDates = datenum(end_dates, formatData);

% save notional of each payment period
data.Notionals = xlsread(filename, 1, 'E2:E61');

end