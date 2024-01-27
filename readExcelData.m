function data = readExcelData(filename,formatData,N_sheet)
% reads data from excel
%
% INPUT
% filename:    excel file name where data are stored
% formatData:  data format in Excel
% N_sheet:     1 -> 24 June 2022, 2 -> 31 January 2023
% 
% OUTPUT
% data:        struct updated with referenceDate, settleDate, OIS data, 3M, Futures, Swaption vols

% reference date
[~, refDate] = xlsread(filename, N_sheet, 'B1');

% date conversion
data.refDate = datenum(refDate, formatData);

% settlement date
[~, settlement] = xlsread(filename, N_sheet, 'B2');

% date conversion
data.settlement = datenum(settlement, formatData);

%% OIS data

% market rates
data.OIS.mktRates = xlsread(filename, N_sheet, 'C6:C41')/100;

%% 3m data

% depos
data.m3.depos.mktRates = xlsread(filename, N_sheet, 'J6')/100;

% fra
data.m3.fra.mktRates = xlsread(filename, N_sheet, 'J7:J9')/100;

% futures
[~, settle_fut] = xlsread(filename, N_sheet, 'Q6:Q13');
[~, expiries_fut] = xlsread(filename, N_sheet, 'R6:R13');
data.m3.futures.settle = datenum(settle_fut, formatData);
data.m3.futures.expiries = datenum(expiries_fut, formatData);
data.m3.futures.futPrices = xlsread(filename, N_sheet, 'J10:J17')/100;

% swaps
data.m3.swaps.mktRates = xlsread(filename, N_sheet, 'J18:J34')/100;

end 
