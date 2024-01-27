function vols = readExcelData_vol_calibration(filename,N_sheet,expiries,tenors)
% reads implied volatilies only for the selected expiries and tenors
%
% INPUT
% filename:  excel file name where data are stored
% N_sheet:   1 -> 24 June 2022, 2 -> 31 January 2023
% expiries:  vector of expiries of the swaptions of interest
% tenors:    vector of tenors of the swaptions of interest
%
% OUTPUT
% vols:      vector of implied volatilies of the swaptions of interest

% implied volatilities
volatilities = xlsread(filename, N_sheet, 'Q17:AE37')*1e-4;

% expiries
[~, exp] = xlsread(filename, N_sheet, 'P17:P37');

% tenors
[~, ten] = xlsread(filename, N_sheet, 'Q16:AE16');

% initialization of volatilities
vols = zeros(length(expiries),1);

% extraction of the volatilities of interest
for i = 1:length(expiries)
    
    for j = 1:length(exp)

        idx_e = find(isletter(exp{j})==1);
        if contains(exp{j},'Y')
            if str2num(exp{j}(1:idx_e-1))==expiries(i)
                idx_expiry = j;
            end
        end

    end

    for j = 1:length(ten)

        idx_e = find(isletter(ten{j})==1);
        if contains(ten{j},'Y')
            if str2num(ten{j}(1:idx_e-1))==tenors(i)
                idx_tenors = j;
            end
        end

    end

    vols(i) = volatilities(idx_expiry, idx_tenors);

end

end
