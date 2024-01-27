function dates = find_busdates_monthly(initial_date,freq,mat)
% finds the requested business dates (monthly)
%
% INPUT
% initial_date:  initial date
% freq:          frequency of dates in a year
% mat:           maturity (number of years)
%
% OUTPUT
% dates:         requested dates

% find the dates
fun = @(i) busdate(addtodate(initial_date-1,12/freq*i,'month'));
I = 1:freq*mat;
dates = arrayfun(fun,I);

end 