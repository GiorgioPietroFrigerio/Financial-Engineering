function [discount_dates, discount_curve,data] = discount_curve_bootstrap(data,settle,Basis_swap)
% computes discount factors curve
%
% INPUT
% data:            struct containing dates and quoted prices of the needed financial instruments 
% settle:          settlement date
% Basis_swap:      30/360 European
%
% OUTPUT
% discount_dates:  dates of the financial instruments used for the curve construction
% discount_curve:  discount factors curve
% data:            updated struct containing dates and quoted prices of the needed financial instruments 

%% <= 1y swaps

% allocation of dates
aux = @(x) change_weekends(addtodate(settle,x,'day'));
days = [1, 7:7:21]';
data.OIS.dates = arrayfun(aux,days);
aux = @(x) change_weekends(addtodate(settle,x,'month'));
months = [1:12, 15 18 21]';
data.OIS.dates = [data.OIS.dates; arrayfun(aux,months)];
aux = @(x) change_weekends(addtodate(settle,x,'year'));
years = [2:12, 15:5:30, 40, 50]';
data.OIS.dates = [data.OIS.dates; arrayfun(aux,years)];
 
% find the index of the expiries prior to 1 year
idx_1y = find(data.OIS.dates<=change_weekends(addtodate(settle,1,'year')));

% creation of discount curve and dates up to 1 year
discount_curve = 1./(1 + yearfrac(settle, data.OIS.dates(idx_1y),Basis_swap).*data.OIS.mktRates(idx_1y));
discount_dates = [data.OIS.dates(idx_1y)];

%% > 1y swaps

% find the index of the expiries after 1 year
idx_over1y = find(data.OIS.dates>change_weekends(addtodate(settle,1,'year')));
aux = @(x) change_weekends(addtodate(settle, x, 'year'));
swap_dates_interp = [data.OIS.dates(idx_over1y(1:3));arrayfun(aux,2:50)'];
swap_rates_interp = interp1(data.OIS.dates(idx_over1y), data.OIS.mktRates(idx_over1y), swap_dates_interp, 'spline');

% discount curve from swaps from 15m up to 21m
BPV_01y = discount_curve(end)*yearfrac(settle,discount_dates(end),Basis_swap);
date_1y = discount_dates(end);
discount_curve = [discount_curve; (1-swap_rates_interp(1:3)*BPV_01y)./(1+swap_rates_interp(1:3).*yearfrac(date_1y, data.OIS.dates(idx_over1y(1:3)), Basis_swap))];

% discount curve from swaps at 2y
partial_sum = BPV_01y;
discount_curve = [discount_curve; (1-swap_rates_interp(4)*partial_sum)./(1+yearfrac(date_1y, swap_dates_interp(4), Basis_swap)*swap_rates_interp(4))];
partial_sum = partial_sum + discount_curve(end)*yearfrac(date_1y, swap_dates_interp(4), Basis_swap);

% discount curve from swaps over 2y
for i = 5:length(swap_dates_interp)
    discount_curve = [discount_curve; (1-swap_rates_interp(i)*partial_sum)./(1+yearfrac(swap_dates_interp(i-1), swap_dates_interp(i), Basis_swap)*swap_rates_interp(i))];
    partial_sum = partial_sum + discount_curve(end)*yearfrac(swap_dates_interp(i-1), swap_dates_interp(i), Basis_swap);
end

% all discount dates
discount_dates = [discount_dates; swap_dates_interp];

end