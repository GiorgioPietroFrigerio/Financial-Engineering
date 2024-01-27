function [B, Bp, Bp_fwd] = swap_discounts_modified(settle,discount_dates,discount_curve,P_dates,P_curve,StartDates,EndDates,Euribor,idx)
% computes discount factors and pseudo-discount factors in each swap payment date;
% must be permormed separately due to the presence of a Euribor settled on
% 28 December, paid on 28 March and being today the 31 January
%
% INPUT
% settle:          settlement date      
% discount_dates:  dates of the financial instruments used for the curve construction
% discount_curve:  discount factors curve
% P_dates:         dates of the financial instruments used for the pseudo-curve construction
% P_curve:         pseudo-discount factors curve
% StartDates:      start dates of each swap payment period
% EndDates:        end dates of each swap payment period
% Euribor:         Euribor established on 28 December and paid on 28 March, taken from web market value
% idx:             vector of indexes of payment dates following the new settlement, namely 31 January
%
% OUTPUT
% B:               discount factors in each payment date
% Bp:              pseudo-discount factors in each payment date
% Bp_fwd:          forward pseudo-discount factors between two consecutive payment dates 

% daycount for swap
swap_daycount = 2;  % act/360

% discount factors for swap
B = interp_ZR(settle,discount_dates,discount_curve,EndDates(idx));

% pseudo-discounts for swap
Bp = interp_ZR(settle,P_dates,P_curve,EndDates(idx));

% first time interval
delta_1 = yearfrac(StartDates(idx(1)),EndDates(idx(1)),swap_daycount);

% pseudo-discount in first time interval
Bp_1 = 1/(1+delta_1*Euribor);

% forward pseudo-discounts for swap
Bp_fwd = [Bp_1; Bp(2:end)./Bp(1:end-1)];

end