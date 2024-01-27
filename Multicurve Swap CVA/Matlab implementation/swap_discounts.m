function [B, Bp, Bp_fwd] = swap_discounts(settle,discount_dates,discount_curve,P_dates,P_curve,Dates)   
% computes discount factors and pseudo-discount factors in each swap payment date
%
% INPUT
% settle:          settlement date      
% discount_dates:  dates of the financial instruments used for the curve construction
% discount_curve:  discount factors curve
% P_dates:         dates of the financial instruments used for the pseudo-curve construction
% P_curve:         pseudo-discount factors curve
% Dates:           payment dates of the swap contract
%
% OUTPUT
% B:               discount factors in each payment date
% Bp:              pseudo-discount factors in each payment date
% Bp_fwd:          forward pseudo-discount factors between two consecutive payment dates 

% discount factors for swap
B = interp_ZR(settle,discount_dates,discount_curve,Dates);

% pseudo-discounts for swap
Bp = interp_ZR(settle,P_dates,P_curve,Dates);

% forward pseudo-discounts for swap
Bp_fwd = [Bp(1); Bp(2:end)./Bp(1:end-1)];

end