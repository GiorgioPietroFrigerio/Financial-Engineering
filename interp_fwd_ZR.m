function new_discounts = interp_fwd_ZR(dates,discounts,new_dates)
% computes the zero-rate linear interpolation on forward discounts
%
% INPUT
% settle:         settlement date
% dates:          dates of the forward discount factors set
% discounts:      forward discount factors set
% new_dates:      dates of interest
%
% OUTPUT
% new_discounts:  forward discount factors in the desidered dates

% act/365 day count for linear ZR interpolation
basis_interpolation = 3; 

% forward expiries
I = (1:length(dates))';
forward_expiries_fun = @(i) change_weekends(addtodate(dates(i), 3, 'month'));
forward_expiries = arrayfun(forward_expiries_fun,I);

% new forward expiries
I = (1:length(new_dates))';
new_forward_expiries_fun = @(i) change_weekends(addtodate(new_dates(i), 3, 'month'));
new_forward_expiries = arrayfun(new_forward_expiries_fun,I);

% find the corresponding zero rates
zr = -log(discounts)./yearfrac(dates,forward_expiries,basis_interpolation);

% interpolate on the zero rates
zr_interp = interp1(dates,zr,new_dates);

% compute new forward discount factors
new_discounts = exp(-zr_interp.*yearfrac(new_dates,new_forward_expiries,basis_interpolation));

end