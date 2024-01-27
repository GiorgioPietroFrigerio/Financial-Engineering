function new_discounts = interp_ZR(settle,dates,discounts,new_dates)
% computes the zero-rate linear interpolation on discounts
%
% INPUT
% settle:         settlement date
% dates:          dates of the discount factors curve
% discounts:      discount factors curve
% new_dates:      dates of interest
%
% OUTPUT
% new_discounts:  discount factors in the desidered dates

% act/365 for linear ZR interpolation
basis_interpolation = 3; 

% find the corresponding zero rates
zr = -log(discounts)./yearfrac(settle, dates, basis_interpolation);

% interpolate on the zero rates
zr_interp = interp1(dates, zr, new_dates,'linear','extrap');

% compute new discount factors
new_discounts = exp(-zr_interp.*yearfrac(settle,new_dates,basis_interpolation));

end

