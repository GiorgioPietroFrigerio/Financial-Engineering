function [P_dates, P_curve] = pseudo_discount_curve_bootstrap(data,settle,daycounts,discount_dates,discount_curve)
% computes pseudo-discount factors curve
%
% INPUT
% data:                  struct containing dates and quoted prices of the needed financial instruments 
% settle:                settlement date
% daycounts:             basis_depos, basis_fra, basis_futures, basis_swap, basis_swap_floating
% discount_dates:        dates of the financial instruments used for the curve construction
% discount_curve:        discount factors curve
%
% OUTPUT
% P_dates:               dates of the financial instruments used for the pseudo-curve construction
% P_curve:               pseudo-discount factors curve

% allocation of dates
data.m3.depos.dates = change_weekends(addtodate(settle,3,'month'));
aux = @(x) change_weekends(addtodate(settle,x,'month'));
months = (1:6)';
data.m3.fra.dates = arrayfun(aux,months);

% 3m depos + fra
P_03 = 1/(1+yearfrac(settle,data.m3.depos.dates,daycounts(1))*data.m3.depos.mktRates);
P_fwd_36 = 1/(1+yearfrac(data.m3.depos.dates,data.m3.fra.dates(end),daycounts(2))*data.m3.fra.mktRates(end));
P_06 = P_03*P_fwd_36;

% interpolation for 4 and 5 months
P_04 = interp_ZR(settle, [data.m3.depos.dates, data.m3.fra.dates(end)], [P_03,P_06], data.m3.fra.dates(4));
P_05 = interp_ZR(settle, [data.m3.depos.dates, data.m3.fra.dates(end)], [P_03,P_06], data.m3.fra.dates(5));

% derivation of pseudo-discounts at 1 and 2 months
P_fwd_14 = 1/(1+yearfrac(data.m3.fra.dates(1),data.m3.fra.dates(4),daycounts(2))*data.m3.fra.mktRates(1));
P_fwd_25 = 1/(1+yearfrac(data.m3.fra.dates(2),data.m3.fra.dates(5),daycounts(2))*data.m3.fra.mktRates(2));
P_01 = P_04/P_fwd_14;
P_02 = P_05/P_fwd_25;

% forward pseudo-discount curve
pseudo_dates = [settle;data.m3.fra.dates(1:3)];
pseudo_curve = [P_03;P_fwd_14;P_fwd_25;P_fwd_36];

% pseudo-discount curve
P_dates = data.m3.fra.dates;
P_curve = [P_01;P_02;P_03;P_04;P_05;P_06];

% forward rates curve
fwd_curve = [data.m3.depos.mktRates;data.m3.fra.mktRates]; % fwd_dates = pseudo_dates

%% futures
nFutures = 7;                                            % first 7 futures
futureRates = 1-data.m3.futures.futPrices(1:nFutures);   % future rates from future prices
fwdRates = futureRates;                                  % no convexity term

% forward pseudo-discounts from futures
P_fwd = 1./(1+yearfrac(data.m3.futures.settle(1:nFutures),data.m3.futures.expiries(1:nFutures),daycounts(3)).*fwdRates);

% update curves (fwd pseudo-discount curve and dates, fwd rates)
pseudo_dates = [pseudo_dates; data.m3.futures.settle(1:nFutures)];
pseudo_curve = [pseudo_curve; P_fwd(1:7)];
fwd_curve = [fwd_curve; fwdRates];

% update pseudo-discount curve and dates
for i = 1:nFutures
    P_0s = interp_ZR(settle, P_dates, P_curve, data.m3.futures.settle(i));
    P_dates = [P_dates; data.m3.futures.expiries(i)];
    P_curve = [P_curve; P_0s.*P_fwd(i)];
end

% sort the curves so that they respect the order of dates
[P_dates, index] = sort(P_dates);
P_curve = P_curve(index);
[pseudo_dates, index] = sort(pseudo_dates);
pseudo_curve = pseudo_curve(index);

%% last forward rate of the second year (Point 5 of the bootstrap)

% construction of the vector w
index_2y = find(data.OIS.dates>=change_weekends(addtodate(settle,2,'year')));
data.m3.swaps.dates = data.OIS.dates(index_2y);
coupon_reset_dates_2y = change_weekends(cfdates(addtodate(data.m3.swaps.dates(1),-2,'year'),data.m3.swaps.dates(1),4));
delta_vec_2ycoupon = yearfrac([settle, coupon_reset_dates_2y(1:end-1)], coupon_reset_dates_2y, daycounts(5));
P_curve_interp = interp_ZR(settle, discount_dates, discount_curve, coupon_reset_dates_2y');
w_vec = delta_vec_2ycoupon'.*P_curve_interp;

% first 7 forward rates
P_fwd_curve_interp = interp_fwd_ZR(pseudo_dates, pseudo_curve, [settle; coupon_reset_dates_2y(1:end-2)']);
fwd_1to7 = (1./P_fwd_curve_interp - 1)./delta_vec_2ycoupon(1:end-1)';

% computation of I_2
swap_rate_2y = data.m3.swaps.mktRates(1);
yearly_reset_dates = change_weekends(cfdates(settle,data.m3.swaps.dates(1),1));
delta_vec_yearly = yearfrac([settle, yearly_reset_dates(1:end-1)], yearly_reset_dates, daycounts(4));
P_curve_interp_1y = interp_ZR(settle, discount_dates, discount_curve, yearly_reset_dates');
I_2 = swap_rate_2y*sum(delta_vec_yearly'.*P_curve_interp_1y);

% fwd rate at 2y and corresponding fwd pseudo-discount
fwd_8 = 1/w_vec(end)*(I_2-sum(w_vec(1:end-1)'*fwd_1to7));
P_fwd_8 = 1./(1+yearfrac(coupon_reset_dates_2y(end-1),coupon_reset_dates_2y(end), daycounts(5)).*fwd_8);

% update curves
pseudo_dates = [pseudo_dates; coupon_reset_dates_2y(end-1)];
pseudo_curve = [pseudo_curve; P_fwd_8];
P_0s = interp_ZR(settle, P_dates, P_curve, coupon_reset_dates_2y(end-1));
P_dates = [P_dates; coupon_reset_dates_2y(end)];
P_curve = [P_curve; P_0s*P_fwd_8];
fwd_curve = [fwd_curve; fwd_8];

%% swaps over 3y
aux = @(x) change_weekends(addtodate(settle, x, 'year'));
swap_dates_interp = arrayfun(aux,3:50)';  % swaps over 3 years
swap_rates_interp = interp1(data.m3.swaps.dates, data.m3.swaps.mktRates, swap_dates_interp, 'spline');
n_swaps = length(swap_rates_interp);      % last year swap

% floating legs
freq = 4; % coupon frequency
floating_dates = change_weekends(cfdates(settle,swap_dates_interp(end),freq));
discounts_floating_dates = interp_ZR(settle, discount_dates, discount_curve, floating_dates');
delta_floating_dates = yearfrac([settle, floating_dates(1:end-1)], floating_dates, daycounts(5));
w_vec = delta_floating_dates'.*discounts_floating_dates;
floating_dates = [settle; floating_dates'];

% fixed legs
fixed_yearly_dates = change_weekends(cfdates(settle,swap_dates_interp(end),1));
discounts_fixed_dates = interp_ZR(settle, discount_dates, discount_curve, fixed_yearly_dates');
delta_fixed_dates = yearfrac([settle, fixed_yearly_dates(1:end-1)], fixed_yearly_dates, daycounts(4));

% allocation of dates and forwards
fwd_dates_new = pseudo_dates;
fwd_curve_new = fwd_curve;

for i = 1:n_swaps
    j = i+2;  % maturity of the swap
    swap_rate = swap_rates_interp(i);
    I_i = swap_rate * sum(delta_fixed_dates(1:j)'.*discounts_fixed_dates(1:j));  % I_i
    
    % derivation of necessary forwards (using interpolation)
    F = @(x) interp1([fwd_dates_new; floating_dates(freq*j)], [fwd_curve_new; x], floating_dates(1:freq*j));
    fun = @(x) sum(w_vec(1:freq*j).*F(x))-I_i;
    last_fwd = fzero(fun, 0.02); % last forward (from 9m to 12 m)
    F_vec = F(last_fwd);        
    fwd_dates_new = floating_dates(1:freq*j);
    fwd_curve_new = F_vec;
    
    % add the yearly pseudo-discount and forward pseudo-discount
    pseudo_dates = [pseudo_dates; floating_dates(freq*j)];
    pseudo_curve = [pseudo_curve; prod(1./(1+delta_floating_dates(freq*j-3:freq*j)'.*F_vec(end-3:end)))];
    P_0s = interp_ZR(settle, P_dates, P_curve, floating_dates(freq*j-3)');
    P_dates = [P_dates; floating_dates(freq*j+1)];
    P_curve = [P_curve;  P_0s*pseudo_curve(end)];
end

end