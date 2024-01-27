function price_put = swaption_price_Jam(settle,dates,discounts,expiry_date,expiry_in_months,tenor_in_months,maturity,sigma,a,k,Notionals,Initial_Notional,beta,freq,flag)
% computes the price of a swaption payer using the Jamshidian formula
%
% INPUT
% settle:             settlement date          
% dates:              dates from bootstrap
% discounts:          vector of discount factors from bootstrap
% expiry_date:        expiry date of the swaption
% expiry_in_months:   expiry in months of the swaption
% tenor_in_months:    tenor of the swaption in months
% sigma:              calibrated sigma for HULL-WHITE model 
% a:                  calibrated a for HULL-WHITE model 
% k:                  strike, namely the fixed rate
% Notionals:          vector of amortized notionals
% Initial_Notional:   initial notional
% beta:               beta between two consecutive payment dates
% freq:               annual payments frequency
% flag:               "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% price_put:          price of the swaption payer

% date corresponding to the expiry of the swaption
Start_date = expiry_date;

% coupon payment dates
PaymentDates = cfdates(Start_date,maturity-1,freq); 

% time intervals (act/365)
deltas_365 = yearfrac(settle,[Start_date PaymentDates],3);

% time intervals (act/360)
deltas_360 = yearfrac([Start_date PaymentDates(1:end-1)], PaymentDates(1:end),2);

% discount factors at coupon payment dates
disc_fact = interp_ZR(settle,dates,discounts,[Start_date PaymentDates]);

% forward discount factors
forw_disc_fact = disc_fact/disc_fact(1);

if flag == "SA"

    % Computing the coupons and adding the face value for the last one
    coupons = k*deltas_360(1:end-1)'.*Notionals(1:end-1)/Initial_Notional + 1 - beta(2:end);
    coupons = [coupons; k*deltas_360(end)*Notionals(end)/Initial_Notional + 1];
    
    % option strike
    option_strike = beta(1);

elseif flag == "DA"

    % Computing the coupons 
    coupons = (1 + k*deltas_360(1:end-1)').*Notionals(1:end-1) - beta(2:end).*Notionals(2:end);
    coupons = [coupons; (1 + k*deltas_360(end))*Notionals(end)];

    % option strike
    option_strike = beta(1)*Notionals(1);

end

% price of the bond - strike (function to be put = 0)
P = @(x) MHW_price(x,a,sigma,coupons',forw_disc_fact(2:end),deltas_365,expiry_in_months,freq) - option_strike;  

% finding x_star
x_star = fzero(P,-1);

% finding the K_i
[~,K_i] = MHW_price(x_star,a,sigma,coupons',forw_disc_fact(2:end),deltas_365,expiry_in_months,freq);

% sigmas for black formula
I = 1:tenor_in_months*freq/12;                           % vector of indexes
Z = sqrt( sigma^2*(1-exp(-2*a*deltas_365(1)))/(2*a) );   % Z
v = Z*(1-exp(-a*(deltas_365(I)-deltas_365(1))))/a;       % v
z = v;                                                   % z (=v if gamma=0)

% prices of options on ZC bonds
puts = disc_fact(1)*K_i.*(1-normcdf(x_star)) - disc_fact(2:end).*(1-normcdf(x_star+z));

% swaption price (payer case)
if flag == "SA"
    price_put = Initial_Notional*puts*coupons;
elseif flag == "DA"
    price_put = puts*coupons;
end

end 