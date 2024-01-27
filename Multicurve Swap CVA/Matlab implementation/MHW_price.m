function [CB_price, ZCB_prices] = MHW_price(x,a,sigma,coupons,forw_disc,deltas_365,expiry,freq)
% computes Coupon Bond price and ZCB prices using Multicurve HW model
%
% INPUT
% x:           value of the ORNSTEIN-UHLENBECK process
% a:           a for HULL-WHITE model
% sigma:       sigma for HULL-WHITE model
% coupons:     coupons of the CB
% forw_disc:   forward discounts
% deltas_365:  time intervals in act/365
% expiry:      expiry in months
% freq:        payment frequency
%
% OUTPUT
% CB_price:    price of the Coupon Bond
% ZCB_prices:  prices of the ZC Bonds

% vector of times
t = deltas_365(2:end);

% Z in the dynamics of B
Z = sqrt(sigma^2*(1-exp(-2*a*deltas_365(1)))/(2*a));

% v in the dynamics of B 
v = Z*(1-exp(-a*(t-deltas_365(1))))/a;

% z in the dynamics of B with gamma=0
z = v;

% ZCB price function 
ZCB_prices = forw_disc(round(freq*t)-expiry*freq/12).*exp(-z.*x - 0.5*z.^2);

% CB price
CB_price = sum(coupons.*ZCB_prices);

end 