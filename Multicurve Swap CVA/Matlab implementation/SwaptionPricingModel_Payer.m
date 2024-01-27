function SP = SwaptionPricingModel_Payer(gamma,r,B,beta,a,sigma,time_interval_Z,deltas,time_intervals_v)
% computes price of the single Payer swaption by using MHW model
% 
% INPUT
% gamma:             fixed value of gamma
% r:                 strike, namely the fixed rate
% B:                 discount factors in each payment date
% beta:              vector of beta between two consecutive payment dates
% a:                 calibrated parameter a of the model
% sigma:             calibrated parameter sigma of the model
% time_interval_Z:   time to swaption expiry
% deltas:            deltas between two consecutive payment dates, but starting from the expiry
% time_intervals_v:  time intervals between expiry and payment dates (with also expiry)
%
% OUTPUT
% SP:                price of the Payer swaption

% Coupon vector
R = r*ones(size(deltas));
c = R.*deltas;
c(end) = c(end) + 1;

% Z 
Z = sigma*sqrt((1-exp(-2*a*time_interval_Z))/(2*a))*(a>0) + sigma*sqrt(time_interval_Z)*(a==0);

% v 
v = Z*(1-exp(-a*time_intervals_v))/a*(a>0) + Z*time_intervals_v*(a==0);

% z 
z = (1-gamma)*v; 

% nu
nu = v(1:end-1) - gamma*v(2:end); 

% forward discount factors B(t0,expiry,Ti)
B_fwd = B/B(1);  

% function f(x) that has to be put equal to zero
f = @(x) sum(c.*B_fwd(2:end).*exp(-z(2:end)*x-z(2:end).^2/2)) + ...
         sum(B_fwd(2:end-1).*exp(-z(2:end-1)*x-z(2:end-1).^2/2)) - ...
         sum(beta.*B_fwd(1:end-1).*exp(-nu*x-nu.^2/2));

% unique value of x such that f(x) = 0
x_star = fzero(f,-1);

% swaption Receiver price 
SR = B(1)*(sum(c.*B_fwd(2:end).*normcdf(x_star+z(2:end))) + ...
           sum(B_fwd(2:end-1).*normcdf(x_star+z(2:end-1))) - ...
           sum(beta.*B_fwd(1:end-1).*normcdf(x_star+nu)));

% Put-Call parity
BPV = sum(deltas.*B_fwd(2:end));
swap_rate = (1 - B_fwd(end) + sum(B_fwd(1:end-1).*(beta-1)))/BPV;

% swaption Payer price
SP = SR + B(1)*BPV*(swap_rate - r);

end