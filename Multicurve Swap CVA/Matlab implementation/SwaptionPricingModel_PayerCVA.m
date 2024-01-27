function SP = SwaptionPricingModel_PayerCVA(gamma,r,B,beta,a,sigma,settle,expiry,payDates,Notionals,Initial_Notional,flag)
% computes price of the single Payer swaption by using the calibrated model
% 
% INPUT
% gamma:             fixed value of gamma
% r:                 strike, namely the fixed rate
% B:                 discount factors in each payment date
% beta:              vector of beta between two consecutive payment dates
% a:                 calibrated parameter a of the model
% sigma:             calibrated parameter sigma of the model
% settle:            settlement date
% expiry:            expiry date Ti
% payDates:          payment dates of the underlying swap
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% PS:                price of the Payer swaption 

% deltas between two consecutive payment dates, but starting from the expiry
delta = yearfrac([expiry; payDates(1:end-1)],payDates,2);

% Coupon vector
R = r*ones(size(delta));
c = R.*delta;

% Z 
time_interval_Z = yearfrac(settle,expiry,3);
Z = sigma*sqrt((1-exp(-2*a*time_interval_Z))/(2*a))*(a>0) + sigma*sqrt(time_interval_Z)*(a==0);

% v 
time_intervals_v = yearfrac(expiry,[expiry; payDates],3);
v = Z*(1-exp(-a*time_intervals_v))/a*(a>0) + Z*time_intervals_v*(a==0);

% z
z = (1-gamma)*v; 

% nu
nu = v(1:end-1) - gamma*v(2:end); 

% forward discount factors B(0,expiry,Ti)
B_fwd = B/B(1); 
BPV = sum(delta.*B_fwd(2:end).*Notionals);

if flag == "SA"

    % take in account the notionals as ratios
    c = c.*Notionals/Initial_Notional;
    
    % updated function f(x) with different notionals
    f = @(x) sum(c.*B_fwd(2:end).*exp(-z(2:end)*x-z(2:end).^2/2)) + ...
             sum(B_fwd(2:end).*exp(-z(2:end)*x-z(2:end).^2/2)) - ...
             sum(beta.*B_fwd(1:end-1).*exp(-nu*x-nu.^2/2));

    % unique value of x such that f(x) = 0
    x_star = fzero(f,-1);

    % swaption Receiver price, updated due to different notionals in single amortizing
    SR = Initial_Notional*B(1)*(sum(c.*B_fwd(2:end).*normcdf(x_star+z(2:end))) + ...
                                sum(B_fwd(2:end).*normcdf(x_star+z(2:end))) - ...
                                sum(beta.*B_fwd(1:end-1).*normcdf(x_star+nu)));

    swap_rate = Initial_Notional*(1 - B_fwd(end) + sum(B_fwd(1:end-1).*(beta-1)))/BPV;

elseif flag == "DA"
    
    % do not consider ratios since Notionals aren't constant in any leg
    c = c.*Notionals;

    % updated function f(x) with different notionals
    f = @(x) sum(c.*B_fwd(2:end).*exp(-z(2:end)*x-z(2:end).^2/2)) + ...
             sum(Notionals.*B_fwd(2:end).*exp(-z(2:end)*x-z(2:end).^2/2)) - ...
             sum(Notionals.*beta.*B_fwd(1:end-1).*exp(-nu*x-nu.^2/2));

    % unique value of x such that f(x) = 0
    x_star = fzero(f,-1);

    % swaption Receiver price, updated due to different notionals in single amortizing
    SR = B(1)*(sum(c.*B_fwd(2:end).*normcdf(x_star+z(2:end))) + ...
               sum(Notionals.*B_fwd(2:end).*normcdf(x_star+z(2:end))) - ...
               sum(Notionals.*beta.*B_fwd(1:end-1).*normcdf(x_star+nu)));

    % forward discount factors B(0,Ti,Ti+1)
    B_fwd_i = B(2:end)./B(1:end-1);

    swap_rate = sum(B_fwd(1:end-1).*(beta - B_fwd_i).*Notionals)/BPV;

end

    % swaption Payer price thanks to Put-Call parity 
    SP = SR + B(1)*BPV*(swap_rate - r);

end