function output_vector = ModelCalibration_Receiver(gamma,r,discount_dates,discount_curve,P_dates,P_curve,settle,expiries,tenors,freq,vol,init_cond)
% calibrates parameters a and sigma for the selected value of gamma by using Receiver swaptions
%
% INPUT
% gamma:           selected value of gamma
% r:               strike, namely the fixed rate
% discount_dates:  dates of the financial instruments used for the curve construction
% discount_curve:  discount factors curve
% P_dates:         dates of the financial instruments used for the pseudo-curve construction
% P_curve:         pseudo-discount factors curve
% settle:          settlement date
% expiries:        vector of expiries
% tenors:          vector of tenors
% freq:            annualy payment frequency
% vol:             vector of implied volatilities for the considered swaptions
% init_cond:       vector of initial conditions for both a and sigma
%
% OUTPUT
% output_vector:   a and sigma calibrated for the model, minimized distance, errors between market prices and model prices 

% expiries in terms of actual dates
upgradeExpiries = @(x) change_weekends(addtodate(settle,expiries(x),'year'));
expiriesDates = arrayfun(upgradeExpiries,1:length(expiries));

% payment dates in terms of actual dates
createDates = @(i) busdate(addtodate(settle-1,12/freq*i,'month'));
payDates = arrayfun(createDates,1:freq*max(expiries+tenors))';

% discounts, pseudo-discounts and betas
B = interp_ZR(settle,discount_dates,discount_curve,payDates); % discount factors in each payment date
B_fwd = B(2:end)./B(1:end-1);                                 % forward discount factors
BP = interp_ZR(settle,P_dates,P_curve,payDates);              % pseudo-discount factors in each payment date
BP_fwd = BP(2:end)./BP(1:end-1);                              % forward pseudo-discount factors 
beta = B_fwd./BP_fwd;                                         % beta between two payment dates

% market prices by using implied volatilies
MktPrice = @(i) SwaptionReceiver(B(expiries(i)*freq:(expiries(i)+tenors(i))*freq),beta(freq*expiries(i):(expiries(i)+tenors(i))*freq-1),payDates(expiries(i)*freq+1:(expiries(i)+tenors(i))*freq),expiriesDates(i),settle,r,vol(i));
Mkt_Prices = arrayfun(MktPrice,1:length(tenors));

% times to expiries of the swaptions
times_to_exp = yearfrac(settle,expiriesDates,3);

% time intervals between expiry and payment dates (with also expiry)
time_intervals_fun = @(i) yearfrac(expiriesDates(i),[expiriesDates(i); payDates(expiries(i)*freq+1:(expiries(i)+tenors(i))*freq)],3);
time_intervals = arrayfun(time_intervals_fun,1:length(expiries),'UniformOutput',false);

% deltas between two consecutive payment dates, but starting from the expiry
deltas_fun = @(i) yearfrac([expiriesDates(i); payDates(expiries(i)*freq+1:((expiries(i)+tenors(i))*freq-1))],payDates(expiries(i)*freq+1:(expiries(i)+tenors(i))*freq),2);
deltas = arrayfun(deltas_fun,1:length(expiries),'UniformOutput',false);

% model price of a single considered swaption
fun = @(i,parameters) SwaptionPricingModel_Receiver(gamma,r,B(expiries(i)*freq:(expiries(i)+tenors(i))*freq),beta(freq*expiries(i):(expiries(i)+tenors(i))*freq-1),parameters(1),parameters(2),times_to_exp(i),cell2mat(deltas(i)),cell2mat(time_intervals(i)));

% model prices of all the considered swaptions
ModelPrices = @(parameters) arrayfun(@(i) fun(i,parameters),1:length(tenors));

% function to minimize
distance = @(parameters) sum((Mkt_Prices-ModelPrices(parameters)).^2);

% % we performed an optimization of the initial conditions for fmincon in order to find the best...
% ... input for this function
% % this is the script that can be runned:
% Initial_conditions_optimization;

% model calibration & minimized distance
[param, dist] = fmincon(distance,init_cond,[],[],[],[],[],[],@(parameters) constraints(parameters));

% model prices computed with the obtained parameters a and sigma
Model_Prices = ModelPrices(param);

% difference between model prices and market prices
err = abs(Model_Prices-Mkt_Prices);

% plot of model and market prices
figure();
plot(expiries, Model_Prices, '-o','LineWidth',2);
hold on
plot(expiries, Mkt_Prices, '-o','LineWidth',2);
legend('Model Prices', 'Mkt Prices');
title_string = string(datestr(busdate(addtodate(settle,-1,'day'),"previous"))) + "   Gamma: " + string(gamma);
title('Receiver', title_string);
xlabel('expiries');
ylabel('prices')

% Generate parameter values
x = linspace(0.0001, 0.003, 25);
y = linspace(0.008, 0.01, 25);
[X, Y] = meshgrid(x, y);

% Evaluate the objective function for each combination of parameter values
Z = zeros(size(X));
parfor i = 1:numel(X)
    Z(i) = distance([X(i), Y(i)]);
end

% Create a 3D plot
figure;
surf(X, Y, Z);

% Set labels and title
xlabel('a');
ylabel('sigma');
zlabel('dist');
title('Distance: Receiver ', title_string)

% output as a, sigma, minimized distance, errors between model prices and market prices
output_vector = [param dist err];

end