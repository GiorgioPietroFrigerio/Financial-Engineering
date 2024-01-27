function NPV_CVA = swap_pricing_tree(Prob,beta,NPV_risk_free,swap_maturity,settle,payment_dates,Notionals,Initial_Notional,K,a,sigma,discount_dates,discount_curve,freq,month_steps,LGD,flag)
% computes NPV of the swap considering CVA, by using trinomial tree
%
% INPUT
% Prob:              matrix of survival probabilities (each column for a different lambda)
% beta:              vector of betas
% NPV_risk_free:     Net Present Value without any risk assumption
% swap_maturity:     maturity of the swap in years
% settle:            settlement date
% payment_dates:     payment dates of the swap contract
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% K:                 strike, namely the fixed rate
% a:                 calibrated parameter a of the model
% sigma:             calibrated parameter sigma of the model
% discount_dates:    dates of the financial instruments used for the curve construction
% discount_curve:    discount factors curve
% freq:              annual payments frequency
% month_steps:       number of steps in a month
% LGD:               Loss Given Default
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% NPV_CVA:           Net Present Values considering CVA with different lambdas

% expiries in months
expiries_in_months = (1/freq:1/freq:(swap_maturity-1/freq))*12; 

% time step
dt_monthly = 1/month_steps;
dt = dt_monthly/12;

% Tree construction
sigma_tree = sigma*sqrt((1-exp(-2*a*dt))/(2*a));   % sigma_tree
mu_tree = 1-exp(-a*dt);                            % mu_tree
delta_x = sqrt(3)*sigma_tree;                      % delta_x
l_max = ceil((1-sqrt(2/3))/mu_tree);               % l_max
l_min = -l_max;                                    % l_min
l = (l_max :-1: l_min)';                           % vector of possible l
x = l*delta_x;                                     % vector of possible x in the tree

% transition probabilities in the cases A B and C
prob_A = [1/2*(1/3-l*mu_tree+(l*mu_tree).^2), 2/3-(l*mu_tree).^2, 1/2*(1/3+l*mu_tree+(l*mu_tree).^2)];                                         % transition probabilities for Case A
prob_B = [1/2*(7/3-3*l_max*mu_tree+(l_max*mu_tree).^2), -1/3+2*l_max*mu_tree-(l_max*mu_tree).^2, 1/2*(1/3-l_max*mu_tree+(l_max*mu_tree).^2)];  % transition probabilities for Case B
prob_C = [1/2*(1/3+l_min*mu_tree+(l_min*mu_tree).^2), -1/3-2*l_min*mu_tree-(l_min*mu_tree).^2, 1/2*(7/3+3*l_min*mu_tree+(l_min*mu_tree).^2)];  % transition probabilities for Case C

% probability matrix
prob_matrix = prob_A;
prob_matrix(end,:) = prob_B;
prob_matrix(1,:) = prob_C;

% sigma hat star
sigma_hat_star = sigma/a*sqrt(dt-2*(1-exp(-a*dt))/a +(1-exp(-2*a.*dt))/(2*a)); 

% auxiliary vectors to create the proper discount D(ti,ti+1)
aux_vec = exp(-0.5*sigma_hat_star^2-sigma_hat_star/sigma_tree*mu_tree*x);      % first vector of exponentials
vec_delta_x = exp(-sigma_hat_star/sigma_tree*exp(-a*dt)*delta_x*(2:-1:-2)');   % second vector of exponentials (with all possible delta x)

% % prices of the payer swaptions needed for the CVA
SwaptionPricing = @(i) swaption_price_tree(settle,discount_dates,discount_curve,payment_dates(i),expiries_in_months(i),payment_dates(end),sigma,a,month_steps,K,Notionals(i+1:end),Initial_Notional,beta(i:end),freq,x,prob_matrix,aux_vec,vec_delta_x,dt_monthly,l_max,flag);
   
PS = zeros(1,length(payment_dates)-1);
parfor i=1:length(payment_dates)-1
    PS(i) = feval(SwaptionPricing,i);
end  

% CVA computations
compute_CVA = @(i) LGD*PS*(Prob(1:end-1,i)-Prob(2:end,i));
CVA = cell2mat(arrayfun(compute_CVA,1:size(Prob,2),'UniformOutput',false));

% final NPVs
NPV_CVA = NPV_risk_free - CVA;

end