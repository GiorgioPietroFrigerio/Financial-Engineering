function swaption_price = swaption_price_tree(settle,dates,discounts,expiry_date,expiry,maturity,sigma,a,month_steps,k,Notionals,Initial_Notional,beta,freq,x,prob_matrix,aux_vec,vec_delta_x,dt,l_max,flag)
% computes the price of a swaption payer with trinomial tree
%
% INPUT
% settle:            settlement date
% dates:             dates from bootstrap
% discounts:         discounts from bootstrap
% expiry_date:       expiry date of the swaption
% expiry:            expiry of the swaption in number of months
% maturity:          maturity date of the swap
% sigma:             calibrated sigma for HULL-WHITE model 
% a:                 calibrated a for HULL-WHITE model 
% month_steps:       number of monthly steps
% k:                 strike, namely the fixed rate
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% beta:              beta between two consecutive payment dates
% freq:              annual payments frequency
% x:                 vector of OU process values
% prob_matrix:       matrix of transition probabilities
% aux_vec:           first vector of exponentials to obtain D(ti,ti+1)
% vec_delta_x:       second vector of exponentials (with all possible delta x) to obtain D(ti,ti+1)
% dt:                monthly time step
% l_max:             l max
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% swaption_price:    price of the swaption 

% tree parameters
Dates_tree = dates_for_tree(settle,expiry_date,expiry,dt,month_steps);   % dates in the tree
dt = dt/12;                                                              % yearly time interval
deltas_365_tree = yearfrac(settle,[settle Dates_tree],3);                % time intervals (act/365)       
disc_fact_tree = interp_ZR(settle,dates,discounts,[settle Dates_tree]);  % discount factors in dates_tree
forward_disc_tree = disc_fact_tree(2:end)./disc_fact_tree(1:end-1);      % forward discount factors in dates_tree   

% swap parameters
Start_date = expiry_date;                                                           % start date of the swap (expiry of swaption)
PaymentDates = cfdates(Start_date,maturity-1,freq);                                 % payment dates of the coupons (swap)
deltas_365 = yearfrac(settle,[Start_date PaymentDates],3);                          % time intervals (act/365)
deltas_360 = yearfrac([Start_date PaymentDates(1:end-1)], PaymentDates(1:end),2);   % time intervals (act/360)
disc_fact = interp_ZR(settle,dates,discounts,[Start_date PaymentDates]);            % discount factors at payment dates
forw_disc_fact = disc_fact/disc_fact(1);                                            % forward discount factors

% Number of steps of the tree
N = month_steps*expiry+1;      

% volatility of x at the expiry of the swaption
Z = sqrt( sigma^2*(1-exp(-2*a*deltas_365(1)))/(2*a) );

if flag == "SA"

    % coupons for the CB
    coupons = k.*deltas_360'.*Notionals/Initial_Notional;                               
    coupons(end) = coupons(end) + 1; 

    % coupond bond in the strike of the put option
    CB_strike_fun = @(x) MHW_price(x,a,sigma,(beta(2:end)-1)',forw_disc_fact(2:end-1),deltas_365(1:end-1),expiry,freq);
    CB_strike = arrayfun(CB_strike_fun,x/Z);
    
    % put option strike
    option_strike = beta(1) + CB_strike;

elseif flag == "DA"

    % coupons for the CB
    coupons = k*deltas_360'.*Notionals;
    coupons(end) = coupons(end) + Notionals(end);
    
    % coupons for option strike
    strike_coupons = beta(2:end).*Notionals(2:end)-Notionals(1:end-1);

    % CB for option strike
    CB_strike_fun = @(x) MHW_price(x,a,sigma,strike_coupons',forw_disc_fact(2:end-1),deltas_365(1:end-1),expiry,freq);
    CB_strike = arrayfun(CB_strike_fun,x/Z);

    % put option strike
    option_strike = beta(1)*Notionals(1) + CB_strike;

end

% CB in T_alpha
fun_CB = @(x) MHW_price(x/Z,a,sigma,coupons',forw_disc_fact(2:end),deltas_365,expiry,freq);
CB = arrayfun(fun_CB,x);

% payoffs of the swaption (payer case)
payoff_payer = max(option_strike-CB,0);    

% allocation of matrices of nodes
tree_matrix_payer = zeros(length(x),N);

% last columns of the matrices are the payoffs
tree_matrix_payer(:,end) = payoff_payer;         

% construction of the matrix of discount factors B(ti,ti+1) depending on x                                   
I = 1:(length(deltas_365_tree)-1);                                          % vector of indexes
Z = sqrt(sigma^2*(1-exp(-2*a*dt*(I-1)))/(2*a));                             % Z
v = Z.*(1-exp(-a*dt))/a;                                                    % v
z = v;                                                                      % z (=v if gamma=0)
matrix_1 = -x.*z./Z;                                                        % first part of the exponential
matrix_1(:,1) = zeros(size(x));                                             % first column of zeros
matrix_2 = -0.5*z.^2;                                                       % second part of the exponential
matrix_exp = exp(matrix_1 + matrix_2);                                      % matrix with exponentials
matrix_disc_fact = forward_disc_tree.*matrix_exp;                           % matrix with discunts B(ti,ti+1)

% backward algorithm (to discount the payoffs of the swaption)
for j = N-1:-1:1
    % l_max
    tree_matrix_payer(1,j) = tree_matrix_payer(1,j+1)*prob_matrix(1,1)*matrix_disc_fact(1,j)*aux_vec(1)*vec_delta_x(3) + tree_matrix_payer(2,j+1)*prob_matrix(1,2)*matrix_disc_fact(1,j)*aux_vec(2)*vec_delta_x(4)... 
                             + tree_matrix_payer(3,j+1)*prob_matrix(1,3)*matrix_disc_fact(1,j)*aux_vec(3)*vec_delta_x(5);
    % l_min
    tree_matrix_payer(end,j) = tree_matrix_payer(end-2,j+1)*prob_matrix(end,1)*matrix_disc_fact(end,j)*aux_vec(end-2)*vec_delta_x(1) + tree_matrix_payer(end-1,j+1)*prob_matrix(end,2)*matrix_disc_fact(end,j)*aux_vec(end-1)*vec_delta_x(2)... 
                               + tree_matrix_payer(end,j+1)*prob_matrix(end,3)*matrix_disc_fact(end,j)*aux_vec(end)*vec_delta_x(3);
    % medium
    tree_matrix_payer(2:end-1,j) = tree_matrix_payer(1:end-2,j+1).*prob_matrix(2:end-1,1).*matrix_disc_fact(2:end-1,j).*aux_vec(1:end-2)*vec_delta_x(2) + tree_matrix_payer(2:end-1,j+1).*prob_matrix(2:end-1,2).*matrix_disc_fact(2:end-1,j).*aux_vec(2:end-1)*vec_delta_x(3)...
                                   + tree_matrix_payer(3:end,j+1).*prob_matrix(2:end-1,3).*matrix_disc_fact(2:end-1,j).*aux_vec(3:end)*vec_delta_x(4);
end

% price in t_0 of the swaption payer
if flag == "SA"
    swaption_price = Initial_Notional*tree_matrix_payer(l_max+1,1);  
elseif flag == "DA"
    swaption_price = tree_matrix_payer(l_max+1,1);
end

end 