function P = SwaptionPayer(B,beta,deltas,N,Initial_Notional,Ti,settle,K,vol,flag)
% computes price of a Payer Swaption with underlying in Ti and ending in Tw
%
% INPUT
% B:                 vector of discount factors from B(0,Ti) to B(0,Tw)
% beta:              vector of spreads beta(0,Ti,Ti+1) with time window 3 months
% deltas:            vector of delta(Ti,Ti+1) with time window 3 months
% N:                 vector of notionals from Ti until Tw-1
% Initial_Notional:  initial notional
% Ti:                start date of the underlying swap
% settle:            settlement date
% K:                 strike, namely the fixed rate
% vol:               implied volatility, namely sigma(Ti,Tw)
% flag:              "SA" -> Single Amortizing, "DA" -> Double Amortizing
%
% OUTPUT
% P:                 price of the Swaption

% forward discount factors B(0,expiry,Ti)
B_fwd_expiry = B/B(1);          

if flag == "SA"  

    % single amortizing case
    N_i_w = Initial_Notional*(1 - B_fwd_expiry(end) + sum(B_fwd_expiry(1:end-1).*(beta-1)));

elseif flag == "DA" 

    % double amortizing case 
    N_i_w = sum(B_fwd_expiry(1:end-1).*(beta.*N(2:end)-N(1:end-1))) + N(1) - B_fwd_expiry(end)*N(end); 

end        

% BPV_i_w (t0)
BPV = sum(N(2:end).*deltas.*B_fwd_expiry(2:end)); 

% S_i_w (t0)
S = N_i_w/BPV; 

% time to maturity of the swaption
time_interval = yearfrac(settle,Ti,3);

% dn
dn = (S - K)/(vol*sqrt(time_interval));

% Swaption price with expiry Ti
P = B(1)*BPV*((S-K)*normcdf(dn) + vol*sqrt(time_interval)*normpdf(dn));

end