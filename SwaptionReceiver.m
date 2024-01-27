function P = SwaptionReceiver(B,beta,payDates,Ti,settle,K,vol)
% computes price of a Receiver Swaption with underlying in Ti and ending in Tw by using Bachelier formula
%
% INPUT
% B:         vector of discount factors from B(0,Ti) to B(0,Tw)
% beta:      vector of spreads beta(0,Ti,Ti+1) with time window 3 months
% payDates:  payment dates of the underlying swap
% Ti:        start date of the underlying swap
% settle:    settlement date
% K:         strike -> fixed rate
% vol:       implied volatility -> sigma(Ti,Tw)
%
% OUTPUT
% P:         price of the Swaption
 
% deltas between two consecutive payment dates, but starting from the expiry Ti
deltas = yearfrac([Ti; payDates(1:end-1)],payDates,2);

% forward discount factors B(0,expiry,Ti)
B_fwd = B/B(1); 

% Niw (t0)
N_i_w = 1 - B_fwd(end) + sum(B_fwd(1:end-1).*(beta-1));
 
% BPViw (t0)
BPV = sum(deltas.*B_fwd(2:end)); 

% Siw (t0)
S = N_i_w/BPV; 

% time to maturity of the swaption
time_interval = yearfrac(settle,Ti,3);

% dn
dn = (S - K)/(vol*sqrt(time_interval));

% price
P = B(1)*BPV*((K - S)*normcdf(-dn) + vol*sqrt(time_interval)*normpdf(dn));

end