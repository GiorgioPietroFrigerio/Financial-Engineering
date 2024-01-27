function [Prob, beta, deltas] = swap_elements(Bp,B,payment_dates,settle,spreads,LGD)
% computes probabilities, betas and time intervals
%
% INPUT
% Bp:             pseudo-discount factors
% B:              discount factors
% payment_dates:  payment dates
% settle:         settlement date
% spreads:        spreads
% LGD:            loss given default
%
% OUTPUT
% Prob:           probabilities
% beta:           betas
% deltas:         time intervals

lambdas = spreads/LGD;           % vector of intensities lambda
Bp_fwd = Bp(2:end)./Bp(1:end-1); % forward pseudo-discount factors
B_fwd = B(2:end)./B(1:end-1);    % forward discount factors
beta = B_fwd./Bp_fwd;            % beta

% time intervals
deltas = yearfrac(payment_dates(1:end-1),payment_dates(2:end),2);

% survival probabilities for each lambda
Prob_daycount = 3; % act/365
Probability = @(lambda) exp(-lambda.*yearfrac(settle,[settle; payment_dates(1:end-1)],Prob_daycount)); 
Prob = cell2mat(arrayfun(Probability,lambdas,'UniformOutput',false));

end