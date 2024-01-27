function NPV_CVA = swap_pricing_CVA(Prob,beta,deltas,NPV_risk_free,vol,settle,payment_dates,B,Notionals,Initial_Notional,K,LGD,flag)
% computes swap pricing taking into account also counterparty risk
%
% INPUT
% Prob:              matrix of survival probabilities (each column for a different lambda)
% beta:              vector of betas
% deltas:            time intervals between payment dates
% NPV_risk_free:     Net Present Value without any risk assumption
% vol:               vector of interpolated implied volatilities for the amortized case
% settle:            settlement date
% payment_dates:     payment dates of the swap contract
% B:                 discount factors in each payment date
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% K:                 strike, namely the fixed rate
% LGD:               Loss Given Default
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% NPV_CVA:           Net Present Values considering CVA with different lambdas

% compute price of each swaption on the residual value of the underlying swap, moving forward on the expiry of the swaption
SwaptionPricing = @(x) SwaptionPayer(B(x:end),beta(x:end),deltas(x:end),Notionals(x:end),Initial_Notional,payment_dates(x),settle,K,vol(x),flag);
PS = arrayfun(SwaptionPricing,1:(length(payment_dates)-1));

% CVAs computation
compute_CVA = @(i) LGD*PS*(Prob(1:end-1,i)-Prob(2:end,i));
CVA = cell2mat(arrayfun(compute_CVA,1:size(Prob,2),'UniformOutput',false));

% final NPVs
NPV_CVA = NPV_risk_free - CVA;

end