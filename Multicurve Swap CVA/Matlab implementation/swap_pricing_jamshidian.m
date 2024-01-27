function NPV_CVA = swap_pricing_jamshidian(Prob,beta,NPV_risk_free,swap_maturity,settle,payment_dates,Notionals,Initial_Notional,K,a,sigma,discount_dates,discount_curve,freq,LGD,flag)
% computes NPV of the swap considering CVA, by using Jamshidian approach
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
% LGD:               Loss Given Default
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% NPV_CVA:           Net Present Values considering CVA with different lambdas

% expiries and tenors in months
expiries_in_months = (1/freq:1/freq:(swap_maturity-1/freq))*12;
tenors_in_months = ((swap_maturity-1/freq):-1/freq:1/freq)*12;

% prices of the payer swaptions needed for the CVA
SwaptionPricing = @(x) swaption_price_Jam(settle,discount_dates,discount_curve,payment_dates(x),expiries_in_months(x),tenors_in_months(x),payment_dates(end),sigma,a,K,Notionals(x+1:end),Initial_Notional,beta(x:end),freq,flag);
PS = arrayfun(SwaptionPricing,1:(length(payment_dates)-1));

% CVAs computation
compute_CVA = @(i) LGD*PS*(Prob(1:end-1,i)-Prob(2:end,i));
CVA = cell2mat(arrayfun(compute_CVA,1:size(Prob,2),'UniformOutput',false));

% final NPVs
NPV_CVA = NPV_risk_free - CVA;

end