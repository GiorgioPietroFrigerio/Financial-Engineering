function NPV_CVA = swap_pricing_CVA_Payer_Model(Prob,beta,NPV_risk_free,settle,payment_dates,B,Notionals,Initial_Notional,K,LGD,a,sigma,flag)
% computes NPV of the swap considering CVA, by using the calibrated parameters of the Model
%
% INPUT
% Prob:              matrix of survival probabilities (each column for a different lambda)
% beta:              vector of betas
% NPV_risk_free:     Net Present Value without any risk assumption
% settle:            settlement date
% payment_dates:     payment dates of the swap contract
% B:                 discount factors in each payment date
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% K:                 strike, namely the fixed rate
% LGD:               Loss Given Default
% a:                 calibrated parameter a of the model
% sigma:             calibrated parameter sigma of the model
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% NPV_CVA:           Net Present Values considering CVA with different lambdas

% gamma in our setting
gamma = 0;

% prices of the payer swaptions needed for the CVA
SwaptionPricing = @(i) SwaptionPricingModel_PayerCVA(gamma,K,B(i:end),beta(i:end),a,sigma,settle,payment_dates(i),payment_dates(i+1:end),Notionals(i+1:end),Initial_Notional,flag);
PS = arrayfun(SwaptionPricing,1:(length(payment_dates)-1));

% CVAs computation
compute_CVA = @(i) LGD*PS*(Prob(1:end-1,i)-Prob(2:end,i));
CVA = cell2mat(arrayfun(compute_CVA,1:size(Prob,2),'UniformOutput',false));

% final NPVs
NPV_CVA = NPV_risk_free - CVA;

end