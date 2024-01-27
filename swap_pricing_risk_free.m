function NPV_risk_free = swap_pricing_risk_free(start_dates,end_dates,B,Bp_fwd,Notionals,Initial_Notional,K,flag)
% computes swap price without any risk assumption
%
% INPUT
% start_dates:       start date of each payment period
% end_dates:         end date of each payment period
% B:                 discount factors in each payment date
% Bp_fwd:            forward pseudo-discount factors between two consecutive payment dates 
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% K:                 strike, namely the fixed rate
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% NPV_risk_free:     NPV of the swap as the difference between between floating and fixed leg

% day count for swap
swap_daycount = 2; % act/360

% time intervals
deltas = yearfrac(start_dates,end_dates,swap_daycount);

% floating NPV
if flag == "SA" 

    NPVA = sum(Initial_Notional*B.*(1./Bp_fwd-1));  % only fixed leg amortizing case

elseif flag == "DA"

    NPVA = sum(B.*Notionals.*(1./Bp_fwd-1));        % both legs amortizing case
    
end

% fixed NPV
NPVB = K * sum(B.*deltas.*Notionals);

% swap NPV risk free
NPV_risk_free = NPVA - NPVB;

end 