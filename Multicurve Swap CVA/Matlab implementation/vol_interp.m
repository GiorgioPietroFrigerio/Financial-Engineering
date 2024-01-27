function vols = vol_interp(settle,Notionals,Initial_Notional,discount_dates,discount_curve,B,payment_dates,expiries,tenors,freq,swap_maturity,vol_matrix,flag)
% creates a complete matrix of volatilities for each expiry and maturity,
% then interpolates in terms of amortized BPV
%
% INPUT
% settle:            settlement date  
% Notionals:         vector of amortized notionals
% Initial_Notional:  initial notional
% discount_dates:    dates of the financial instruments used for the curve construction
% discount_curve:    discount factors curve
% B:                 discount factors in each payment date
% payment_dates:     payment dates of the underlying swap
% expiries:          vector of expiries
% tenors:            vector of tenors
% freq:              yearly payment frequency
% swap_maturity:     maturity of the underlying swap
% vol_matrix:        matrix of market implied volatilities
% flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg
%
% OUTPUT
% vols:              needed volatilities for the corresponding swaptions computation

exp_val = payment_dates(1:end-1);               % expiries
ten_val = 1/freq:1/freq:(swap_maturity-1/freq); % tenors

% initialization of auxiliary matrix
mat_aux = [];

% construction of the first interpolation
for i = 1:size(vol_matrix,1)
    mat_aux = [mat_aux; interp1(tenors,vol_matrix(i,:),ten_val,'linear','extrap')];
end

% initialization of the matrix of interpolated volatilities
vol_interp = [];

% creation of volatility matrix for each expiry and tenor of plain vanilla swaptions
for i = 1:size(mat_aux,2)
    vol_interp = [vol_interp, interp1(expiries,mat_aux(:,i),exp_val,'linear','extrap')];
end

if flag == "SA" % single amortizing case

    % needed setting for regular and amortized BPVs computation
    fun = @(i) busdate(addtodate(payment_dates(end)-1,12/freq*i,"month"));
    I = (1:(length(payment_dates)-2))';
    additive_dates = arrayfun(fun,I);
    payment_dates = [payment_dates; additive_dates];
    deltas = yearfrac(payment_dates(1:end-1),payment_dates(2:end),2);
    B_add_dates = interp_ZR(settle,discount_dates,discount_curve,additive_dates);
    discounts = [B; B_add_dates];
    len = length(additive_dates)+2;
    vols = [];

    % interpolate in terms of amortized BPV 
    for i = 1:length(exp_val)
        disc_fwd = discounts((i+1):(len+i-1))./discounts(i);
        BPVs = cumsum(deltas(i:len+i-2).*disc_fwd);
        ratios = Notionals(i+1:end)./Initial_Notional;
        BPV_hat = sum(ratios.*deltas(i:len-1).*disc_fwd(1:len-i));
        vols = [vols; interp1(BPVs,vol_interp(i,:),BPV_hat,'linear','extrap')];
    end  

elseif flag == "DA" % double amortizing case 

    vols = diag(fliplr(vol_interp)); % we just need the main diagonal of the flipped volatility matrix

end

end 