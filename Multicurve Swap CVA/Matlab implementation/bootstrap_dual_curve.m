function [B_dates, B_curve, P_dates, P_curve] = bootstrap_dual_curve(data)
% computes discount factors and pseudo-discount factors in the multicurve framework
% 
% INPUT
% data:     struct containing settlement and quoted prices of the financial instruments 
%           needed for the bootstrap of both curves
% 
% OUTPUT
% B_dates:  dates of the financial instruments used for the curve construction
% B_curve:  discount factors curve
% P_dates:  dates of the financial instruments used for the pseudo-curve construction
% P_curve:  pseudo-discount factors curve

% parameters
settle = data.settlement;   % settlement date
Basis_swap = 6;             % 30/360 European
Basis_swap_floating = 2;    % act/360
Basis_depos = 2;            % depos act/360
Basis_fra = 2;              % act/360
Basis_futures = 2;          % futures act/360

% discount curve
[B_dates, B_curve,data] = discount_curve_bootstrap(data,settle,Basis_swap);

% Pseudo discount curve
daycounts = [Basis_depos, Basis_fra, Basis_futures, Basis_swap, Basis_swap_floating];
[P_dates, P_curve] = pseudo_discount_curve_bootstrap(data,settle,daycounts,B_dates,B_curve);

% plot of discounts and pseudo-discounts
figure();
plot(B_dates, B_curve,P_dates, P_curve,'LineWidth',2);
legend('Discount','Pseudo-Discount');
title('Discounts and Pseudo-discounts on ',datestr(data.refDate));
xlabel('years');
ylim([min([B_curve;P_curve])-0.05 max([B_curve;P_curve])+0.05]);
datetick

end