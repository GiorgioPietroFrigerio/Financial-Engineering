% Project 1B: Artioli Fiammetta, Cantù Davide, Frigerio Giorgio Pietro
clc;
close all;
clear all;

%% EX1 
%% Reading data from Excel (24 Jun 2022)
N_sheets = [1 2];
formatData = 'dd/mm/yyyy';
data_24_jun = readExcelData('CVAProjectMktData',formatData,N_sheets(1));

%% Bootstrap 24 Jun 2022
[discount_dates_24_jun, discount_curve_24_jun, P_dates_24_jun, P_curve_24_jun] = bootstrap_dual_curve(data_24_jun);

%% EX 2 
%% Reading swap data from excel 24 Jun 2022
SwapData = readExcelData_SwapAnnex('SwapAmortizingPlan',formatData);

%% Swap pricing risk free 24 Jun 2022
K = 0.0221;                                % strike
settle_1 = data_24_jun.settlement;         % settlement date for 24 jun 2022 case
Initial_Notional = SwapData.Notionals(1);  % initial notional

% discounts, pseudo-discounts and forward pseudo-discounts for swap
[B_24_jun, Bp_24_jun, Bp_fwd_24_jun] = swap_discounts(settle_1,discount_dates_24_jun,discount_curve_24_jun,P_dates_24_jun,P_curve_24_jun,SwapData.EndDates);

% swap NPV risk free (only fixed leg amortizing case)
NPV_SA_risk_free_24_jun = swap_pricing_risk_free(SwapData.StartDates,SwapData.EndDates,B_24_jun,Bp_fwd_24_jun,SwapData.Notionals,Initial_Notional,K,"SA");

% swap NPV risk free (both legs amortizing case)
NPV_DA_risk_free_24_jun = swap_pricing_risk_free(SwapData.StartDates,SwapData.EndDates,B_24_jun,Bp_fwd_24_jun,SwapData.Notionals,Initial_Notional,K,"DA");

%% EX 3
%% Reading swaption data from excel 24 Jun 2022
VolData_1 = readExcelData_SwaptionVol('CVAProjectMktData',formatData,N_sheets(1));

%% Swap pricing with CVA 24 Jun 2022
LGD = 0.40;                          % loss given default
spreads = [300, 500]*1e-4;           % vector of spreads
freq = 4;                            % yearly payment frequency
swap_maturity = 15;                  % swap maturity in years

% probabilities, betas and deltas for swap pricing with CVA
[Prob_24_jun, beta_24_jun, deltas_24_jun] = swap_elements(Bp_24_jun,B_24_jun,SwapData.PayDates,settle_1,spreads,LGD);

% interpolated volatilities
vol_SA_24_jun = vol_interp(settle_1,SwapData.Notionals,Initial_Notional,discount_dates_24_jun,discount_curve_24_jun,B_24_jun,SwapData.PayDates,VolData_1.expiry,VolData_1.tenor,freq,swap_maturity,VolData_1.vol,"SA");
vol_DA_24_jun = vol_interp(settle_1,SwapData.Notionals,Initial_Notional,discount_dates_24_jun,discount_curve_24_jun,B_24_jun,SwapData.PayDates,VolData_1.expiry,VolData_1.tenor,freq,swap_maturity,VolData_1.vol,"DA");

% swap NPVs with CVA (only fixed leg amortizing case)
NPVs_CVA_SA_24_jun = swap_pricing_CVA(Prob_24_jun,beta_24_jun,deltas_24_jun,NPV_SA_risk_free_24_jun,vol_SA_24_jun,settle_1,SwapData.PayDates,B_24_jun,SwapData.Notionals,Initial_Notional,K,LGD,"SA");  

% swap NPVs with CVA (both legs amortizing case)
NPVs_CVA_DA_24_jun = swap_pricing_CVA(Prob_24_jun,beta_24_jun,deltas_24_jun,NPV_DA_risk_free_24_jun,vol_DA_24_jun,settle_1,SwapData.PayDates,B_24_jun,SwapData.Notionals,Initial_Notional,K,LGD,"DA");  

%% EX 4
%% Reading data from Excel (31 Jan 2023)
data_31_jan = readExcelData('CVAProjectMktData',formatData,N_sheets(2));

%% Bootstrap 31 Jan 2023
[discount_dates_31_jan, discount_curve_31_jan, P_dates_31_jan, P_curve_31_jan] = bootstrap_dual_curve(data_31_jan);

%% Reading swaption data from excel 31 Jan 2023
VolData_2 = readExcelData_SwaptionVol('CVAProjectMktData',formatData,N_sheets(2));

%% Swap pricing risk free (unwinding at 31 Jan 2023) 
settle_2 = data_31_jan.settlement;                   % settlement date for 31 jan 2023 case
Euribor_28_12 = 0.02202;                             % Euribor 3m on 28 dec 2022 from market data
idx = find(settle_2<SwapData.EndDates);              % indexes of dates of interest

% discounts, pseudo-discounts and forward pseudo-discounts for swap
[B_31_jan, Bp_31_jan, Bp_fwd_31_jan] = swap_discounts_modified(settle_2,discount_dates_31_jan,discount_curve_31_jan,P_dates_31_jan,P_curve_31_jan,SwapData.StartDates,SwapData.EndDates,Euribor_28_12,idx);

% swap NPV risk free (only fixed leg amortizing case)
NPV_SA_risk_free_31_jan = swap_pricing_risk_free(SwapData.StartDates(idx),SwapData.EndDates(idx),B_31_jan,Bp_fwd_31_jan,SwapData.Notionals(idx),Initial_Notional,K,"SA");

% swap NPV risk free (both legs amortizing case)
NPV_DA_risk_free_31_jan = swap_pricing_risk_free(SwapData.StartDates(idx),SwapData.EndDates(idx),B_31_jan,Bp_fwd_31_jan,SwapData.Notionals(idx),Initial_Notional,K,"DA");

%% Swap pricing with CVA (unwinding at 31 Jan 2023)
swap_maturity_new = swap_maturity - (idx(1) - 1)/freq;   % new swap maturity in years

% probabilities, betas and deltas for swap pricing with CVA
[Prob_31_jan, beta_31_jan, deltas_31_jan] = swap_elements(Bp_31_jan,B_31_jan,SwapData.PayDates(idx),settle_2,spreads,LGD);

% interpolated volatilities
vol_SA_31_jan = vol_interp(settle_2,SwapData.Notionals(idx),Initial_Notional,discount_dates_31_jan,discount_curve_31_jan,B_31_jan,SwapData.PayDates(idx),VolData_2.expiry,VolData_2.tenor,freq,swap_maturity_new,VolData_2.vol,"SA");
vol_DA_31_jan = vol_interp(settle_2,SwapData.Notionals(idx),Initial_Notional,discount_dates_31_jan,discount_curve_31_jan,B_31_jan,SwapData.PayDates(idx),VolData_2.expiry,VolData_2.tenor,freq,swap_maturity_new,VolData_2.vol,"DA");

% swap NPVs with CVA (only fixed leg amortizing case)
NPVs_CVA_SA_31_jan = swap_pricing_CVA(Prob_31_jan,beta_31_jan,deltas_31_jan,NPV_SA_risk_free_31_jan,vol_SA_31_jan,settle_2,SwapData.PayDates(idx),B_31_jan,SwapData.Notionals(idx),Initial_Notional,K,LGD,"SA");  

% swap NPVs with CVA (both legs amortizing case)
NPVs_CVA_DA_31_jan = swap_pricing_CVA(Prob_31_jan,beta_31_jan,deltas_31_jan,NPV_DA_risk_free_31_jan,vol_DA_31_jan,settle_2,SwapData.PayDates(idx),B_31_jan,SwapData.Notionals(idx),Initial_Notional,K,LGD,"DA");  

%% EX 5 
%% settings for model calibration
gamma = [0; 0.5; 1];                % fixed values of gamma
expiries = [1 3 5 8 10 12 15];      % expiries of the considered swaptions for the calibration
tenors = [15 12 10 7 5 3 1];        % tenors of the considered swaptions for the calibration
init_cond_24_jun = [0.001, 0.01];   % initial conditions for 24 June 2022
init_cond_31_jan = [0.005, 0.005];  % initial conditions for 31 January 2023

% volatilities for calibration
vol_cal_24_jun = readExcelData_vol_calibration('CVAProjectMktData',N_sheets(1),expiries,tenors);
vol_cal_31_jan = readExcelData_vol_calibration('CVAProjectMktData',N_sheets(2),expiries,tenors);

% each of the following matrix of parameters is structured in the following way :
% calibrated a, calibrated sigma, minimized distance, differences between model and market prices
% where the i-th row corresponds to the i-th value of gamma

%% Model Calibration considering Receiver swaptions at 24 June 2022
calibration_24_jun_Receiver = @(x)  ModelCalibration_Receiver(x,K,discount_dates_24_jun, discount_curve_24_jun,P_dates_24_jun, P_curve_24_jun,settle_1,expiries,tenors,freq,vol_cal_24_jun, init_cond_24_jun);
parameters_24_jun_Receiver = cell2mat(arrayfun(calibration_24_jun_Receiver,gamma,'UniformOutput',false));

%% Model Calibration considering Receiver swaptions at 31 January 2023
calibration_31_jan_Receiver = @(x)  ModelCalibration_Receiver(x,K,discount_dates_31_jan, discount_curve_31_jan,P_dates_31_jan, P_curve_31_jan,settle_2,expiries,tenors,freq,vol_cal_31_jan, init_cond_31_jan);
parameters_31_jan_Receiver = cell2mat(arrayfun(calibration_31_jan_Receiver,gamma,'UniformOutput',false));

%% Model Calibration considering Payer swaptions at 24 June 2022
calibration_24_jun_Payer = @(x)  ModelCalibration_Payer(x,K,discount_dates_24_jun, discount_curve_24_jun,P_dates_24_jun, P_curve_24_jun,settle_1,expiries,tenors,freq,vol_cal_24_jun,init_cond_24_jun);
parameters_24_jun_Payer = cell2mat(arrayfun(calibration_24_jun_Payer,gamma,'UniformOutput',false));

%% Model Calibration considering Payer swaptions at 31 January 2023
calibration_31_jan_Payer = @(x)  ModelCalibration_Payer(x,K,discount_dates_31_jan, discount_curve_31_jan,P_dates_31_jan, P_curve_31_jan,settle_2,expiries,tenors,freq,vol_cal_31_jan,init_cond_31_jan);
parameters_31_jan_Payer = cell2mat(arrayfun(calibration_31_jan_Payer,gamma,'UniformOutput',false));

%% EX 6 
%% NPV with CVA with alternative approaches assuming Gamma = 0
month_steps = 6;
a = parameters_24_jun_Receiver(1,1); 
sigma = parameters_24_jun_Receiver(1,2);

%% Model used for the calibration
NPVs_CVA_SA_Model = swap_pricing_CVA_Payer_Model(Prob_24_jun,beta_24_jun,NPV_SA_risk_free_24_jun,settle_1,SwapData.PayDates,B_24_jun,SwapData.Notionals,Initial_Notional,K,LGD,a,sigma,"SA");
NPVs_CVA_DA_Model = swap_pricing_CVA_Payer_Model(Prob_24_jun,beta_24_jun,NPV_DA_risk_free_24_jun,settle_1,SwapData.PayDates,B_24_jun,SwapData.Notionals,Initial_Notional,K,LGD,a,sigma,"DA");

%% Generalized trinomial Tree
NPVs_CVA_SA_tree = swap_pricing_tree(Prob_24_jun,beta_24_jun,NPV_SA_risk_free_24_jun,swap_maturity,settle_1,SwapData.PayDates,SwapData.Notionals,Initial_Notional,K,a,sigma,discount_dates_24_jun,discount_curve_24_jun,freq,month_steps,LGD,"SA");
NPVs_CVA_DA_tree = swap_pricing_tree(Prob_24_jun,beta_24_jun,NPV_DA_risk_free_24_jun,swap_maturity,settle_1,SwapData.PayDates,SwapData.Notionals,Initial_Notional,K,a,sigma,discount_dates_24_jun,discount_curve_24_jun,freq,month_steps,LGD,"DA");

%% Jamshidian approach
NPVs_CVA_SA_Jam = swap_pricing_jamshidian(Prob_24_jun,beta_24_jun,NPV_SA_risk_free_24_jun,swap_maturity,settle_1,SwapData.PayDates,SwapData.Notionals,Initial_Notional,K,a,sigma,discount_dates_24_jun,discount_curve_24_jun,freq,LGD,"SA");
NPVs_CVA_DA_Jam = swap_pricing_jamshidian(Prob_24_jun,beta_24_jun,NPV_DA_risk_free_24_jun,swap_maturity,settle_1,SwapData.PayDates,SwapData.Notionals,Initial_Notional,K,a,sigma,discount_dates_24_jun,discount_curve_24_jun,freq,LGD,"DA");

