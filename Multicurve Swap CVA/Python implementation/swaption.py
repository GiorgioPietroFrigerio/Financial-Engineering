import numpy as np
from scipy.stats import norm
from scipy.optimize import fsolve, minimize, LinearConstraint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.model_selection import ParameterGrid
from utilities import *

def SwaptionPayer(B,beta,deltas,N,Initial_Notional,Ti,settle,K,vol,flag):
    """
    Function that computes price of a Payer Swaption with underlying in Ti and ending in Tw

    INPUT
    B:                 vector of discount factors from B(0,Ti) to B(0,Tw)
    beta:              vector of spreads beta(0,Ti,Ti+1) with time window 3 months
    deltas:            vector of delta(Ti,Ti+1) with time window 3 months
    N:                 vector of notionals from Ti until Tw-1
    Initial_Notional:  initial notional
    Ti:                start date of the underlying swap
    settle:            settlement date
    K:                 strike, namely the fixed rate
    vol:               implied volatility, namely sigma(Ti,Tw)
    flag:              "SA" -> Single Amortizing, "DA" -> Double Amortizing

    OUTPUT
    P:                 price of the Swaption
    """

    B_fwd_expiry = B / B[0]         # forward discount factors B(0,expiry,Ti)
    B_fwd = B[1:]/B[:-1]            # forward discount factors B(0,Ti,Ti+1)
    Bp_fwd = B_fwd/beta             # forward pseudo discounts

    if flag == "SA": #single amortizing case
        N_i_w = Initial_Notional * (1-B_fwd_expiry[-1] + sum(B_fwd_expiry[:-1]*(beta - 1)))
    elif flag == "DA":
        N_i_w = sum(B_fwd_expiry[1:] * N[1:] * (1/Bp_fwd-1))

    # BPV_i_w(t0)
    BPV = sum(N[1:]*deltas * B_fwd_expiry[1:])

    # S_i_w(t0)
    S = N_i_w / BPV

    # dn
    dn = (S - K) / (vol * np.sqrt(yearfrac(settle, np.array([Ti]), 3)))

    # Swaption price with expiry Ti
    P = B[0] * BPV * ((S - K) * norm.cdf(dn) + vol * np.sqrt(yearfrac(settle, np.array([Ti]), 3)) * norm.pdf(dn))

    return P

def SwaptionPayer_calib(B,beta,payDates,Ti,settle,K,vol):
    """
    Function that computes price of a Payer Swaption with underlying in Ti and ending in Tw
    by using Bachelier formula

    INPUT
    B:         vector of discount factors from B(0,Ti) to B(0,Tw)
    beta:      vector of spreads beta(0,Ti,Ti+1) with time window 3 months
    payDates:  payment dates of the underlying swap
    Ti:        start date of the underlying swap
    settle:    settlement date
    K:         strike -> fixed rate
    vol:       implied volatility -> sigma(Ti,Tw)

    OUTPUT
    P:         price of the Swaption
    """
    # time intervals
    deltas = yearfrac(np.hstack((Ti, payDates[:-1])),payDates,2)

    # forward discount factors B(0,expiry,Ti)
    B_fwd = B/B[0]

    # Niw (t0)
    N_i_w = 1 - B_fwd[-1] + sum(B_fwd[:-1]*(beta-1))

    # BPViw (t0)
    BPV = sum(deltas*B_fwd[1:])

    # Siw (t0)
    S = N_i_w/BPV

    # dn
    dn = (S - K)/(vol*np.sqrt(yearfrac(settle,np.array([Ti]),3)))

    # price
    P = B[0]*BPV*((S-K)*norm.cdf(dn) + vol*np.sqrt(yearfrac(settle,np.array([Ti]),3))*norm.pdf(dn))

    return P

def SwaptionReceiver(B,beta,payDates,Ti,settle,K,vol):
    """
    Function that computes price of a Receiver Swaption with underlying in Ti and ending in Tw by using Bachelier formula

    INPUT
    B:         vector of discount factors from B(0,Ti) to B(0,Tw)
    beta:      vector of spreads beta(0,Ti,Ti+1) with time window 3 months
    payDates:  payment dates of the underlying swap
    Ti:        start date of the underlying swap
    settle:    settlement date
    K:         strike -> fixed rate
    vol:       implied volatility -> sigma(Ti,Tw)

    OUTPUT
    P:         price of the Swaption
    """
    # time intervals
    deltas = yearfrac(np.hstack((Ti, payDates[:-1])),payDates,2)

    # forward discount factors B(0,expiry,Ti)
    B_fwd = B/B[0]

    # Niw (t0)
    N_i_w = 1 - B_fwd[-1] + sum(B_fwd[:-1]*(beta-1))

    # BPViw (t0)
    BPV = sum(deltas*B_fwd[1:])

    # Siw (t0)
    S = N_i_w/BPV

    # dn
    dn = (S - K)/(vol*np.sqrt(yearfrac(settle,np.array([Ti]),3)))

    # price
    P = B[0]*BPV*((K-S)*norm.cdf(-dn) + vol*np.sqrt(yearfrac(settle,np.array([Ti]),3))*norm.pdf(dn))

    return P

def SwaptionPricingModel_Payer(gamma,r,B,beta,a,sigma,settle,expiry,payDates):
    """
    Function that computes price of the single Payer swaption by using MHW model

    INPUT
    gamma:             fixed value of gamma
    r:                 strike, namely the fixed rate
    B:                 discount factors in each payment date
    beta:              vector of beta between two consecutive payment dates
    a:                 calibrated parameter a of the model
    sigma:             calibrated parameter sigma of the model
    time_interval_Z:   time to swaption expiry
    deltas:            deltas between two consecutive payment dates, but starting from the expiry
    time_intervals_v:  time intervals between expiry and payment dates (with also expiry)

    OUTPUT
    SP:                price of the Payer swaption
    """

    # deltas between two consecutive payment dates, but starting from the expiry
    delta = yearfrac(np.hstack((expiry, payDates[:-1])),payDates,2)

    # Coupon vector
    R = r*np.ones(np.shape(delta))
    c = R*delta
    c[-1] = c[-1] + 1

    # Z
    Z = sigma*np.sqrt((1-np.exp(-2*a*yearfrac(settle,np.array([expiry]),3)))/(2*a))+ sigma*np.sqrt(yearfrac(settle,np.array([expiry]),3))*(a==0)

    # v
    v = Z*(1-np.exp(-a*yearfrac(expiry,np.hstack((expiry, payDates)),3)))/a*(a>0) + Z*yearfrac(expiry,np.hstack((expiry, payDates)),3)*(a==0)

    # z
    z = (1-gamma)*v

    # nu
    nu = v[:-1] - gamma*v[1:]

    # forward discount factors B(t0,expiry,Ti)
    B_fwd = B/B[0]

    # function f(x) that has to be put equal to zero
    f = lambda x: sum(c*B_fwd[1:]*np.exp(-z[1:]*x-z[1:]**2/2))+\
                  sum(B_fwd[1:-1]*np.exp(-z[1:-1]*x-z[1:-1]**2/2))-\
                  sum(beta*B_fwd[:-1]*np.exp(-nu*x-nu**2/2))

    # unique value of x such that f(x) = 0
    x_star = fsolve(f,-1)

    # swaption Receiver price
    SR = B[0]*(sum(c*B_fwd[1:]*norm.cdf(x_star+z[1:])) +\
               sum(B_fwd[1:-1]*norm.cdf(x_star+z[1:-1])) -\
               sum(beta*B_fwd[:-1]*norm.cdf(x_star+nu)))

    # Put-Call parity
    BPV = sum(delta*B_fwd[1:])
    swap_rate = (1 - B_fwd[-1] + sum(B_fwd[:-1]*(beta-1)))/BPV

    # swaption Payer price
    SP = SR + B[0]*BPV*(swap_rate - r)

    return SP

def SwaptionPricingModel_Receiver(gamma,r,B,beta,a,sigma,settle,expiry,payDates):
    """
    Function that computes price of the single Receiver swaption by using MHW model

    INPUT
    gamma:             fixed value of gamma
    r:                 strike, namely the fixed rate
    B:                 discount factors in each payment date
    beta:              vector of beta between two consecutive payment dates
    a:                 calibrated parameter a of the model
    sigma:             calibrated parameter sigma of the model
    time_interval_Z:   time to swaption expiry
    deltas:            deltas between two consecutive payment dates, but starting from the expiry
    time_intervals_v:  time intervals between expiry and payment dates (with also expiry)

    OUTPUT
    SR:                price of the Receiver swaption
    """

    # deltas between two consecutive payment dates, but starting from the expiry
    delta = yearfrac(np.hstack((expiry, payDates[:-1])),payDates,2)

    # Coupon vector
    R = r*np.ones(np.shape(delta))
    c = R*delta
    c[-1] = c[-1] + 1

    # Z
    Z = sigma*np.sqrt((1-np.exp(-2*a*yearfrac(settle,np.array([expiry]),3)))/(2*a))+ sigma*np.sqrt(yearfrac(settle,np.array([expiry]),3))*(a==0)

    # v
    v = Z*(1-np.exp(-a*yearfrac(expiry,np.hstack((expiry, payDates)),3)))/a*(a>0) + Z*yearfrac(expiry,np.hstack((expiry, payDates)),3)*(a==0)

    # z
    z = (1-gamma)*v

    # nu
    nu = v[:-1] - gamma*v[1:]

    # forward discount factors B(t0,expiry,Ti)
    B_fwd = B/B[0]

    # function f(x) that has to be put equal to zero
    f = lambda x: sum(c*B_fwd[1:]*np.exp(-z[1:]*x-z[1:]**2/2))+\
                  sum(B_fwd[1:-1]*np.exp(-z[1:-1]*x-z[1:-1]**2/2))-\
                  sum(beta*B_fwd[:-1]*np.exp(-nu*x-nu**2/2))

    # unique value of x such that f(x) = 0
    x_star = fsolve(f,-1)

    # swaption Receiver price
    SR = B[0]*(sum(c*B_fwd[1:]*norm.cdf(x_star+z[1:])) +\
               sum(B_fwd[1:-1]*norm.cdf(x_star+z[1:-1])) -\
               sum(beta*B_fwd[:-1]*norm.cdf(x_star+nu)))

    return SR

def ModelCalibration_Payer(gamma,r,discount_dates,discount_curve,P_dates,P_curve,settle,expiries,tenors,freq,vol,init_cond):
    """
    Function that calibrates parameters a and sigma for the selected value of gamma by using Payer swaptions

    INPUT
    gamma:           selected value of gamma
    r:               strike, namely the fixed rate
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    P_dates:         dates of the financial instruments used for the pseudo-curve construction
    P_curve:         pseudo-discount factors curve
    settle:          settlement date
    expiries:        vector of expiries
    tenors:          vector of tenors
    freq:            annualy payment frequency
    vol:             vector of implied volatilities for the considered swaptions
    init_cond:       vector of initial conditions for both a and sigma

    OUTPUT
    output_vector:   a and sigma calibrated for the model, minimized distance, model and market prices, errors between market prices and model prices
    """

    # expiries in terms of actual dates
    expiriesDates = change_weekends(addtodate(settle,expiries,'year'))

    # payment dates in terms of actual dates
    createDates = np.array([12/freq*i for i in range(1,freq*max(expiries+tenors)+1)])
    payDates = change_weekends(addtodate(settle,createDates,'month'))

    # discounts, pseudo-discounts and betas
    B = interp_ZR(settle,discount_dates,discount_curve,payDates)    # discount factors in each payment date
    B_fwd = B[1:]/B[:-1]                                            # forward discount factors
    BP = interp_ZR(settle,P_dates,P_curve,payDates)                 # pseudo-discount factors in each payment date
    BP_fwd = BP[1:]/BP[:-1]                                         # forward pseudo-discount factors
    beta = B_fwd/BP_fwd                                             # beta between two payment dates

    # market prices by using implied volatilies
    MktPrice = lambda i: SwaptionPayer_calib(B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                            payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq],expiriesDates[i],settle,r,vol[i])
    Mkt_Prices = np.array([MktPrice(i) for i in range(len(tenors))]).squeeze()

    # model price of a single swaption
    def fun(i,parameters):
        return SwaptionPricingModel_Payer(gamma,r,B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                parameters[0],parameters[1],settle,expiriesDates[i],payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq])

    # model prices of all the considered swaptions
    ModelPrices = lambda parameters: np.array([fun(i,parameters) for i in range(len(tenors))])

    # function to minimize
    distance = lambda parameters: np.sum((Mkt_Prices-ModelPrices(parameters))**2)

    # parameters must be positive (linear constraints)
    linear_constraint = LinearConstraint([[1,0],[0,1]], [0, 0], [np.inf, np.inf])

    # parameter constrained optimization
    res = minimize(distance, init_cond, constraints=linear_constraint)#bounds=[(0,np.inf),(0,np.inf)])

    # optimal parameters a and sigma
    params = res.x

    # model prices computed with the obtained parameters a and sigma
    Model_Prices = ModelPrices(params)

    # difference between model prices and market prices
    err = abs(Model_Prices-Mkt_Prices)

    output_vector = {}
    output_vector["Params"] = res.x
    output_vector["dist"] = res.fun
    output_vector["model"] = Model_Prices
    output_vector["mkt"] = Mkt_Prices
    output_vector["err"] = err

    return output_vector

def ModelCalibration_Receiver(gamma,r,discount_dates,discount_curve,P_dates,P_curve,settle,expiries,tenors,freq,vol,init_cond):
    """
    Function that calibrates parameters a and sigma for the selected value of gamma by using Receiver swaptions

    INPUT
    gamma:           selected value of gamma
    r:               strike, namely the fixed rate
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    P_dates:         dates of the financial instruments used for the pseudo-curve construction
    P_curve:         pseudo-discount factors curve
    settle:          settlement date
    expiries:        vector of expiries
    tenors:          vector of tenors
    freq:            annualy payment frequency
    vol:             vector of implied volatilities for the considered swaptions
    init_cond:       vector of initial conditions for both a and sigma

    OUTPUT
    output_vector:   a and sigma calibrated for the model, minimized distance, model and market prices, errors between market prices and model prices

    """

    # expiries in terms of actual dates
    expiriesDates = change_weekends(addtodate(settle,expiries,'year'))

    # payment dates in terms of actual dates
    createDates = np.array([12/freq*i for i in range(1,freq*max(expiries+tenors)+1)])
    payDates = change_weekends(addtodate(settle,createDates,'month'))

    # discounts, pseudo-discounts and betas
    B = interp_ZR(settle,discount_dates,discount_curve,payDates)    # discount factors in each payment date
    B_fwd = B[1:]/B[:-1]                                            # forward discount factors
    BP = interp_ZR(settle,P_dates,P_curve,payDates)                 # pseudo-discount factors in each payment date
    BP_fwd = BP[1:]/BP[:-1]                                         # forward pseudo-discount factors
    beta = B_fwd/BP_fwd                                             # beta between two payment dates

    # market prices by using implied volatilies
    MktPrice = lambda i: SwaptionReceiver(B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                            payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq],expiriesDates[i],settle,r,vol[i])
    Mkt_Prices = np.array([MktPrice(i) for i in range(len(tenors))]).squeeze()

    # model price of a single swaption
    def fun(i,parameters):
        return SwaptionPricingModel_Receiver(gamma,r,B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                parameters[0],parameters[1],settle,expiriesDates[i],payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq])

    # model prices of all the considered swaptions
    ModelPrices = lambda parameters: np.array([fun(i,parameters) for i in range(len(tenors))])

    # function to minimize
    distance = lambda parameters: np.sum((Mkt_Prices-ModelPrices(parameters))**2)

    # parameters must be positive (linear constraints)
    linear_constraint = LinearConstraint([[1,0],[0,1]], [0, 0], [np.inf, np.inf])

    # parameter constrained optimization
    res = minimize(distance, init_cond,constraints=linear_constraint)#bounds=[(0,np.inf),(0,np.inf)])

    # optimal parameters a and sigma
    params = res.x

    # model prices computed with the obtained parameters a and sigma
    Model_Prices = ModelPrices(params)

    # difference between model prices and market prices
    err = abs(Model_Prices-Mkt_Prices)

    output_vector = {}
    output_vector["Params"] = res.x
    output_vector["dist"] = res.fun
    output_vector["model"] = Model_Prices
    output_vector["mkt"] = Mkt_Prices
    output_vector["err"] = err

    return output_vector

def SwaptionPricingModel_PayerCVA(gamma,r,B,beta,a,sigma,settle,expiry,payDates,Notionals,Initial_Notional,flag):
    """
    Function that computes price of the single Payer swaption by using the calibrated model

    INPUT
    gamma:             fixed value of gamma
    r:                 strike, namely the fixed rate
    B:                 discount factors in each payment date
    beta:              vector of beta between two consecutive payment dates
    a:                 calibrated parameter a of the model
    sigma:             calibrated parameter sigma of the model
    settle:            settlement date
    expiry:            expiry date Ti
    payDates:          payment dates of the underlying swap
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    SP:                price of the Payer swaption

    """

    # deltas between two consecutive payment dates, but starting from the expiry
    delta = yearfrac(np.hstack((expiry, payDates[:-1])),payDates,2)

    # Coupon vector
    R = r*np.ones(np.shape(delta))
    c = R*delta

    # Z
    Z = sigma*np.sqrt((1-np.exp(-2*a*yearfrac(settle,np.array([expiry]),3)))/(2*a))*(a>0) +\
        sigma*np.sqrt(yearfrac(settle,np.array([expiry]),3))*(a==0)

    # v
    v = Z*(1-np.exp(-a*yearfrac(expiry,np.hstack((expiry, payDates)),3)))/a*(a>0) +\
        Z*yearfrac(expiry,np.hstack((expiry, payDates)),3)*(a==0)

    # z
    z = (1-gamma)*v

    # nu
    nu = v[:-1] - gamma*v[1:]

    # forward discount factors B(0,expiry,Ti)
    B_fwd = B/B[0]
    BPV = sum(delta*B_fwd[1:]*Notionals)

    if flag == "SA":
        # take in account the notionals as ratios
        c = c*Notionals/Initial_Notional

        # updated function f(x) with different notionals
        f = lambda x: sum(c * B_fwd[1:] * np.exp(-z[1:] * x - z[1:] ** 2 / 2)) + \
                      sum(B_fwd[1:] * np.exp(-z[1:] * x - z[1:] ** 2 / 2)) - \
                      sum(beta * B_fwd[:-1] * np.exp(-nu * x - nu ** 2 / 2))

        # unique value of x such that f(x) = 0
        x_star = fsolve(f, -1)

        # swaption Receiver price
        SR = Initial_Notional * B[0] * (sum(c * B_fwd[1:] * norm.cdf(x_star + z[1:])) + \
                                        sum(B_fwd[1:] * norm.cdf(x_star + z[1:])) - \
                                        sum(beta * B_fwd[:-1] * norm.cdf(x_star + nu)))

        swap_rate = Initial_Notional*(1 - B_fwd[-1] + sum(B_fwd[:-1]*(beta-1)))/BPV

    elif flag == "DA":

        #do not consider ratios since Notionals aren't constant in any leg
        c = c*Notionals

        # updated function f(x) with different notionals
        f = lambda x: sum(c * B_fwd[1:] * np.exp(-z[1:] * x - z[1:] ** 2 / 2)) + \
                  sum(Notionals*B_fwd[1:] * np.exp(-z[1:] * x - z[1:] ** 2 / 2)) - \
                  sum(Notionals*beta * B_fwd[:-1] * np.exp(-nu * x - nu ** 2 / 2))

        # unique value of x such that f(x) = 0
        x_star = fsolve(f,-1)

        # swaption Receiver price, updated due to different notionals in single amortizing
        SR = B[0] * (sum(c * B_fwd[1:] * norm.cdf(x_star + z[1:])) + \
                     sum(Notionals * B_fwd[1:] * norm.cdf(x_star + z[1:])) - \
                     sum(Notionals * beta * B_fwd[:-1] * norm.cdf(x_star + nu)))


        # forward discount factors B(0,Ti,Ti+1)
        B_fwd_i = B[1:]/B[:-1]
        swap_rate = sum(B_fwd[:-1]*(beta - B_fwd_i)*Notionals)/BPV

    # swaption Payer price thanks to Put-Call parity
    SP = SR + B[0]*BPV*(swap_rate - r)

    return SP


def swap_pricing_CVA_Payer_Model(Prob,beta,NPV_risk_free,settle,payment_dates,B,Notionals,Initial_Notional,K,LGD,a,sigma,flag):
    """Function that computes NPV of the swap considering CVA, by using the calibrated parameters of the Model

    INPUT
    Prob:              matrix of survival probabilities (each column for a different lambda)
    beta:              vector of betas
    NPV_risk_free:     Net Present Value without any risk assumption
    settle:            settlement date
    payment_dates:     payment dates of the swap contract
    B:                 discount factors in each payment date
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    K:                 strike, namely the fixed rate
    LGD:               Loss Given Default
    a:                 calibrated parameter a of the model
    sigma:             calibrated parameter sigma of the model
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    NPV_CVA:           Net Present Values considering CVA with different lambdas
    """

    # gamma in our setting
    gamma = 0

    # prices of the payer swaptions needed for the CVA
    SwaptionPricing = lambda i: SwaptionPricingModel_PayerCVA(gamma,K,B[i:],beta[i:],a,sigma,settle,payment_dates[i],
                                                              payment_dates[i+1:],Notionals[i+1:],Initial_Notional,flag)
    PS = np.array([SwaptionPricing(i) for i in range(len(payment_dates)-1)])

    # CVAs computation
    compute_CVA = lambda i: LGD * sum(PS * (Prob[i][:-1] - Prob[i][1:]))
    CVA = np.hstack([compute_CVA(i) for i in range(Prob.shape[0])])

    # final NPVs
    NPV_CVA = NPV_risk_free - CVA

    return NPV_CVA


def surf_plot(gamma,r,discount_dates,discount_curve,P_dates,P_curve,settle,expiries,tenors,freq,vol,init_cond):
    """
    Function that computes the vector Z of distances as function of the parameters a and sigma in the Payer calibration case,
    used for the surface plot

    INPUT
    gamma:           selected value of gamma
    r:               strike, namely the fixed rate
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    P_dates:         dates of the financial instruments used for the pseudo-curve construction
    P_curve:         pseudo-discount factors curve
    settle:          settlement date
    expiries:        vector of expiries
    tenors:          vector of tenors
    freq:            annualy payment frequency
    vol:             vector of implied volatilities for the considered swaptions
    init_cond:       vector of initial conditions for both a and sigma

    OUTPUT
    X:               vector of parameters a
    Y:               vector of parameters sigma
    Z:               vector Z of distances as function of the parameters a and sigma
    """

    # expiries in terms of actual dates
    expiriesDates = change_weekends(addtodate(settle,expiries,'year'))

    # payment dates in terms of actual dates
    createDates = np.array([12/freq*i for i in range(1,freq*max(expiries+tenors)+1)])
    payDates = change_weekends(addtodate(settle,createDates,'month'))

    # discounts, pseudo-discounts and betas
    B = interp_ZR(settle,discount_dates,discount_curve,payDates)    # discount factors in each payment date
    B_fwd = B[1:]/B[:-1]                                            # forward discount factors
    BP = interp_ZR(settle,P_dates,P_curve,payDates)                 # pseudo-discount factors in each payment date
    BP_fwd = BP[1:]/BP[:-1]                                         # forward pseudo-discount factors
    beta = B_fwd/BP_fwd                                             # beta between two payment dates

    # market prices by using implied volatilies
    MktPrice = lambda i: SwaptionPayer_calib(B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                            payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq],expiriesDates[i],settle,r,vol[i])
    Mkt_Prices = np.array([MktPrice(i) for i in range(len(tenors))]).squeeze()

    # model price of a single swaption
    def fun(i,parameters):
        return SwaptionPricingModel_Payer(gamma,r,B[expiries[i]*freq-1:(expiries[i]+tenors[i])*freq],beta[freq*expiries[i]-1:(expiries[i]+tenors[i])*freq-1],
                parameters[0],parameters[1],settle,expiriesDates[i],payDates[expiries[i]*freq:(expiries[i]+tenors[i])*freq])

    # model prices of all the considered swaptions
    ModelPrices = lambda parameters: np.array([fun(i,parameters) for i in range(len(tenors))])

    # function to minimize
    distance = lambda parameters: np.sum((Mkt_Prices-ModelPrices(parameters))**2)

    x = np.linspace(0.00001, 0.002, 40)
    y = np.linspace(0.006, 0.012, 40)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    for i in range(len(x)):
        for j in range(len(y)):
            Z[i, j] = distance([X[i, j], Y[i, j]])

    return X, Y, Z