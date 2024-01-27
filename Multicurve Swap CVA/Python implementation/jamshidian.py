import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm
from HW_tree import MHW_price
from utilities import *

def swap_pricing_jamshidian(Prob,beta,NPV_risk_free,swap_maturity,settle,payment_dates,Notionals,Initial_Notional,K,a,
                            sigma,discount_dates,discount_curve,freq,LGD,flag):
    """ Function that computes NPV of the swap considering CVA, by using Jamshidian approach

    INPUT
    Prob:              matrix of survival probabilities (each column for a different lambda)
    beta:              vector of betas
    NPV_risk_free:     Net Present Value without any risk assumption
    swap_maturity:     maturity of the swap in years
    settle:            settlement date
    payment_dates:     payment dates of the swap contract
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    K:                 strike, namely the fixed rate
    a:                 calibrated parameter a of the model
    sigma:             calibrated parameter sigma of the model
    discount_dates:    dates of the financial instruments used for the curve construction
    discount_curve:    discount factors curve
    freq:              annual payments frequency
    LGD:               Loss Given Default
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    NPV_CVA:           Net Present Values considering CVA with different lambdas
    """

    # expiries and tenors in months
    expiries_in_months = np.arange(1 / freq, swap_maturity, 1 / freq) * 12
    tenors_in_months = np.arange((swap_maturity - 1/freq),0,-1 / freq) * 12

    # prices of the payer swaptions needed for the CVA
    SwaptionPricing = lambda x: swaption_price_Jam(settle,discount_dates,discount_curve,payment_dates[x],expiries_in_months[x],
                            tenors_in_months[x],payment_dates[-1],sigma,a,K,Notionals[x+1:],Initial_Notional,beta[x:],freq,flag)
    PS = np.array([SwaptionPricing(i) for i in range(len(payment_dates) - 1)])

    # CVAs computation
    compute_CVA = lambda i: LGD * sum(PS * (Prob[i][:-1] - Prob[i][1:]))
    CVA = np.hstack([compute_CVA(i) for i in range(Prob.shape[0])])

    # final NPVs
    NPV_CVA = NPV_risk_free - CVA

    return NPV_CVA

def swaption_price_Jam(settle,dates,discounts,expiry_date,expiry_in_months,tenor_in_months,maturity,sigma,a,k,Notionals,Initial_Notional,beta,freq,flag):
    """
    Function that computes the price of a swaption payer using the Jamshidian formula

    INPUT
    settle:             settlement date
    dates:              dates from bootstrap
    discounts:          vector of discount factors from bootstrap
    expiry_date:        expiry date of the swaption
    expiry_in_months:   expiry in months of the swaption
    tenor_in_months:    tenor of the swaption in months
    sigma:              calibrated sigma for HULL-WHITE model
    a:                  calibrated a for HULL-WHITE model
    k:                  strike, namely the fixed rate
    Notionals:          vector of amortized notionals
    Initial_Notional:   initial notional
    beta:               beta between two consecutive payment dates
    freq:               annual payments frequency
    flag:               "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    price_put:          price of the swaption payer
    """

    # date corresponding to the expiry of the swaption
    Start_date = expiry_date

    # coupon payment dates
    PaymentDates = cfdates(Start_date, maturity + relativedelta(days=-1), freq)

    # time intervals (act/365)
    deltas_365 = yearfrac(settle,np.hstack((Start_date, PaymentDates)),3)

    # time intervals (act/360)
    deltas_360 = yearfrac(np.hstack((Start_date, PaymentDates[:-1])), PaymentDates,2)

    # discount factors at coupon payment dates
    disc_fact = interp_ZR(settle,dates,discounts,np.hstack((Start_date, PaymentDates)))

    # forward discount factors
    forw_disc_fact = disc_fact/disc_fact[0]

    if flag == "SA":

        # Computing the coupons and adding the face value for the last one
        coupons = np.hstack((k*deltas_360[:-1]*Notionals[:-1]/Initial_Notional + 1 - beta[1:], k*deltas_360[-1]*Notionals[-1]/Initial_Notional + 1))

        # option strike
        option_strike = beta[0]

    elif flag == "DA":

        # Computing the coupons
        coupons = np.hstack(((1 + k*deltas_360[:-1])*Notionals[:-1] - beta[1:]*Notionals[1:],(1 + k*deltas_360[-1])*Notionals[-1]))

        # option strike
        option_strike = beta[0]*Notionals[0]


    # price of the bond - strike (function to be put = 0)
    P = lambda x: MHW_price(x,a,sigma,coupons,forw_disc_fact[1:],deltas_365,expiry_in_months,freq)[0].astype(np.float64) - option_strike

    # finding x_star
    x_star = fsolve(P,-1.0)

    # finding the K_i
    K_i = MHW_price(x_star, a, sigma, coupons, forw_disc_fact[1:], deltas_365, expiry_in_months, freq)[1]

    # sigmas for black formula
    I = np.array([i for i in range(int(tenor_in_months*freq/12))])
    Z = np.sqrt(sigma**2*(1-np.exp(-2*a*deltas_365[0]))/(2*a))
    v = Z*(1-np.exp(-a*(deltas_365[I]-deltas_365[0])))/a
    z = v

    # prices of options on ZC bonds
    puts = disc_fact[0]*K_i*(1-norm.cdf(x_star)) - disc_fact[1:]*(1-norm.cdf(x_star+z))

    # swaption price (payer case)
    if flag == "SA":
        price_put = Initial_Notional*sum(puts*coupons)
    elif flag == "DA":
        price_put = sum(puts*coupons)

    return price_put


