import numpy as np
from scipy.interpolate import interp1d
from utilities import *
from swaption import *


def swap_pricing_risk_free(start_dates,end_dates,B,Bp_fwd,Notionals,Initial_Notional,K,flag):
    """
    Function that computes swap price without any risk assumption

    INPUT
    start_dates:       start date of each payment period
    end_dates:         end date of each payment period
    B:                 discount factors in each payment date
    Bp_fwd:            forward pseudo-discount factors between two consecutive payment dates
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    K:                 strike, namely the fixed rate
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    NPV_risk_free:     NPV of the swap as the difference between between floating and fixed leg
    """

    # day count for swap
    swap_daycount = 2 # act/360

    # time intervals
    deltas = yearfrac(start_dates, end_dates, swap_daycount).squeeze()

    # floating NPV
    if flag == "SA":
        NPVA = sum(Initial_Notional * B * (1 / Bp_fwd - 1))   # only fixed leg amortizing
    elif flag == "DA":
        NPVA = sum(B * Notionals.squeeze() * (1 / Bp_fwd - 1))         # both leg amortizing leg

    # fixed NPV
    NPVB = K * sum(B * deltas * Notionals.squeeze())

    # swap NPV risk free
    NPV_risk_free = NPVA - NPVB

    return NPV_risk_free

def swap_pricing_CVA(Prob,beta,deltas,NPV_risk_free,vol,settle,payment_dates,B,Notionals,Initial_Notional,K,LGD,flag):
    """
    Function that computes swap pricing taking into account also counterparty risk

    INPUT
    Prob:              matrix of survival probabilities (each column for a different lambda)
    beta:              vector of betas
    deltas:            time intervals between payment dates
    NPV_risk_free:     Net Present Value without any risk assumption
    vol:               vector of interpolated implied volatilities for the amortized case
    settle:            settlement date
    payment_dates:     payment dates of the swap contract
    B:                 discount factors in each payment date
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    K:                 strike, namely the fixed rate
    LGD:               Loss Given Default
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    NPV_CVA:           Net Present Values considering CVA with different lambdas
    """

    # compute price of each swaption on the residual value of the underlying swap, moving forward on the expiry of the swaption
    SwaptionPricing = lambda x: SwaptionPayer(B[x:],beta[x:],deltas[x:],Notionals[x:],Initial_Notional,payment_dates[x],settle,K,vol[x],flag)
    PS = np.array([SwaptionPricing(i) for i in range(len(payment_dates)-1)]).squeeze()

    # CVAs computation
    compute_CVA = lambda i: LGD * sum(PS * (Prob[i][:-1] - Prob[i][1:]))
    CVA = np.hstack([compute_CVA(i) for i in range(Prob.shape[0])])

    # final NPVs
    NPV_CVA = NPV_risk_free - CVA

    return NPV_CVA

def swap_discounts(settle,discount_dates,discount_curve,P_dates,P_curve,Dates):
    """
    Function that computes discount factors and pseudo-discount factors in each swap payment date

    INPUT
    settle:          settlement date
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    P_dates:         dates of the financial instruments used for the pseudo-curve construction
    P_curve:         pseudo-discount factors curve
    Dates:           payment dates of the swap contract

    OUTPUT
    B:               discount factors in each payment date
    Bp:              pseudo-discount factors in each payment date
    Bp_fwd:          forward pseudo-discount factors between two consecutive payment dates
    """

    # Discount factors for swap
    B = interp_ZR(settle, discount_dates, discount_curve, Dates)

    # pseudo-discounts for swap
    Bp = interp_ZR(settle, P_dates, P_curve, Dates)

    # forward pseudo-discounts for swap
    Bp_fwd = np.hstack((Bp[0],Bp[1:]/Bp[:-1]))

    return B,Bp,Bp_fwd

def swap_discounts_modified(settle,discount_dates,discount_curve,P_dates,P_curve,StartDates,EndDates,Euribor,idx):
    """
    Function that computes discount factors and pseudo-discount factors in each swap payment date;
    must be performed separately due to the presence of a Euribor settled on
    28 December, paid on 28 March and since today is the 31 January

    INPUT
    settle:          settlement date
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    P_dates:         dates of the financial instruments used for the pseudo-curve construction
    P_curve:         pseudo-discount factors curve
    StartDates:      start dates of each swap payment period
    EndDates:        end dates of each swap payment period
    Euribor:         Euribor established on 28 December and paid on 28 March, taken from web market value
    idx:             vector of indexes of payment dates following the new settlement, namely 31 January

    OUTPUT
    B:               discount factors in each payment date
    Bp:              pseudo-discount factors in each payment date
    Bp_fwd:          forward pseudo-discount factors between two consecutive payment dates
    """

    swap_daycount = 2   # act/360

    # discount factors for swap
    B = interp_ZR(settle, discount_dates, discount_curve, EndDates[idx])

    # pseudo-discounts for swap
    Bp = interp_ZR(settle,P_dates,P_curve,EndDates[idx])

    # first time interval
    delta_1 = yearfrac(np.array([StartDates[idx[0]]]),np.array([EndDates[idx[0]]]),swap_daycount)

    # pseudo-discount in first time interval
    Bp_1 = 1/(1+delta_1*Euribor)

    # forward pseudo-discounts for swap
    Bp_fwd = np.hstack((Bp_1, Bp[1:]/Bp[:-1]))

    return B, Bp, Bp_fwd

def swap_elements(Bp,B,payment_dates,settle,spreads,LGD):
    """
    Function that computes probabilities, betas and time intervals

    INPUT
    Bp:             pseudo-discount factors
    B:              discount factors
    payment_dates:  payment dates
    settle:         settlement date
    spreads:        spreads
    LGD:            loss given default

    OUTPUT
    Prob:           probabilities
    beta:           betas
    deltas:         time intervals
    """

    lambdas = spreads / LGD                 # vector of intensities lambdas
    Bp_fwd = Bp[1:]/Bp[:-1]                 # forward pseudo-discount factors
    B_fwd = B[1:]/B[:-1]                    # forward discount factors
    beta = B_fwd/Bp_fwd                     # beta

    # time intervals
    deltas = yearfrac(payment_dates[:-1], payment_dates[1:], 2)

    # survival probabilities for each lambda
    Prob_daycount = 3                       # act/365
    Probability = lambda lam: np.exp(- lam*yearfrac(settle,np.hstack((settle, payment_dates[:-1])), Prob_daycount))

    Prob = np.vstack([Probability(l) for l in lambdas])

    return Prob, beta, deltas

def vol_interp(settle,Notionals,Initial_Notional,discount_dates,discount_curve,B,payment_dates,expiries,tenors,freq,swap_maturity,vol_matrix,flag):
    """

    Function that creates a complete matrix of volatilities for each expiry and maturity,
    then interpolates in terms of amortized BPV

    INPUT
    settle:            settlement date
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    discount_dates:    dates of the financial instruments used for the curve construction
    discount_curve:    discount factors curve
    B:                 discount factors in each payment date
    payment_dates:     payment dates of the underlying swap
    expiries:          vector of expiries
    tenors:            vector of tenors
    freq:              yearly payment frequency
    swap_maturity:     maturity of the underlying swap
    vol_matrix:        matrix of market implied volatilities
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    vols:              needed volatilities for the corresponding swaptions computation
    """

    exp_val = payment_dates[:-1]                           # expiries
    ten_val = np.arange(1/freq,swap_maturity,1/freq)       # tenors

    # construction of the first interpolation
    mat_aux = np.vstack([interp1d(tenors, vol_matrix[i][:],fill_value="extrapolate")(ten_val) for i in range(vol_matrix.shape[0])])

    # creation of volatility matrix for each expiry and tenor of plain vanilla swaptions
    datenum_exp = [date.toordinal() for date in expiries]
    datenum_newexp = [date.toordinal() for date in exp_val]
    mat_aux = mat_aux.T
    vol_interp = np.vstack([interp1d(datenum_exp, mat_aux[i], fill_value="extrapolate")(datenum_newexp) for i in range(mat_aux.shape[0])]).T

    if flag == "SA":
        # needed setting for regular and amortized BPVs computation
        additive_dates = change_weekends(addtodate(payment_dates[-1],12/freq*np.array(range(1,len(payment_dates)-1)),"month"))
        payment_dates = np.hstack((payment_dates, additive_dates))
        deltas = yearfrac(payment_dates[:-1], payment_dates[1:], 2)
        B_add_dates = interp_ZR(settle, discount_dates, discount_curve, additive_dates)
        discounts = np.hstack((B,B_add_dates))
        length = len(additive_dates) + 2

        #interpolate in terms of amortized BPV

        vols = np.empty(1)
        for i in range(len(exp_val)):
            disc_fwd = discounts[i+1:(length+i)]/discounts[i]
            BPVs = np.cumsum(deltas[i:length+i-1]*disc_fwd)
            ratios = Notionals[i+1:]/Initial_Notional
            BPV_hat = sum(ratios*deltas[i:length-1]*disc_fwd[:length-i-1])
            vols = np.vstack((vols,interp1d(BPVs, vol_interp[i], fill_value="extrapolate")(BPV_hat)))
        vols = vols[1:]

    elif flag == "DA": # double amortizing case

        vols = np.diag(np.fliplr(vol_interp)) # we just need the main diagonal of the flipped volatility matrix

    return vols

