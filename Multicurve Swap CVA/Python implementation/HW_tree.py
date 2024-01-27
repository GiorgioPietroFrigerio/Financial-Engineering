import numpy as np
from utilities import *


def swap_pricing_tree(Prob,beta,NPV_risk_free,swap_maturity,settle,payment_dates,Notionals,Initial_Notional,K,a,sigma,discount_dates,discount_curve,freq,month_steps,LGD,flag):
    """
    Function that computes NPV of the swap considering CVA, by using trinomial tree

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
    month_steps:       number of steps in a month
    LGD:               Loss Given Default
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    NPV_CVA:           Net Present Values considering CVA with different lambdas
    """

    # expiries in months
    expiries_in_months = np.arange(1/freq,swap_maturity,1/freq)*12

    # time step
    dt_monthly = 1 / month_steps
    dt = dt_monthly / 12

    # Tree construction
    sigma_tree = sigma*np.sqrt((1-np.exp(-2*a*dt))/(2*a))   # sigma_tree
    mu_tree = 1-np.exp(-a*dt)                               # mu_tree
    delta_x = np.sqrt(3)*sigma_tree                         # delta_x
    l_max = int(np.ceil((1-np.sqrt(2/3))/mu_tree))          # l_max
    l_min = -l_max                                          # l_min
    l = np.array([i for i in range(l_max,l_min-1,-1)])      # vector of possible l
    x = l*delta_x                                           # vector of possible x in the tree


    # transition probabilities in the cases A B and C
    prob_A = np.vstack((1/2*(1/3-l*mu_tree+(l*mu_tree)**2), 2/3-(l*mu_tree)**2, 1/2*(1/3+l*mu_tree+(l*mu_tree)**2)))                                         # transition probabilities for Case A
    prob_B = np.hstack((1/2*(7/3-3*l_max*mu_tree+(l_max*mu_tree)**2), -1/3+2*l_max*mu_tree-(l_max*mu_tree)**2, 1/2*(1/3-l_max*mu_tree+(l_max*mu_tree)**2)))  # transition probabilities for Case B
    prob_C = np.hstack((1/2*(1/3+l_min*mu_tree+(l_min*mu_tree)**2), -1/3-2*l_min*mu_tree-(l_min*mu_tree)**2, 1/2*(7/3+3*l_min*mu_tree+(l_min*mu_tree)**2)))  # transition probabilities for Case C

    # probability matrix
    prob_matrix = prob_A
    prob_matrix[:,-1] = prob_B
    prob_matrix[:,0] = prob_C

    # sigma hat star
    sigma_hat_star = sigma/a*np.sqrt(dt-2*(1-np.exp(-a*dt))/a +(1-np.exp(-2*a*dt))/(2*a))

    # auxiliary vectors to create the proper discount D(ti,ti+1)
    aux_vec = np.exp(-0.5*sigma_hat_star**2-sigma_hat_star/sigma_tree*mu_tree*x)                    # first vector of exponentials
    vec_delta_x = np.exp(-sigma_hat_star/sigma_tree*np.exp(-a*dt)*delta_x*np.array([2,1,0,-1,-2]))  #second vector of exponentials (with all possible delta x)

    # prices of the payer swaptions needed for the CVA
    SwaptionPricing = lambda i: swaption_price_tree(settle,discount_dates,discount_curve,payment_dates[i],expiries_in_months[i],
                                            payment_dates[-1],sigma,a,month_steps,K,Notionals[i+1:],Initial_Notional,beta[i:],freq,
                                            x,prob_matrix,aux_vec,vec_delta_x,dt_monthly,l_max,flag)
    PS = np.array([SwaptionPricing(i) for i in range(len(payment_dates)-1)])

    # CVAs computation
    compute_CVA = lambda i: LGD * sum(PS * (Prob[i][:-1] - Prob[i][1:]))
    CVA = np.hstack([compute_CVA(i) for i in range(Prob.shape[0])])

    # final NPVs
    NPV_CVA = NPV_risk_free - CVA

    return NPV_CVA

def swaption_price_tree(settle,dates,discounts,expiry_date,expiry,maturity,sigma,a,month_steps,k,Notionals,
                        Initial_Notional,beta,freq,x,prob_matrix,aux_vec,vec_delta_x,dt,l_max,flag):
    """
    Function that computes the price of a swaption payer with trinomial tree

    INPUT
    settle:            settlement date
    dates:             dates from bootstrap
    discounts:         discounts from bootstrap
    expiry_date:       expiry date of the swaption
    expiry:            expiry of the swaption in number of months
    maturity:          maturity date of the swap
    sigma:             calibrated sigma for HULL-WHITE model
    a:                 calibrated a for HULL-WHITE model
    month_steps:       number of monthly steps
    k:                 strike, namely the fixed rate
    Notionals:         vector of amortized notionals
    Initial_Notional:  initial notional
    beta:              beta between two consecutive payment dates
    freq:              annual payments frequency
    x:                 vector of OU process values
    prob_matrix:       matrix of transition probabilities
    aux_vec:           first vector of exponentials to obtain D(ti,ti+1)
    vec_delta_x:       second vector of exponentials (with all possible delta x) to obtain D(ti,ti+1)
    dt:                monthly time step
    l_max:             l max
    flag:              "SA" -> single amortizing only on the fixed leg, "DA" -> double amortizing also on the floating leg

    OUTPUT
    swaption_price:    price of the swaption """

    # tree parameters
    Dates_tree = dates_for_tree(settle,expiry_date,expiry,dt,month_steps)                   # dates in the tree
    dt = dt/12                                                                              # yearly time interval
    deltas_365_tree = yearfrac(settle, np.hstack((settle,Dates_tree)),3)                    # time intervals (act/365)
    disc_fact_tree = interp_ZR(settle,dates,discounts,np.hstack((settle,Dates_tree)) )      # discount factors in dates_tree
    forward_disc_tree = disc_fact_tree[1:]/disc_fact_tree[:-1]                              # forward discount factors in dates_tree

    # swap parameters
    Start_date = expiry_date                                                                # start date of the swap (expiry of swaption)
    PaymentDates = cfdates(Start_date,maturity+relativedelta(days=-1),freq)                 # payment dates of the coupons (swap)
    deltas_365 = yearfrac(settle,np.hstack((Start_date,PaymentDates)),3)                    # time intervals (act/365)
    deltas_360 = yearfrac(np.hstack((Start_date,PaymentDates[:-1])),PaymentDates,2)    # time intervals (act/360)
    disc_fact = interp_ZR(settle,dates,discounts,np.hstack((Start_date,PaymentDates)))      # discount factors at payment dates
    forw_disc_fact = disc_fact/disc_fact[0]                                                 # forward discount factors

    # Number of steps of the tree
    N = month_steps*int(expiry)+1

    # volatility of x at the expiry of the swaption
    Z = np.sqrt(sigma**2*(1-np.exp(-2*a*deltas_365[0]))/(2*a))

    if flag == "SA":
        #coupons for the CB
        coupons = k*deltas_360*Notionals/Initial_Notional
        coupons[-1] = coupons[-1] + 1

        # coupond bond in the strike of the put option
        CB_strike, _ = MHW_price(x/Z,a,sigma,beta[1:]-1,forw_disc_fact[1:-1],deltas_365[:-1],expiry,freq)

        # put option strike
        option_strike = beta[0] + CB_strike

    elif flag == "DA":

        # coupons for the CB
        coupons = k*deltas_360*Notionals
        coupons[-1] = coupons[-1] + Notionals[-1]

        # coupons for option strike
        strike_coupons = beta[1:]*Notionals[1:]-Notionals[:-1]

        # CB for option strike
        CB_strike, _ = MHW_price(x/Z,a,sigma,strike_coupons,forw_disc_fact[1:-1],deltas_365[:-1],expiry,freq)

        # put option strike
        option_strike = beta[0]*Notionals[0] + CB_strike

    # CB in T_alpha
    CB, _ = MHW_price(x/Z,a,sigma,coupons,forw_disc_fact[1:],deltas_365,expiry,freq)

    # payoffs of the swaption (payer case)
    payoff_payer = np.maximum(option_strike-CB,0).squeeze()

    # allocation of matrices of nodes
    tree_matrix_payer = np.zeros((np.shape(x)[0],N))

    # last columns of the matrices are the payoffs
    tree_matrix_payer[:,-1] = payoff_payer


    # construction of the matrix of discount factors B(ti,ti+1) depending on x
    Z = lambda i: np.sqrt(sigma**2*(1-np.exp(-2*a*dt*(i-1)))/(2*a))                      # Z
    v = lambda i: Z(i)*(1-np.exp(-a*dt))/a                                               # v
    z = v                                                                                # z (=v in gamma=0)
    vec_Z = np.array([Z(i) for i in range(1,len(deltas_365_tree))])                      # Z vector
    vec_z = np.array([z(i) for i in range(1,len(deltas_365_tree))])                      # z vector
    matrix_1 = np.outer(-x,np.hstack((0,vec_z[1:]/vec_Z[1:])))                           # first part of the exponential
    matrix_1[:,0] = 0                                                                    # first column of zeros
    matrix_2 = -0.5*vec_z**2                                                             # second part of the exponential
    matrix_exp = np.exp(matrix_1 + matrix_2)                                             # matrix with exponentials
    matrix_disc_fact = forward_disc_tree*matrix_exp                                      # matrix with discunts B(ti,ti+1)


    # backward algorithm (to discount the payoffs of the swaption)
    for j in range(N-2,-1,-1):
        # l_max
        tree_matrix_payer[0,j] = tree_matrix_payer[0,j+1]*prob_matrix[0,0]*matrix_disc_fact[0,j]*aux_vec[0]*vec_delta_x[2] + tree_matrix_payer[1,j+1]*prob_matrix[1,0]*matrix_disc_fact[0,j]*aux_vec[1]*vec_delta_x[3]\
                                 + tree_matrix_payer[2,j+1]*prob_matrix[2,0]*matrix_disc_fact[0,j]*aux_vec[2]*vec_delta_x[4]
        # l_min
        tree_matrix_payer[-1,j] = tree_matrix_payer[-3,j+1]*prob_matrix[0,-1]*matrix_disc_fact[-1,j]*aux_vec[-3]*vec_delta_x[0] + tree_matrix_payer[-1,j+1]*prob_matrix[1,-1]*matrix_disc_fact[-1,j]*aux_vec[-2]*vec_delta_x[1]\
                                + tree_matrix_payer[-1,j+1]*prob_matrix[2,-1]*matrix_disc_fact[-1,j]*aux_vec[-1]*vec_delta_x[2]
        # medium
        tree_matrix_payer[1:-1,j] = tree_matrix_payer[:-2,j+1]*prob_matrix[0,1:-1]*matrix_disc_fact[1:-1,j]*aux_vec[:-2]*vec_delta_x[1] + tree_matrix_payer[1:-1,j+1]*prob_matrix[1,1:-1]*matrix_disc_fact[1:-1,j]*aux_vec[1:-1]*vec_delta_x[2]\
                                    + tree_matrix_payer[2:,j+1]*prob_matrix[2,1:-1]*matrix_disc_fact[1:-1,j]*aux_vec[2:]*vec_delta_x[3]

    # price in t_0 of the swaption payer
    if flag == "SA":
        swaption_price = Initial_Notional*tree_matrix_payer[l_max,0]
    elif flag == "DA":
        swaption_price = tree_matrix_payer[l_max,0]

    return swaption_price

def MHW_price(x,a,sigma,coupons,forw_disc,deltas_365,expiry,freq):
    """
    Function that computes Coupon Bond price and ZCB prices using Multicurve HW model

    INPUT
    x:           value of the ORNSTEIN-UHLENBECK process
    a:           a for HULL-WHITE model
    sigma:       sigma for HULL-WHITE model
    coupons:     coupons of the CB
    forw_disc:   forward discounts
    deltas_365:  time intervals in act/365
    expiry:      expiry in months
    freq:        payment frequency

    OUTPUT
    CB_price:    price of the Coupon Bond
    ZCB_prices:  prices of the ZC Bonds
    """

    # Z in the dynamics of B
    Z = np.sqrt(sigma**2*(1-np.exp(-2*a*deltas_365[0]))/(2*a))

    # v in the dynamics of B
    v = lambda tau: Z*(1-np.exp(-a*(tau-deltas_365[0])))/a

    # z in the dynamics of B with gamma=0
    z = v

    # ZCB price function
    fun_ZCB_price = lambda tau: forw_disc[int(np.round(freq*tau)-expiry*freq/12-1)]*np.exp(-z(tau)*x - 0.5*z(tau)**2)

    # ZCB prices
    ZCB_prices = np.array([fun_ZCB_price(i) for i in deltas_365[1:]])

    # CB price
    CB_price = (ZCB_prices.T@coupons[:,None]).flatten()

    return CB_price, ZCB_prices.flatten()

