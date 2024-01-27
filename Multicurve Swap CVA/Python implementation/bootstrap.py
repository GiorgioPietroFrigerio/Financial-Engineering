import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from scipy.optimize import fsolve
from readExcelData import *
from utilities import *


def discount_curve_bootstrap(data,settle,Basis_swap):
    """
    Function that computes discount factors curve

    INPUT
    data:            struct containing dates and quoted prices of the needed financial instruments
    settle:          settlement date
    Basis_swap:      30/360 European

    OUTPUT
    discount_dates:  dates of the financial instruments used for the curve construction
    discount_curve:  discount factors curve
    data:            updated struct containing dates and quoted prices of the needed financial instruments
    """

    # Swaps with dates < 1y
    days = np.array([1,7,14,21])

    data['OIS']['dates'] = change_weekends(addtodate(settle, days, 'day'))

    months = np.array([1,2,3,4,5,6,7,8,9,10,11,12,15,18,21])

    data['OIS']['dates'] = np.concatenate((data['OIS']['dates'], change_weekends(addtodate(settle, months, 'month'))))

    years = np.array([2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,40,50])

    data['OIS']['dates'] = np.concatenate((data['OIS']['dates'], change_weekends(addtodate(settle, years, 'year'))))

    # find the index of the expiries prior to 1y
    idx_1y = np.argwhere(data['OIS']['dates'] <= change_weekends(addtodate(settle, np.array([1]), 'year'))).flatten()

    discount_curve = 1/(1 + yearfrac(settle, data['OIS']['dates'][idx_1y], Basis_swap)*data['OIS']['mktRates'][idx_1y.T])
    discount_dates = data['OIS']['dates'][idx_1y]

    # Swaps after 1 year
    #find the index of the expiries after 1 year
    idx_over1y = np.argwhere(data['OIS']['dates'] > change_weekends(addtodate(settle, np.array([1]), 'year'))).flatten()
    swap_dates_interp = np.concatenate((data['OIS']['dates'][idx_over1y[:3]],change_weekends(addtodate(settle, np.array(range(2,51)), 'year'))))
    # conversion to ordinal dates just for interpolation
    datenums_OIS = [date.toordinal() for date in data['OIS']['dates'][idx_over1y]]
    datenums = [date.toordinal() for date in swap_dates_interp]
    swap_rates_interp = interp1d(datenums_OIS, data['OIS']['mktRates'][idx_over1y], kind='cubic')(datenums)

    # discount curve from swaps from 15m to 21m
    BPV_01y = discount_curve[-1] * yearfrac(settle, np.array(discount_dates[-1]), Basis_swap)
    date_1y = discount_dates[-1]
    discount_curve = np.concatenate((discount_curve,(1 - swap_rates_interp[:3] * BPV_01y)/(1 + swap_rates_interp[:3] *
                            yearfrac(date_1y, data['OIS']['dates'][idx_over1y[:3]], Basis_swap))))

    # discount curve from swaps at 2y
    partial_sum = BPV_01y
    discount_curve = np.concatenate((discount_curve, (1-swap_rates_interp[3] * partial_sum)/
                                     (1 + yearfrac(date_1y, np.array(swap_dates_interp[3]), Basis_swap) * swap_rates_interp[3])))
    partial_sum = partial_sum + discount_curve[-1] * yearfrac(date_1y, np.array(swap_dates_interp[3]), Basis_swap)

    # discount curve from swaps over 2y
    for i in range(4, len(swap_dates_interp)):
        discount_curve = np.concatenate((discount_curve,(1 - swap_rates_interp[i] * partial_sum)/
                                (1 + yearfrac(swap_dates_interp[i - 1], np.array(swap_dates_interp[i]), Basis_swap) * swap_rates_interp[i])))
        partial_sum = partial_sum + discount_curve[-1] * yearfrac(swap_dates_interp[i - 1], np.array(swap_dates_interp[i]),Basis_swap)

    # all discount dates
    discount_dates = np.concatenate((discount_dates, swap_dates_interp))
    return discount_dates, discount_curve, data

def pseudo_discount_curve_bootstrap(data, settle, Basis_depos, Basis_fra, Basis_futures,
                                                         Basis_swap, Basis_swap_floating, discount_dates,discount_curve):
    """
    Function that computes pseudo-discount factors curve
    INPUT
    data:                  struct containing dates and quoted prices of the needed financial instruments
    settle:                settlement date
    daycounts:             basis_depos, basis_fra, basis_futures, basis_swap, basis_swap_floating
    discount_dates:        dates of the financial instruments used for the curve construction
    discount_curve:        discount factors curve
    OUTPUT
    P_dates:               dates of the financial instruments used for the pseudo-curve construction
    P_curve:               pseudo-discount factors curve
    """

    # Allocation of dates
    data['m3']['depos']['dates'] = change_weekends(addtodate(settle, np.array([3]), 'month'))
    data['m3']['fra']['dates'] = change_weekends(addtodate(settle, range(1,7), 'month'))

    # 3m depos + fra
    P_03 = 1 / (1 + yearfrac(settle, data['m3']['depos']['dates'], Basis_depos) * data['m3']['depos']['mktRates'])
    P_fwd_36 = 1 / (1 + yearfrac(data['m3']['depos']['dates'], data['m3']['fra']['dates'][-1], Basis_fra) * data['m3']['fra']['mktRates'][-1])
    P_06 = P_03 * P_fwd_36

    # interpolation for 4 and 5 months
    P_04 = interp_ZR(settle, np.concatenate((data['m3']['depos']['dates'], np.array([data['m3']['fra']['dates'][-1]]))),
                     np.concatenate((P_03, P_06)), np.array([data['m3']['fra']['dates'][3]]))
    P_05 = interp_ZR(settle, np.concatenate((data['m3']['depos']['dates'], np.array([data['m3']['fra']['dates'][-1]]))),
                     np.concatenate((P_03, P_06)), np.array([data['m3']['fra']['dates'][4]]))

    # derivation of pseudo-discounts at 1 and 2 months
    P_fwd_14 = 1 / (1 + yearfrac(data['m3']['fra']['dates'][0], np.array([data['m3']['fra']['dates'][3]]), Basis_fra) * data['m3']['fra']['mktRates'][0])
    P_fwd_25 = 1 / (1 + yearfrac(data['m3']['fra']['dates'][1], np.array([data['m3']['fra']['dates'][4]]), Basis_fra) * data['m3']['fra']['mktRates'][1])
    P_01 = P_04 / P_fwd_14
    P_02 = P_05 / P_fwd_25

    # forward pseudo-discount curve
    pseudo_dates = np.concatenate((np.array([settle]), data['m3']['fra']['dates'][:3]))
    pseudo_curve = np.array([P_03,P_fwd_14,P_fwd_25,P_fwd_36])

    # pseudo-discount curve
    P_dates = data['m3']['fra']['dates']
    P_curve = np.array([P_01,P_02,P_03,P_04,P_05,P_06]).flatten()

    # forward rates curve
    fwd_curve = np.concatenate((np.array([data['m3']['depos']['mktRates']]), data['m3']['fra']['mktRates']))

    # futures
    nFutures = 7
    futureRates = 1 - data['m3']['futures']['futPrices'][:nFutures]
    fwdRates = futureRates #no convexity

    # forward pseudo-discounts
    P_fwd = 1/(1 + yearfrac(data['m3']['futures']['settle'][:nFutures], data['m3']['futures']['expiries'][:nFutures], Basis_futures)*fwdRates)

    # update curves (fwd pseudo-discount curve and dates, fwd rates)
    pseudo_dates = np.concatenate((pseudo_dates, data['m3']['futures']['settle'][:nFutures]))
    pseudo_curve = np.concatenate((pseudo_curve.flatten(),P_fwd[:7]))
    fwd_curve = np.concatenate((fwd_curve, fwdRates))


    # update pseudo-discount curve and dates
    for i in range(nFutures):
        P_0s = interp_ZR(settle, P_dates, P_curve, np.array([data['m3']['futures']['settle'][i]]))
        P_dates = np.concatenate((P_dates,np.array([data['m3']['futures']['expiries'][i]])))
        P_curve = np.concatenate((P_curve, P_0s*P_fwd[i]))

    # sort the curves to respect the order of dates
    idx = np.argsort(P_dates)
    P_dates.sort()
    P_curve = P_curve[idx]
    idx = np.argsort(pseudo_dates)
    pseudo_dates.sort()
    pseudo_curve = pseudo_curve[idx]

    # last forward rate of the second year (Point 5 of the bootstrap)
    # construction of the vector w
    index_2y = np.argwhere(data['OIS']['dates'] >= change_weekends(addtodate(settle, 2, 'year'))).flatten()
    data['m3']['swaps']['dates'] = data['OIS']['dates'][index_2y]
    coupon_reset_dates_2y = change_weekends(cfdates(addtodate(data['m3']['swaps']['dates'][0], -2, 'year'), data['m3']['swaps']['dates'][0], 4))
    delta_vec_2ycoupon = yearfrac(np.concatenate((np.array([settle]), coupon_reset_dates_2y[:-1])), coupon_reset_dates_2y, Basis_swap_floating)
    P_curve_interp = interp_ZR(settle, discount_dates, discount_curve, coupon_reset_dates_2y)
    w_vec = delta_vec_2ycoupon*P_curve_interp

    # first 7 forward rates
    P_fwd_curve_interp = interp_fwd_ZR(pseudo_dates, pseudo_curve, np.concatenate((np.array([settle]),coupon_reset_dates_2y[:-2])))
    fwd_1to7 = (1/P_fwd_curve_interp - 1)/delta_vec_2ycoupon[:-1]

    # computation of I_2
    swap_rate_2y = data['m3']['swaps']['mktRates'][0]
    yearly_reset_dates = change_weekends(cfdates(settle, data['m3']['swaps']['dates'][0], 1))
    delta_vec_yearly = yearfrac(np.concatenate((np.array([settle]), yearly_reset_dates[:-1])), yearly_reset_dates, Basis_swap)
    P_curve_interp_1y = interp_ZR(settle, discount_dates, discount_curve, yearly_reset_dates)
    I_2 = swap_rate_2y * sum(delta_vec_yearly*P_curve_interp_1y)

    # fwd rate at 2y and corresponding fwd pseudo-discount
    fwd_8 = 1 / w_vec[-1] * (I_2 - sum(w_vec[:-1]*fwd_1to7))
    P_fwd_8 = 1 / (1 + yearfrac(np.array([coupon_reset_dates_2y[-2]]), np.array([coupon_reset_dates_2y[-1]]), Basis_swap_floating) * fwd_8)

    # update curves
    pseudo_dates = np.concatenate((pseudo_dates,np.array([coupon_reset_dates_2y[-2]])))
    pseudo_curve = np.concatenate((pseudo_curve, P_fwd_8))
    P_0s = interp_ZR(settle, P_dates, P_curve, np.array([coupon_reset_dates_2y[-2]]))
    P_dates = np.concatenate((P_dates,np.array([coupon_reset_dates_2y[-1]])))
    P_curve = np.concatenate((P_curve,P_0s * P_fwd_8))
    fwd_curve = np.concatenate((fwd_curve,np.array([fwd_8])))

    # swaps over 3y
    swap_dates_interp = change_weekends(addtodate(settle, range(3,51), 'year'))
    datenums_swapDates = [date.toordinal() for date in data['m3']['swaps']['dates']]
    datenums_swapDatesInterp = [date.toordinal() for date in swap_dates_interp]
    swap_rates_interp = interp1d(datenums_swapDates, data['m3']['swaps']['mktRates'],kind='cubic')(datenums_swapDatesInterp)
    n_swaps = np.shape(swap_rates_interp)[0]

    # floating legs
    freq = 4 #coupon frequency
    floating_dates = change_weekends(cfdates(settle, swap_dates_interp[-1], freq))
    discounts_floating_dates = interp_ZR(settle, discount_dates, discount_curve, floating_dates)
    delta_floating_dates = yearfrac(np.concatenate((np.array([settle]), floating_dates[:-1])), floating_dates, Basis_swap_floating)
    w_vec = delta_floating_dates*discounts_floating_dates
    floating_dates = np.concatenate((np.array([settle]),floating_dates))

    #fixed legs
    fixed_yearly_dates = change_weekends(cfdates(settle, swap_dates_interp[-1], 1))
    discounts_fixed_dates = interp_ZR(settle, discount_dates, discount_curve, fixed_yearly_dates)
    delta_fixed_dates = yearfrac(np.concatenate((np.array([settle]),  fixed_yearly_dates[:-1])), fixed_yearly_dates, Basis_swap)

    #allocation of dates and forwards
    fwd_dates_new = pseudo_dates
    fwd_curve_new = fwd_curve

    for i in range(n_swaps):
        j = i + 2 +1 #swap maturity
        swap_rate = swap_rates_interp[i]
        I_i = swap_rate * sum(delta_fixed_dates[:j]*discounts_fixed_dates[:j])  #I_i

        # derivation of necessary forwards (using interpolation)
        datesInterp = [date.toordinal() for date in np.concatenate((fwd_dates_new, np.array([floating_dates[freq*j-1]])))]
        newDatesInterp = [date.toordinal() for date in floating_dates[:freq*j]]
        F = lambda x: np.interp(newDatesInterp, datesInterp,np.hstack((fwd_curve_new.flatten(),x)))#np.concatenate((fwd_curve_new, np.array([x]).reshape(-1))))
        fun = lambda x: sum(w_vec[:freq*j]*F(x))-I_i
        last_fwd = fsolve(fun, 0.02)   #last forward (from 9m to 12 m)
        F_vec = F(last_fwd)
        fwd_dates_new = floating_dates[:freq*j]
        fwd_curve_new = F_vec

        # add the yearly pseudo-discount and forward pseudo-discount
        pseudo_dates = np.hstack((pseudo_dates,floating_dates[freq*j-1]))
        pseudo_curve = [pseudo_curve, np.prod(1 / (1 + delta_floating_dates[freq*j-4:freq * j]*F_vec[-4:]))]
        P_0s = interp_ZR(settle, P_dates, P_curve, np.array([floating_dates[freq*j - 4]]))
        P_dates = np.hstack((P_dates,floating_dates[freq*j]))
        P_curve = np.concatenate((P_curve,P_0s * pseudo_curve[-1]))

    return P_dates, P_curve

def bootstrap_dual_curve(data):
    """
    Function that computes discount factors and pseudo-discount factors in the multicurve framework

    INPUT
    data:     struct containing settlement and quoted prices of the financial instruments
              needed for the bootstrap of both curves

    OUTPUT
    B_dates:  dates of the financial instruments used for the curve construction
    B_curve:  discount factors curve
    P_dates:  dates of the financial instruments used for the pseudo-curve construction
    P_curve:  pseudo-discount factors curve
    """

    # Parameters
    settle = data['settlement']
    Basis_swap = 6
    Basis_swap_floating = 2
    Basis_depos = 2
    Basis_fra = 2
    Basis_futures = 2

    # Discount curve
    [B_dates, B_curve, data] = discount_curve_bootstrap(data, settle, Basis_swap)

    # Pseudo discount curve
    #daycounts = [Basis_depos, Basis_fra, Basis_futures, Basis_swap, Basis_swap_floating]
    [P_dates, P_curve] = pseudo_discount_curve_bootstrap(data, settle, Basis_depos, Basis_fra, Basis_futures, Basis_swap,
                                                         Basis_swap_floating,B_dates, B_curve)

    '''
    #plot of discounts and pseudo-discounts
    plt.plot_date(B_dates, B_curve, linestyle='-', marker=None)
    plt.plot_date(P_dates, P_curve, linestyle='-', marker=None)

    plt.xlabel('Years')
    plt.legend(['Discounts', 'Pseudo-Discounts'])
    plt.title("Discounts and Pseudo Discounts on {}".format(data['settlement']))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    plt.show()'''

    return B_dates, B_curve, P_dates, P_curve


