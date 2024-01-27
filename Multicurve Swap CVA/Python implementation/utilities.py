import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

def day30_360(startDate, endDate):
    """
    This function computes the year fraction with Act 30/360 between the given dates
    INPUT:
    startDate           Initial dates
    endDate             Final dates

    OUTPUT:
    delta               Year fraction vector between given dates
    """
    if not isinstance(startDate,np.ndarray):
        startDate = np.array([startDate])


    n_days = []
    if startDate.flatten().shape[0] == 1:
        if startDate[0].day == 31:
            startDate = startDate + relativedelta(days=-1)
        for endD in endDate.flatten():
            if endD.day == 31:
                endD = endD + relativedelta(days=-1)
            n_days.append(360 * (endD.year - startDate[0].year) + 30 * (endD.month - startDate[0].month) + endD.day - startDate[0].day)
    else:
        for startD,endD in zip(startDate.flatten(),endDate.flatten()):
            if endD.day == 31:
                endD = endD + relativedelta(days=-1)
            if startD.day == 31:
                startD = startD + relativedelta(days=-1)
            n_days.append(360 * (endD.year - startD.year) + 30 * (endD.month - startD.month) + endD.day - startD.day)

    delta = np.array(n_days)/360
    return delta

def yearfrac(startDate, endDate, basis):
    """
    This function returns an array of year fractions between the given initial and final dates

    INPUT:
    startDate:       Initial dates
    endDate:         Final dates
    basis:           Basis for year fraction: 2 (Act/360), 3 (Act/365), 6 (Act/30/360)

    OUTPUT:
    result:     Converted dates

    """

    if basis == 2: #Act 360
        return (endDate - startDate).astype('timedelta64[D]') / np.timedelta64(360, 'D')
    elif basis == 3: #Act 365
        return (endDate - startDate).astype('timedelta64[D]') / np.timedelta64(365, 'D')
    elif basis == 6: #Act 30/360 European
        return day30_360(startDate, endDate)
    else:
        raise Exception("Unknown basis for yearfrac function")

def change_weekends(date):
    """
    This function returns an array of dates with the following business day convention

    INPUT:
    date:       Dates that need to be converted to following business day convention

    OUTPUT:
    result:     Converted dates

    """
    if not isinstance(date,np.ndarray):
        date_out = np.array(date)
    else:
        date_out = date.copy()

    for i in range(np.shape(date_out)[0]):
        if date_out[i].weekday() == 5:
            date_out[i] = date_out[i] + timedelta(days=2)
        if date_out[i].weekday() == 6:
            date_out[i] = date_out[i] + timedelta(days=1)
    return date_out

def addtodate(date, num, unit):
    """
    This function returns an array with the requested dates

    INPUT:
    date:           Initial date
    num:            Number of units to add to the initial date
    unit:           Unit of num (day, month or year)

    OUTPUT:
    result:         Array of required dates
    """
    if not isinstance(num,np.ndarray):
        num = np.array(num)

    num = num.flatten()
    date_out = []
    if unit == 'day':
        for i in range(np.shape(num)[0]):
            date_out.append(date + timedelta(days=int(num[i])))
        return np.array(date_out).flatten()
    elif unit == 'month':
        for i in range(np.shape(num)[0]):
            date_out.append(date + relativedelta(months=int(num[i])))
        return np.array(date_out).flatten()
    elif unit == 'year':
        for i in range(np.shape(num)[0]):
            date_out.append(date + relativedelta(years=int(num[i])))
        return np.array(date_out).flatten()
    else:
        raise Exception("Only day, month or year unit allowed")


def cfdates(startDate, endDate, freq):
    """
    Function that returns an array of dates between startDate, endDate with frequency freq

    INPUT:
    startDate:              Initial date
    endDate:                Final date
    freq:                   Frequency

    OUTPUT:
    date_out:               Dates output vector
    """
    n_months = 12/freq
    if endDate + relativedelta(months=-n_months) > startDate:
        date_out = np.array([endDate + relativedelta(months=-n_months), endDate])
    else:
        date_out = np.array([endDate])
    while date_out[0] + relativedelta(months=-n_months) > startDate:
        date_out = np.concatenate((np.array([date_out[0] + relativedelta(months=-n_months)]), date_out))
    return date_out

def interp_ZR(settle, dates, discounts, new_dates):
    """
    Function that computes the zero-rate linear interpolation on discounts

    INPUT
    settle:         settlement date
    dates:          dates of the discount factors curve
    discounts:      discount factors curve
    new_dates:      dates of interest

    OUTPUT
    new_discounts:  discount factors in the desidered dates
    """
    basis_interpolation = 3 #Act / 365
    # find the corresponding zero rates
    zr = -np.log(discounts)/yearfrac(settle, dates, basis_interpolation)
    b = yearfrac(settle, dates, basis_interpolation)
    # interpolate on the zero rates
    datenum_dates = [date.toordinal() for date in dates]
    datenum_newdates = [date.toordinal() for date in new_dates]
    zr_interp = interp1d(datenum_dates, zr,fill_value="extrapolate")(datenum_newdates) #linearly interpolates and extrapolates

    # compute new discount factors
    new_discounts = np.exp(-zr_interp*yearfrac(settle, new_dates, basis_interpolation))

    return new_discounts

def interp_fwd_ZR(dates, discounts, new_dates):
    """
    Function that computes the zero-rate linear interpolation on forward discounts

    INPUT
    settle:         settlement date
    dates:          dates of the forward discount factors set
    discounts:      forward discount factors set
    new_dates:      dates of interest

    OUTPUT
    new_discounts:  forward discount factors in the desidered dates
    """
    # act/365 daycount for linear ZR interpolation
    basis_interpolation = 3
    # forward expiries
    forward_expiries = np.array([change_weekends(addtodate(dates[i], 3, 'month')) for i in range(np.shape(dates.flatten())[0])]).flatten()
    # new forward expiries
    new_forward_expiries = np.array([change_weekends(addtodate(new_dates[i], 3, 'month')) for i in range(np.shape(new_dates.flatten())[0])]).flatten()
    # find the corresponding zero rates
    zr = -np.log(discounts)/yearfrac(dates, forward_expiries, basis_interpolation)
    # interpolate on the zero rates
    datenum_dates = [date.toordinal() for date in dates]
    datenum_newdates = [date.toordinal() for date in new_dates]
    zr_interp = interp1d(datenum_dates, zr)(datenum_newdates)
    # compute new forward discount factors
    new_discounts = np.exp(-zr_interp*yearfrac(new_dates, new_forward_expiries, basis_interpolation))
    return new_discounts

def dates_for_tree(initial_date,expiry_date,expiry,dt,month_steps):
    """
    Function that finds the dates for the trinomial tree

    INPUT
    initial_date:  initial date
    expiry_date:   expiry date
    expiry:        expiry in number of months
    dt:            time interval for the tree
    month_steps:   number of monthly steps for the tree

    OUTPUT
    dates_tree:    dates in the tree """

    # mean number of days in a month (between settlement date and expiry date)
    mean_number_of_days = (expiry_date - initial_date).days/expiry

    # dates for the tree
    Dates_tree = addtodate(initial_date,np.array([round(dt*i*mean_number_of_days) for i in range(1,int(expiry)*month_steps+1)]),'day')

    return Dates_tree



