import pandas as pd
import numpy as np
from datetime import datetime
import re
from utilities import change_weekends, addtodate


def readExcelData(filename, formatData, n_sheet):
    """
    Function that reads data from excel

    INPUT
    filename:    excel file name where data are stored
    formatData:  data format in Excel
    N_sheet:     1 -> 24 June 2022, 2 -> 31 January 2023

    OUTPUT
    data:        struct updated with referenceDate, settleDate, OIS data, 3M, Futures, Swaption vols
    """

    data = {}

    # Read settlement date
    settlement = pd.read_excel(filename, sheet_name=n_sheet, usecols="B", nrows=1, dtype=str).values[0][0]

    # Convert settlement date to datetime
    data['settlement'] = datetime.strptime(settlement, formatData)

    # Read OIS data
    ois_mkt_rates = pd.read_excel(filename, sheet_name=n_sheet, usecols="C", skiprows=4,
                                  nrows=36).values.flatten() / 100
    data['OIS'] = {'mktRates': np.array(ois_mkt_rates)}

    # Read 3m data
    data['m3'] = {}

    # Read depos
    depos_mkt_rate = pd.read_excel(filename, sheet_name=n_sheet, usecols="J", skiprows=4, nrows=1).values[0][0] / 100
    data['m3']['depos'] = {'mktRates': np.array(depos_mkt_rate)}

    # Read fra
    fra_mkt_rates = pd.read_excel(filename, sheet_name=n_sheet, usecols="J", skiprows=5, nrows=3).values.flatten() / 100
    data['m3']['fra'] = {'mktRates': np.array(fra_mkt_rates)}

    # Read futures
    settle_fut = pd.read_excel(filename, sheet_name=n_sheet, usecols="Q", skiprows=4, nrows=8,dtype=str).values.flatten()
    expiries_fut = pd.read_excel(filename, sheet_name=n_sheet, usecols="R", skiprows=4, nrows=8,dtype=str).values.flatten()
    fut_prices = pd.read_excel(filename, sheet_name=n_sheet, usecols="J", skiprows=8, nrows=8).values.flatten() / 100
    data['m3']['futures'] = {'settle': [datetime.strptime(date, formatData) for date in settle_fut],
                             'expiries': [datetime.strptime(date, formatData) for date in expiries_fut],
                             'futPrices': np.array(fut_prices)}

    data['m3']['futures']['settle'] = np.array(data['m3']['futures']['settle'])
    data['m3']['futures']['expiries'] = np.array(data['m3']['futures']['expiries'])

    # Read swaps
    swap_mkt_rates = pd.read_excel(filename, sheet_name=n_sheet, usecols="J", skiprows=16,
                                   nrows=17).values.flatten() / 100
    data['m3']['swaps'] = {'mktRates': np.array(swap_mkt_rates)}

    return data


def readExcelData_SwapAnnex(filename, formatData):
    """
    Function that reads swap annex for the amortized contract

    INPUT
    filename:    excel file name where data are stored
    formatData:  data format in Excel

    OUTPUT
    data:        struct updated with start and end dates of each payment period, with the corresponding notional
    """
    data = {}

    # Save payment dates
    pay_dates = pd.read_excel(filename, usecols="A", nrows=61, header=None).values[1:]
    data["PayDates"] = pay_dates.flatten()

    # Save start date of each payment period
    start_dates = pd.read_excel(filename, usecols="B", nrows=61, header=None).values[1:]
    data["StartDates"] = start_dates.flatten()

    # Save end date of each payment period
    end_dates = pd.read_excel(filename, usecols="C", nrows=61, header=None).values[1:]
    data["EndDates"] = end_dates.flatten()

    # Save notional of each payment period
    data["Notionals"] = pd.read_excel(filename, usecols="E", nrows=61, header=None).values[1:].flatten()

    return data

def readExcelData_SwaptionVol(filename,formatData,N_sheet):
     """ Function that reads swaption volatilities from market values

        INPUT
        filename:    excel file name where data are stored
        formatData:  data format in Excel
        N_sheet:     1 -> 24 June 2022, 2 -> 31 January 2023

        OUTPUT
        data:        updated struct containing swaption volatilities, with the relative expiry and tenor
     """

     data = {}

     # Read implied volatilities
     data['vol'] = pd.read_excel(filename, sheet_name=N_sheet, usecols="Q:AE", skiprows=15, nrows=21).values * 1e-4

     # Read expiries
     expiries = pd.read_excel(filename, sheet_name=N_sheet, usecols="P", skiprows=15, nrows=21).values.flatten()
     data['expiry'] = [None]*len(expiries)

     # Read settlement date
     settlement = pd.read_excel(filename, sheet_name=N_sheet, usecols="B", nrows=1, dtype=str).values[0][0]

     # Date conversion
     settle = datetime.strptime(settlement, formatData)

     # Convert expiries from years and months to number
     for i, expiry in enumerate(expiries):
         time_unit = np.array([int(re.findall(r'\d+',expiry)[0])])
         if 'M' in expiry:
             data['expiry'][i] = change_weekends(addtodate(settle,time_unit, 'month'))
         else:
             data['expiry'][i] = change_weekends(addtodate(settle,time_unit, 'year'))
     data['expiry'] = np.array(data['expiry']).flatten()

     # Read tenors
     tenors = pd.read_excel(filename, sheet_name=N_sheet, usecols="Q:AE", skiprows=14, nrows=1).values.flatten()
     data['tenor'] = [None]*len(tenors)

     # Convert tenors from years to number
     for i, tenor in enumerate(tenors):
         data['tenor'][i] = np.array([int(re.findall(r'\d+', tenor)[0])])
     data['tenor'] = np.array(data['tenor']).flatten()

     return data

def readExcelData_vol_calibration(filename,N_sheet,expiries,tenors):
    """
        Function that reads implied volatilies only for the selected expiries and tenors

        INPUT
        filename:  excel file name where data are stored
        N_sheet:   1 -> 24 June 2022, 2 -> 31 January 2023
        expiries:  vector of expiries of the swaptions of interest
        tenors:    vector of tenors of the swaptions of interest

        OUTPUT
        vols:      vector of implied volatilies of the swaptions of interest
    """

    # implied volatilities
    volatilities = pd.read_excel(filename, sheet_name=N_sheet, usecols="Q:AE", skiprows=15, nrows=21).values * 1e-4

    # Read expiries
    exp = pd.read_excel(filename, sheet_name=N_sheet, usecols="P", skiprows=15, nrows=21).values.flatten()

    # Read tenors
    ten = pd.read_excel(filename, sheet_name=N_sheet, usecols="Q:AE", skiprows=14, nrows=1).values.flatten()

    # initialization of volatilities
    vols = np.zeros((len(expiries), 1))

    # extraction of volatilities of interest
    for i in range(len(expiries)):
        for j in range(len(exp)):
            if 'Y' in exp[j]:
                if int(re.findall(r'\d+',exp[j])[0]) == expiries[i]:
                    idx_exp = j

        for j in range(len(ten)):
            if 'Y' in ten[j]:
                if int(re.findall(r'\d+',ten[j])[0]) == tenors[i]:
                    idx_ten = j

        vols[i] = volatilities[idx_exp][idx_ten]

    return vols

