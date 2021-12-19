import pandas as pd
import numpy as np
from datetime import datetime, timedelta


def load_convbond_list_from_excel(path, client, upsert=True):
    """Load convbond list to rmi_convbond_info from excels,
    which are retrieved from Wind.

    Parameters
    ----------
    path : list of strings
        A list of strings, which contain the path of all the excel forms.
    client :
        A mongodb client.
    Returns
    -------
        None
    """
    for p in path:
        excel = pd.ExcelFile(p)
        df = excel.parse("file")

        convbond_info = client.convbond.rmi_convbond_info

        for i in range(0, len(df)):
            convbond_info.update_one({'windcode': df['Code'][i]},
                                     {'$set': {'windcode': df['Code'][i]}},
                                     upsert)


def get_windcodes(client):
    convbond_info = client.convbond.rmi_convbond_info
    cursor = convbond_info.find({}, {'windcode': 1})
    df = pd.DataFrame(list(cursor))
    wind_ids = df['windcode'].values
    return wind_ids


def upsert_wssdata(filter, wssdata, w, collection, upsert=True):
    key_list = [item.lower() for item in wssdata.Fields]
    value_list = [item for sublist in wssdata.Data for item in sublist]

    data_dict = dict(zip(key_list, value_list))

    collection.update_one(filter,
                         {'$set': data_dict},
                         upsert)
    return data_dict


def update_one_to_info(windcode, w, client, upsert=True):
    def update_info_paymentdate():
        paymentdate_list = []
        cols = "paymentdate"
        for i in range(0, interest_freq):
            params = "N={}".format(i)
            data = w.wss(windcode, cols, params)
            paymentdate_list += data.Data[0]
            bond_info = client.convbond.rmi_convbond_info
            bond_info.update_one({'windcode': windcode},
                                 {'$set': {'paymentdate': paymentdate_list}},
                                 upsert)

    cols = "windcode,trade_code,sec_name,fullname,par,issueamount,carrydate,carryenddate,maturitydate,term,interesttype,coupontxt,paymenttype,actualbenchmark,coupon,form,interestfrequency,coupondatetxt,taxfree,taxrate,mktpricetype,redemption_beginning,redemption_feeration,repaymentmethod,paymentorder"
    params = "unit=1"
    data = w.wss(windcode, cols, params)
    dict_inserted = upsert_wssdata({'windcode': windcode}, data, w, client.convbond.rmi_convbond_info, upsert)

    interest_freq = dict_inserted['interestfrequency']
    update_info_paymentdate()

    cols = "ipo_date,clause_interest_5,clause_interest_8,clause_interest_6,clause_interest_compensationinterest,underlyingcode,clause_conversion_2_relativeswapsharemonth,clause_conversion_2_swapsharestartdate,clause_conversion_2_swapshareenddate,clause_conversion_code,clause_conversion_2_conversionproportion,clause_conversion2_tosharepriceadjustitem,clause_reset_item,clause_reset_isexitreset,clause_reset_resetstartdate,clause_reset_resetperiodenddate,clause_reset_resetmaxtimespan,clause_reset_resettimespan,clause_reset_resettriggerratio,clause_reset_stockpricelowestlimit,clause_calloption_relativecalloptionperiod,clause_calloption_redemptiontimesperyear,clause_calloption_conditionalredeemstartdate,clause_calloption_conditionalredeemenddate,clause_calloption_redeemmaxspan,clause_calloption_redeemspan,clause_calloption_redemptionprice,clause_calloption_redemptionmemo,clause_calloption_interestdisposal,clause_calloption_redeemitem,clause_putoption_putbackperiodobs,clause_putoption_conditionalputbackstartenddate,clause_putoption_conditionalputbackenddate,clause_putoption_putbacktriggermaxspan,clause_putoption_putbacktriggerspan,clause_putoption_redeem_triggerproportion,clause_putoption_resellingprice,clause_putoption_resellingpriceexplaination,clause_putoption_putbacktimesperyear,clause_putoption_interestdisposing,clause_putoption_sellbackitem,clause_putoption_timeputbacktimes,clause_putoption_timeputbackclause"
    params = ""
    data = w.wss(windcode, cols, params)
    upsert_wssdata({'windcode': windcode}, data, w, client.convbond.rmi_convbond_info, upsert)

    print('Inserted {} to rmi_convbond_info'.format(windcode))


def update_coupon_rates_hardcode(client, upsert=True):
    hardcoded_coupon_rates = {
        '110030.SH': {
            '1': 0.6, '2': 0.8, '3': 1., '4': 1.5, '5': 2.0
        },
        '110031.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 1.6
        },
        '110032.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.6, '6': 2.0
        },
        '110033.SH': {
            '1': 0.3, '2': 0.5, '3': 0.9, '4': 1.4, '5': 1.7, '6': 2.0
        },
        '110034.SH': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 2.0
        },
        '110038.SH': {
            '1': 0.2, '2': 0.5, '3': 0.8, '4': 1.5, '5': 6.0
        },
        '110039.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '110040.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '110041.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '110042.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '110043.SH': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 1.8
        },
        '113008.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 1.6
        },
        '113009.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 1.6
        },
        '113010.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113011.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113012.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '113013.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113014.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113015.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '113016.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113017.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '113018.SH': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 1.8
        },
        '113019.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 2.0
        },
        '113502.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113503.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113504.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '113505.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '120001.SZ': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0
        },
        '123001.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '123002.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '123003.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '123004.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '123005.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '123006.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 1.0, '5': 1.5, '6': 2.0
        },
        '123007.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '123008.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '123009.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '127003.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '127004.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.7, '5': 1.8, '6': 2.0
        },
        '127005.SZ': {
            '1': 0.2, '2': 0.4, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '127006.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 2.0
        },
        '128010.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.6, '5': 1.6, '6': 1.6
        },
        '128012.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.3, '5': 1.3, '6': 1.6
        },
        '128013.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128014.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128015.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128016.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128017.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128018.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128019.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128020.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 1.8
        },
        '128021.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128022.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128023.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128024.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.8, '4': 1.2, '5': 1.6, '6': 2.0
        },
        '128025.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 1.8
        },
        '128026.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128027.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128028.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.5, '6': 1.8
        },
        '128029.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.5
        },
        '128030.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128032.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128033.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128034.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 1.8
        },
        '128035.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 2.0
        },
        '128036.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128037.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '128038.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 2.0
        },
        '132002.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0
        },
        '132003.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0
        },
        '132004.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0, '6': 1.0
        },
        '132005.SH': {
            '1': 1.7, '2': 1.7, '3': 1.7, '4': 1.7, '5': 1.7
        },
        '132006.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0
        },
        '132007.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0
        },
        '132008.SH': {
            '1': 1.7, '2': 1.7, '3': 1.7, '4': 1.7, '5': 1.7
        },
        '132009.SH': {
            '1': 1., '2': 1., '3': 1., '4': 1., '5': 1.
        },
        '132010.SH': {
            '1': 1., '2': 1., '3': 1.
        },
        '132011.SH': {
            '1': 1., '2': 1., '3': 1., '4': 1., '5': 1.
        },
        '132012.SH': {
            '1': 1., '2': 1., '3': 1.
        },
        '132013.SH': {
            '1': 1., '2': 1., '3': 1.
        },
        '132015.SH': {
            '1': 1.4, '2': 1.4, '3': 1.4, '4': 1.4, '5': 1.4
        }
    }

    for key, value in hardcoded_coupon_rates.items():
        bond_info = client.convbond.rmi_convbond_info
        bond_info.update_one({'windcode': key},
                             {'$set': {'coupon_rate': value}},
                             upsert)


def update_adjusted_coupon_rates_hardcode(client, upsert=True):
    hardcoded_coupon_rates = {
        '110030.SH': {
            '1': 0.6, '2': 0.8, '3': 1., '4': 1.5, '5': 6.0
        },
        '110031.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 7
        },
        '110032.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.6, '6': 6.0
        },
        '110033.SH': {
            '1': 0.3, '2': 0.5, '3': 0.9, '4': 1.4, '5': 1.7, '6': 8.0
        },
        '110034.SH': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 8.0
        },
        '110038.SH': {
            '1': 0.2, '2': 0.5, '3': 0.8, '4': 1.5, '5': 6.0
        },
        '110039.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '110040.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '110041.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '110042.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 5.0
        },
        '110043.SH': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 6.0
        },
        '113008.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 6.6
        },
        '113009.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.5, '6': 6.0
        },
        '113010.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8
        },
        '113011.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 5.0
        },
        '113012.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 5.0
        },
        '113013.SH': {
            '1': 0.2, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 5.0
        },
        '113014.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '113015.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '113016.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '113017.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '113018.SH': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 6.0
        },
        '113019.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 10.0
        },
        '113502.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '113503.SH': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '113504.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '113505.SH': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '120001.SZ': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0+7
        },
        '123001.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '123002.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '123003.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '123004.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '123005.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 7.0
        },
        '123006.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 1.0, '5': 1.5, '6': 7.0
        },
        '123007.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 7.0
        },
        '123008.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '123009.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '127003.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 10.0
        },
        '127004.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.7, '5': 1.8, '6': 10.0
        },
        '127005.SZ': {
            '1': 0.2, '2': 0.4, '3': 1.0, '4': 1.5, '5': 1.8, '6': 5.0
        },
        '127006.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 5.0
        },
        '128010.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.6, '5': 1.6, '6': 8.0
        },
        '128012.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.3, '5': 1.3, '6': 3
        },
        '128013.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '128014.SZ': {
            '1': 0.5, '2': 0.7, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '128015.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 8.0
        },
        '128016.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128017.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128018.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128019.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128020.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 8.0
        },
        '128021.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128022.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128023.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '128024.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.8, '4': 1.2, '5': 1.6, '6': 5.0
        },
        '128025.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.3, '5': 1.5, '6': 6.0
        },
        '128026.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128027.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128028.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.5, '6': 6.0
        },
        '128029.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 6.0
        },
        '128030.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '128032.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128033.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128034.SZ': {
            '1': 0.3, '2': 0.5, '3': 0.8, '4': 1.0, '5': 1.3, '6': 6.0
        },
        '128035.SZ': {
            '1': 0.2, '2': 0.4, '3': 0.6, '4': 0.8, '5': 1.6, '6': 5.0
        },
        '128036.SZ': {
            '1': 0.4, '2': 0.6, '3': 1.0, '4': 1.5, '5': 1.8, '6': 6.0
        },
        '128037.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 8.0
        },
        '128038.SZ': {
            '1': 0.3, '2': 0.5, '3': 1.0, '4': 1.5, '5': 1.8, '6': 9.0
        },
        '132002.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0
        },
        '132003.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0
        },
        '132004.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0, '6': 1.0+3
        },
        '132005.SH': {
            '1': 1.7, '2': 1.7, '3': 1.7, '4': 1.7, '5': 1.7+7.5
        },
        '132006.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0+12-1
        },
        '132007.SH': {
            '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0+5
        },
        '132008.SH': {
            '1': 1.7, '2': 1.7, '3': 1.7, '4': 1.7, '5': 1.7+6
        },
        '132009.SH': {
            '1': 1., '2': 1., '3': 1., '4': 1., '5': 5.0+1.0
        },
        '132010.SH': {
            '1': 1., '2': 1., '3': 3.
        },
        '132011.SH': {
            '1': 1., '2': 1., '3': 1., '4': 1., '5': 3.
        },
        '132012.SH': {
            '1': 1., '2': 1., '3': 1.0+2.0
        },
        '132013.SH': {
            '1': 1., '2': 1., '3': 1.0+3.0
        },
        '132015.SH': {
            '1': 1.4, '2': 1.4, '3': 1.4, '4': 1.4, '5': 1.4+5.0
        }
    }

    for key, value in hardcoded_coupon_rates.items():
        bond_info = client.convbond.rmi_convbond_info
        bond_info.update_one({'windcode': key},
                             {'$set': {'adjusted_coupon_rate': value}},
                             upsert)


def update_call_put_clause_hardcode(client, upsert=True):
    hardcoded_call_put_clause = {
        '110030.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '110031.SH': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '110032.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '110033.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '110034.SH': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '110038.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 125, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 50, 'put_price': None
        },
        '110039.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '110040.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '110041.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '110042.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '110043.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '113008.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '113009.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113010.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 80, 'put_price': 103
        },
        '113011.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '113012.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113013.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '113014.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113015.SH': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113016.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113017.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113018.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '113019.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113502.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113503.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113504.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '113505.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '120001.SZ': {
            'call_n': 10, 'call_m': 20, 'call_percentage': 130, 'call_price': None,
            'put_n': 10, 'put_m': 20, 'put_percentage': 80, 'put_price': 103
        },
        '123001.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123002.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123003.SZ': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123004.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123005.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123006.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123007.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123008.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '123009.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '127003.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '127004.SZ': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '127005.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '127006.SZ': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128010.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128012.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': 103,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '128013.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128014.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128015.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128016.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128017.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128018.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128019.SZ': {
            'call_n': 20, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128020.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128021.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128022.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128023.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128024.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '128025.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128026.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128027.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128028.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128029.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128030.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128032.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128033.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128034.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '128035.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': 103
        },
        '128036.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128037.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '128038.SZ': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132002.SH': {
            'call_n': 10, 'call_m': 20, 'call_percentage': 135, 'call_price': 107,
            'put_n': 10, 'put_m': 20, 'put_percentage': 80, 'put_price': 107
        },
        '132003.SH': {
            'call_n': 10, 'call_m': 30, 'call_percentage': 120, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 80, 'put_price': None
        },
        '132004.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132005.SH': {
            'call_n': None, 'call_m': None, 'call_percentage': None, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '132006.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 120, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132007.SH': {
            'call_n': 10, 'call_m': 20, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 80, 'put_price': None
        },
        '132008.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132009.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132010.SH': {
            'call_n': 10, 'call_m': 20, 'call_percentage': 120, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132011.SH': {
            'call_n': 10, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132012.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        },
        '132013.SH': {
            'call_n': None, 'call_m': None, 'call_percentage': None, 'call_price': None,
            'put_n': None, 'put_m': None, 'put_percentage': None, 'put_price': None
        },
        '132015.SH': {
            'call_n': 15, 'call_m': 30, 'call_percentage': 130, 'call_price': None,
            'put_n': 30, 'put_m': 30, 'put_percentage': 70, 'put_price': None
        }
    }

    for key, value in hardcoded_call_put_clause.items():
        bond_info = client.convbond.rmi_convbond_info
        bond_info.update_one({'windcode': key},
                             {'$set': {'call_put_clause': value}},
                             upsert)


def update_one_to_daily(windcode, date, w, client, upsert=True):
    cols = "clause_calloption_triggerproportion,clause_conversion2_swapshareprice,clause_conversion2_conversionproportion,clause_conversion2_bondlot,clause_conversion2_bondproportion,close,high,low,volume,latestissurercreditrating,amount,rate_latest,rate_changesofrating,rate_style"
    params = "tradeDate={};priceAdj=CP;cycle=D;unit=1".format(date.strftime("%Y%m%d"))
    data = w.wss(windcode, cols, params)
    upsert_wssdata({'windcode': windcode, 'date': date.replace(hour=0, minute=0, second=0, microsecond=0)}, data, w,
                   client.convbond.rmi_convbond_daily, upsert)

    print('Inserted {} to rmi_convbond_daily on {}'.format(windcode, date))


def get_hist_vol_f(days, client, wind_id, date):
    cursor = client.convbond.rmi_convbond_underlying.find({'windcode': wind_id,
                                                           'date': {'$lte': date},
                                                           'volume': {'$nin': [float('nan'), 0.0]}},
                                                          {'close_f': 1, 'date': 1, 'volume': 1, '_id': 0}).sort('date', -1).limit(days)
    ls = [x['close_f'] for x in list(cursor)]
    S_a = ls[0:-1]
    S_b = ls[1:]
    S_a = np.array(S_a, dtype=float)
    S_b = np.array(S_b, dtype=float)
    mean = np.sum(np.log(np.divide(S_b, S_a)))/(len(ls)-1)
    return np.sqrt(np.sum(np.square(np.log(np.divide(S_b, S_a)) - mean))*252/(len(ls)-2)) * 100


def get_hist_vol(days, client, wind_id, date):
    cursor = client.convbond.rmi_convbond_underlying.find({'windcode': wind_id,
                                                           'date': {'$lte': date},
                                                           'volume': {'$nin': [float('nan'), 0.0]}},
                                                          {'close': 1, 'date': 1, 'volume': 1, '_id': 0}).sort('date', -1).limit(days)
    ls = [x['close'] for x in list(cursor)]
    S_a = ls[0:-1]
    S_b = ls[1:]
    S_a = np.array(S_a, dtype=float)
    S_b = np.array(S_b, dtype=float)
    mean = np.sum(np.log(np.divide(S_b, S_a)))/(len(ls)-1)
    return np.sqrt(np.sum(np.square(np.log(np.divide(S_b, S_a)) - mean))*252/(len(ls)-2)) * 100


def update_one_to_underlying(windcode, date, w, client, upsert=True):
    underlyingcode = client.convbond.rmi_convbond_info.find_one({'windcode': windcode}, {'windcode': 1, 'underlyingcode': 1})['underlyingcode']
    cols = "open,high,low,close,volume,amt"
    params = "tradeDate={};cycle=D;priceAdj=U".format(date.strftime("%Y%m%d"))
    data = w.wss(underlyingcode, cols, params)
    upsert_wssdata({'windcode': underlyingcode, 'date': date.replace(hour=0, minute=0, second=0, microsecond=0)}, data, w,
                   client.convbond.rmi_convbond_underlying, upsert)

    data = w.wss(underlyingcode, "open,high,low,close", "tradeDate={};priceAdj=F;cycle=D".format(date.strftime("%Y%m%d"))).Data
    client.convbond.rmi_convbond_underlying.update_one(
                                                        {'windcode': underlyingcode, 'date': date.replace(hour=0, minute=0, second=0, microsecond=0)},
                                                        {'$set': {'open_f': data[0][0],
                                                                  'high_f': data[1][0],
                                                                  'low_f': data[2][0],
                                                                  'close_f': data[3][0]}},
                                                        upsert)

    hist_vol_15D_f = get_hist_vol_f(15+1, client, underlyingcode, date)
    hist_vol_1M_f = get_hist_vol_f(30+1, client, underlyingcode, date)
    hist_vol_45D_f = get_hist_vol_f(45+1, client, underlyingcode, date)
    hist_vol_2M_f = get_hist_vol_f(60+1, client, underlyingcode, date)
    hist_vol_3M_f = get_hist_vol_f(90+1, client, underlyingcode, date)
    hist_vol_15D = get_hist_vol(15 + 1, client, underlyingcode, date)
    hist_vol_1M = get_hist_vol(30 + 1, client, underlyingcode, date)
    hist_vol_45D = get_hist_vol(45 + 1, client, underlyingcode, date)
    hist_vol_2M = get_hist_vol(60 + 1, client, underlyingcode, date)
    hist_vol_3M = get_hist_vol(90 + 1, client, underlyingcode, date)
    client.convbond.rmi_convbond_underlying.update_one({'windcode': underlyingcode, 'date': date.replace(hour=0, minute=0, second=0, microsecond=0)},
                         {'$set': {'hist_vol_15D_f': hist_vol_15D_f,
                                   'hist_vol_1M_f': hist_vol_1M_f,
                                   'hist_vol_45D_f': hist_vol_45D_f,
                                   'hist_vol_2M_f': hist_vol_2M_f,
                                   'hist_vol_3M_f': hist_vol_3M_f,
                                   'hist_vol_15D': hist_vol_15D,
                                   'hist_vol_1M': hist_vol_1M,
                                   'hist_vol_45D': hist_vol_45D,
                                   'hist_vol_2M': hist_vol_2M,
                                   'hist_vol_3M': hist_vol_3M
                                   }},
                         upsert)
    print('Inserted {} to rmi_convbond_underlying on {}'.format(underlyingcode, date))


def patch_one_to_daily(windcode, beg_date, end_date, w, client, upsert=True):
    beg_date = beg_date.replace(hour=0, minute=0, second=0, microsecond=0)
    date = client.convbond.rmi_convbond_info.find_one({'windcode': windcode}, {'windcode': 1, 'carrydate': 1})['carrydate']
    date = date if beg_date < date else beg_date
    while date <= end_date:
        update_one_to_daily(windcode, date, w, client, upsert)
        date += timedelta(days=1)


def patch_one_to_underlying(windcode, beg_date, end_date, w, client, upsert=True):
    beg_date = beg_date.replace(hour=0, minute=0, second=0, microsecond=0)
    date = client.convbond.rmi_convbond_info.find_one({'windcode': windcode}, {'windcode': 1, 'carrydate': 1})['carrydate']
    date = date-timedelta(180) if beg_date < date else beg_date
    # end_date = client.convbond.rmi_convbond_info.find_one({'windcode': windcode}, {'windcode': 1, 'carrydate': 1})['carrydate']
    while date <= end_date:
        update_one_to_underlying(windcode, date, w, client, upsert)
        date += timedelta(days=1)


def update_to_info(w, client, upsert=True):
    for windcode in get_windcodes(client):
        update_one_to_info(windcode, w, client, upsert)

    update_coupon_rates_hardcode(client)
    update_adjusted_coupon_rates_hardcode(client)
    update_call_put_clause_hardcode(client)


def update_to_daily(date, w, client, upsert=True):
    for windcode in get_windcodes(client):
        update_one_to_daily(windcode, date, w, client, upsert)


def update_to_underlying(date, w, client, upsert=True):
    for windcode in get_windcodes(client):
        update_one_to_underlying(windcode, date, w, client, upsert)


def patch_to_daily(beg_date, end_date, w, client, upsert=True, only_new=False):
    for windcode in get_windcodes(client):
        if only_new:
            if client.convbond.rmi_convbond_daily.find_one({'windcode': windcode}) is not None:
                print('skipping {}'.format(windcode))
                continue
        print('patching {} to daily'.format(windcode))
        patch_one_to_daily(windcode, beg_date, end_date, w, client, upsert)


def patch_to_underlying(beg_date, end_date, w, client, upsert=True, only_new=False):
    for windcode in get_windcodes(client):
        if only_new:
            underlyingcode = \
                client.convbond.rmi_convbond_info.find_one({'windcode': windcode})['underlyingcode']
            if client.convbond.rmi_convbond_underlying.find_one({'windcode': underlyingcode}) is not None:
                print('skipping {}'.format(windcode))
                continue

        print('patching {} to underlying'.format(windcode))
        patch_one_to_underlying(windcode, beg_date, end_date, w, client, upsert)


if __name__ == '__main__':
    import pymongo
    client = pymongo.MongoClient("mongodb://siteRootAdmin:rmi%40119613@172.18.101.218:27017,172.18.101.218:27010,172.18.101.218:27012/admin?replicaSet=rs0")

    path=[r'Bonds Convertible to Stocks.xls',]
    load_convbond_list_from_excel(path, client)

    from libRMI.wind import w
    # update_to_info(w, client, True)

    date = datetime(2018, 7, 18)

    # update_to_daily(date, w, client)
    # update_to_underlying(date, w, client)

    patch_to_daily(datetime.min, date, w, client, upsert=True, only_new=True)
    patch_to_underlying(datetime.min, date, w, client, upsert=True)

    # patch_to_underlying(datetime(2017, 9, 5), date, w, client, upsert=True)
