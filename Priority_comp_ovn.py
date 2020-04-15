import requests
import pandas as pd
import matplotlib.pyplot as plt
import ast
import numpy as np
from  more_itertools import unique_everseen
import difflib
import datetime
import numpy as np
from astropy.time import Time
now = datetime.datetime.now().isoformat()
t = Time(now)
w_time = Time('2020-02-26', format = 'iso')
priority_old = pd.read_csv("/Users/bhagyasubrayan/Desktop/REFITT/022620_priority_old.csv", sep=",",
                       names = ['Locus_id','Band','Magnitude','Mag_error'])
priority_new = pd.read_csv("/Users/bhagyasubrayan/Desktop/REFITT/priority_new.csv", sep=",",
                       names = ['Locus_id','Band','Magnitude','Mag_error'])
def duplicates(x):
     x = x.drop_duplicates(subset = ['Locus_id'], keep = 'first')
     x = x.reset_index(drop = True)
     return list(x['Locus_id'])
po_locus = duplicates(priority_old)
pn_locus = duplicates(priority_new)
def con_floa(x):
    new_lis = [float(i) for i in x]
    return new_lis
lis_names = {}
def choose(x,y):
    if (len(x) != 0 and len(y) == 0):
        return x
    elif(len(x)== 0 and len(y) != 0):
        return y
    else:
        return x
def priority(x):
    for m in range(0,len(x)):
        url = 'https://antares.noao.edu/api/alerts/?locus_id=' + str(x[m])
        r = requests.get(url, allow_redirects=True)
        open( '/Users/bhagyasubrayan/Desktop/REFITT/'+ str(x[m])+'.txt', 'wb').write(r.content)
        data = []
        with open('/Users/bhagyasubrayan/Desktop/REFITT/'+str(x[m])+".txt", "r") as inFile:
            data = ast.literal_eval(inFile.read())
        a = data['result'][0]
        #print(a)
        #object_id = a['properties']['ztf_object_id']
        locus_id = a['locus_id']
        #print(locus_id)
        ra = a['ra']
        dec = a['dec']
        #print('Locus_id:',locus_id)
        #print('RA:',ra)
        #print('Dec:',dec)
        data_g = pd.DataFrame(columns= ['Date','Magnitude','Mag_error','Band'])
        data_r = pd.DataFrame(columns= ['Date','Magnitude','Mag_error','Band'])
        date_g =[]
        date_r =[]
        magn_g =[]
        magn_r = []
        band_g = []
        band_r = []
        mag_err_g =[]
        mag_err_r = []
        for i in range(len(data['result'])):
            val = data['result'][i]['properties']
            if 'ztf_object_id' in a:
                object_id_m = a['properties']['ztf_object_id']
            elif 'ztf_object_id' in val:
                object_id_m = val['ztf_object_id']
            #val = data['result'][i]['properties']
            if ('ztf_magpsf'in val and 'passband' in val):
                #print(data['result'][i])
                if(val['passband'] == 'g'):
                    mjd_g = data['result'][i]['mjd']
                    mag_g = val['ztf_magpsf']
                    date_g.append(mjd_g)
                    magn_g.append(mag_g)
                    band_g.append(val['passband'])
                    mag_err_g.append(val['ztf_sigmapsf'])
                else:
                    mjd_r = data['result'][i]['mjd']
                    mag_r = val['ztf_magpsf']
                    date_r.append(mjd_r)
                    magn_r.append(mag_r)
                    band_r.append(val['passband'])
                    mag_err_r.append(val['ztf_sigmapsf'])
        #print('Object_id:',object_id_m)
        magn_g = con_floa(magn_g)
        magn_r = con_floa(magn_r)
        non_empty = choose(magn_g,magn_r)
        if(min(non_empty) < 19):
            lis_names.update({locus_id :[object_id_m,min(non_empty),m+1]})
        #print(magn_g)
        #print(locus_id)
        mag_err_g = con_floa(mag_err_g)
        mag_err_r = con_floa(mag_err_r)
        data_g['Date'] = date_g
        data_g['Magnitude'] = magn_g
        data_g['Band'] = band_g
        data_g['Mag_error'] = mag_err_g
        data_g.sort_values('Date',inplace =True)
        data_r['Date'] = date_r
        data_r['Magnitude'] = magn_r
        data_r['Band'] = band_r
        data_r['Mag_error'] = mag_err_r
        data_r.sort_values(['Date'],inplace =True)
        #print(data_g)
        #print(data_r)
        final = pd.concat([data_g,data_r]).reset_index(drop = True)
        plt.figure(figsize=(5,5))
        plt.axvline(w_time.mjd)
        plt.legend(['Feb 26 2020'])
        plt.errorbar(date_g, magn_g,mag_err_g,fmt= 'go')
        plt.errorbar(date_r,magn_r,mag_err_r,fmt = 'ro')
        plt.gca().invert_yaxis()
        #plt.xticks(np.arange(58000, 59500, 100.0))
        plt.title(object_id_m + ':Current mjd:'+ str(t.mjd))
        plt.xlabel('MJD')
        plt.ylabel('magpsf')
        plt.ylim(21.5,17.5)
        plt.xlim(58870,t.mjd)
        plt.show()
    return lis_names
old = {k: v for k, v in sorted(priority(po_locus).items(), key=lambda x: x[1][1])}
print('Starting with the new priority')
new = {k: v for k, v in sorted(priority(pn_locus).items(), key=lambda x: x[1][1])}
old
new
