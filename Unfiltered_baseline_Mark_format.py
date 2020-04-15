import requests
import pandas as pd
import matplotlib.pyplot as plt
import ast
import numpy as np
from  more_itertools import unique_everseen
import difflib
from matplotlib.backends.backend_pdf import PdfPages
#day = input('Enter the date(yy-mm-dd):')
baseline = pd.read_csv("/Users/bhagyasubrayan/Desktop/REFITT/Unfiltered/20-03-19/targets.baseline.2020-03-19", sep=",",
                       names = ['Locus','RA','Dec','Un1','Un2','Band','Un3','Un4','Un5','Un6'])
priority = pd.read_csv("/Users/bhagyasubrayan/Desktop/REFITT/Unfiltered/20-03-19/targets.priority.2020-03-19", sep=",",
                       names = ['Locus','RA','Dec','Un1','Un2','Band','Un3','Un4','Un5','Un6'])
baseline.sort_values('Locus',inplace =True)
baseline.drop_duplicates(subset = ['Locus'], keep = 'first')
priority.sort_values('Locus',inplace =True)
priority = priority.drop_duplicates(subset = ['Locus'], keep = 'first')
priority = priority.reset_index(drop = True)
pr_locus = list(priority['Locus'])
for i in range(len(pr_locus)):
    baseline.drop(baseline.loc[baseline['Locus']==pr_locus[i]].index, inplace=True)
baseline = baseline.reset_index(drop = True)
baseline_locus = baseline['Locus']
def con_floa(x):
    new_lis = [float(i) for i in x]
    return new_lis
lis_names = {}
def min_mag(x):
    for r in range(len(x)):
        if min(x) < 19:
            lis_names.update({locus_id :[object_id_m,min(x)]})
            break
    return lis_names
for m in range(0,len(baseline_locus)):
    url = 'https://antares.noao.edu/api/alerts/?locus_id=' + str(baseline_locus[m])
    r = requests.get(url, allow_redirects=True)
    open( '/Users/bhagyasubrayan/Desktop/REFITT/Unfiltered/20-03-17/'+ str(baseline_locus[m])+'.txt', 'wb').write(r.content)
    data = []
    with open('/Users/bhagyasubrayan/Desktop/REFITT/Unfiltered/20-03-17/'+str(baseline_locus[m])+".txt", "r") as inFile:
        data = ast.literal_eval(inFile.read())
    a = data['result'][0]
    #print(a)
    #object_id = a['properties']['ztf_object_id']
    locus_id = a['locus_id']
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
    lis_names = min_mag(magn_g)
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
    if (max(final['Date'])- min(final['Date']) < 50) and (len(final) > 5) and (final['Magnitude'].mean( ) <= 18.5) :
        print('Locus_id:',locus_id)
        print('Object_id:',object_id_m)
        print('RA:',ra)
        print('Dec:',dec)
        print(final)
        plt.figure(figsize=(5,5))
        plt.errorbar(date_g, magn_g,mag_err_g,fmt= 'go')
        plt.errorbar(date_r,magn_r,mag_err_r,fmt = 'ro')
        plt.gca().invert_yaxis()
        #plt.xticks(np.arange(58000, 59500, 100.0))
        plt.title(object_id_m)
        plt.xlabel('MJD')
        plt.ylabel('magpsf')
        #plt.yticks(np.arange(22,15,0.5))
        #pp = PdfPages( '/Users/bhagyasubrayan/Desktop/REFITT/Unfiltered/'+ day + '/'+ day+'.pdf')
        #pp.savefig(f)
        plt.show()
    else:
        continue
