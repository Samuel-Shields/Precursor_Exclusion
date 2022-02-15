#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #UV-PUFF_v1.py
    #Samuel W. Shields
    #2022/02/15
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#%% package import
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks
#from scipy.interpolate import interp1d
from pyteomics import mzml
from csv import writer
import time

#excecution time
start_time = time.time()

#%%
#----------------------------------User defined variables----------------------------------#
#output file name
output_name = 'AA_BLE_PCs_PExUVPD_5E5agc_128msmax'

#peak picking params
height_cutoff = 1000 #height threshold (0 for full profile data)
prom = 1000
    # mass tolerances
ppm1 = 14.99 #ppm error for Fatty acids and premz
ppm2 = 14.99 #ppm error for DB frags
ppm3 = 9.99 #ppm error for precursor !!!! Need to fix error when no premz is found 

#phospholipids to search
PA = 0
PC = 1
PE = 1
PG = 0
PI = 0
PS = 0
SM = 0

#adducts
    #neutral
mass = 0
    #POS mode
MpH = 1
MpNH4 = 0
MpNa = 0
    #NEG mode
MmH = 0
MmHpCHO2 = 0


#iterate through scans
scans_it = 0

#%% preamble
#find which GPLS and which adducts to search from YES/NO table
Lclass = [PA, PC, PE, PG, PI, PS, SM]
searchL = ['PAs', 'PCs', 'PEs', 'PGs', 'PIs', 'PSs', 'SMs']

Adduct = [mass, MpH, MpNH4, MpNa, MmH, MmHpCHO2]
searchA = ['mass','MpH', 'MpNH4', 'MpNa', 'MmH', 'MmHpCHO2']

idxC = [i for i, x in enumerate(Lclass) if x == 1]
idxA = [i for i, x in enumerate(Adduct) if x == 1]

#import DB table
DB = np.genfromtxt('z_DB_NLs.csv', dtype= None, delimiter=',', encoding='UTF-8')
DBname = DB[0]
#read CSV and change into array
FA = np.genfromtxt('z_fattyacidloss.csv', dtype= None, delimiter=',', encoding='UTF-8')
FAname = FA[0]
# mass of water
h2o = 18.01056

hg = np.array([183.0660, 141.0191]).astype('float')

#newcsv
with open ('%s.csv' %(output_name), 'a+', newline='') as file:
    output = writer(file, dialect='excel')
    
#%% functions

def main():
    names = list_file()
    for z in names:
        fileName = z
        spectral_data, idx_scan = read_mzml(fileName)
        
        for e in idx_scan:
            premz, a, b = process_mzml(spectral_data, e, fileName)
            premz_t, found = premz_search(premz)
            mz, ints, noiseLevel = peaks(premz, a, b)
            find_HG(premz_t, mz, ints)
            find_FA(premz_t, mz, ints)
            find_DB(premz_t, mz, ints)
    return

def list_file():
    files = os.listdir()
    names=[]
    for name in files:
        if name.endswith(".mzML"):
            names += [name]       
    return names

def read_mzml(name):
    reader = mzml.read(name)
    spectral_data = [s for s in reader]
    reader.close()
    idx_scan = np.arange(0, int(spectral_data[-1]['index'])+1, 1)
    return spectral_data, idx_scan

def process_mzml(spectral_data, e, fileName):    
    if spectral_data[e]['ms level'] == 2:
        premz = spectral_data[e]['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
    else:
        premz = spectral_data[e]['precursorList']['precursor'][1]['isolationWindow']['isolation window target m/z']
    
    x = spectral_data[e]['m/z array']
    y = spectral_data[e]['intensity array']
    
    # cubic spline interpolation
    #a = np.linspace(x.min(), x.max(), len(x)*2)
    #f = interp1d(x, y, kind='cubic', bounds_error=False)
    #b = f(a)
    
    scan = spectral_data[e]['id'][-6:]
    Rt = round(float(spectral_data[e]['scanList']['scan'][0]['scan start time']), 3)
    
    with open ('%s.csv' %(output_name), 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel')
        output.writerow([fileName, 'exp_premz', premz, scan, 'RT', Rt])
        output.writerow(['found precursor', '----'])
    return premz, x, y

def premz_search(premz):
    found = pd.DataFrame(columns=['ID', 'mz', 'ppm', 'adduct'])
    for e in idxC:
        #GPL table variable
        cl = searchL[e]
        #read in the appropriate theoretical table
        t = pd.read_csv('%s.csv' %(cl))
        # loop for matching GPLS in one table
        for el in idxA:
            #which adduct
            adds = searchA[el]
            # define ppm
            ppm = ppm3/1E6
            #find list of matching precursors
            y = t[(t[adds] > premz-premz*ppm) & (t[adds] < premz+premz*ppm)]
            
            # dataframe formatting and appending
            y = y[['ID', adds]]
            y = y.rename(columns={adds: 'mz'})
            y['ppm'] = (premz - y.mz)/y.mz*1E6
            y['adduct'] = adds
            found = found.append(y)
    #list of matched precursors
    found = found.sort_values('ppm', ascending=True).reset_index().drop('index', axis=1)
    # first match with lowers ppm error
    premz_t = found.mz[0]            
    #convert to nparray
    found = found.to_numpy()
    #write to csv
    with open ('%s.csv' %(output_name), 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel')
        output.writerows(found)     
    return premz_t, found 

def peaks(premz, a, b):
    if height_cutoff > 0:
        noiseLevel = height_cutoff
    else:
        noiseLevel = (np.std((b[np.where(a > premz + 10)]))*3)*3
    
    peaks, _ = find_peaks(b, height = noiseLevel, prominence=prom)
    mz = a[peaks]
    ints = b[peaks]
    
    return mz, ints, noiseLevel

def find_FA(premz_t, mz, ints):
    FAmz = premz_t - FA[1].astype('float')   
    with open ('%s.csv' %(output_name), 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel')
        output.writerow(['--', '--', '--', '--'])
        output.writerow(['FA ID', 'exp mz', 'ppm', 'Int', '', 'exp mz', 'ppm', 'Int'])
    
        for e in mz:
                ppm = e*(ppm1/1E6)
                i = np.where((e-ppm < FAmz) & (e+ppm > FAmz) & (e < premz_t-2))
                if len(FAname[i]) > 0: #write to x only if peak is found
                   idx = np.where(e == mz)
                   w = str(FAname[i]).strip("',[,]"), str(round(e,4)).strip("',[,]"),\
                            str(np.around((((e-FAmz[i])/FAmz[i])*1E6), decimals=2)).strip("',[,]"), \
                            str(np.around(ints[idx], decimals=1)).strip("',[,]")
                   ppmk = (e+h2o)*(ppm1/1E6)
                   ik = np.where((e+h2o-ppmk < mz) & (e+h2o+ppmk > mz))
                   x = str('ketene'), str(np.around(mz[ik], decimals=4)).strip("',[,]"),\
                                str(np.around((((mz[ik]-(FAmz[i]+h2o))/(FAmz[i]+h2o))*1E6), decimals=2)).strip("',[,]"),\
                                str(np.around(ints[ik], decimals=1)).strip("',[,]")
                   if len(mz[ik]) > 0:
                       output.writerows([w+x])
    return

def find_DB(premz_t, mz, ints):
    DBlow = premz_t - DB[1].astype('float')
    with open ('%s.csv' %(output_name), 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel')
        output.writerow(['--', '--', '--', '--'])
        output.writerow(['DB ID', 'DB1', 'DB2', 'ppm_DB1', 'ppm_DB2', 'Int_DB1', 'Int_DB2'])
        
        for e in mz:
           ppm = e*(ppm2/1E6) 
           i = np.where((e-ppm < DBlow) & (e+ppm > DBlow) & (e < premz_t-2))
           if len(DBname[i]) > 0: #go to next step only if e matched in DBlow
               ppmDB = (e+24)*(ppm2/1E6)
               #look for the other double bond fragment
               j = np.where((e+24-ppmDB < mz) & (e+24+ppmDB > mz))
               if len(mz[j])>0: # go to next step only if a DB2 is matched
                   #find index where e is in peak list to get intensity 
                   idx = np.where(e == mz)
                   
                   output.writerow([str(DBname[i]).strip("',[]").replace('"', ''), str(round(e,4)).strip("',[,]"), 
                       str(np.around(mz[j], decimals=4)).strip("',[,]"), 
                       str(np.around((((e-DBlow[i])/DBlow[i])*1E6), decimals=2)).strip("',[,]"),
                       str(np.around((((mz[j]-(DBlow[i]+24))/(DBlow[i]+24))*1E6), decimals=2)).strip("',[,]"),
                       str(np.around(ints[idx], decimals=1)).strip("',[,]"),
                       str(np.around(ints[j], decimals=1)).strip("',[,]")
                       ])
        output.writerow(['--', '--', '--', '--'])
        output.writerow(['--', '--', '--', '--'])
    return

def find_HG(premz_t, mz, ints):
    hg_l = premz_t-hg
    with open ('%s.csv' %(output_name), 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel')
        output.writerow(['--', '--', '--', '--'])
        output.writerow(['HG_loss', 'mz', 'ppm', 'Int'])
    
        for e in hg_l:
            ppm = e*(ppm1/1E6)
            n = np.where((mz > e-ppm) & (mz < e+ppm))
            if len(n) > 0:
                err = np.around(((mz[n]-e)/e)*1E6, 1)
                output.writerow([str(premz_t-e).strip("',[,]"), str(mz[n]).strip("',[,]"),
                                 str(err).strip("',[,]"), str(ints[n]).strip("',[,]")])
    return

main()

#prints excecution time
print("--- %s seconds ---" % (time.time() - start_time))



