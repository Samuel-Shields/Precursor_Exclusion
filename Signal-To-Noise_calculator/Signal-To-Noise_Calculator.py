#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Signal-To-Noise_Calculator.py
    #Samuel W. Shields
    #2022/02/15
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.signal import find_peaks, peak_widths
from pyteomics import mzml
from csv import writer

#excecution time
start_time = time.time()


tPeaks = np.genfromtxt('z_targeted_fragment_list.csv', delimiter=',')

#user defined varables
ppms = 30
hgt = 10000
prom = 750

mzWindo = 1.5 #width of window (i.e. peak+/-mzWindo)
FWHMx = 3 #number of FWHM below and above max in window
sigma = 3 #number of standard deviations for noise level calculation

plot = 0 #plotting? yes=1 no=0
report = 1 #generate report? yes=1 no=0

#functions
def file_lst():
    files = os.listdir() 
    x = []
    files = os.listdir()
    for name in files:
        if name.endswith(".mzML"):
            x += [name]
    return x

def read_mzML():
    names = file_lst()
    datas = []    
    for z in names:
        fileName = z  
        reader = mzml.read(fileName)
        spectral_data = [s for s in reader]
        idx_scan = np.arange(0, int(spectral_data[-1]['index'])+1, 1)
        
        scans = []
        data = []
        for e in idx_scan:
            mz = spectral_data[e]['m/z array']
            inten = spectral_data[e]['intensity array']
            peak, _ = find_peaks(inten, height = hgt, prominence=prom)
            scans += [fileName, idx_scan, spectral_data[e]['id'][-6:], 
                      mz[peak], inten[peak], mz, inten] 
            data += [scans]
            scans = []

        reader.close()
        datas += [data]
    return datas

def find():
    names = file_lst()
    datas = read_mzML()
    y = np.arange(0, len(names), 1)
    line_all = []
    for e in y: #loop though files
        w = np.arange(0, len(datas[e][0][1]), 1)
        for el in w: #loop through scans
            mz_x = datas[e][el][3]
            ints_x = datas[e][el][4]
            name_x = datas[e][el][0]
            mz_full = datas[e][el][5]
            ints_full = datas[e][el][6]
                        
            
            for ele in tPeaks:
                ppm = ele*(ppms/1E6) #define ppm as (10/1E6)
                idx_f = np.where((ele-ppm < mz_x) & (ele+ppm > mz_x))
                
                mz_f = np.around(mz_x[idx_f],4)
                ints_f = np.around(ints_x[idx_f],0)
                
                if len(ints_f) > 0:
                    if len(ints_f) < 2:
                        y =  np.where((mz_f-(mzWindo/2) < mz_full) & 
                                      (mz_f+(mzWindo/2) > mz_full))
                    else:
                        a = np.where(ints_f == max(ints_f))
                        mz_f = mz_f[a]
                        ints_f = ints_f[a]
                        y =  np.where((mz_f-(mzWindo/2) < mz_full) & 
                                      (mz_f+(mzWindo/2) > mz_full))

                    ints_w = ints_full[y]
                    mz_w = mz_full[y]
                    
                    max_i = np.asarray(np.where(ints_w == max(ints_w))).flatten()
                    width_w = peak_widths(ints_w, max_i, rel_height=0.5)
                    mzStep = (max(mz_w)-min(mz_w)) / len(mz_w)
                    peakwidth = mzStep*width_w[0]
                    
                    noise_i = np.asarray(np.where((mz_w < mz_w[max_i] - FWHMx*peakwidth))).flatten()
                    topnoise_i, _ = find_peaks(ints_w[noise_i])
                    noise_j = np.asarray(np.where((mz_w > mz_w[max_i] + FWHMx*peakwidth))).flatten()                     
                    topnoise_j, _ = find_peaks(ints_w[noise_j])
                    
                    noise_w = np.concatenate((ints_w[noise_i[topnoise_i]], ints_w[noise_j[topnoise_j]]))
                    noise_mz_w = np.concatenate((mz_w[noise_i[topnoise_i]], mz_w[noise_j[topnoise_j]]))
                    noiseLevel = np.around(np.mean(noise_w)+np.std(noise_w)*sigma,0)
                    
                    SNR = np.around(ints_f/noiseLevel,1)
                    
                    err = np.around(((mz_f-ele)/ele)*1E6,2)

                    line = [str(name_x), 
                            str(datas[e][el][2]), 
                            float(mz_f), float(ints_f), 
                            float(err), float(SNR), 
                            noiseLevel, len(noise_w), float(np.around(peakwidth, 4))]
                    
                    line_all += [line] 
                    
                else:
                    line = [str(name_x), str(datas[e][el][2]), ele, 'none']
                    line_all += [line]
                    
                    if plot == 1:
                        plt.figure()
                        plt.title('$\it{m/z}$ ' + str(np.around(mz_w[max_i], 4)).strip("[]"))
                        plt.xlabel('$\it{m/z}$')
                        plt.ylabel('Intensity')
                        plt.plot(mz_w, ints_w, 'k')
                        plt.hlines(noiseLevel, xmin=min(mz_w), xmax=max(mz_w), color="C3")        
                        
                        plt.ylim(-(max(ints_w)*0.01), max(ints_w)/10)
                        plt.text(min(mz_w), max(ints_w)/11, "10x zoom")
                        plt.text(min(mz_w), noiseLevel*1.2, 'Noise = ' + str(round(noiseLevel)))
                        
                        plt.plot(noise_mz_w, noise_w, 'bo')
                        plt.plot(max(mz_w)-0.28, max(ints_w)/11, 'bo')
                        plt.text(max(mz_w)-0.25, max(ints_w)/11.2, "noise peak")
                        plt.text(max(mz_w)-0.3, noiseLevel*1.2, 'S/N = ' + str(SNR).strip("[]"))
    return line_all

run = find()

header = ['Filename', 'Scan', 'm/z', 'Intensity', 'ppm', 
          'S/N', 'noise', 'noise points']
params = ['ppm', float(ppms), 'height', float(hgt), 
          'prominence', float(prom), 'window', float(mzWindo), 
          'FWHMx', float(FWHMx), 'sigma', float(sigma)]

if report == 1:
    with open('z_results.csv', 'a+', newline='') as file:
        output = writer(file, delimiter=',', dialect='excel') 
        output.writerow(params)
        output.writerow(header)
        output.writerows(run)

print("--- %s seconds ---" % (time.time() - start_time)) 

