#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import linecache
import sys
sys.path.append("/Users/tiger/Desktop/chandra/code")
# import GL_algorithm.funcs as funcs
from GL_algorithm import funcs

font1 = {'weight': 'normal',
         'size': 20, }
font2 = {'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=18 )

def get_LS(time, flux, freq, outpath=None, outname=None, save=False, show=True):
    x = time
    y = flux
    LS = LombScargle(x, y, normalization='standard')
    power = LS.power(freq)
    max_NormLSP = np.max(power)
    FP = LS.false_alarm_probability(power.max(), minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')
    FP_99 = LS.false_alarm_level(0.0027, minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    fig1=plt.figure(1,figsize=(12, 7.5))
    plt.plot(freq, power)
    plt.semilogx()
    out_period = 1. / freq[np.where(power == np.max(power))][0]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0], FP_99, '1-FAP 99.73%', font1)
    plt.text(freq[0], FP_95, '95%', font1)
    plt.text(freq[0], FP_68, '68%', font1)
    plt.xlabel('Frequency (1/s)', font1)
    plt.ylabel('Normalized LS Periodogram', font1)
    plt.tick_params(labelsize=16)
    if save:plt.savefig(outpath + outname + '_LS.pdf',dpi=1200)
    if show:plt.show()
    plt.close()
    return [FP, out_period, max_NormLSP]

def main_process():
    dataname='93'
    bin_len = 500
    path = '/Volumes/SSDISK/chandra/CentaurusA_ACISI/merge_data/txt/txt_all_obs_p90/'
    data_file = path + f'{dataname}.txt'
    epoch_file = path + f'epoch_src_{dataname}.txt'
    bkg_file = path + f'{dataname}_bkg.txt'
    epoch_info = np.loadtxt(epoch_file) 
    if epoch_info.ndim == 1:
        epoch_info = np.array([epoch_info])
    src_evt = np.loadtxt(data_file)
    bkg_evt = np.loadtxt(bkg_file)

    # obsid=[18052,18054,18875,18047,18817,19685,18049,19688,18048,19981,19982,18053,18051,18050,19993,19992,19991]
    # obsid=[18047]
    choose=0
    if choose:
        CR=funcs.count_rate(epoch_info=epoch_info,src_evt=src_evt)
        (useid, epoch_info_use) = funcs.choose_obs(epoch_info, flux_info=CR,
                                                        flux_filter=1e-2, expT_filter=0,          
                                                        if_flux_high=1, if_expT_high=1, obsID=None)
        [src_evt_use,bkg_evt_use] = funcs.filter_obs(src_evt, useid, bkg_evt)
    else:
        epoch_info_use=epoch_info
        src_evt_use=src_evt
    print(len(src_evt_use))
    # print(useid)
    
    path_fig = '/Volumes/SSDISK/chandra/CentaurusA_ACISI/merge_data/fig/'
    time = src_evt_use[:, 0]
    lc=funcs.get_hist(time,len_bin=bin_len,tstart=epoch_info[:,0][0],tstop=epoch_info[:,1][-1])
    T_tot=epoch_info_use[:,1][-1]-epoch_info_use[:,0][0]
    print(T_tot)
    # freq = np.arange(1 / 100000, 1 / 30000, 1/T_tot)
    freq = np.arange(1 / 100000, 1 / 30000, 1e-8)
    (FP, out_period, max_NormLSP)=get_LS(lc.time,lc.counts,freq=freq,outpath=path_fig, outname=str(dataname)+'_1',save=1,show=1)
    print('Period=',format(out_period))

main_process()
