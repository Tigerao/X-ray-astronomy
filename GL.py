#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 20:37:05 2022

@author: tiger
Created on Sun Sep 20 18:13:40 2015
@author: Felix Darvas
@author: baotong
compute the Gregory-Laredo algorithm on arival times
This function computes the likelihood of a set of arrival times originating
from a periodic system rather than constant rate (poisson) background noise
based on
Gregory, P. C. and Thomas. J. Loredo, 1992,
"A New Method For The Detection Of A Periodic Signal Of Unknown Shape And Period"
in
The Astrophysical Journal, Astrophysical J., 398, p.146
inputs:
Tlist    -  list of arrival times, numpy int array
m_max    -  max number of bins typically 12-15, we use 12 as default
w_range  -  frequency range to scan numpy float array of frequency values
           default is  w_lo=20*pi/T at delta w = pi/T to w_hi=pi*N/T
           where N=#arrival times, T=observation time
ni       - number of integration steps, default ni=10
parallel - use parallel execution - default is off
outut:
O_period - Odds ratio for a periodic process vs. constant rate process
p_period - probability of a periodic process 0<=p_period<=1
m_opt    - optimal bin size 1<= m_opt <=m_max
S        - The probability spectrum
w        - The frequency range for S
w_peak   - the expected frequency of the process
w_conf   - 95% confidence interval of w_peak
"""
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import sys
sys.path.append("/Users/tiger/Desktop/chandra/code")
import GL_algorithm.funcs as funcs
import GL_algorithm.GL as GL

starttime = datetime.datetime.now()

def compute_bin(Tlist, m, w, fi):
    n = np.zeros(m, 'int')
    j = np.floor(m * np.mod(w * Tlist + fi, 2 * np.pi) / (2 * np.pi))
    j.astype(int)
    for u in range(0, m):
        n[u] = np.size(np.extract(j == u, j))
    return n

def compute_factor(N, m, v):
    # compute m^N *(N+m-1)! / (2pi*v*(m-1)!
    # which is used to scale the multiplicity function (see eqn.5.25 of the paper )
    # we return the log of the integral scale here to avoid numerical overflow
    f1 = N * np.log(m)  # log of m^N
    f2 = np.sum(np.log(np.arange(1, N + m)))  # log of (N+m-1)!
    #f3 = np.sum(np.log(np.arange(1, m + 1)))  # log of m!
    f3 = np.sum(np.log(np.arange(1, m )))
    #should be (m-1)! --Tong
    f = f1 + f3 - f2 - np.log(2 * np.pi * v)
    return f

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return T

def precompute_binmult(N):
    #return [lg(n!) for n in range(1,N)], 第一项忽略
    # precompute all potential bin factorials
    fbin = np.zeros(int(N) + 1)
    for i in range(2, int(N) + 1):  # n=0 -> define log(n)=0, n=1, log(n)=0, so no need to compute n=0,1
        fbin[i] = fbin[i - 1] + np.log(i)
    return fbin

def compute_S(epoch_info,Tlist,w,m,fi):
    n = compute_bin(Tlist, m, w, fi)
    tao=get_T_in_mbins(epoch_info,w,m,fi)
    s_wfi=tao/(sum(tao)/m)
    ln_S_wfi=-sum(n*np.log(s_wfi))
    #if no photons in some bin, none result returned
    return ln_S_wfi

def compute_W_scaled(Tlist, m, w, fi, factor, fbin,epoch_info):
    # compute the scaled multiplicity 1/Wm(w,fi) (see eqn. 5.25)
    # actually, it compute only n1!n2!...nm!,cause N! is been reduced in compute_factor --Tong
    # note that for large arrival time numbers the multiplicity becomes
    # excessively large. Since W_m(fi,w) is never needed explicitly,
    # we use the scaled version.
    # input factor is the log of the actual factor
    n = compute_bin(Tlist, m, w, fi)  # find bin histogram for a given bin number m, frequency w and fi
    f = 0
    for i in range(0, m):
        # f=f+np.sum(np.log(np.arange(2,n[i]+1)))
        f = f + fbin[n[i]]
    if epoch_info.ndim==2:
        ln_S_wfi=compute_S(epoch_info,Tlist,w,m,fi)
    elif epoch_info==None:ln_S_wfi=0

    y = np.exp(f + factor+ln_S_wfi)

    return y

def compute_Om1(w, Tlist, m, factor, fbin, ni,epoch_info):
    # compute  specific odds-ratios (eqn. 5.25)
    p = np.arange(0, ni) / float(ni) * 2 * np.pi / m
    #p 即为输入的fi
    # intgration range, only integrate over a single bin, as values of the integral repeat over bins
    y = np.zeros(np.size(p), 'float')  # array to hold values of W_scaled over the integration range
    for i in range(0, np.size(y)):
        y[i] = compute_W_scaled(Tlist, m, w, p[i], factor, fbin,epoch_info)
    return np.trapz(y, p) * m  # return intregrated W_Scaled


def compute_Om1wPar(Tlist, m_max, w, fa, fbin, ni,epoch_info):
    # compute odds-ratios for bins and frequencies
    # parallel version
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    pool = mp.Pool()  # use all workers
    for m in range(0, m_max):
        Om1w[m, :] = pool.map(functools.partial(compute_Om1, Tlist=Tlist, m=(m + 1), factor=fa[m], fbin=fbin, ni=ni,epoch_info=epoch_info), w)
    return Om1w


def compute_Om1w(Tlist, m_max, w, fa, fbin, ni,epoch_info):  # compute odds-ratios for bins and frequencies
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    for m in range(0, m_max):
        for wi in range(0, np.size(w)):
            Om1w[m, wi] = compute_Om1(w[wi], Tlist, m + 1, fa[m], fbin, ni,epoch_info)
    return Om1w


def compute_GL(Tlist,epoch_info, m_max=20, w_range=None, ni=10, parallel=False):
    # initialize output values
    O_period = None
    p_period = None
    m_opt = None
    S = None
    w = None
    w_peak = None
    w_mean = None
    w_conf = None
    N = float(np.size(Tlist))  # need float value to avoid int/int
    if N > 0:
        # compute GL algorithm
        fbin = precompute_binmult(N)
        v = m_max - 1
        T = float(np.max(Tlist)-np.min(Tlist))  # duration of the observation
        if w_range is None:  # use default frequencies
            w_hi = 5*np.pi * N / T  # max default frequency

            w_lo = np.minimum(20, N / 10) * np.pi / T  # min default frequency
            dw = np.pi / T / 10  # step size
            w = np.arange(w_lo, w_hi, dw)
            if np.size(w) < 2:
                print
                "error "
                raise ValueError('bad arrival time list')
        else:  # use user supplied frequency vector
            w = w_range
            w_hi = np.max(w_range)
            w_lo = np.min(w_range)
            dw=w[2]-w[1]
        if w_lo == 0:
            # cannot have w_lo =0
            print("minimum frequency cannot be 0!\n")
            return

        fa = np.zeros(m_max)
        for m in range(0, m_max):  # precompute factors for each m
            fa[m] = compute_factor(N, m + 1, v)
        if parallel:
            Om1w = compute_Om1wPar(Tlist, m_max, w, fa, fbin, ni,epoch_info)
        else:
            Om1w = compute_Om1w(Tlist, m_max, w, fa, fbin, ni,epoch_info)
        pw = 1./ (w*np.log(w_hi / w_lo))
        O1m = np.zeros(m_max)
        for i in range(0, m_max):  # intgreate odd ratios across frequencies
            O1m[i] = np.trapz(pw * Om1w[i], w)

        m_opt = np.argmax(O1m) # find optimum bin number, i.e largest odds-ratio
        S = Om1w[m_opt] / w  # compute Spectral probability
        m_opt = m_opt + 1  # start bin index with 1
        C = np.trapz(S, w)  # compute normalization
        S = S / C  # normalized probability
        O_period = np.sum(O1m[1:])  # integrated odds ratio
        p_period = O_period / (1 + O_period)  # likelihood of periodic event
        cdf = np.array(S)
        # plt.plot(w,cdf)
        # plt.show()
        for i in range(0, np.size(S)):
            cdf[i] = np.trapz(S[0:i], w[0:i])

        wr = np.extract(np.logical_and(cdf > .0025, cdf < .9975), w)
        w_peak = w[np.argmax(S)]
        w_mean = np.trapz(S * w, w)
        if np.size(wr) > 0:
            #print('exist')
            w_conf = [np.min(wr), np.max(wr)]
            if np.min(wr)==np.max(wr):
                print('asymmetric')
                w_conf=[np.min(wr)-dw, np.max(wr)+dw]
        else:
            w_conf = [w_peak-dw, w_peak+dw]
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf,cdf
    else:
        # throw an error
        print("No valid arrival time array provided!\n")
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf


def get_T_in_mbins(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    t_start = epoch_info[:, 0];
    t_end = epoch_info[:, 1]

    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)

    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

def filter_obs(src_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    return src_evt_use

def result(dataname,src_evt,epoch_info,w_range=2 * np.pi * np.arange(1/5000., 1 / 1000., 1e-6),filter=False):
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    if filter:src_evt_use = funcs.filter_obs(src_evt,epoch_info[:,-2])
    else: src_evt_use=src_evt

    time=src_evt_use[:,0]
    energy=src_evt_use[:,1]
    counts=len(time)
    starttime = datetime.datetime.now()
    GL_R=compute_GL(time,epoch_info=epoch_info,w_range=w_range,m_max=12,parallel=False)
    endtime = datetime.datetime.now()
    srcid=dataname
    runtime=(endtime - starttime).seconds
    Prob=GL_R[1]
    wpeak=GL_R[5]
    wmean=GL_R[6]
    mopt=GL_R[2]
    wconf_lo=GL_R[7][0]
    wconf_hi=GL_R[7][1]
    result=[int(srcid),runtime,Prob,2*np.pi/wpeak,wpeak,wmean,mopt,wconf_lo,wconf_hi,counts]
    #print(result)
    # get result directly with src_evt and epoch_info
    return result


def get_result(data_file,epoch_file,w_range,dataname='1',if_filter=False):
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    src_evt=np.loadtxt(data_file)
    if if_filter:
        CR = funcs.count_rate(src_evt=src_evt, epoch_info=epoch_info)
        plt.close()
        choseobs = [12112]
        (useid, epoch_info_use) = funcs.choose_obs(epoch_info,flux_info=CR,
                                                flux_filter=0,expT_filter=0,
                                                if_flux_high=True, if_expT_high=True,obsID=choseobs)
        if len(useid)==0:
            return [int(dataname), 0, 0, 0, 0, 0, 0, 0, 0, 0]
        epoch_info = epoch_info_use
        src_evt_use = funcs.filter_obs(src_evt, useid)
        src_evt = src_evt_use
        # select the epoch wanted
    
    counts=len(src_evt)
    if counts <200:
         return [int(dataname), 0, 0, 0, 0, 0, 0, 0, 0, counts]
    time=src_evt[:,0]

    starttime = datetime.datetime.now()
    GL_R=compute_GL(time,epoch_info=epoch_info,w_range=w_range,m_max=20,parallel=True)
    endtime = datetime.datetime.now()
    srcid=dataname
    runtime=(endtime - starttime).seconds
    Prob=GL_R[1]
    wpeak=GL_R[5]
    wmean=GL_R[6]
    mopt=GL_R[2]
    wconf_lo=GL_R[7][0]
    wconf_hi=GL_R[7][1]
    # get result from file name
    return [int(srcid),runtime,Prob,2*np.pi/wpeak,wpeak,wmean,mopt,wconf_lo,wconf_hi,counts]


def get_result_fromid(id_range,path,path_result,w_range=2 * np.pi * np.arange(1/10000., 1 / 5000., 1e-7)):
    w_range = 2 * np.pi * np.arange(1 / 100000., 1 / 50000., 2e-9) # 5000
    # w_range = 2 * np.pi * np.arange(1 / 50000., 1 / 20000., 2e-9) # 6000
    # w_range = 2 * np.pi * np.arange(1 / 20000., 1 / 10000., 1e-8) # 5000
    # w_range=2 * np.pi * np.arange(1/10000., 1 / 5000., 2e-8) # 5000
    # w_range=2 * np.pi * np.arange(1/5000., 1 / 1000., 1e-7) # 8000
    # w_range = 2 * np.pi * np.arange(1 / 1000., 1 / 500., 2e-7) # 5000     
    # w_range = 2 * np.pi * np.arange(1 / 500., 1 / 100., 1e-6) # 8000
    # w_range = 2 * np.pi * np.arange(1 / 500., 1 / 200., 5e-7) # 6000

    for id in id_range :
        id=int(id)
        dataname=str(id)
        data_file = path + f'txt_all_obs_p90/{dataname}.txt'
        epoch_file = path + f'txt_all_obs_p90/epoch_src_{dataname}.txt'

        result = get_result(data_file, epoch_file,w_range, dataname=dataname,if_filter=False)
        print(result)
        result_new=np.reshape(np.array(result),(1,10))
        np.savetxt(path_result+'result_100000_50000_{0}.txt'.format(dataname), result_new,
                   fmt='%10s %10.5f %10.5f %10.16f %10.16f %10.16f %10d %10.16f %10.16f %10d')
        
if __name__ == '__main__':

    # id_range = np.arange(1, 324)


    # path='/data/home/tiger/chandra/M81/txt_HETG/'
    # path_out=path+'result/'

    # obs_ID_all=[6174, 6346, 6347, 5601, 5600, 6892, 6893, 6894, 6895, 6896, 6897, 6898, 6899, 6900, 6901, 20624] # use all obervations for period searching below 5ks
    
    # path='/data/home/tiger/chandra/M81/txt_ACISI/'
    # path_out=path+'result/'
    # obs_ID_all=[18052,18054,18875,18047,18817,19685,18049,19688,18048,19981,19982,18053,18051,18050,19993,19992,19991,21384]
    # get_result_fromid(id_range, path, path_out)
    id_range=[436]
    path='/Volumes/SSDISK/chandra/CentaurusA_ACISI/merge_data/txt/'
    path_out=path+'result/'
    get_result_fromid(id_range, path, path_out)