#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 16:13:09 2022

@author: tiger
"""

#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket


path_in='/Volumes/WDDISK/chandra/M31HRC_new/merge_data/xdata/'
path_out='/Volumes/WDDISK/chandra/M31HRC_new/merge_data/timing/'

obs_ID_all=[267,    268,    269,    270,     271,   272,    273,    275,    276,    277,    278,
            1569,   1570,  1912,    2904,   2905,   2906,   5925,   5926,   5927, 5928,
            6177,   6202,  7283,    7284,   7285,   7286,   8526,   8527, 8528, 8529, 8530,
            9825,   9826,  9827,    9828,   9829,
            10683,  10684,  10838,  10882,  10883,  10884,  10885,  10886,
            11808,  11809,  12110,  12111,  12112,  12113,  12114, 13178, 13179, 13180,
		    13227,  13228,  13229,  13230,  13231,  13278,  13279,  13280, 13281]

# path_in='/media/tiger/WDDISK/graduate/M31ACIS/merge_data/xdata/'
# #path_in='/media/tiger/TIGERDISK/graduate/M31star/merge_data/timing/'
# path_out_reg='/media/tiger/WDDISK/graduate/M31ACIS/merge_data/timing/reg/'
# path_out_txt='/media/tiger/WDDISK/graduate/M31ACIS/merge_data/timing/txt/'
# obs_ID_all=[    10551,  11277,  12163,  13302,  14196,  15326,  16298,  308,   4719,  7140,  9520,
# 10552,  11278,  12164,  13825,  14197,  15327,  17443,  309  , 4720  ,8183 , 9521,
# 10553,  11279,  12970,  13826,  14198,  15328,  17444,  310  , 4721  ,8184 , 9522,
# 10554,  11838,  12971,  13827,  14927,  1575 ,  17445,  311  , 4722  ,8185 , 9523,
# 10555,  11839,  12972,  13828,  14928,  1581 ,  17446,  312  , 4723  ,8186 , 9524,
# 10715,  11840,  12973,  13833,  14929,  1582 ,  17447,  4360 , 7064  ,8187 , 9529,
# 10716,  11841,  12974,  13834,  14930,  1583 ,  1854 ,  4678 , 7068  ,8191,
# 10717,  11842,  13298,  13835,  14931,  16294,  303  ,  4679 , 7136  ,8192,
# 10719,  12160,  13299,  13836,  15267,  16295,  305  ,  4680 , 7137  ,8193,
# 11275,  12161,  13300,  13837,  15324,  16296,  306  ,  4681 , 7138  ,8194,
# 11276,  12162,  13301,  14195,  15325,  16297,  307  ,  4682 , 7139  ,8195,
# ]

def input_srcinfo():
    cat=fits.open(path_in+'M31hrc_catalog_final.fits')
    ra = cat[1].data['RA'][0]
    dec = cat[1].data['DEC'][0]
    srcID_list=np.arange(1,len(ra)+1,1)
    return (srcID_list,ra,dec)
#%%
def main_process():
    makeregion=0
    makefk5=0
    gettxt=0
    mergedata=1
    if makeregion: 
        rocket.make_region_each_obs(path_in,path_out+'reg/',ra=ra,dec=dec,
                                    wcsimage=wcsimage,obs_ID_all=obs_ID_all,
                                    bkg=1,ecf=90,srcid=srcID_list,multiple_src=1)
        rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,
                               outpath=path_out,outname='M31HRC_epoch')
    if makefk5:
        rocket.make_fk5region_each_obs(path_in,path_out=path_out+'fk5reg/'
                                ,ra=ra,dec=dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,srcid=srcID_list)
    if gettxt:
        epoch_info=np.loadtxt(path_out+'txt/M31HRC_epoch.txt')
        for id in obs_ID_all:
            rocket.get_txt(path_in=path_in,path_in_reg=path_out+'reg/',path_out=path_out+'txt/',srcid_list=srcID_list,obs_id=id,ecf=90,suffix='_p90')
            rocket.extract_evtlist_bkg(path_in=path_in,path_in_reg=path_out+'reg/', path_out=path_out+'txt/',
                                       obs_id=id, srcid_list=srcID_list, 
                                       ecf=90,suffix='_p90')
    if mergedata:
        epoch_info=np.loadtxt(path_out+'txt/M31HRC_epoch.txt')
        for srcid in srcID_list:
            rocket.merge_txt(srcid,epoch_info,inpath=path_out+'txt/',
                             outpath=path_out+'txt/',bkg=1,suffix='_p90',outname='txt_all_obs_p90')

if __name__=='__main__':
    (srcID_list, ra, dec) = input_srcinfo()
    useid=obs_ID_all
    wcsimage ='M31HRC_100_10000.fits'
    main_process()
    # rocket.select_src_bypos(srcID_list,ra,dec,ra_c=ra_center,dec_c=dec_center,
    #                         inter_radius=inter_radius,outpath=path_out_txt+'txt_all_obs_p50/',
    #                         outname='inter_src_60arcsec.txt')