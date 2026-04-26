#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:48:10 2022

@author: tiger
"""

#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits
import sys,os
import math
from tkinter import _flatten
from astropy.wcs import WCS
import rocket.make_region as reg_func
import rocket.position as pos
def delete_photon_ID(time, energy, ID,e_low,e_high):
    i = 0
    while i < len(energy):
        if energy[i] > e_high or energy[i] < e_low:
            energy = np.delete(energy, i)
            time = np.delete(time, i)
            ID = np.delete(ID, i)
            i = i - 1
        i = i + 1
    return [time, energy, ID]

def where_region(x,y,reg):
    r=np.array((x-reg[0],y-reg[1]))
    len_r=np.sqrt(r[0]**2+r[1]**2)
    temp=len_r-reg[2]
    return np.where(temp<=0)

def trans_radec2xy(imagename,ra,dec):
    w = WCS(imagename)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    binx,biny=np.shape(fits.open(imagename)[0].data)
    phy_x = src_x + 4096-binx/2
    phy_y = src_y + 4096-biny/2
    return (phy_x,phy_y)

def make_fk5region_each_obs(path_in,path_out,ra,dec,wcsimage,obs_ID_all,srcid=None):
    if srcid is None:srcid=np.arange(1,len(ra)+1,1)
    os.chdir(path_in)
    fitsname = wcsimage
    w = WCS(path_in + fitsname)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    binx,biny=np.shape(fits.open(path_in+fitsname)[0].data)
    print(src_x,src_y)
    print(ra,dec)

    out_index=np.union1d(np.where(src_x>binx-1),np.where(src_y>biny-1))
    src_x[out_index]=binx-1;src_y[out_index]=biny-1
    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)

    for i in range(len(obs_ID_all)):
        os.chdir(path_out)
        if not os.path.exists('region_{0}'.format(obs_ID_all[i])):
            os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        # p90_list='reproj_psf90_{0}_b50_300.fits'.format(obs_ID_all[i])
        p90_list='reproj_psf90_{0}_300_8000.fits'.format(obs_ID_all[i])
        hdul_p90 = fits.open(path_in + p90_list)
        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        print(obs_ID_all[i])
        reg_func.make_fk5_reg(srcid=srcid,ra=ra,dec=dec,psfradii=src_radius,outpath=
                              path_out+'region_{0}/'.format(obs_ID_all[i]))

    return None

def make_region_each_obs(path_in,path_out,ra,dec,wcsimage,obs_ID_all,bkg=1,ecf=90,srcid=None,multiple_src=1,single_name=0,single_srcradius=0):
    if srcid is None:srcid=np.arange(1,len(ra)+1,1)
    os.chdir(path_in)
    fitsname = wcsimage
    w = WCS(path_in + fitsname)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    binx,biny=np.shape(fits.open(path_in+fitsname)[0].data)
    print(src_x,src_y)
    print(ra,dec)
    phy_x = src_x + 4096-binx/2
    phy_y = src_y + 4096-biny/2 #acis
    # phy_x = src_x + 16384-binx/2
    # phy_y = src_y + 16384-biny/2 #hrc
    # phy_x = 8*(src_x) + 8384.5
    # phy_y = 8*(src_y) + 8384.5 #hrc_new bin2000
    print(phy_x,phy_y)
    out_index=np.union1d(np.where(src_x>binx-1),np.where(src_y>biny-1))
    src_x[out_index]=binx-1;src_y[out_index]=biny-1

    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)


    for i in range(len(obs_ID_all)):
        os.chdir(path_out)
        if not os.path.exists('region_{0}'.format(obs_ID_all[i])):
            os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        #p90_list='reproj_psf{0}_{1}_b4.fits'.format(ecf,obs_ID_all[i])
        # p90_list='reproj_psf90_{0}_b50_300.fits'.format(obs_ID_all[i])
        p90_list='reproj_psf{0}_{1}_300_8000.fits'.format(ecf,obs_ID_all[i])
        hdul_p90 = fits.open(path_in + p90_list)
        
        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        
        if single_srcradius:src_radius=np.array([single_srcradius])
        src_radius =src_radius* 2.032521 #acis
        # src_radius =src_radius*7.5872534 #hrc
        
        os.chdir(path_out+'region_{0}'.format(obs_ID_all[i]))
        os.system('mkdir region_{0}'.format(ecf))
        if multiple_src:singlename=None
        elif single_name:singlename=single_name
        else:
            print("Must specify multiple or single source")
            return None

        reg_func.make_phy_reg(srcid=srcid,x=phy_x,y=phy_y,psfradii=src_radius,
                              outpath=path_out+'region_{0}/'.format(obs_ID_all[i])+'region_{0}/'.format(ecf),
                              coordtype='physical',singlename=singlename)
        if bkg:
            srcid_list, x_list, y_list, psf_list = srcid, phy_x, phy_y, src_radius

            for n in range(srcid_list.size):
                f = open(path_out + 'region_{}/region_{}/{}_bkg.reg'.format(obs_ID_all[i],ecf,srcid_list[n]), 'w+')
                x = x_list[n]
                y = y_list[n]
                psf = psf_list[n]
                distance_square = (x_list - x) ** 2 + (y_list - y) ** 2
                idx_list = np.where((distance_square - (4 * psf + psf_list) ** 2) < 0)[0]
                f.write('annulus({},{},{},{})\n'.format(x, y, psf * 2, psf * 4))
                for idx in idx_list:
                    f.write('-circle({:.6f},{:.6f},{:.6f})\n'.format(x_list[idx], y_list[idx], 1.5 * psf_list[idx]))
                f.close()

    return None
# path_in,path_in_reg,path_out,reg_name,obs_id,ecf=90,suffix=''

def get_txt(path_in,path_in_reg,path_out,srcid_list,obs_id,ecf=90,suffix=''):
    ##对单个obs_id##
    evt_list='all_bcc_{0}_reproj_evt.fits'.format(obs_id)
    hdul_evt= fits.open(path_in+evt_list)
    x=hdul_evt[1].data['x']
    y=hdul_evt[1].data['y']
    energy=hdul_evt[1].data['energy']
    # pi=hdul_evt[1].data['pi']
    time=hdul_evt[1].data['time']
    obs_ID=np.array([obs_id for i in range(len(time))])
    os.chdir(path_out)
    if not os.path.exists('txt_{0}{1}'.format(obs_id,suffix)):
        print('making dir')
        os.system('mkdir txt_{0}{1}'.format(obs_id,suffix))

    for item in srcid_list:
        reg = reg_func.read_region(path_in_reg+'region_{0}/'.format(obs_id)+'region_{0}/'.format(ecf)+str(item)+'.reg')
        src_index=where_region(x,y,reg)
        src_x=x[src_index]
        src_y=y[src_index]
        src_t=time[src_index]
        src_E=energy[src_index]
        src_ID=obs_ID[src_index]

        # [src_t,src_E,src_ID]=delete_photon_ID(src_t,src_E,src_ID,e_low=50,e_high=300)
        [src_t,src_E,src_ID]=delete_photon_ID(src_t,src_E,src_ID,e_low=500,e_high=8000)
        src_t=src_t.astype('float')
        src_E =src_E.astype('float')
        src_ID=src_ID.astype('int')
        src_txt=np.column_stack((src_t,src_E,src_ID))
        src_txt = src_txt[src_txt[:,0].argsort()]

        np.savetxt(path_out+'txt_{0}{1}/'.format(obs_id,suffix)+str(item)+'.txt',src_txt,fmt="%.7f  %5.3f  %d")

    return None


def extract_evtlist_bkg(path_in, path_in_reg, path_out, obs_id,srcid_list,ecf=90,suffix=''):
    print('processing ' + str(obs_id))
    for srcid in srcid_list:
        f = open(path_in_reg + 'region_{}/region_{}/{}_bkg.reg'.format(obs_id, ecf,srcid)).readlines()
        anu = f[0][8:-2]
        rad = np.array(f[1:])
        x, y, inner_radius, outer_radius = np.array(anu.split(',')).astype('float')
        info = np.zeros((rad.size, 3))
        for i in range(rad.size):
            info[i, :] = rad[i][8:-2].split(',')
        x1 = info.astype('float')[:, 0]
        y1 = info.astype('float')[:, 1]
        r1 = info.astype('float')[:, 2]
        x1 = np.append(x1, x)
        y1 = np.append(y1, y)
        r1 = np.append(r1, inner_radius)
        evtfile = path_in + 'all_bcc_{0}_reproj_evt.fits'.format(obs_id)
        hdu = fits.open(evtfile)
        arrival_time_list = hdu[1].data['time']
        x_list = hdu[1].data['x']
        y_list = hdu[1].data['y']
        # energy_list = hdu[1].data['PI']
        energy_list = hdu[1].data['energy']
        idx1 = np.where((x_list - x) ** 2 + (y_list - y) ** 2 < outer_radius ** 2)[0]
        arrival_time_list = arrival_time_list[idx1]
        x_list = x_list[idx1]
        y_list = y_list[idx1]
        energy_list = energy_list[idx1]
        idx = []
        for i in range(x_list.size):
            xxx = x_list[i]
            yyy = y_list[i]
            distant_square = (x1 - xxx) ** 2 + (y1 - yyy) ** 2 - r1 ** 2
            if (distant_square > 0).all():
                idx.append(i)
        idx = np.array(idx)
        if idx.size != 0:
            x_list = x_list[idx]
            y_list = y_list[idx]
            energy_list = energy_list[idx]
            arrival_time_list = arrival_time_list[idx]
            arrival_time_list=np.array(arrival_time_list);energy_list=np.array(energy_list)
            src_ID=np.array([obs_id for i in range(len(arrival_time_list))])
            # [arrival_time_list, energy_list, src_ID] = delete_photon_ID(arrival_time_list, energy_list, src_ID, e_low=50, e_high=300)
            [arrival_time_list, energy_list, src_ID] = delete_photon_ID(arrival_time_list, energy_list, src_ID, e_low=500, e_high=8000)
            f = open(path_out + 'txt_{}{}/{}_bkg.txt'.format(obs_id, suffix,srcid), 'w+')
            #print(energy_list.size,src_ID.size)
            for i in range(energy_list.size):
                f.write("{} {} {}\n".format(arrival_time_list[i], energy_list[i], src_ID[i]))
            f.close()
        else:
            f = open(path_out + 'txt_{}{}/{}_bkg.txt'.format(obs_id, suffix,srcid), 'w+')
            f.close()
    return None

def make_epoch_file(obsid,inpath,outpath,outname):
    TSTART=[]
    TSTOP=[]
    for id in obsid:
        evtfits=inpath+'all_bcc_{0}_reproj_evt.fits'.format(str(id))
        hdul=fits.open(evtfits)
        # TSTART.append(hdul[1].header['TSTART'])
        # TSTOP.append(hdul[1].header['TSTOP'])
        TSTART.append(hdul[1].data['time'][0])
        TSTOP.append(hdul[1].data['time'][-1])

    TSTART=np.array(TSTART)
    TSTOP=np.array(TSTOP)
    obs_id=np.array(obsid)
    T_exp=TSTOP-TSTART

    epoch=np.column_stack((TSTART,TSTOP,obs_id,T_exp))
    epoch = epoch[epoch[:,0].argsort()]
    np.savetxt(outpath+outname+'.txt',epoch,fmt='%15.2f %15.2f %10d %20.2f')

def merge_txt(src_id,epoch_info,inpath,outpath,suffix,bkg=1,outname='txt_all_obs_50_300'):
    ## 对单个点源 ##
    ## Watch out the name of directory  ##
    ## Must follow the style given by get_txt##
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    obs_tstart=epoch_info[:,0]
    obs_tstop = epoch_info[:,1]
    obs_ID_all=epoch_info[:,2]
    obs_expt=epoch_info[:,-1]

    obs_ID_all=obs_ID_all.astype(int)

    res_t=[];res_E=[];res_ID=[];res_bkg_t=[];res_bkg_E=[];res_bkg_ID=[]
    epoch_ID=[];epoch_start=[];epoch_stop=[];epoch_expt=[]

    for i in range(len(obs_ID_all)):
        res_temp=np.loadtxt(inpath+'txt_{0}{1}/'.format(obs_ID_all[i],suffix)+str(src_id)+'.txt')
        if len(res_temp)==0:continue
        if res_temp.ndim==1:res_temp=np.array(res_temp)
        else:
            epoch_ID.append(obs_ID_all[i])
            epoch_start.append(obs_tstart[i])
            epoch_stop.append(obs_tstop[i])
            epoch_expt.append(obs_expt[i])

            res_t.append(list(res_temp[:, 0]))
            res_E.append(list(res_temp[:, 1]))
            res_ID.append(list((res_temp[:, 2])))
        if bkg:
            res_bkg_temp=np.loadtxt(inpath+'txt_{0}{1}/'.format(obs_ID_all[i],suffix)+str(src_id)+'_bkg.txt')
            if len(res_bkg_temp) == 0: continue
            if res_bkg_temp.ndim == 1: res_bkg_temp = np.array([res_bkg_temp])
            res_bkg_t.append(list(res_bkg_temp[:, 0]))
            res_bkg_E.append(list(res_bkg_temp[:, 1]))
            res_bkg_ID.append(list((res_bkg_temp[:, 2])))

    res_t = list(_flatten(res_t))
    res_E = list(_flatten(res_E))
    res_ID = list(_flatten(res_ID))
    result = np.column_stack((res_t, res_E, res_ID))
    if bkg:
        res_bkg_t = list(_flatten(res_bkg_t))
        #print(res_bkg_t)
        res_bkg_E = list(_flatten(res_bkg_E))
        res_bkg_ID = list(_flatten(res_bkg_ID))
        result_bkg= np.column_stack((res_bkg_t, res_bkg_E, res_bkg_ID))
        result_bkg=result_bkg.squeeze()
    epoch_info = np.column_stack((epoch_start, epoch_stop, epoch_ID, epoch_expt))
    if not os.path.exists(outpath+'{0}/'.format(outname)):
        os.mkdir(outpath+'{0}/'.format(outname))
    np.savetxt(outpath + '{0}/'.format(outname) + 'epoch_src_' + str(src_id) + '.txt', epoch_info,
               fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(outpath + '{0}/'.format(outname) + str(src_id) + '.txt', result, fmt="%.7f  %5.3f  %d")
    if bkg:
        #print('bkg:',result_bkg)
        np.savetxt(outpath + '{0}/'.format(outname) + str(src_id) + '_bkg.txt', result_bkg)
        #, fmt="%.7f  %5.3f  %d")