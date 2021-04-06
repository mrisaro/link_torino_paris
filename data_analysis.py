#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 12:02:44 2021
Traces analysis.
@author: guest1
"""
#%% Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os
#%% Function definitions
def fun_pks_beat(f,sp):
    fh = 36e6
    sh = 2*fh
    beat = 150e6
    
    kk_fh = np.where((f>fh-5e6)&(f<fh+5e6))
    pk_fh = np.max(sp[kk_fh[0]])
    
    kk_sh = np.where((f>sh-5e6)&(f<sh+5e6))
    pk_sh = np.max(sp[kk_sh[0]])
    
    kk_beat = np.where((f>beat-5e6)&(f<beat+5e6))
    pk_beat = np.max(sp[kk_beat[0]])
    
    return np.array([pk_fh, pk_sh, pk_beat])

def annot_max(xmax,ymax, ax=None):
    text= "Pk:{:.1f}dBm".format(ymax)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="data",
              arrowprops=arrowprops, bbox=bbox_props, ha="left", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(xmax+.5,ymax+10), **kw)
#%% find files in folder

path = 'data_susa/'
files = os.listdir(path)    
files = list(filter(lambda f: f.endswith('.csv'), files))

#%%
data = []

for f in files:    
    my_data = np.genfromtxt(path+f, delimiter=',',skip_header=45)
    data.append(my_data)

props = dict(boxstyle='round', facecolor='white', alpha=0.6)
    
#%% Plots of Susa 2/4/21

aa = 6
title = 'EDFA Susa Port 1 + Filters'
pks = fun_pks_beat(data[aa][:,0],data[aa][:,1])
pk_3a = fun_pks_beat(data[1][:,0],data[1][:,1])
pk_3b = fun_pks_beat(data[1][:,0],data[1][:,1])

fig = plt.figure(1,figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(data[aa][:,0]*1e-6,data[aa][:,1],label='Susa+EDFA')
ax.plot(data[1][:,0]*1e-6,data[1][:,1],label='Susa Ref')
ax.set_xlim(-5,190) 
ax.set_ylim(-72,2)
ax.set_xlabel(r'Freq (MHz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title(title,fontsize=14)
annot_max(150, round(pks[2],1))
annot_max(150, round(pk_3[2],1))
ax.legend(fontsize=14,loc=2)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%% Some plots
props = dict(boxstyle='round', facecolor='white', alpha=0.6)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(data[0][:,0],data[0][:,1])
ax.set_ylim(-92,2) 
ax.set_xlabel(r'Freq (Hz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title('Trace reference F3',fontsize=14)

ax.text(0.10, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.30, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.80, 0.85, 'Pk: -6 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%%
fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(data[1][:,0],data[1][:,1]) 
ax.set_xlabel(r'Freq (Hz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title('First Measurement F1',fontsize=14)

ax.text(0.10, 0.85, 'Pk: -25 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.30, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.80, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%%
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(data[2][:,0],data[2][:,1]) 
ax.set_xlabel(r'Freq (Hz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title('First Measurement F3',fontsize=14)

ax.text(0.10, 0.85, 'Pk: -25 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.30, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.80, 0.85, 'Pk: -24 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%%
fig = plt.figure(4)
ax = fig.add_subplot(111)
ax.plot(data[3][:,0],data[3][:,1])
ax.set_ylim(-92,2) 
ax.set_xlabel(r'Freq (Hz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title('EDFA Measurement F1',fontsize=14)

ax.text(0.10, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.30, 0.85, 'Pk: -20 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.80, 0.85, 'Pk: -6 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%%
fig = plt.figure(5)
ax = fig.add_subplot(111)
ax.plot(data[4][:,0],data[4][:,1])
ax.set_ylim(-92,2) 
ax.set_xlabel(r'Freq (Hz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.set_title('EDFA Measurement F3',fontsize=14)

ax.text(0.10, 0.85, 'Pk: -26 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.30, 0.85, 'Pk: -20 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.text(0.80, 0.85, 'Pk: -6 dBm', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top', 
        bbox=props)
ax.grid(which='both',linestyle='--')
fig.tight_layout()

#%%
fig = plt.figure(3,figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(data[0][:,0]*1e-6,data[0][:,1],label='BPF-310')
ax.plot(data[1][:,0]*1e-6,data[1][:,1],label='BPF-328')
ax.set_xlabel(r'Freq (MHz)', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)
ax.legend(fontsize=14)    
ax.set_title('BPF Mini-circuits',fontsize=14)
ax.grid(which='both',linestyle='--')
fig.tight_layout()