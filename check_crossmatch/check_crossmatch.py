import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table, unique
from astropy.visualization import astropy_mpl_style
from astropy.visualization import simple_norm

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3

import csv
import matplotlib.ticker as ticker
from scipy.stats import kde
from math import pi
import pandas as pd
#import match 

from astropy.table import QTable
from astropy.table import hstack

# 2MASS
files = open('2mass/2mass_files.txt', 'r')
lst_fl = []
for line in files:
    archivos = line.strip()
    lst_fl.append(archivos)

data_number = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','15.0']
non_duplicated_data = []
total_data = []
for i in range(len(lst_fl)):
    dat_ = Table.read(lst_fl[i], format='csv')
    df_ = dat_.to_pandas()
    #print(f'total datos en {lst_fl[i]}: {len(df_)}')
    dfg = df_.drop_duplicates(subset='srcid', keep=False)
    #print(f'total unicos en{lst_fl[i]}: {len(dfg)}')
    non_duplicated_data.append(len(dfg))
    total_data.append(len(df_))
duplicated_data = [total_data[x] - non_duplicated_data[x] for x in range(len(total_data))]
crt_2mass = Table([data_number,total_data, non_duplicated_data, duplicated_data],names=('radius','Total data', 'Non duplicated', 'duplicated'))


z_2mass = Table.read('zang_x_2mass_galaxies.csv', format='csv')
t_2mass = z_2mass[z_2mass['Px']>0.7]
z_non_star = z_2mass[z_2mass['X']!=0]
t_non_star = t_2mass[t_2mass['X']!=0]

fig, (ax0, ax1) = plt.subplots(2,1, figsize=(10,7))
ax0.hist(z_2mass['Px'], bins=20,histtype='step', edgecolor="darkcyan", linewidth=1.5)
ax0.hist(z_non_star['Px'], bins=5,histtype='step', edgecolor="red", linewidth=1.5)
ax0.text(0.8, 0.8, f'Total={len(z_2mass)}', fontsize=13, color='darkcyan', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax0.text(0.8, 0.7, f'Non-star={len(z_non_star)}', fontsize=13, color='red', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax0.text(0.10, 0.85, f'2MASS', fontsize=14, color='black', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax0.set_yscale('log')

ax1.hist(t_2mass['Px'], bins=15,histtype='step', edgecolor="darkcyan", linewidth=1.5)
ax1.hist(t_non_star['Px'], bins=5,histtype='step', edgecolor="red", linewidth=1.5)
ax1.set_xlim(0.6,1)
ax1.text(0.8, 0.8, f'Total={len(t_2mass)}', fontsize=13, color='darkcyan',horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.8, 0.7, f'Non-star={len(t_non_star)}', fontsize=13, color='red', horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)
ax1.set_xlabel('Px', fontsize=13)
ax1.set_ylabel('n sources', fontsize=13)
ax0.set_ylabel('n sources', fontsize=13)
ax1.set_yscale('log')
plt.savefig('2mass_hist.jpg')

t_2mass.write('zang_galaxies_70qual.csv', format='csv')
t_non_star.write('2mass_zang_galaxies_70qual.csv', format='csv')


# GAIA EDR3
files = open('gaia/gaia_files.txt', 'r')
lst_fl = []
for line in files:
    archivos = line.strip()
    lst_fl.append(archivos)

non_duplicated_data = []
total_data = []
for i in range(len(lst_fl)):
    dat_ = Table.read(lst_fl[i], format='csv')
    df_ = dat_.to_pandas()
    #print(f'total datos en {lst_fl[i]}: {len(df_)}')
    dfg = df_.drop_duplicates(subset='srcid', keep=False)
    #print(f'total unicos en{lst_fl[i]}: {len(dfg)}')
    non_duplicated_data.append(len(dfg))
    total_data.append(len(df_))
duplicated_data = [total_data[x] - non_duplicated_data[x] for x in range(len(total_data))]
crt_gaia = Table([data_number,total_data, non_duplicated_data, duplicated_data], names=('radius','Total data', 'Non duplicated', 'duplicated'))

plt.figure(figsize=(8,5))
#plt.bar(crt_gaia['radius'],crt_gaia['Total data'], color='darkblue')
plt.bar(crt_gaia['radius'][:-1],crt_gaia['duplicated'][:-1], color='gold')
plt.xlabel('Radius separation (arcsec)', fontsize=15)
plt.ylabel('Zang x GAIA Duplicated Data', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=14, rotation=90)
plt.yticks(fontsize=14)
plt.savefig('duplicated_gaia.jpg')


# VVV DR2
files = open('vvv/vvv_files.txt', 'r')
lst_fl = []
for line in files:
    archivos = line.strip()
    lst_fl.append(archivos)

non_duplicated_data = []
total_data = []
for i in range(len(lst_fl)):
    dat_ = Table.read(lst_fl[i], format='csv')
    df_ = dat_.to_pandas()
    #print(f'total datos en {lst_fl[i]}: {len(df_)}')
    dfg = df_.drop_duplicates(subset='srcid', keep=False)
    #print(f'total unicos en{lst_fl[i]}: {len(dfg)}')
    non_duplicated_data.append(len(dfg))
    total_data.append(len(df_))
duplicated_data = [total_data[x] - non_duplicated_data[x] for x in range(len(total_data))]
crt_vvv = Table([data_number,total_data, non_duplicated_data, duplicated_data], names=('radius','Total data', 'Non duplicated', 'duplicated'))

plt.figure(figsize=(8,5))
#plt.bar(crt_vvv['radius'],crt_vvv['Total data'], color='darkblue')
plt.bar(crt_vvv['radius'][:-1],crt_vvv['duplicated'][:-1], color='lime')
plt.xlabel('Radius separation (arcsec)', fontsize=15)
plt.ylabel('Zang x VVV Duplicated Data', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=14, rotation=90)
plt.yticks(fontsize=14)
plt.savefig('duplicated_vvv.jpg')


# AllWISE
data_number = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5',
               '2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','10.0','15.0']

files = open('allwise/allwise_files.txt', 'r')
lst_fl = []
for line in files:
    archivos = line.strip()
    lst_fl.append(archivos)

non_duplicated_data = []
total_data = []
for i in range(len(lst_fl)):
    dat_ = Table.read(lst_fl[i], format='csv')
    df_ = dat_.to_pandas()
    #print(f'total datos en {lst_fl[i]}: {len(df_)}')
    dfg = df_.drop_duplicates(subset='srcid', keep=False)
    #print(f'total unicos en{lst_fl[i]}: {len(dfg)}')
    non_duplicated_data.append(len(dfg))
    total_data.append(len(df_))
duplicated_data = [total_data[x] - non_duplicated_data[x] for x in range(len(total_data))]
crt_allwise = Table([data_number,total_data, non_duplicated_data, duplicated_data],
                 names=('radius','Total data', 'Non duplicated', 'duplicated'))

z_allwise = Table.read('zang_x_2mass_galaxies.csv', format='csv')
t_allwise = z_allwise[z_allwise['Px']>0.7]
len(t_allwise)


# GLIMPSE
data_number = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5',
               '2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','15.0']

files = open('glimpse/glimpse_files.txt', 'r')
lst_fl = []
for line in files:
    archivos = line.strip()
    lst_fl.append(archivos)

non_duplicated_data = []
total_data = []
for i in range(len(lst_fl)):
    dat_ = Table.read(lst_fl[i], format='csv')
    df_ = dat_.to_pandas()
    #print(f'total datos en {lst_fl[i]}: {len(df_)}')
    dfg = df_.drop_duplicates(subset='srcid', keep=False)
    #print(f'total unicos en{lst_fl[i]}: {len(dfg)}')
    non_duplicated_data.append(len(dfg))
    total_data.append(len(df_))
duplicated_data = [total_data[x] - non_duplicated_data[x] for x in range(len(total_data))]
crt_glimpse = Table([data_number,total_data, non_duplicated_data, duplicated_data],
                 names=('radius','Total data', 'Non duplicated', 'duplicated'))


# PLOTS
d_2mass = Table.read('2mass/Zang_x_2MASS_15.0.csv', format='csv')
d_gaia = Table.read('gaia/Zang_x_GAIA_15.0.csv', format='csv')
d_vvv = Table.read('vvv/Zang_x_VVV_15.0.csv', format='csv')
d_glimpse = Table.read('glimpse/Zang_x_GLIMPSE_15.0.csv', format='csv')
d_sdss = Table.read('Zang_x_SDSS_DR12_15.0.csv', format='csv')
d_allwise = Table.read('allwise/Zang_x_AllWISE_15.0.csv', format='csv')
d_wise = Table.read('Zang_x_WISE_15.0.csv', format='csv')

plt.figure(figsize=(10,8))
#hx1, hy1, _ = plt.hist(d_2mass['angDist'], bins=20,cumulative=True, density=1, fill=False, edgecolor="blue")
hx1, hy1, _ = plt.hist(d_2mass['angDist'], bins=20,cumulative=True, density=1, histtype='step',
                       edgecolor="darkcyan", linewidth=1.3, label='2MASS')
hx2, hy2, _ = plt.hist(d_gaia['angDist'], bins=20,cumulative=True, density=1, histtype='step' ,
                       edgecolor="gold", linewidth=1.3, label='GAIA EDR3')
hx3, hy3, _ = plt.hist(d_vvv['angDist'], bins=20,cumulative=True, density=1, histtype='step' ,
                       edgecolor="lime", linewidth=1.3, label='VVV DR2')
hx4, hy4, _ = plt.hist(d_glimpse['angDist'], bins=20,cumulative=True, density=1, histtype='step' ,
                       edgecolor="darkturquoise", linewidth=1.3, label='GLIMPSE')
#hx5, hy5, _ = plt.hist(d_sdss['angDist'], bins=20,cumulative=True, density=1, histtype='step' ,
#                       edgecolor="fuchsia", linewidth=1.3, label='SDSS')
hx6, hy6, _ = plt.hist(d_allwise['angDist'], bins=20,cumulative=True, density=1, histtype='step' ,
                       edgecolor="blue", linewidth=1.3, label='AllWISE')

#plt.ylim(0.0,max(hx1)+0.05)
plt.ylabel('Cumulative Fraction', fontsize=16)
plt.xlabel('all Data separation (arcsec)', fontsize=16)
plt.tick_params(labelsize=13)
plt.legend(loc=(0.75,0.15))
#plt.set_y
plt.grid()
plt.savefig("all_data.jpg", bbox_inches='tight')
