import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.cm import ScalarMappable
from matplotlib.patheffects import withStroke
import matplotlib.ticker as ticker

import math
import csv
from math import pi
import pandas as pd
import glob
#import match 

import seaborn as sns
import statistics as stat

from scipy.stats import kde #for the plot very soft
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity

from astropy.io import fits
from astropy.table import Table, unique, join, vstack, QTable, hstack
from astropy.visualization import astropy_mpl_style
from astropy.visualization import simple_norm
from astropy.visualization import make_lupton_rgb
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy import constants as Cons
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Galactic, Angle

import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection

#from astroquery.vizier import Vizier
#from astroquery.xmatch import XMatch

plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3

import nilo.loaders as ld

#------------------------------------------------------------------------------
# TRANSFORM COORDINATES -------------------------------------------------------
def galactic(ra,dec):
    c = SkyCoord(ra*u.degree, dec*u.degree)
    l = c.galactic.l*u.degree
    b = c.galactic.b*u.degree
    return l, b

def equatorial(l,b):
    c = SkyCoord(l*u.degree, b*u.degree, frame='galactic')
    ra = c.fk5.ra
    dec = c.fk5.dec
    return ra, dec





def angle_coordinates(ra, dec):
    RA = []
    DEC = []
    ra_h = Angle(ra, u.deg).hms[0][i]
    ra_m = Angle(ra, u.deg).hms[1][i]
    ra_s = Angle(ra, u.deg).hms[2][i]
    RA.append(str(str(int(ra_h))+' '+str(int(ra_m))+' '+str(round(ra_s,4))))
    dec_d = Angle(dec, u.deg).dms[0][i]
    dec_m = abs(Angle(dec, u.deg).dms[1][i])
    dec_s = abs(Angle(dec, u.deg).dms[2][i])
    DEC.append(str(str(int(dec_d))+' '+str(int(dec_m))+' '+str(round(dec_s,4))))
    return RA, DEC


# CATALOGO DE ZHANG+21 ---------------------------------------------------------
zang = Table.read('zang_cat.csv', format='csv')
print('Zang+21 completo: ',+len(zang))
print('Galaxias en Zang+21 completo: ',+len(zang[zang['Class_x']=='GALAXY']))

## Hacer el corte entre 0<l<10
zang['l_LONG'], zang['b_LAT'] = galactic(zang['sc_ra'],zang['sc_dec'])

#disc_all = zang[zang['GAL_LONG']>250]#
disc_all = zang[zang['l_LONG']>295]#
disc_all = disc_all[disc_all['l_LONG']<350]
disc_all = disc_all[disc_all['b_LAT']<2]
disc_all = disc_all[disc_all['b_LAT']>-2]
print('Total objetos en el disco del VVV: ',+len(disc_all))

#disc_all.write('disc_all.csv', format='csv', overwrite=True)

disc = disc_all[disc_all['Class_x']=='GALAXY']
disc_QSO = disc_all[disc_all['Class_x']=='QSO']
disc_S = disc_all[disc_all['Class_x']=='STAR']

#disc.write('disc_gal.csv', format='csv')
#disc_QSO.write('disc_QSO.csv', format='csv')
#disc_S.write('disc_S.csv', format='csv')

print()
print('Total de galaxias en el disco del VVV: ',+len(disc))
print('Total de estrellas en el disco del VVV: ',+len(disc_S))
print('Total de QSO en el disco del VVV: ',+len(disc_QSO))
disc.write('zang_vvv_galaxies.csv', format='csv')
disc.write('zang_vvv_galaxies.cat', format='ascii')

all_cross = Table.read('new_cross_nirgc.cat', format='ascii')
#false = np.zeros(10)
print('GALX: ',+len(zang[zang['Class_x']=='GALAXY']))
print('QSOs: ',+len(zang[zang['Class_x']=='QSO']))
print('STAR: ',+len(zang[zang['Class_x']=='STAR']))

#### estadistica de todo zang?
z_ = zang.to_pandas()

print('Total de datos con WISE: ',+ len(z_[~z_['Pxi'].isna()]))
print('Total de datos con SDSS: ',+ len(z_[~z_['Pxo'].isna()]))
print('Total de datos con WISE+SDSS: ',+ len(z_[~z_['Pxio'].isna()]))
dzc = disc_all.to_pandas()
print()
print('Total de disc con WISE: ',+ len(dzc[~dzc['Pxi'].isna()]))
print('Total de disc con SDSS: ',+ len(dzc[~dzc['Pxo'].isna()]))
print('Total de disc con WISE+SDSS: ',+ len(dzc[~dzc['Pxio'].isna()]))

#rad = Table([radi], names=['size'])
ext = Table([disc_all['sc_ra'], disc_all['sc_dec']])
radi = np.ones(len(ext))*0.1
ext.add_column(radi)
ext.rename_column('col2', 'size')
#ext.write('disc_all_forEXT.cat', format='ascii')

extin = Table.read('ext_disc_all.tbl', format='ascii')
disc_df = disc.to_pandas()
print('Total de galaxias con información en el IR',+len(disc_df[~disc_df['Pxi'].isna()]))

gal_zang = Table.read('zang_x_NIRGC_1.3.csv', format='csv')
gal = Table.read('zang_x_NIRGC_gals.csv', format='csv')


# CRUCE CON NIRGC Y ZHANG+21 -------- 1.3 arcsecs -------------------
all_dat = Table.read('zhang_disc_gal_cross_VVV_1.3arcsec.csv', format='csv')
all_dat['l_LONG'], all_dat['b_LAT'] = galactic(all_dat['RA_SEARCH'],all_dat['DEC_SEARCH'])
df = all_dat.to_pandas()
print(f'total datos con 1.3arcsec de separacion : {len(df)}')
dfg = df.drop_duplicates(subset='ID_SEARCH')
print(f'Total de objetos únicos en "dat": {len(dfg)}')
nan = dfg[dfg['ALPHA_J2000'].isna()]
t_nan = Table.from_pandas(nan)
#t_nan.write('nan_rows_1.3.csv', format='csv', overwrite=True)
print('Total objetos "NaN": ',+ len(nan))
no_nan = dfg[~dfg['ALPHA_J2000'].isna()]
t_no_nan = Table.from_pandas(no_nan)
print('Total objetos sin valores "NaN": ',+ len(no_nan))
print()
star = dfg[dfg['EXTENDED']==0]
extended = dfg[dfg['EXTENDED']==1]

print('Cantidad de fuentes puntuales: ',+len(star))
print('Cantidad de fuentes extendidas: ',+len(extended)-4)

fig, (ax0, ax1, ax2) = plt.subplots(1,3, sharey=True, figsize=(20,6))
ax0.scatter(dfg['SPREAD_MODEL'], dfg['MAG_AUTO'], color='black', s=7)
ax0.scatter(gal['SPREAD_MODEL'], gal['MAG_AUTO_Ks'], color='red', s=14, marker='D')
ax1.scatter(dfg['CLASS_STAR'], dfg['MAG_AUTO'], color='black', s=7)
ax1.scatter(gal['CLASS_STAR'], gal['MAG_AUTO'], color='red', s=14, marker='D')
ax2.scatter(dfg['ELLIPTICITY'], dfg['MAG_AUTO'], color='black', s=7)
ax2.scatter(gal['ELLIPTICITY_1'], gal['MAG_AUTO'], color='red', s=14, marker='D')
ax0.set_xlim(-0.025,0.025)
ax0.set_ylim(9,19)
ax0.invert_yaxis()
ax0.set_ylabel(r'MAG AUTO', size=15)
ax0.set_xlabel(r'SPREAD MODEL', size=15)
ax0.minorticks_on()
ax0.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax0.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
ax0.axvline(x=0, c='gray', linestyle='--')
ax1.set_xlabel(r'CLASS STAR', size=15)
ax1.minorticks_on()
ax1.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax1.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
ax2.set_xlabel(r'ELIPTICITY', size=15)
ax2.minorticks_on()
ax2.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax2.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
plt.tight_layout()


# SIN CONTRAPARTE REVISED ----------------------------------------------
_crowded = new_nan[new_nan['GENERAL']=='a']
_cstar = new_nan[new_nan['GENERAL']=='b']
_empty = new_nan[new_nan['GENERAL']=='c']
_saturated = new_nan[new_nan['GENERAL']=='d']
_group = new_nan[new_nan['GENERAL']=='e']
_spikes = new_nan[new_nan['GENERAL']=='f']
_SFR = new_nan[new_nan['GENERAL']=='g']
_PUVS = new_nan[new_nan['GENERAL']=='h']
_assoc = new_nan[new_nan['GENERAL']=='j']
_PUVSFR = new_nan[new_nan['GENERAL']=='I']
_error = new_nan[new_nan['GENERAL']=='error']

_spikes_empty = _spikes[_spikes['GENERAL']=='c']
_saturated_empty = _saturated[_saturated['GENERAL']=='c']


print('Total sources with normal crowdded region: ',len(_crowded))
print('Total sources with a central star: ',len(_cstar))
print('Total sources with empty central region: ',len(_empty))
print('Total sources with Saturated Star: ',len(_saturated))
print('Total sources with Stellar group: ',len(_group))
print('Total sources with Spikes: ',len(_spikes))
print('Total sources with Normal SFR: ',len(_SFR))
print('Total sources with PUVS: ',len(_PUVS))
print('Total sources with Massive star association: ',len(_assoc))
print('Total sources with PUVSFR: ',len(_PUVSFR))
print()
print('Total sources with Spikes+empty region: ',len(_spikes_empty))
print('Total sources with Saturated+empty region: ',len(_saturated_empty))

### Porcentajes
_p_crowded = ((len(_crowded)+1)/len(new_nan))*100
_p_cstar = (len(_cstar)/len(new_nan))*100
_p_empty = (len(_empty)/len(new_nan))*100
_p_saturated = (len(_saturated)/len(new_nan))*100
_p_group = (len(_group)/len(new_nan))*100
_p_spikes = (len(_spikes)/len(new_nan))*100
_p_SFR = (len(_SFR)/len(new_nan))*100
_p_PUVS = (len(_PUVS)/len(new_nan))*100
_p_assoc = (len(_assoc)/len(new_nan))*100
_p_PUVSFR = (len(_PUVSFR)/len(new_nan))*100

_p_spikes_empty = len(_spikes_empty)/len(new_nan)
_p_saturated_empty = len(_saturated_empty)/len(new_nan)

print('Total sources with normal crowdded region: ''%.2f'%_p_crowded)
print('Total sources with a central star: ''%.2f'%_p_cstar)
print('Total sources with empty central region: ''%.2f'%_p_empty)
print('Total sources with Saturated Star: ''%.2f'%_p_saturated)
print('Total sources with Stellar group: ''%.2f'%_p_group)
print('Total sources with Spikes: ''%.2f'%_p_spikes)
print('Total sources with Normal SFR: ''%.2f'%_p_SFR)
print('Total sources with PUVS: ''%.2f'%_p_PUVS)
print('Total sources with Massive star association: ''%.2f'%_p_assoc)
print('Total sources with PUVSFR: ''%.2f'%_p_PUVSFR)


print('Total sources with Spikes+empty region: ''%.3f'%_p_spikes_empty)
print('Total sources with Saturated+empty region: ''%.3f'%_p_saturated_empty)

nnnnn = new_nan.to_pandas()
pd.unique(nnnnn['GENERAL'])
T['ID_n'] = new_nan['ID_n']
nn_nan = join(T,new_nan, keys='ID_n')
zang['ID_SEARCH'] = zang['srcid']
_nan_ = join(nn_nan,zang, keys=['ID_SEARCH'])
_df_nan_ = _nan_.to_pandas()

nan_rev_gals = _nan_[_nan_['Galaxy']==1]
nan_revisar = Table([nan_rev_gals['RA_SEARCH'], nan_rev_gals['DEC_SEARCH']])
nan_revisar.write('nan_revisar.cat', format='ascii', overwrite=True)

# EJEMPLOS PARA TESIS ------------------------------------------------------------
ejeee = _nan_[_nan_['Notes']=='ejemplo']

ejemplos = Table(ejeee[0])
ejemplos.remove_row(0)

ejemplos.add_row(ejeee[ejeee['GENERAL']=='a'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='b'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='c'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='d'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='e'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='f'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='g'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='h'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='I'][0])
ejemplos.add_row(ejeee[ejeee['GENERAL']=='j'][0])
ejemplos.add_row(_nan_[4446])
ejemplos.add_row(_nan_[4650])

#ejeee[ejeee['GENERAL']=='I']

ejemplos_coord = Table([ejemplos['RA_SEARCH'], ejemplos['DEC_SEARCH']])
ejemplos_coord.write('ejemplos_coord.cat', format='ascii', overwrite=True)

# ---> HIGH, LOW
te_ex = Table(_nan_[0])
te_ex.remove_row(0)
te_ex.add_row(_nan_[518-1])   #NORMAL CROWDED HIGH
te_ex.add_row(_nan_[5488-1])
te_ex.add_row(_nan_[1813-1])  #STAR HIGH
te_ex.add_row(_nan_[3571-1])
te_ex.add_row(_nan_[655-1])  #EMPTY CENTRAL HIGH
te_ex.add_row(_nan_[2025-1])
te_ex.add_row(_nan_[2927-1]) #STAR ASSOC HIGH
te_ex.add_row(_nan_[2646-1])
te_ex.add_row(_nan_[4484-1]) #SATURATED HIGH
te_ex.add_row(_nan_[4601-1])
te_ex.add_row(_nan_[2555-1]) #SPIKES HIGH
te_ex.add_row(_nan_[149-1])
te_ex.add_row(_nan_[5194-1]) #NORMAL SFR HIGH
te_ex.add_row(_nan_[2930-1])
te_ex.add_row(_nan_[2869-1]) #PUVS HIGH
te_ex.add_row(_nan_[1958-1])
te_ex.add_row(_nan_[2189-1]) #MASSIVE ASSOC HIGH
te_ex.add_row(_nan_[2188-1])
te_ex.add_row(_nan_[1361-1]) #PUVSFR HIGH
te_ex.add_row(_nan_[520-1])
te_ex.add_row(_nan_[4446]) #PROBLEM
te_ex.add_row(_nan_[4650]) #PROBLEM

te_ex_coord = Table([te_ex['RA_SEARCH'], te_ex['DEC_SEARCH']])
te_ex_coord.write('ejemplos_coord_tesis.cat', format='ascii')

names = ['1_normal_crowded','2_star','3_empty_central','4_saturated','5_stellar_group','6_spikes','7_normal_sfr',
        '8_PUVS','9_massive_star','10_PUVSFR','11_special_1','12_special_2']

data_path = '/home/pmc/Documents/astro/astro.master/tesis/Zang_etal/ejemplos/fits/12/'

files = glob.glob(data_path+'/*.fits.gz')
files.sort()

filetrs = ['Z','Y','J','H','Ks']
wavelengths = [0.87*u.um, 1.02*u.um, 1.25*u.um, 1.64*u.um, 2.14*u.um]

fig = plt.figure(figsize=(20, 5))

for i,file in zip(range(len(files)),files):
    hdu = fits.open(file)[1]
    wcs = WCS(hdu.header)

    ax = fig.add_subplot(1,5, i+1)
    norm = simple_norm(hdu.data, 'asinh', min_percent=0.5, max_percent=99.5)
    ax.imshow(hdu.data, origin='lower', norm=norm, cmap=plt.cm.gray_r, interpolation=None)
    ax.set_xticks([ ])
    ax.set_yticks([ ])
    ax.set_xlabel(' ')
    ax.set_title(r'${}\,\,({})$'.format(filetrs[i], wavelengths[i]), size=23)
    ax.plot(hdu.data.shape[1] // 2, hdu.data.shape[0] // 2, color='red', marker='+', markersize=16, mew=1)

fig.text(-1.6,-0.15, r'$RA,DEC:\,{},\,{}\, ;\,\,\quad P_x:{}$'.format(ejemplos['RA_SEARCH'][11],ejemplos['DEC_SEARCH'][11],ejemplos['Px'][11]),va='bottom', ha='center', transform=ax.transAxes, color='k', fontsize=23)

fig.tight_layout()
fname = 'ejemplos/fits/{}.jpg'.format(names[11])
#plt.savefig(fname, bbox_inches='tight')

names = ['1_normal_crowded_H','1_normal_crowded_L','2_star_H','2_star_L','3_empty_central_H','3_empty_central_L',
         '4_stellar_group_H','4_stellar_group_L','5_saturated_H','5_saturated_L','6_spikes_H','6_spikes_L',
         '7_normal_sfr_H','7_normal_sfr_L','8_PUVS_H','8_PUVS_L','9_massive_star_H','9_massive_star_L',
         '10_PUVSFR_H','10_PUVSFR_L','11_special_1','12_special_2']


data_path = '/home/pmc/Documents/astro/astro.master/tesis/Zang_etal/ejemplos/tesis/1_'

files = glob.glob(data_path+'/*.fits.gz')
files.sort()

filetrs = ['Z','Y','J','H','Ks']
wavelengths = [0.87*u.um, 1.02*u.um, 1.25*u.um, 1.64*u.um, 2.14*u.um]

fig = plt.figure(figsize=(20, 5))

for i,file in zip(range(len(files)),files):
    hdu = fits.open(file)[1]
    wcs = WCS(hdu.header)

    ax = fig.add_subplot(1,5, i+1)
    norm = simple_norm(hdu.data, 'asinh', min_percent=0.5, max_percent=99.5)
    ax.imshow(hdu.data, origin='lower', norm=norm, cmap=plt.cm.gray_r, interpolation=None)
    ax.set_xticks([ ])
    ax.set_yticks([ ])
    ax.set_xlabel(' ')
    ax.set_title(r'${}\,\,({})$'.format(filetrs[i], wavelengths[i]), size=23)
    ax.plot(hdu.data.shape[1] // 2, hdu.data.shape[0] // 2, color='red', marker='+', markersize=16, mew=1)

fig.text(-1.6,-0.15, r'$RA,DEC:\,{},\,{}\, ;\,\,\quad P_x:{}$'.format(te_ex['RA_SEARCH'][21],te_ex['DEC_SEARCH'][21],te_ex['Px'][21])
         ,va='bottom', ha='center', transform=ax.transAxes, color='k', fontsize=23)

fig.tight_layout()
fname = 'ejemplos/tesis/{}.jpg'.format(names[21])
#plt.savefig(fname, bbox_inches='tight')

# 4XMM COMPLETE
xmm = Table.read('4XMM-DR9/dr9-4xmm.fit', format='fits')
xmm['srcid'] = xmm['Source']
xmm_nan = Table.read('_nan_flux.csv', format='csv')
xmm_no_nan_ = Table.read('xmm_no_nan.csv', format='csv')
xmm_no_nan = join(xmm_no_nan_, disc, keys=['srcid'])
vphas_gals = join(xmm, vphasp, keys=['srcid'])

gal_xmm = join(xmm, disc, keys=['srcid'])
print('Total de fuentes: ', len(xmm_nan))
fg_ = xmm_nan[np.logical_or(xmm_nan['Flux4']>0,xmm_nan['Flux5']>0)]
print('Total de fuentes con emisión sobre los > 2.0 keV (Hard Xray):', len(fg_))
print('Corresponde a un porcentaje de: ', len(fg_)*100/len(xmm_nan))

print('Total de fuentes sin emisión en Soft Xray:', len( fg_[np.logical_or(fg_['Flux1']==0,fg_['Flux2']==0)]))

print('Total de fuentes: ', len(gal_xmm))
fg = gal_xmm[np.logical_or(gal_xmm['Flux4']>0, gal_xmm['Flux5']>0)]
print('Total de fuentes con emisión sobre los > 2.0 keV (Hard Xray):', len(fg))
print('Corresponde a un porcentaje de: ', len(fg)*100/len(gal_xmm))

print('Total de fuentes sin emisión en Soft Xray:', len( fg[np.logical_or(fg['Flux1']==0,fg['Flux2']==0)] ))


# CORRELATIONS OF nIR, XRAY LUMINOSITY OF EARLY TYPE GALAXIES.....
