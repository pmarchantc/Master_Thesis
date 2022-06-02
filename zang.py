''' Zang Catalog'''
zang = Table.read('zang_cat.csv', format='csv')
print('Zang+21 completo: ',+len(zang))
print('Galaxis en Zang+21 completo: ',+len(zang[zang['Class_x']=='GALAXY']))

zang['l_LONG'], zang['b_LAT'] = galactic(zang['sc_ra'],zang['sc_dec'])
disc = zang[zang['l_LONG']>295]
disc = disc[disc['l_LONG']<350]
disc = disc[disc['b_LAT']<2]
disc = disc[disc['b_LAT']>-2]
disc = disc[disc['Class_x']=='GALAXY']
print('Total de galaxias en el disco del VVV: ',+len(disc))
disc.write('zang_vvv_galaxies.csv', format='csv')
disc.write('zang_vvv_galaxies.cat', format='ascii')

all_cross = Table.read('new_cross_nirgc.cat', format='ascii')
good_gal = disc[disc['Px']>0.9]
print('cantidad de "Galaxias" con Px > 0.9: ',+len(good_gal))

thin_disc = disc[disc['b_LAT']>-1]
thin_disc = thin_disc[thin_disc['b_LAT']<1]
print('Total de galaxias en el disco fino del VVV: ',+len(thin_disc))

gal_zang = Table.read('zang_x_NIRGC_1.3.csv', format='csv')
gal = Table.read('zang_x_NIRGC_gals.csv', format='csv')
print(gal_zang)

## CRUCE CON NIRGC Y Zhang21 -------------------------

print('Cruce con radio de 3arcsec y trabajar con columnas con valores NaN')
dat = Table.read('zhang_disc_gal_cross_VVV_3arcsec.csv', format='csv')
dt = dat.to_pandas()
print(f'total datos con 3arcsec de separacion : {len(dt)}')
dft = dt.drop_duplicates(subset='ID_SEARCH')
print(f'Total de objetos únicos en "dat": {len(dft)}')
### Guardar tablas NaN ###
nan = dft[dft['ALPHA_J2000'].isna()]
t_nan = Table.from_pandas(nan)
#t_nan.write('nan_rows_3.csv', format='csv', overwrite=True)
print('Total objetos "NaN": ',+ len(nan))
print('Total objetos sin valores "NaN": ',+ len(dft[~dft['ALPHA_J2000'].isna()]))
print()
s_ = dft[dft['EXTENDED']==0]
e_ = dft[dft['EXTENDED']==1]
print('Cantidad de fuentes puntuales: ',+len(s_))
print('Cantidad de fuentes extendidas: ',+len(e_))

### DEPURACION DE DATOS
dt_ = dt[~dt['ALPHA_J2000'].isnull()]    #mantener sólo solumnas con datos en RA
dt__ = dt_.drop_duplicates(subset='ID_SEARCH')  #eliminar duplicados
print('-Total con datos en todas las columnas y sin duplicados:  ',+len(dt__))
dt__ = dt__[dt__['DISTANCE']<=1.3]          #cortar en r<=1.3
print('-Datos con radio menor a 1.3, con datos en todas las columnas y sin duplicados:  ',+len(dt__))

print()
print()
print('Cruce con radio de 1.3 arcsec')
all_dat = Table.read('zhang_disc_gal_cross_VVV_1.3arcsec.csv', format='csv')
df = all_dat.to_pandas()
print(f'total datos con 1.3arcsec de separacion : {len(df)}')
dfg = df.drop_duplicates(subset='ID_SEARCH')
print(f'Total de objetos únicos en "dat": {len(dfg)}')
nan = dfg[dfg['ALPHA_J2000'].isna()]
t_nan = Table.from_pandas(nan)
#t_nan.write('nan_rows_1.3.csv', format='csv', overwrite=True)
print('Total objetos "NaN": ',+ len(nan))
dfg = dfg[~dfg['ALPHA_J2000'].isna()]
print('Total objetos sin valores "NaN": ',+ len(dfg))
print()
s_ = dfg[dfg['EXTENDED']==0]
e_ = dfg[dfg['EXTENDED']==1]
print('Cantidad de fuentes puntuales: ',+len(s_))
print('Cantidad de fuentes extendidas: ',+len(e_))
dfg['l_LONG'], dfg['b_LAT'] = galactic(dfg['RA_SEARCH'],dfg['DEC_SEARCH'])
### DEPURACION DE DATOS
df_ = df[~df['ALPHA_J2000'].isnull()]     #mantener sólo solumnas con datos en RA
df__ = df_.drop_duplicates(subset='ID_SEARCH')   #eliminar duplicados
print('-Datos con radio menor a 1.3, con datos en todas las columnas y sin duplicados:  ',+len(df__))
print('CONCLUSION: De las dos formac, devolvemos el mismo número de datos')

#------
dg = df_.loc[df_.groupby('ID').DISTANCE.idxmin()]
print(('Total datos, luego del slice: ',+len(dg)))
print('rows nan',+len(dfg[dfg['ALPHA_J2000'].isna()]))
#dg.drop_duplicates(subset='ID_SEARCH')

### RANDOM OBJECTS ###
T = Table.read('non_dr2_1.3.cat', format='ascii')
print('Cantidad de objetos sin contraparte por estudiar: ',+len(T))
#T.sort('ID_SEARCH')
#an = list(T['ID_SEARCH'])
#an[0:30]
#X1 = np.random.randint(low=0, high=len(T), size=(1000,))

### transform RA from deg to hourangle
from astropy.coordinates import Angle

ra = []
dec = []
for i in range(len(T)):
    ra_h = Angle(T['RA_SEARCH'], u.deg).hms[0][i]
    ra_m = Angle(T['RA_SEARCH'], u.deg).hms[1][i]
    ra_s = Angle(T['RA_SEARCH'], u.deg).hms[2][i]
    ra.append(str(str(int(ra_h))+' '+str(int(ra_m))+' '+str(round(ra_s,4))))
    dec_d = Angle(T['DEC_SEARCH'], u.deg).dms[0][i]
    dec_m = abs(Angle(T['DEC_SEARCH'], u.deg).dms[1][i])
    dec_s = abs(Angle(T['DEC_SEARCH'], u.deg).dms[2][i])
    dec.append(str(str(int(dec_d))+' '+str(int(dec_m))+' '+str(round(dec_s,4))))
T['ra'] = ra
T['dec']  = dec

#print('ID: ',+T['ID_SEARCH'][1])
#print('ra: ',ra)
#print('dec: ',dec)
#(T['RA_SEARCH'][3]*u.deg).to(u.hourangle)
#(T['RA_SEARCH'][3]*u.deg).to_string(unit=u.hour)

t = Table([T['ra'][:100],T['dec'][:100]])
t.write('4XMM-DR9/buscar.csv',format='csv', overwrite=True)
t1 = Table([T['ID_SEARCH'][:100],T['ra'][:100],T['dec'][:100]])
t1.write('4XMM-DR9/new_talbe.csv', format='csv', overwrite=True)

star = dfg[dfg['EXTENDED']==0]
extended = dfg[dfg['EXTENDED']==1]
datos = Table([extended['ALPHA_J2000'], extended['DELTA_J2000']], names=['RA', 'DEC'])
datos.write('extended_.cat', format='ascii', overwrite=True)

#dfg['DISTANCE'].unique()
#print(dfg.columns)
print('\n "FLUX_RADIUS" corresponde a f20,f50 y f80, respectivamente')

data = dfg[~dfg['ALPHA_J2000'].isna()]
dfg_ = Table.from_pandas(data)
#dfg_['l_LONG'], dfg_['b_LAT'] = galactic(dfg_['RA_SEARCH'],dfg_['DEC_SEARCH'])

dfg_e = dfg_[dfg_['EXTENDED']==1]
len(dfg_e)

dfg_e_coord = Table([dfg_e['ALPHA_J2000'], dfg_e['DELTA_J2000']])
dfg_e_coord.write('extended_buenas.cat', format='ascii', overwrite=True)

### PARA PODER EXTRAER EL FLUJO f50 DE LA COLUMNA ###
lst = []
for i in range(len(dfg_)):
    if len(str(dfg_['FLUX_RADIUS'][i])) > 5:
        if len(dfg_['FLUX_RADIUS'][i]) < 7 :
            kalam = dfg_['FLUX_RADIUS'][i].split(" ")[1:]
            cont = 0
            while True:
                try:
                    float(kalam[cont])
                    lst += [kalam[cont]]
                    break
                except:
                    cont += 1

        if len(dfg_['FLUX_RADIUS'][i]) >= 7 :
            kalam = dfg_['FLUX_RADIUS'][i].split(" ")[1:]
            cont = 0
            while True:
                try:
                    float(kalam[cont])
                    lst += [kalam[cont]]
                    break
                except:
                    cont += 1
    else:
        lst += [-10]
lst = np.array(lst)


print(len(dfg_) == len(lst))
print()
#for i in range(len(lst)):
    #print(lst[i])
flux = [float(x) for x in lst]
dfg_['F50'] = flux

## SCATTER PLOTS   -----
fig, (ax0, ax1, ax2) = plt.subplots(1,3, sharey=True, figsize=(20,6))
#fig, (ax0, ax1, ax2) = plt.subplots(1,3)
##fig.suptitle('2mass')
ax0.scatter(dfg['SPREAD_MODEL'], dfg['MAG_AUTO'], color='black', s=7)
ax0.scatter(gal['SPREAD_MODEL'], gal['MAG_AUTO_Ks'], color='red', s=14, marker='D')
ax1.scatter(dfg['CLASS_STAR'], dfg['MAG_AUTO'], color='black', s=7)
ax1.scatter(gal['CLASS_STAR'], gal['MAG_AUTO'], color='red', s=14, marker='D')
ax2.scatter(dfg['ELLIPTICITY'], dfg['MAG_AUTO'], color='black', s=7)
ax2.scatter(gal['ELLIPTICITY_1'], gal['MAG_AUTO'], color='red', s=14, marker='D')

ax0.set_xlim(-0.025,0.025)
#ax0.set_ylim(bottom=0)
ax0.set_ylim(9,19)
ax0.invert_yaxis()
ax0.set_ylabel(r'MAG AUTO', size=15)
ax0.set_xlabel(r'SPREAD MODEL', size=15)
ax0.minorticks_on()
ax0.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax0.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
ax0.axvline(x=0, c='gray', linestyle='--')
#ax1.set_ylim(top=30)
ax1.set_xlabel(r'CLASS STAR', size=15)
ax1.minorticks_on()
ax1.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax1.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
#ax2.set_ylim(top=30)
ax2.set_xlabel(r'ELIPTICITY', size=15)
ax2.minorticks_on()
ax2.tick_params(which='major', size=12, labelsize=15, width=1, left=True, right=True)
ax2.tick_params(which='minor', size=6, labelsize=15, width=1, left=True, right=True)
plt.tight_layout()
plt.savefig('plots.jpg')


## EXRINCTION MAP   -----
dust, interpolator = ld.load_dust(return_interpolator=True)
fig = plt.figure(figsize=(30,10))
ax = fig.add_subplot(211,aspect='equal')
ax.contour(dust['xi'], dust['yi'], dust['zi'],levels=[11,15,20,25],origin='lower',alpha=0.6,cmap='gist_yarg')
ax.scatter(disc['l_LONG'], disc['b_LAT'], s=3, c='darkred', label='Zang+21')
ax.scatter(good_gal['l_LONG'], good_gal['b_LAT'], s=9, c='green', marker='D', label=r'Px$>0.9$')
ax.scatter(dfg_e['l_LONG'], dfg_e['b_LAT'], s=50, c='deepskyblue', marker='D',alpha=0.7, label='VVV extended')
ax.scatter(all_cross['GAL_LONG'], all_cross['GAL_LAT'], s=50, c='yellow', marker='D', label='all galaxies NIRGC')
ax.scatter(gal_zang['GAL_LONG_1'], gal_zang['GAL_LAT_1'], s=50, c='blue', marker='D', label='galaxies NIRGC')
ax.set_xlabel('Galactic Longitude [deg]', fontsize=15)
ax.set_ylabel('Galactic Latitude [deg]', fontsize=15)
plt.grid(True, alpha=0.5)
ax.set_xlim(350,322)
ax.set_ylim(-2,2)
ax.minorticks_on()
ax.tick_params(which='major', size=12, labelsize=15, width=2)
ax.tick_params(which='minor', size=8, width=1)
ax.tick_params(labelsize=15, width=3)

ax.text(348, 0.05, 'Z1', fontsize=18, color='black',horizontalalignment='center',verticalalignment='center')
ax.text(343.5, 1.55, 'Z2', fontsize=18, color='black',horizontalalignment='center',verticalalignment='center')
ax.text(326.5, -1.1, 'Z3', fontsize=18, color='black',horizontalalignment='center',verticalalignment='center')
ax.add_patch(Rectangle((346.3, -1.1), 1.7, 1.25, edgecolor = 'black', facecolor = 'none', fill=False,lw=1))
ax.add_patch(Rectangle((343.2, 0.9), 0.63, 0.5, edgecolor = 'black', facecolor = 'none', fill=False,lw=1))
ax.add_patch(Rectangle((325.9, -2.0), 0.9, 0.6, edgecolor = 'black', facecolor = 'none', fill=False,lw=1))

ax = fig.add_subplot(212,aspect='equal')
ax.contour(dust['xi'], dust['yi'], dust['zi'],levels=[11,15,20,25],origin='lower',alpha=0.6,cmap='gist_yarg')
ax.scatter(disc['l_LONG'], disc['b_LAT'], s=3, c='darkred', label='Zang+21')
ax.scatter(good_gal['l_LONG'], good_gal['b_LAT'], s=9, c='green', marker='D', label=r'Px$>0.9$')
ax.scatter(dfg_e['l_LONG'], dfg_e['b_LAT'], s=50, c='deepskyblue', marker='D',alpha=0.7, label='VVV extended')
ax.scatter(all_cross['GAL_LONG'], all_cross['GAL_LAT'], s=50, c='yellow', marker='D', label='all galaxies NIRGC')
ax.scatter(gal_zang['GAL_LONG_1'], gal_zang['GAL_LAT_1'], s=50, c='blue', marker='D', label='galaxies NIRGC')
ax.set_xlabel('Galactic Longitude [deg]', fontsize=15)
ax.set_ylabel('Galactic Latitude [deg]', fontsize=15)
plt.grid(True, alpha=0.5)
ax.set_xlim(322,295)
ax.set_ylim(-2,2)
ax.minorticks_on()
ax.tick_params(which='major', size=12, labelsize=15, width=2)
ax.tick_params(which='minor', size=8, width=1)
plt.legend(loc=1, framealpha=0.5, prop={"size":16})

ax.text(320, -1, 'Z4', fontsize=18, color='black',horizontalalignment='center',verticalalignment='center')
ax.text(296, 0, 'Z5', fontsize=18, color='black',horizontalalignment='center',verticalalignment='center')
ax.add_patch(Rectangle((320.1, -1.6), 0.6, 0.7, edgecolor = 'black', facecolor = 'none', fill=False,lw=1))
ax.add_patch(Rectangle((295.8, -0.9), 0.7, 0.9, edgecolor = 'black', facecolor = 'none', fill=False,lw=1))

plt.tight_layout()
plt.savefig('map.jpg')

### ZONAS
## z1
z1 = good_gal[good_gal['l_LONG']>346.3]
z1 = z1[z1['l_LONG']<348]
z1 = z1[z1['b_LAT']>-1.1]
z1 = z1[z1['b_LAT']<0.15]
print('Z1: ',+len(z1))
## z2
z2 = good_gal[good_gal['l_LONG']>343.2]
z2 = z2[z2['l_LONG']<343.83]
z2 = z2[z2['b_LAT']>0.9]
z2 = z2[z2['b_LAT']<1.5]
print('Z2: ',+len(z2))
## z3
z3 = good_gal[good_gal['l_LONG']>325.9]
z3 = z3[z3['l_LONG']<326.8]
z3 = z3[z3['b_LAT']>-2.0]
z3 = z3[z3['b_LAT']<-1.4]
print('Z3: ',+len(z3))
## z4
z4 = good_gal[good_gal['l_LONG']>320.1]
z4 = z4[z4['l_LONG']<320.7]
z4 = z4[z4['b_LAT']>-1.6]
z4 = z4[z4['b_LAT']<-0.9]
print('Z4: ',+len(z4))
## z5
z5 = good_gal[good_gal['l_LONG']>295.8]
z5 = z5[z5['l_LONG']<296.5]
z5 = z5[z5['b_LAT']>-0.9]
z5 = z5[z5['b_LAT']<0]
print('Z5: ',+len(z5))

### DIAGRAMAS DE VIOLIN -- Zang21
d_25,d_50,d_75 = np.percentile(disc['Px'], [25,50,75])

fig = plt.subplots(1,1, figsize=(9,8))
ax = sns.violinplot(y=disc["Px"], width=0.8, linewidth=3,scale="count")

ax.tick_params(labelsize=13)
#ax.set_xlim(0.4, 1.6)
ax.set_ylabel('Px',fontsize=15)
ax.set_xlabel('Galaxies', fontsize=15)
ax.grid(axis='y', alpha=0.5)
ax.legend((d_25,d_50,d_75),title='Values for sample, for \n"First Quartile", \n "median",\n "Third Quartile"',
          fancybox=True, framealpha=1, shadow=True, borderpad=1, bbox_to_anchor=(1, 1),prop={"size":15})
ax.tick_params(size=12, labelsize=15, width=1)
plt.tight_layout()
#plt.savefig('sns_violin_plot.jpg')
