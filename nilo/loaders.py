import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
from matplotlib import tri

import warnings
warnings.filterwarnings('ignore')

def load_xmass():
    dtypes=np.dtype([
    ('name',np.str),
    ('RAJ2000',np.float32),
    ('DEJ2000',np.float32),
    ('l',np.float32),
    ('b',np.float32),
    ('vel',np.float32),
    ('velerr',np.float32),
    ('redshift',np.float32)
    ])
    columns=['name','RAJ2000','DEJ2000','l','b','vel','velerr','redshift']
    
    data = np.empty(0,dtype=dtypes)
    xmas = pd.DataFrame(data)
    xmas = xmas.append({'name':'J11494881-6400073','RAJ2000':11+49/60.0+48.81/3600.0,'DEJ2000':-(64+ 0/60.0+ 7.3/3600.0),'vel':2119,'velerr':7},ignore_index=True)
    xmas = xmas.append({'name':'J11523245-5950248','RAJ2000':11+52/60.0+32.46/3600.0,'DEJ2000':-(59+50/60.0+24.9/3600.0),'vel':5400,'velerr':15},ignore_index=True)
    xmas = xmas.append({'name':'J12092077-6229125','RAJ2000':12+ 9/60.0+20.78/3600.0,'DEJ2000':-(62+29/60.0+12.6/3600.0),'vel':1449,'velerr':14},ignore_index=True)
    xmas = xmas.append({'name':'J12280968-6054558','RAJ2000':12+28/60.0+09.68/3600.0,'DEJ2000':-(60+54/60.0+55.8/3600.0),'vel':5011,'velerr':51},ignore_index=True)
    xmas = xmas.append({'name':'J12294880-6107045','RAJ2000':12+29/60.0+48.80/3600.0,'DEJ2000':-(61+ 7/60.0+ 4.5/3600.0),'vel':5391,'velerr':319},ignore_index=True)
    xmas = xmas.append({'name':'J13464910-6024299','RAJ2000':13+46/60.0+49.11/3600.0,'DEJ2000':-(60+24/60.0+29.9/3600.0),'vel':3872,'velerr':20},ignore_index=True)
    xmas = xmas.append({'name':'J14155209-5855348','RAJ2000':14+15/60.0+52.09/3600.0,'DEJ2000':-(58+55/60.0+34.8/3600.0),'vel':5399,'velerr':4},ignore_index=True)
    c = SkyCoord(ra=xmas['RAJ2000']*15.0, dec=xmas['DEJ2000'], frame='icrs',unit='degree')
    xmas['l'] = c.galactic.l
    xmas['b'] = c.galactic.b
    xmas['redshift'] = xmas['vel']/300000.0
    
    return(xmas)

def load_cluster_d115():
    filename = "Catalogos/gxy_cluster_d015.txt"
    dtype = {'ID':str,'RA':str,'DEC':str,'l':np.float32,'b':np.float32,
         'A_Ks':np.float32,'MAG_AUTO_J':np.float32,'MAG_AUTO_H':np.float32,'MAG_AUTO_Ks':np.float32,
         'MAG_APER_J':np.float32,'MAG_APER_H':np.float32,'MAG_APER_Ks':np.float32,
         'R':np.float32,'C':np.float32,'ELLIPTICITY':np.float32,'n':np.float32,'WISE':np.int32,'z':np.float32}
    parse_dates = ['RA', 'DEC']
    cluster_d115 = pd.read_csv(filename,delim_whitespace=True,skiprows=1,names=dtype.keys(),dtype=dtype)
    return(cluster_d115)


def load_catalogue(filename):
    
    df = pd.read_csv(filename,delim_whitespace=True)
    df.replace('********','99.99999',inplace=True)
    
    a = df.dtypes == 'object'
    for column in a.index:
        if a[column]:
            df[column] = df[column].astype(np.float)
            
    c = SkyCoord(ra=df['ALPHA_J2000'], dec=df['DELTA_J2000'], frame='icrs',unit='degree')
    df['l'] = c.galactic.l
    df['b'] = c.galactic.b
    df['colorJKs'] = df['MAG_APER_J'] - df['MAG_APER_Ks']
    df['colorJH'] = df['MAG_APER_J'] - df['MAG_APER_H']
    df['colorHKs'] = df['MAG_APER_H'] - df['MAG_APER_Ks']
    
    return(df)

def load_catalogue_soto(filename):
    
    df = pd.read_csv(filename,delim_whitespace=True)
    c = SkyCoord(ra=df['ALPHA_J2000'], dec=df['DELTA_J2000'], frame='icrs',unit='degree')
    df['l'] = c.galactic.l
    df['b'] = c.galactic.b
    df['colorJKs'] = df['MAG_APER_J'] - df['MAG_APER_Ks']
    df['colorJH'] = df['MAG_APER_J'] - df['MAG_APER_H']
    df['colorHKs'] = df['MAG_APER_H'] - df['MAG_APER_Ks']
    
    return(df)

def load_voro(file_group,file_glxs):
    
    dtype = {'nmem':np.int32,'l': np.float32,'b': np.float32,'radius': np.float32,'prob': np.float32,'id':np.int32}
    voro = pd.read_csv(file_group,delim_whitespace=True,skiprows=1,names=dtype.keys(),dtype=dtype)

    dtype = {'l':np.float32,'b':np.float32,'jmag':np.float32,'hmag':np.float32,'kmag':np.float32,'id_obj':np.int32,'id_clu':np.int32}
    galxs_voro = pd.read_csv(file_glxs,delim_whitespace=True,skiprows=1,names=dtype.keys(),dtype=dtype)

    galxs_voro['colorJKs'] = galxs_voro['jmag'] - galxs_voro['kmag']

    voro['B1'] = pd.Series(index=voro.index)
    voro['B3'] = pd.Series(index=voro.index)
    voro['BF'] = pd.Series(index=voro.index)

    for iclu, k in enumerate(voro['id'].values):
        
        if k == -1:
            continue
            
        class_member_mask = (galxs_voro['id_clu'] == k)
        nmem = class_member_mask.sum()
        assert(nmem == voro['nmem'][iclu])
        
        if nmem >= 2:
            
            xy = galxs_voro[class_member_mask]
            
            Kmag = xy['kmag']
            Kmag = Kmag.sort_values()      
            if nmem < 3:
                B1 = Kmag.iloc[0]
                B3 = 0.0
                BF = Kmag.iloc[-1]      
            else:
                B1 = Kmag.iloc[0]
                B3 = Kmag.iloc[2]
                BF = Kmag.iloc[-1]
            
            lmean = xy['l'].mean() 
            dl = np.ptp(xy['l'])*60
            bmean = xy['b'].mean() 
            db = np.ptp(xy['b'])*60
        
            voro['B1'][iclu] = B1
            voro['B3'][iclu] = B3
            voro['BF'][iclu] = BF
            
    return(voro,galxs_voro)


def load_mst(file_group,file_glxs):
    
    dtype = {'nmem':np.int32,'l': np.float32,'b': np.float32,'radius': np.float32,'prob': np.float32,'id':np.int32}
    voro = pd.read_csv(file_group,delim_whitespace=True,skiprows=1,names=dtype.keys(),dtype=dtype)

    dtype = {'l':np.float32,'b':np.float32,'jmag':np.float32,'hmag':np.float32,'kmag':np.float32,'id_obj':np.int32,'id_clu':np.int32}
    galxs_voro = pd.read_csv(file_glxs,delim_whitespace=True,skiprows=1,names=dtype.keys(),dtype=dtype)

    galxs_voro['colorJKs'] = galxs_voro['jmag'] - galxs_voro['kmag']

    voro['B1'] = pd.Series(index=voro.index)
    voro['B3'] = pd.Series(index=voro.index)
    voro['BF'] = pd.Series(index=voro.index)

    for iclu, k in enumerate(voro['id'].values):
        
        if k == -1:
            continue
            
        class_member_mask = (galxs_voro['id_clu'] == k)
        nmem = class_member_mask.sum()
        assert(nmem == voro['nmem'][iclu])
        
        if nmem >= 2:
            
            xy = galxs_voro[class_member_mask]
            
            Kmag = xy['kmag']
            Kmag = Kmag.sort_values()      
            if nmem < 3:
                B1 = Kmag.iloc[0]
                B3 = 0.0
                BF = Kmag.iloc[-1]      
            else:
                B1 = Kmag.iloc[0]
                B3 = Kmag.iloc[2]
                BF = Kmag.iloc[-1]
            
            lmean = xy['l'].mean() 
            dl = np.ptp(xy['l'])*60
            bmean = xy['b'].mean() 
            db = np.ptp(xy['b'])*60
        
            voro['B1'][iclu] = B1
            voro['B3'][iclu] = B3
            voro['BF'][iclu] = BF
            
    return(voro,galxs_voro)

def load_optics(file_groups,file_galxs):
    fields ={'GroupID':np.int32,'Nmem':np.int32,'l':np.float32,'b':np.float32,'delta_l':np.float32,
             'delta_b':np.float32,'Ks1':np.float32,'Ks3':np.float32,'KsF':np.float32,'radius':np.float32,
             'lB1':np.float32,'bB1':np.float32,'radiusB1':np.float32}
    grupos = pd.read_csv(file_groups,delim_whitespace=True,skiprows=1,names=fields.keys(),dtype=fields)
    fields = {'ID':np.int32,'gID':np.int32,'l':np.float32,'b':np.float32}
    glxs = pd.read_csv(file_galxs,delim_whitespace=True,skiprows=1,names=fields.keys(),dtype=fields)
    return(grupos,glxs)


def load_compact_groups():
    #grupos compactos
    file_CGs = 'CGs/new_centre_CG_VVV.dat'
    CGs = pd.read_csv(file_CGs,sep=' ')
    CGs_coords = SkyCoord(ra=CGs['ra'] , dec=CGs['dec'], frame='icrs',unit='degree')
    CGs['l'] = CGs_coords.galactic.l
    CGs['b'] = CGs_coords.galactic.b

    nCGs = len(CGs)

    Galxs_CGs = pd.read_csv("CGs/new_gal_in_CGVVV.dat",sep=' ')

    Galxs_CGs_coords = SkyCoord(ra=Galxs_CGs['ra'] , dec=Galxs_CGs['dec'], frame='icrs',unit='degree')
    Galxs_CGs['l'] = Galxs_CGs_coords.galactic.l
    Galxs_CGs['b'] = Galxs_CGs_coords.galactic.b
    return(CGs,Galxs_CGs)

def load_dust(return_interpolator=False):
    dd = pd.read_csv("extinction12.tbl",names=['ra','dec','Av_SandF','AV_SFD'],header=None,skiprows=16,delim_whitespace=True,skipinitialspace=True,usecols=[0,1,8,14])
    c = SkyCoord(ra=dd['ra'], dec=dd['dec'], frame='icrs',unit='degree')
    dd['l'] = c.galactic.l
    dd['b'] = c.galactic.b
    xi = np.linspace(dd['l'].min(), dd['l'].max(), 1000)
    yi = np.linspace(dd['b'].min(), dd['b'].max(), 100)

    triang = tri.Triangulation(dd['l'], dd['b'])
    interpolator = tri.LinearTriInterpolator(triang, dd['Av_SandF'])
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = interpolator(Xi, Yi)
    
    dust = {'xi': Xi, 'yi': Yi, 'zi': Zi}
    if return_interpolator:
        return(dust,interpolator)
    else:
        return(dust)