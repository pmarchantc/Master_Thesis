#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  loaders.py
#
#  -----------------------------------------------------------------------------------------------------------
#  Copyright 2021 M.Sgro, I.Daza, M.Lares, ... <vanessa.daza@unc.edu.ar>, RESCRIBIR UNA VEZ TERMINADO
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
# -----------------------------------------------------------------------------------------------------------
#
# Module with simple tools for load specific VVV and VVVx data.

import os
import pandas as pd
import glob
from astropy.io import fits
import numpy as np
from astropy.table import Table
import timeit
import time
    
def __init__(self):

    self.Nx = Nx
    self.Ny = Ny
    self.Nb = Nb
    self.root = root
    self.tile = tile
    self.tiles = os.listdir(self.root)

    if self.Nb != 3:
        print('The bands number must be 3')
        return

    if len(self.tiles) == 0:
        print('%s is empty' % self.tile)
        return

    if not os.path.isdir(self.root):
        print('%s doesnt exists' % self.root)
        return

    self.nImages = len(glob.glob(self.root+'*.fits'))
    if self.nImages == 0:
        print('%s is empty' % self.root)
    else:
        print('Number of files in %s: %d' % (self.root,self.nImages))

#@timeit 
def load_images(self, tile=None, Nx=44, Ny=44, Nb=3):
    """This method generates a list of dataframes for each tile 
    from the information contained in the images of galaxies and
    other objects.
    input: path to the folder containing the image folders for each tile.
    output: List of dataframes, each element corresponds to a tile, which 
    contains information about the images, label (galaxy = 1, other object = 0),
    ra, dec and ID. 
    """
    
    if not os.path.isdir(self):
        print('%s doesnt exists' % self)
        return
    
    if Nb != 3:
        print('The bands number must be 3')
        return
     
    tiles = os.listdir(f'{self}Images/')
    if len(tiles) == 0:
        print('%s is empty' % tiles)
        return
    
    if tile is not None and type(tile)==list and type(tile[0]):
        lista_tiles = tile
        print('generating dataframes of the tile(s):', tile)
    else: 
        #print('tile is not a list or the list element is not str')
        print('%s - generating dataframes'% tiles)
        lista_tiles = tiles
        
    df_list = []
    for til in lista_tiles:
        print(til)
        path_tile_gx = f'{self}Images/{til}/individual_15arcsec/'
        path_labels = f'{self}Galaxies/{til}_galaxies.csv'
        
        
        img = os.listdir(path_tile_gx)
        #if len(img) > 500:
        #    continue
        df_labels = pd.read_csv(path_labels)
        df_labels.set_index('ID',inplace=True)

        Npx = Nx*Ny*Nb
        all_gx = np.zeros((len(img),Npx))

        ra = np.zeros(len(img),dtype=float)
        dec = np.zeros_like(ra)
        ID = np.zeros(len(img),dtype='<U24')
        TILE = np.zeros(len(img),dtype='<U5')
        AIRM = np.zeros(len(img),dtype='<f4')
        FWHM = np.zeros(len(img),dtype='<f4')
        CRPIX1 = np.zeros(len(img),dtype='<f8')
        CRPIX2 = np.zeros(len(img),dtype='<f8')
        
        label = np.zeros(len(img),dtype=int) 

        t = 0
        for i,archive in enumerate(img):
            hdul =  fits.open(path_tile_gx+str(archive))
            
            ra[i] = hdul[0].header['RA']
            dec[i] = hdul[0].header['DEC']
            TILE[i] = hdul[0].header['TILE']
            ID[i] = hdul[0].header['ID']
            AIRM[i] = hdul[0].header['HIERARCH ESO TEL AIRM END']
            FWHM[i] = hdul[0].header['HIERARCH ESO TEL AMBI FWHM END']
            CRPIX1[i] = hdul[0].header['CRPIX1']
            CRPIX2[i] = hdul[0].header['CRPIX2']
            
            all_gx[i,:] = np.concatenate([np.array(hdul[0].data).flatten(order='F'),np.array(hdul[1].data).flatten(order='F'),np.array(hdul[2].data).flatten(order='F')])
            # Encuentro la etiqueta
            label[i] = df_labels.loc[ID[i]]

        # Dataset de gx como DataFrame  
        df_gx = pd.DataFrame(all_gx)

        df_gx['Target'] = label
        df_gx['ra'] = ra
        df_gx['dec'] = dec
        df_gx['TILE'] = TILE
        df_gx['ID'] = ID
        df_gx['AIRM'] = AIRM
        df_gx['FWHM'] = FWHM
        df_gx['CRPIX1'] = CRPIX1
        df_gx['CRPIX2'] = CRPIX2
            
        # Lista de DataFrames por cada tile
        df_list.append(df_gx)

    if len(df_list) == 1: 
        df_output = df_gx
        print(len(df_list))
    else:
        print('concatenando', len(df_list), 'tiles')
        df_output = pd.concat(df_list)
    
    df_output.set_index('ID', inplace=True)
    print('number of duplicate values', df_output[df_output.astype(str).duplicated()].shape)
    df_output = df_output.loc[~(df_output.index.astype(str).duplicated(keep="first"))]
    df_output = df_output.reset_index()

    return(df_output)

#@timeit 
def load_extragalactic_table_VVVx(filename):
    """This method generates a list of dataframes for each tile 
    from the information contained in the tables of galaxies and
    other objects.
    input: path to the folder containing the .fits folders for each tile.
    output: List of dataframes, each element corresponds to a tile, which 
    contains information about output of SExtractor.
    """
    data = Table.read(filename, format='fits')

    one_value_columns = [name for name in data.colnames if len(data[name].shape) <= 1]
    multi_values_coulumns = [name for name in data.colnames if len(data[name].shape) > 1]

    _data = data[one_value_columns].to_pandas()

    _data[['MAG_APER_KS_0','MAG_APER_KS_1','MAG_APER_KS_2','MAG_APER_KS_3']] = pd.DataFrame(data['MAG_APER'].tolist())
    _data[['MAGERR_APER_KS_0','MAGERR_APER_KS_1','MAGERR_APER_KS_2','MAGERR_APER_KS_3']] = pd.DataFrame(data['MAGERR_APER'].tolist())
    _data[['FLUX_RADIUS_0','FLUX_RADIUS_1','FLUX_RADIUS_2']] = pd.DataFrame(data['FLUX_RADIUS'].tolist())

    _data.rename(columns = {'MAG_APER_KS_0': 'MAG_APER_KS','MAG_APER_KS_0_CORR': 'MAG_APER_KS_CORR', 'MAGERR_APER_KS_0': 'MAGERR_APER_KS'}, inplace = True)
    _data.rename(columns = {'MAG_APER_KS': 'MAG_APER'}, inplace = True)
    _data.rename(columns = {'MAG_APER_J_0': 'MAG_APER_J', 'MAG_APER_J_0_CORR': 'MAG_APER_J_CORR', 'MAGERR_APER_J_0': 'MAGERR_APER_J'}, inplace = True)
    _data.rename(columns = {'MAG_APER_H_0': 'MAG_APER_H', 'MAG_APER_H_0_CORR': 'MAG_APER_H_CORR', 'MAGERR_APER_H_0': 'MAGERR_APER_H'}, inplace = True)

    _data['ID'] = _data['ID'].str.decode(encoding = 'ASCII')
    _data.set_index('ID',inplace=True)

    #create_ID(_data)
    return(_data)

def load_tables(self, tile=None):
    """This method generates a list of dataframes for each tile 
    from the information contained in the tables of galaxies and
    other objects.
    input: path to the folder containing the .fits folders for each tile.
    output: List of dataframes, each element corresponds to a tile, which 
    contains information about output of SExtractor.
    """
    
    if not os.path.isdir(f'{self}Extragalactic/'):
        print('%s doesnt exists' % self)
        return
     
    tiles = os.listdir(f'{self}Extragalactic/')
    if len(tiles) == 0:
        print('%s is empty' % tiles)
        return
    
    if tile is not None and type(tile)==list and type(tile[0]):
        lista_tiles = tile
        print('generating dataframes of the tile:', tile)
    else: 
        #print('tile is not a list or the list element is not str')
        print('%s - generating dataframes')
        lista_tiles = tiles
        
    df_list = []
    for til in lista_tiles:
        path_tables = f'{self}Extragalactic/{til}'
        data = Table.read(path_tables, format='fits')
        names = [name for name in data.colnames if len(data[name].shape) <= 1]
        data_d = data[names].to_pandas()

        data_d[['FLUX_RADIUS_0','FLUX_RADIUS_1','FLUX_RADIUS_2']] = pd.DataFrame(data['FLUX_RADIUS'].tolist())
        data_d['colorJH']  = data_d['MAG_APER_J_CORR'] - data_d['MAG_APER_H_CORR']
        data_d['colorHKs'] = data_d['MAG_APER_H_CORR'] - data_d['MAG_APER_KS_CORR']
        data_d['colorJKs'] = data_d['MAG_APER_J_CORR'] - data_d['MAG_APER_KS_CORR']

        columns = ['TILE', 'ID', 'ALPHA_J2000', 'DELTA_J2000',
                   'MAG_PSF','MAG_AUTO','MAG_APER','MAG_MODEL','SPREAD_MODEL',
                   'AMODEL_IMAGE','BMODEL_IMAGE','ELONGATION',
                   'ELLIPTICITY','A_IMAGE','B_IMAGE','SPHEROID_SERSICN','MU_MAX',
                   'CLASS_STAR','MAG_APER_H','MAG_PSF_H',
                   'MAG_APER_J', 'MAG_PSF_J', 'A_v', 'A_Ks', 'A_H', 'A_J',
                   'MAG_APER_KS_CORR','MAG_PSF_KS_CORR','MAG_APER_H_CORR',
                   'MAG_PSF_H_CORR','MAG_APER_J_CORR', 'MAG_PSF_J_CORR', 
                   'FLUX_RADIUS_0', 'FLUX_RADIUS_1', 'FLUX_RADIUS_2',
                   'colorJH','colorHKs','colorJKs']

        df_text = data_d[columns]
        df_text.TILE = df_text.TILE.str.decode(encoding = 'UTF-8')
        df_text.ID = df_text.ID.str.decode(encoding = 'UTF-8')
        # Lista de DataFrames por cada tile
        df_list.append(df_text)

    if len(df_list) == 1: 
        df_output = df_text
        print(len(df_list))
    else:
        print('concatenando', len(df_list), 'tiles')
        df_output = pd.concat(df_list)
    
    df_output.set_index('ID', inplace=True)
    print('number of duplicate values', df_output[df_output.astype(str).duplicated()].shape)
    df_output = df_output.loc[~(df_output.index.astype(str).duplicated(keep="first"))]
    df_output = df_output.reset_index()

    return(df_output)


def match(df_images, df_tables):
    """This function performs the union between the DataFrames
    with information of the images and SExtractor.
    
    input: Two Dataframes with information of images and SExtractor
    
    output: A DataFrame.
    """
    
    if df_images.shape[0] != df_tables.shape[0]:
        print('Dataframes have different lengths.')
    
    if set(df_images.TILE) != set(df_tables.TILE):
        print('The list of tiles is different')
        print(list(set(df_images.TILE)^set(df_tables.TILE)))
    else:
        df_images = df_images.drop('TILE', axis=1)
    
    df_img_tab = df_images.set_index(['ID']).join(df_tables.set_index(['ID']))
    df_img_tab = df_img_tab.reset_index()
    return(df_img_tab)
    
