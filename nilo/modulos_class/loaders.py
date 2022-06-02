!/usr/bin/env python
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

class Load:
    
    def __init__(self, tile=['d115'], root='../../tile/', Nx=44, Ny = 44, Nb=3):
        
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
            
    def load_images(self):
        """This method generates a list of dataframes for each tile 
        from the information contained in the images of galaxies and
        other objects.
        input: path to the folder containing the image folders for each tile.
        output: List of dataframes, each element corresponds to a tile, which 
        contains information about the images, label (galaxy = 1, other object = 0),
        ra, dec and ID. 
        """
        if self.tile is None:
            print('generating list of dataframes:', self.tiles)
            lista_tiles = self.tiles
        else:
            print('generating dataframes of the tile:', self.tile)
            lista_tiles = self.tile
            
        
        for til in lista_tiles:
            ## DataFrame galaxies 
            path_tile_gx = f'../../tile/{til}/individual_15arcsec/'
            path_tablita_etiquetas = f'../../tile/{til}/{tile}_galaxies.csv/'
            GX = os.listdir(path_tile_gx)

            df_list = []

            all_gx = []

            ra = []
            dec = []
            ID = []

            for archive in GX:

                data_H = fits.open(path_tile_gx+str(archive))[0].data
                data_J = fits.open(path_tile_gx+str(archive))[1].data
                data_Ks = fits.open(path_tile_gx+str(archive))[2].data

                coord_ra = fits.open(path_tile_gx+str(archive))[0].header['RA']
                coord_dec = fits.open(path_tile_gx+str(archive))[0].header['DEC']
                id_gx = fits.open(path_tile_gx+str(archive))[0].header['ID']
                tile = fits.open(path_tile_gx+str(archive))[0].header['TILE']

                # Cada Ima como matrix por si hace falta
                M_gx =  np.zeros((self.Nx, self.Ny, self.Nb))
                M_gx[:,:,0] =  data_H
                M_gx[:,:,1] =  data_J
                M_gx[:,:,2] =  data_Ks

                # Cada Ima como un vector de dim 5808
                Npx = self.Nx*self.Ny*self.Nb
                gx = np.reshape(M_gx,[Npx])
                all_gx.append(gx)

                ra.append(coord_ra)
                dec.append(coord_dec)
                ID.append(id_gx)


                # Dataset de gx como DataFrame  
                df_gx = pd.DataFrame(all_gx)

            target = np.ones(df_gx.shape[0])
            df_gx['Target'] = target  
            df_gx['ra'] = ra
            df_gx['dec'] = dec
            df_gx['ID'] = ID

            ## DataFrame Other objects
            path_tile_other = f'../../tile/{til}/other_objects_{til}/'
            other_obj = os.listdir(path_tile_other)

            all_oo = []

            ra = []
            dec = []
            ID = []

            for archive in other_obj: 

                data_H = fits.open(path_tile_other+str(archive))[0].data
                data_J = fits.open(path_tile_other+str(archive))[1].data
                data_Ks = fits.open(path_tile_other+str(archive))[2].data

                coord_ra = fits.open(path_tile_other+str(archive))[0].header['RA']
                coord_dec = fits.open(path_tile_other+str(archive))[0].header['DEC']
                id_oo = fits.open(path_tile_other+str(archive))[0].header['ID']

                # Cada Ima como matrix por si hace falta
                M_oo =  np.zeros((self.Nx, self.Ny, self.Nb))
                M_oo[:,:,0] =  data_H
                M_oo[:,:,1] =  data_J
                M_oo[:,:,2] =  data_Ks

                # Cada Ima como un vector de dim Nx*Ny*Nb
                Npx = self.Nx*self.Ny*self.Nb
                oo = np.reshape(M_oo,[Npx])
                all_oo.append(oo)

                ra.append(coord_ra)
                dec.append(coord_dec)
                ID.append(id_oo)

                # Dataset de gx como DataFrame  
                df_oo = pd.DataFrame(all_oo)

            target = np.zeros(df_oo.shape[0])
            df_oo['Target'] = target    
            df_oo['ra'] = ra
            df_oo['dec'] = dec
            df_oo['ID'] = ID

            # Dataframe de galaxias y otros objetos
            df_row = pd.concat([df_gx, df_oo])

            # Lista de DataFrames por cada tile
            df_list.append(df_row)
        
        if len(df_list) == 1: 
            df_output = df_row
        else:
            df_output = df_list
            
        return(df_output)
