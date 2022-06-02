#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  normalization.py
#
#  -----------------------------------------------------------------------------------------------------------
#  Copyright 2021 I.Daza, M.Sgro, M.Lares, ... <vanessa.daza@unc.edu.ar>, RESCRIBIR UNA VEZ TERMINADO
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
# This module contains the normalization functions of the images.

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler


def diff(list1, list2):
    c = set(list1).union(set(list2))  # or c = set(list1) | set(list2)
    d = set(list1).intersection(set(list2))  # or d = set(list1) & set(list2)
    return list(c - d)

#def min_max_one_img(self, Nx=44, Ny=44, Nb=3):
def min_max_one_img(self, Nx=11, Ny=11, Nb=3):
    """This function applies the min max method to each row,
    i.e. to each image of (44, 44, 3).
    input: DataFrame, each row is an image, the first columns are the 
    pixel values and the last 4 columns son: Target, ra, dec and ID.
    output: DataFrame with the same structure that dataframe of input. 
    """

    dim_img = Nx*Ny*Nb
    M_mm = np.zeros((self.shape[0], dim_img))
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))
    M_row = self[px].to_numpy()
    scaler_m = MinMaxScaler()

    for i in range(self.shape[0]):
        scaler_m.fit(M_row[i].reshape(-1,1))
        mm = scaler_m.transform(M_row[i].reshape(-1,1))
        M_mm[i,:] = np.reshape(mm, dim_img)

    df_mm = pd.DataFrame(M_mm)
    
    columns = list(self.columns)
    difer = diff(columns, px)
    df_mm[difer] = self[difer]

    return(df_mm)

def min_max_all_img(self, Nx=44, Ny=44, Nb=3):
    """This method normalizations all the pixel values for a given filter 
    over the full set of images, the method uses min max method, getting
    the min. max and sigma from the pixel values distributions of all images
    in each band.
    input: DataFrame with rows equal to each image and columns equal the 
    pixels more four columns, which specify Target, ra, dec and ID.
    output: DataFrame with the same structure of the DataFrame of input.
    """

    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))
    #px = list(np.arange(dim_img))
    
    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))
    banda_0 = np.float32(M_4[:,:,:,0])
    banda_1 = np.float32(M_4[:,:,:,1])
    banda_2 = np.float32(M_4[:,:,:,2])

    dim = self.shape[0]*Nx*Ny

    oo_0 = np.reshape(banda_0,[dim])
    oo_1 = np.reshape(banda_1,[dim])
    oo_2 = np.reshape(banda_2,[dim])

    scaler_m = MinMaxScaler()

    scaler_m.fit(oo_0.reshape(-1,1))
    M_mm_0 = scaler_m.transform(oo_0.reshape(-1,1))

    scaler_m.fit(oo_1.reshape(-1,1))
    M_mm_1 = scaler_m.transform(oo_1.reshape(-1,1))

    scaler_m.fit(oo_2.reshape(-1,1))
    M_mm_2 = scaler_m.transform(oo_2.reshape(-1,1))

    df_min_max = pd.DataFrame(M_mm_0, columns=['banda_0'])
    df_min_max['banda_1'] = M_mm_1
    df_min_max['banda_2'] = M_mm_2

    df_prueb = np.reshape(df_min_max.to_numpy(), (self.shape[0], dim_img))
    df_pp = pd.DataFrame(df_prueb)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_pp[difer] = self[difer]

    return(df_pp)

def z_score_all(self, Nx=44, Ny=44, Nb=3):
    """This method normalizations all the pixel values for a given filter 
    over the full set of  images,  the method uses z-score method, getting
    the mean and sigma from the pixel values distributions of all images
    in each band.
    input: DataFrame with rows equal to each image and columns equal the 
    pixels more four columns, which specify Target, ra, dec and ID.
    output: DataFrame with the same structure of the DataFrame of input.
    """
    
    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))

    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))
    banda_0 = np.float32(M_4[:,:,:,0])
    banda_1 = np.float32(M_4[:,:,:,1])
    banda_2 = np.float32(M_4[:,:,:,2])

    dim = self.shape[0]*Nx*Ny

    oo_0 = np.reshape(banda_0, [dim])
    oo_1 = np.reshape(banda_1, [dim])
    oo_2 = np.reshape(banda_2, [dim])

    scaler = StandardScaler()

    scaler.fit(np.array(oo_0).reshape(-1,1))
    M_z_0 = scaler.transform(oo_0.reshape(-1,1))

    scaler.fit(oo_1.reshape(-1,1))
    M_z_1 = scaler.transform(oo_1.reshape(-1,1))

    scaler.fit(oo_2.reshape(-1,1))
    M_z_2 = scaler.transform(oo_2.reshape(-1,1))

    df_z = pd.DataFrame(M_z_0, columns=['banda_0'])
    df_z['banda_1'] = M_z_1
    df_z['banda_2'] = M_z_2


    df_prue = np.reshape(df_z.to_numpy(), (self.shape[0], dim_img))
    df_p = pd.DataFrame(df_prue)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_p[difer] = self[difer]

    return(df_p)
