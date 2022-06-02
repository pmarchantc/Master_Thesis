#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  util.py
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
# Module with simple tools for specific VVV and VVVx data handling routines.

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def diff(list1, list2):
    c = set(list1).union(set(list2))  # or c = set(list1) | set(list2)
    d = set(list1).intersection(set(list2))  # or d = set(list1) & set(list2)
    return list(c - d)

def df_3bands(self, Nx=44, Ny=44, Nb=3):
    """This method transforms a DataFrame conformed for several images of 
    (44, 44, 3) in each row and columns whit all pixel values in each filter in a 
    DataFrame conform by unique image, which has three columns, each 
    column is a band and each row is a px. 
    input: DataFrame with rows equal to each images and columns equal to
    pixel values in a image for each filter more five columns with values
    of Target, ra, dec, TILE and ID. 
    output: DataFrame with three columns, each column represent a band and each row
    represents a pixel.
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

    df_bandas = pd.DataFrame(oo_0, columns=['banda_0'])
    df_bandas['banda_1'] = oo_1
    df_bandas['banda_2'] = oo_2

    return(df_bandas)


def color(self, Nx=44, Ny=44, Nb=3):   
    """Create a new DataFrame with the difference of the filters for each image: H-K and J-K, 
    the size of the images must be (44, 44, 3).
    input: DataFrame with rows = images and columns = value in each px more five 
    columns which specify Target, ra, dec. TILE and ID.
    output: Same DataFrame of input more columns with the difference between filters
    """
    
    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))

    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))

    colour = np.zeros((self.shape[0], Nx*Ny*2))

    for i in range(self.shape[0]):
        banda_j = np.float32(M_4[i,:,:,0])
        banda_h = np.float32(M_4[i,:,:,1])
        banda_k = np.float32(M_4[i,:,:,2])
        h_k = np.reshape(banda_h-banda_k, Nx*Ny)
        j_k = np.reshape(banda_j-banda_k, Nx*Ny)
        c = np.concatenate((h_k, j_k))
        colour[i,:] = c


    df_colour = pd.DataFrame(colour)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_colour[difer] = self[difer]

    return(df_colour)

def colorHK_JK(self, Nx=44, Ny=44, Nb=3):   
    """Create a new DataFrame with the difference of the filters for each image: H-K and J-K, 
    the size of the images must be (44, 44, 3).
    input: DataFrame with rows = images and columns = value in each px more five 
    columns which specify Target, ra, dec. TILE and ID.
    output: Same DataFrame of input more columns with the difference between filters
    """
    
    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))

    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))

    colour = np.zeros((self.shape[0], Nx*Ny*2))

    for i in range(self.shape[0]):
        banda_j = np.float32(M_4[i,:,:,0])
        banda_h = np.float32(M_4[i,:,:,1])
        banda_k = np.float32(M_4[i,:,:,2])
        HK = banda_h-banda_k
        JK = banda_j-banda_k
        h_k = np.reshape(HK, Nx*Ny)
        j_k = np.reshape(JK, Nx*Ny)
        c = np.concatenate((h_k, j_k), axis=None)
        colour[i,:] = c
        
    df_colour = pd.DataFrame(colour)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_colour[difer] = self[difer]

    return(df_colour)

def colorJK_J(self, Nx=44, Ny=44, Nb=3):   
    """Create a new DataFrame with the difference of the filters for each image: J-K y J, 
    the size of the images must be (44, 44, 3).
    input: DataFrame with rows = images and columns = value in each px more five 
    columns which specify Target, ra, dec. TILE and ID.
    output: Same DataFrame of input more columns with the difference between filters
    """
    
    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))

    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))

    colour = np.zeros((self.shape[0], Nx*Ny*2))

    for i in range(self.shape[0]):
        banda_j = np.float32(M_4[i,:,:,0])
        banda_h = np.float32(M_4[i,:,:,1])
        banda_k = np.float32(M_4[i,:,:,2])
        #h_k = np.reshape(banda_h-banda_k, Nx*Ny)
        j_k = np.reshape(banda_j-banda_k, Nx*Ny)
        j = np.reshape(banda_j, Nx*Ny)
        
        c = np.concatenate((j_k, j))
        colour[i,:] = c


    df_colour = pd.DataFrame(colour)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_colour[difer] = self[difer]

    return(df_colour)

def colorJK_K(self, Nx=44, Ny=44, Nb=3):   
    """Create a new DataFrame with the difference of the filters for each image: J-K y K, 
    the size of the images must be (44, 44, 3).
    input: DataFrame with rows = images and columns = value in each px more five 
    columns which specify Target, ra, dec. TILE and ID.
    output: Same DataFrame of input more columns with the difference between filters
    """
    
    dim_img = Nx*Ny*Nb
    px = list(map(lambda x: str(x),  list(np.arange(dim_img))))

    M = self[px].to_numpy()
    M_4 = np.reshape(M, (self.shape[0], Nx, Ny, Nb))

    colour = np.zeros((self.shape[0], Nx*Ny*2))

    for i in range(self.shape[0]):
        banda_j = np.float32(M_4[i,:,:,0])
        banda_h = np.float32(M_4[i,:,:,1])
        banda_k = np.float32(M_4[i,:,:,2])
        #h_k = np.reshape(banda_h-banda_k, Nx*Ny)
        j_k = np.reshape(banda_j-banda_k, Nx*Ny)
        k = np.reshape(banda_k, Nx*Ny)
        
        c = np.concatenate((j_k, k))
        colour[i,:] = c


    df_colour = pd.DataFrame(colour)
    columns = list(self.columns)
    difer = diff(columns, px)
    df_colour[difer] = self[difer]

    return(df_colour)





def pca(self):
    if self.shape[0] > 10:
        n = self.shape[0] - 10
        pca = PCA(n_components=450) #450 es un buen número
        pca_m = pca.fit(self)
        X_train_z_all = pca_m.transform(self)

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
        ax1.plot(pca_m.explained_variance_ratio_, 'o-',color='lightgrey')
        ax1.set_xlim(-0.5, 10)
        ax1.set_title('variance ratio')
        ax1.set_xlabel('PCA')
        ax1.set_ylabel('Autovalere')

        ax2.plot(pca_m.explained_variance_, 'o-',color='purple')
        ax2.set_title('variance ratio')
        ax2.set_xlabel('PCA')
    else:
        print('The componets number is <= 10')
