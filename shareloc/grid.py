#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of Shareloc
# (see https://github.com/CNES/shareloc).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np
from shareloc.readwrite import  read_bsq_hd
from shareloc.math_utils import interpol_bilin, interpol_bilin_vectorized

#-------------------------------------------------------------------------------
class grid:
    """ multi H direct localization grid handling class 
    """
    def __init__(self, grid_filename, grid_format ='bsq'):
        """
        Constructor
        :param grid_filename: grid filename
        :type grid_filename: string
        :param grid_format: grid format (by default bsq)
        :type grid_format: string
        """
        self.filename = grid_filename
        self.format    = grid_format
        self.row0      = None
        self.col0      = None
        self.nbrow     = None
        self.nbcol     = None
        self.steprow    = None
        self.stepcol    = None
        self.repter    = None
        self.nbalt     = None
        self.index_alt = {}
        self.gld_lon   = None
        self.gld_lat   = None
        self.alts_down = []
        self.rowmax    = None
        self.colmax    = None
        self.load()
        self.type = 'multi H grid'

    def load(self):
        """
        header and grid loading function
        2 data cubes are defined :
        - gld_lon : [alt,row,col]
        - gld_lat : [alt,row,col]
        bsq grids are stored by increasing altitude H0 ... Hx
        internal structure is decreasing one
        """
        if self.format=='bsq':

            dico_a_lire = {'nbrow':('LINES',int),'nbcol':('COLUMNS',int),'bpp':('BITS PER PIXEL',int),\
            'nbalt':('NB ALT',int),'stepcol':('PAS COL',float),'steprow':('PAS LIG',float),\
            'col0':('COL0',float),'row0':('LIG0',float),\
            'repter':('REFERENTIEL TERRESTRE',str)}

            nom_hd      = self.filename[:-4] + '1.hd'
            dico_hd = read_bsq_hd(nom_hd, dico_a_lire)

            for var in dico_hd:
                setattr(self,var,dico_hd[var])

            """renvoie une structure 3D [i_alt][l,c]"""
            gld_lon = np.zeros((self.nbalt,self.nbrow,self.nbcol))
            gld_lat = np.zeros((self.nbalt,self.nbrow,self.nbcol))

            codage = float

            for i in range(self.nbalt):
                k = self.nbalt-i
                nom_gri_lon = self.filename[:-4] + str(k) + '.c1'
                nom_gri_lat = self.filename[:-4] + str(k) + '.c2'
                nom_hd      = self.filename[:-4] + str(k) + '.hd'

                gld_lon[i,:,:] = np.fromfile(nom_gri_lon,dtype=codage).reshape((self.nbrow,self.nbcol))
                gld_lat[i,:,:] = np.fromfile(nom_gri_lat,dtype=codage).reshape((self.nbrow,self.nbcol))

                dico_hd = read_bsq_hd(nom_hd, {'index':('ALT INDEX', int), 'alt':('ALTITUDE', float)})
                self.index_alt[dico_hd['index']] = dico_hd['alt']

            self.gld_lon = gld_lon
            self.gld_lat = gld_lat
            self.alts_down = [self.index_alt[_] for _ in range(int(self.nbalt-1),-1,-1)]
            self.rowmax = self.row0 + self.steprow*(self.nbrow-1)
            self.colmax = self.col0 + self.stepcol*(self.nbcol-1)
        else:

            print("dtm format is not handled")

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers
        :return alt_min,lat_max
        :rtype list
        """
        return [self.alts_down[-1], self.alts_down[0]]

    def direct_loc_h(self, row, col, alt):
        """
        direct localization at constant altitude
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param alt :  altitude
        :type alt : float
        :return ground position (lon,lat,h)
        :rtype numpy.array
        """
        #faire une controle sur row / col !!!!
        # 0.5 < row < rowmax
        (kh,kb) = self.return_grid_index(alt)
        altbas  = self.alts_down[kb]
        althaut = self.alts_down[kh]
        dh = (alt - altbas)/(althaut - altbas)
        mats =  [self.gld_lon[kh:kb+1,:,:],self.gld_lat[kh:kb+1,:,:]]
        P     = np.zeros(3)
        P[2] = alt
        dl = (row - self.row0)/self.steprow
        dc = (col - self.col0)/self.stepcol
        [vlon,vlat] = interpol_bilin(mats,self.nbrow,self.nbcol,dl,dc)
        P[0] = (dh*vlon[0] +(1-dh)*vlon[1])
        P[1] = (dh*vlat[0] +(1-dh)*vlat[1])
        return P

    def direct_loc_h_vectorized(self, row, col, alt):
        """
        direct localization at constant altitude on multiple points
        :param row :  line sensor position
        :type row : 1D numpy array, dtype=float64
        :param col :  column sensor position
        :type col : 1D numpy array, dtype=float64
        :param alt :  altitude
        :type alt : float
        :return ground position (number of points, lon,lat,h)
        :rtype 3D numpy.array
        """
        #faire une controle sur row / col !!!!
        # 0.5 < row < rowmax
        (kh,kb) = self.return_grid_index(alt)
        altbas  = self.alts_down[kb]
        althaut = self.alts_down[kh]
        dh = (alt - altbas)/(althaut - altbas)
        mats =  [self.gld_lon[kh:kb+1,:,:],self.gld_lat[kh:kb+1,:,:]]

        P = np.zeros((col.size, 3))
        P[:, 2] = alt
        dl = (row - self.row0)/self.steprow
        dc = (col - self.col0)/self.stepcol
        print('dl', dl.shape)
        [vlon,vlat] = interpol_bilin_vectorized(mats,self.nbrow,self.nbcol,dl,dc)

        P[:, 0] = (dh*vlon[0, :] +(1-dh)*vlon[1, :])
        P[:, 1]  = (dh*vlat[0, :] +(1-dh)*vlat[1, :])
        return P

    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param dtm : dtm model
        :type dtm  : shareloc.dtm
        :return ground position (lon,lat,h)
        :rtype numpy.array
        """
        visee = np.zeros((3,self.nbalt))
        vislonlat = self.interpolate_grid_in_plani(row, col)
        visee[0,:] = vislonlat[0]
        visee[1,:] = vislonlat[1]
        visee[2,:] = self.alts_down
        v = visee.T
        (code1, code2, PointB, dH3D) = dtm.checkCubeDTM(v)
        (code3,code4,Point_dtm) = dtm.intersection(v, PointB, dH3D)
        return Point_dtm

    def los_extrema(self,row,col,alt_min,alt_max):
        """
        compute los extrema
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param alt_min : los alt min
        :type alt_min  : float
        :param alt_max : los alt max
        :type alt_max : float
        :return los extrema
        :rtype numpy.array (2x3)
        """
        los_edges = np.zeros([2,3])
        los_edges[0,:] = self.direct_loc_h(row,col,alt_max)
        los_edges[1,:] = self.direct_loc_h(row,col,alt_min)
        return los_edges


    def fct_locdir_dtmopt(self,lig,col, dtm):
        """
        direct localization on 3D cube dtm
        :param lig :  line sensor position
        :type lig : float
        :param col :  column sensor position
        :type col : float
        :param dtm : dtm model
        :type dtm  : shareloc.dtm
        :return boolean true
        :rtype bool
        """
        visee = np.zeros((3,self.nbalt))
        vislonlat = self.fct_interp_visee_unitaire_gld(lig,col)
        visee[0,:] = vislonlat[0]
        visee[1,:] = vislonlat[1]
        visee[2,:] = self.alts_down
        v = visee.T
        (code, code2, PointB, dH3D) = dtm.checkCubeDTM(v)
        #(code,code4,Point_dtm) = self.intersection(v, PointB, dH3D,dtm)
        return code


    def interpolate_grid_in_plani(self, row, col):
        """
        interpolate positions on multi h grid
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :return interpolated positions
        :rtype list
        """
        dl = (row - self.row0)/self.steprow
        dc = (col - self.col0)/self.stepcol
        mats =  [self.gld_lon,self.gld_lat]
        res = interpol_bilin(mats,self.nbrow,self.nbcol,dl,dc)
        return res

    def interpolate_grid_in_altitude(self, nbrow, nbcol, nbalt=None):
        """
        interpolate equally spaced grid (in altitude)
        :param nbrow :  grid nb row
        :type nbrow : int
        :param nbcol :  grid nb col
        :type nbcol : int
        :param nbalt :  grid nb alt, of None self.nbalt is used instead
        :type nbalt : int
        :return equally spaced grid
        :rtype numpy.array
        """
        if not nbalt:
            nbalt = self.nbalt
            list_alts = self.alts_down
        else:
            list_alts = np.linspace(self.alts_down[0],self.alts_down[-1],nbalt)

        gld_lon = np.zeros((nbalt,nbrow,nbcol))
        gld_lat = np.zeros((nbalt,nbrow,nbcol))
        """genere un cube de visee interpole de nrow/ncol visee"""
        #row_max = self.row0 + self.steprow * (self.nbrow-1)
        #col_max = self.col0 + self.stepcol * (self.nbcol-1)

        steprow = (self.rowmax - self.row0)/(nbrow-1)
        stepcol = (self.colmax - self.col0)/(nbcol-1)

        for k,alt in enumerate(list_alts):

            res = self.direct_loc_grid_h(self.row0, self.col0, steprow, stepcol, nbrow, nbcol, alt)
            gld_lon[k] = res[0]
            gld_lat[k] = res[1]
        return gld_lon,gld_lat


    def direct_loc_grid_dtm(self, row0, col0, steprow, stepcol, nbrow, nbcol, dtm):
        """
         direct localization  grid on dtm
         :param row0 :  grid origin (row)
         :type row0 : int
         :param col0 :  grid origin (col)
         :type col0 : int
         :param steprow :  grid step (row)
         :type steprow : int
         :param stepcol :  grid step (col)
         :type stepcol : int
         :param nbrow :  grid nb row
         :type nbrow : int
         :param nbcol :  grid nb col
         :type nbcol : int
         :param dtm : dtm model
         :type dtm  : shareloc.dtm
         :return direct localization grid
         :rtype numpy.array
        """
        glddtm = np.zeros((3,nbrow,nbcol))
        visee = np.zeros((3,self.nbalt))
        for i in range(nbrow):
            for j in range(nbcol):
                col = col0 + stepcol * j
                row = row0 + steprow * i
                vislonlat = self.interpolate_grid_in_plani(row, col)
                visee[0,:] = vislonlat[0]
                visee[1,:] = vislonlat[1]
                visee[2,:] = self.alts_down
                v = visee.T
                (code1, code2, PointB, dH3D) = dtm.checkCubeDTM(v)
                (code3,code4,PointR) = dtm.intersection(v, PointB, dH3D)
                glddtm[:,i,j] = PointR
        return glddtm

    def return_grid_index(self, alt):
        """
         return layer index enclosing a given altitude
         :param alt :  altitude
         :type alt : float
         :return grid index (up,down)
         :rtype tuple
        """
        if alt > self.alts_down[0] :
            (indicehaut,indicebas) = (0,0)
        elif alt < self.alts_down[-1]:
            (indicehaut,indicebas) = (self.nbalt-1,self.nbalt-1)
        else:
            i = 0
            while i < self.nbalt and self.alts_down[i] >= alt:
                i+=1
            if i== self.nbalt: #pour gerer alt min
                i = self.nbalt - 1
            indicebas = i     #indice grille bas
            indicehaut = i-1   #indice grille haut
        return (indicehaut,indicebas)


    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """
         direct localization  grid at constant altitude
         :param row0 :  grid origin (row)
         :type row0 : int
         :param col0 :  grid origin (col)
         :type col0 : int
         :param steprow :  grid step (row)
         :type steprow : int
         :param stepcol :  grid step (col)
         :type stepcol : int
         :param nbrow :  grid nb row
         :type nbrow : int
         :param nbcol :  grid nb col
         :type nbcol : int
         :param alt : altitude of the grid
         :type alt  : float
         :return direct localization grid
         :rtype numpy.array
        """
        "TODO: check extrapolations"
        (kh,kb) = self.return_grid_index(alt)
        gldalt  = np.zeros((3,nbrow,nbcol))
        altbas  = self.alts_down[kb]
        althaut = self.alts_down[kh]

        dh = (alt - altbas)/(althaut - altbas)
        mats =  [self.gld_lon[kh:kb+1,:,:],self.gld_lat[kh:kb+1,:,:]]
        P     = np.zeros(3)
        P[2] = alt
        for i in range(nbrow):
            row = row0 + steprow*i
            dl = (row - self.row0)/self.steprow
            for j in range(nbcol):
                col = col0 + stepcol*j
                dc = (col - self.col0)/self.stepcol
                [vlon,vlat] = interpol_bilin(mats,self.nbrow,self.nbcol,dl,dc)
                P[0] = (dh*vlon[0] +(1-dh)*vlon[1])
                P[1] = (dh*vlat[0] +(1-dh)*vlat[1])
                gldalt[:,i,j] = P
        return gldalt

    def estimate_inverse_loc_predictor(self, nbrow_pred=3, nbcol_pred=3):
        """
        initialize inverse localization polynomial predictor
        it composed of 4 polynoms estimated on 5x5 grid at hmin and hmax :

        col_min = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_min = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        col_max = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_max = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        least squarred method is used to calculate coefficients, which are noramlized in [-1,1]
        :param nbrow_pred :  predictor nb row (3 by default)
        :type nbrow_pred : int
        :param nbcol_pred :  predictor nb col (3 by default)
        :type nbcol_pred : int
        """
        "why 5x5 grid ?"
        nb_alt     = 2
        nb_coeff   = 6
        nb_mes     = nbcol_pred*nbrow_pred

        col_norm    = np.linspace(-1.0,1.0,nbcol_pred)
        row_norm    = np.linspace(-1.0,1.0,nbrow_pred)
        gcol_norm, grow_norm = np.meshgrid(col_norm,row_norm)
        glon,glat = self.interpolate_grid_in_altitude(nbrow_pred, nbcol_pred, nb_alt)

        #normalisation des variables
        (glon_min,glon_max) = (glon.min(),glon.max())
        (glat_min,glat_max) = (glat.min(),glat.max())
        lon_ofset = (glon_max + glon_min)/2.0
        lon_scale = (glon_max - glon_min)/2.0
        lat_ofset = (glat_max + glat_min)/2.0
        lat_scale = (glat_max - glat_min)/2.0

        glon_norm = (glon - lon_ofset) / lon_scale
        glat_norm = (glat - lat_ofset) / lat_scale


        col_ofset = (self.colmax + self.col0)/2.0
        col_scale = (self.colmax - self.col0)/2.0
        row_ofset = (self.rowmax + self.row0)/2.0
        row_scale = (self.rowmax - self.row0)/2.0

        glon2     = glon_norm*glon_norm
        glonlat   = glon_norm*glat_norm
        glat2     = glat_norm*glat_norm

        Amin    = np.zeros((nb_mes,nb_coeff))
        Amax    = np.zeros((nb_mes,nb_coeff))
        Bcol    = np.zeros((nb_mes,1))
        Brow    = np.zeros((nb_mes,1))

        #resolution des moindres carres
        imes = 0
        for irow in range(nbrow_pred):
            for icol in range(nbcol_pred):

                Bcol[imes]       = gcol_norm[irow,icol]
                Brow[imes]       = grow_norm[irow,icol]
                Amin[imes,0]     = 1.0
                Amax[imes,0]     = 1.0
                Amin[imes,1]     = glon_norm[1,irow,icol]
                Amax[imes,1]     = glon_norm[0,irow,icol]
                Amin[imes,2]     = glat_norm[1,irow,icol]
                Amax[imes,2]     = glat_norm[0,irow,icol]
                Amin[imes,3]     = glon2[1,irow,icol]
                Amax[imes,3]     = glon2[0,irow,icol]
                Amin[imes,4]     = glat2[1,irow,icol]
                Amax[imes,4]     = glat2[0,irow,icol]
                Amin[imes,5]     = glonlat[1,irow,icol]
                Amax[imes,5]     = glonlat[0,irow,icol]
                imes +=1

        #Calcul des coeffcients
        matAmin       = np.array(Amin)
        matAmax       = np.array(Amax)

        tAAmin        = matAmin.T@matAmin
        tAAmax        = matAmax.T@matAmax
        tAAmin_inv    = np.linalg.inv(tAAmin)
        tAAmax_inv    = np.linalg.inv(tAAmax)

        coef_col_min = tAAmin_inv@matAmin.T@Bcol
        coef_row_min = tAAmin_inv@matAmin.T@Brow
        coef_col_max = tAAmax_inv@matAmax.T@Bcol
        coef_row_max = tAAmax_inv@matAmax.T@Brow

        setattr(self,'pred_col_min',coef_col_min.flatten())
        setattr(self,'pred_row_min',coef_row_min.flatten())
        setattr(self,'pred_col_max',coef_col_max.flatten())
        setattr(self,'pred_row_max',coef_row_max.flatten())
        setattr(self,'pred_ofset_scale_lon',  [lon_ofset , lon_scale])
        setattr(self,'pred_ofset_scale_lat',  [lat_ofset , lat_scale] )
        setattr(self,'pred_ofset_scale_row',  [row_ofset , row_scale] )
        setattr(self,'pred_ofset_scale_col',  [col_ofset , col_scale] )

    #-------------------------------------------------------------------------------------------------------------------------
    def inverse_loc_predictor(self, lon, lat, alt=0.0):
        """
        evaluate inverse localization predictor at a given geographic position
        :param lon : longitude
        :type lon : float
        :param lat : latitude
        :type lat : float
        :param alt : altitude (0.0 by default)
        :type alt : float
        :return sensor position (row,col, is extrapolated)
        :rtype tuple (float,float,boolean)
        """
        seuil_extrapol = 20.0
        extrapol = False
        altmin = self.alts_down[-1]
        altmax = self.alts_down[0]

        #normalisation
        lon_n = (lon - self.pred_ofset_scale_lon[0]) / self.pred_ofset_scale_lon[1]
        lat_n = (lat - self.pred_ofset_scale_lat[0]) / self.pred_ofset_scale_lat[1]
        if abs(lon_n) > (1+seuil_extrapol/100.0):
            #print "Attention, en extrapolation de 20% en longitude:",lon_n
            extrapol = True
        if abs(lat_n) > (1+seuil_extrapol/100.0):
            #print "Attention, en extrapolation de 20% en latitude:",lat_n
            extrapol = True

        #application polynome
        vect_sol = np.array([1,lon_n,lat_n,lon_n**2,lat_n**2,lon_n*lat_n])
        col_min = ((self.pred_col_min*vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_min = ((self.pred_row_min*vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]
        col_max = ((self.pred_col_max*vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_max = ((self.pred_row_max*vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]

        hx = (alt-altmin)/(altmax-altmin)
        col = (1-hx)*col_min + hx*col_max
        row = (1-hx)*row_min + hx*row_max

        if row > self.rowmax: row = self.rowmax
        if row < self.row0:   row = self.row0
        if col > self.colmax: col = self.colmax
        if col < self.col0:   col = self.col0
        return (row,col,extrapol)

    #-------------------------------------------------------------------------------------------------------------------------
    def inverse_partial_derivative(self, row, col, alt=0):
        """
        calculate partial derivative at a given geographic position
        it gives the matrix to apply to get sensor shifts from geographic ones
        it returns M matrix :
        [dcol,drow]T = M x [dlon,dlat]T
        dlon/dlat in microrad
        M is calculated on each node of grid
        M is necessary for direct localization inversion in iterative inverse loc
        :param lon : longitude
        :type lon : float
        :param lat : latitude
        :type lat : float
        :param alt : altitude (0.0 by default)
        :type alt : float
        :return matrix
        :rtype numpy.array
        """
        dl = (row - self.row0)/self.steprow
        dc = (col - self.col0)/self.stepcol
        il = int(np.floor(dl))
        ic = int(np.floor(dc))
        (kh,kb) = self.return_grid_index(alt)

        lon_h00 = self.gld_lon[kh,il  ,ic]
        lon_h01 = self.gld_lon[kh,il  ,ic+1]
        lon_h10 = self.gld_lon[kh,il+1,ic]

        lon_b00 = self.gld_lon[kb,il  ,ic]
        lon_b01 = self.gld_lon[kb,il  ,ic+1]
        lon_b10 = self.gld_lon[kb,il+1,ic]

        lat_h00 = self.gld_lat[kh,il  ,ic]
        lat_h01 = self.gld_lat[kh,il  ,ic+1]
        lat_h10 = self.gld_lat[kh,il+1,ic]

        lat_b00 = self.gld_lat[kb,il  ,ic]
        lat_b01 = self.gld_lat[kb,il  ,ic+1]
        lat_b10 = self.gld_lat[kb,il+1,ic]

        dlon_ch = np.deg2rad(lon_h01 - lon_h00)/self.stepcol
        dlon_cb = np.deg2rad(lon_b01 - lon_b00)/self.stepcol
        dlon_lh = np.deg2rad(lon_h10 - lon_h00)/self.steprow
        dlon_lb = np.deg2rad(lon_b10 - lon_b00)/self.steprow

        dlat_ch = np.deg2rad(lat_h01 - lat_h00)/self.stepcol
        dlat_cb = np.deg2rad(lat_b01 - lat_b00)/self.stepcol
        dlat_lh = np.deg2rad(lat_h10 - lat_h00)/self.steprow
        dlat_lb = np.deg2rad(lat_b10 - lat_b00)/self.steprow

        hx = (alt - self.alts_down[kb])/(self.alts_down[kh]-self.alts_down[kb])

        dlon_c = ((1-hx)*dlon_cb + (hx)*dlon_ch)*1e6
        dlat_c = ((1-hx)*dlat_cb + (hx)*dlat_ch)*1e6
        dlon_l = ((1-hx)*dlon_lb + (hx)*dlon_lh)*1e6
        dlat_l = ((1-hx)*dlat_lb + (hx)*dlat_lh)*1e6
        det = dlon_c*dlat_l - dlon_l*dlat_c
        if abs(det) > 0.000000000001:
            Matdp = np.array([[dlat_l,-dlon_l],[-dlat_c,dlon_c]])/det
        else:
            print("nul determinant")
        return Matdp
    #-------------------------------------------------------------------------------------------------------------------------
    def inverse_loc(self, lon,lat,alt = 0.0, nb_iterations = 15):
        """
        inverse localization at a given geographic position
        First initialize position,
        * apply inverse predictor lon,lat,at ->  col_0,row_0
        * direct loc col_0,row_0 -> lon_0, lat_0
        Then iterative process:
        * calculate geographic error dlon,dlat
        * calculate senor correction dlon,dlat -> dcol,drow
        * apply direct localization  -> lon_i,lat_i
        :param lon : longitude
        :type lon: float
        :param lat : latitude
        :type lat: float
        :param alt : altitude
        :type alt: float
        :param nb_iterations : max number of iterations (15 by default)
        :type nb_iterations : int
        :return sensor position (row,col, is_valid)
        :rtype tuple (float,float,boolean)
        """


        deg2mrad = np.deg2rad(1.0)*1e6
        k=0
        coslon = np.cos(np.deg2rad(lat))
        Rtx = 1e-12*6378000**2
        (row_i,col_i,extrapol) = self.inverse_loc_predictor(lon, lat, alt)
        erreur_m2    = 10.0
        point_valide = 0
        if not extrapol:
            #Processus iteratif
            #while erreur > seuil:1mm
            while (erreur_m2 > 1e-6) and (k < nb_iterations):
                #print k,row_i,col_i
                P = self.direct_loc_h(row_i, col_i, alt)
                dlon_microrad = (P[0] - lon)*deg2mrad
                dlat_microrad = (P[1] - lat)*deg2mrad
                erreur_m2 = Rtx*(dlat_microrad**2+(dlon_microrad*coslon)**2)
                dsol = np.array([dlon_microrad, dlat_microrad])
                mat_dp = self.inverse_partial_derivative(row_i, col_i, alt)
                dimg = mat_dp@dsol
                col_i += -dimg[0]
                row_i += -dimg[1]
                k +=1
                point_valide = 1
        return (row_i,col_i,point_valide)

    #-------------------------------------------------------------------------------

def coloc(gld_xH_src, gld_xH_dst, dtm, \
          l0_src, c0_src, steprow_src, stepcol_src, nbrow_src, nbcol_src):
    """
    colocalization grid on dtm
    localization on dtm from src grid, then inverse localization in right grid
    :param gld_xH_src : source grid
    :type gld_xH_src : shareloc.grid
    :param gld_xH_dst : destination grid
    :type gld_xH_dst : shareloc.grid
     :param l0_src :  grid origin (row)
     :type l0_src : int
     :param c0_src :  grid origin (col)
     :type c0_src : int
     :param steprow_src :  grid step (row)
     :type steprow_src : int
     :param stepcol_src :  grid step (col)
     :type stepcol_src : int
     :param nbrow_src :  grid nb row
     :type nbrow_src : int
     :param nbcol_src :  grid nb col
     :type nbcol_src : int
     :return colocalization grid
     :rtype numpy.array
    """
    gricoloc = np.zeros((3,nbrow_src,nbcol_src))
    for l in range(nbrow_src):
        row = l0_src + steprow_src * l
        for c in range(nbcol_src):
            col = c0_src + stepcol_src * c
            (lon,lat,alt) = gld_xH_src.direct_loc_dtm(row,col, dtm)
            Pdst = gld_xH_dst.inverse_loc(lon,lat,alt)
            gricoloc[:,l,c] = Pdst
    return gricoloc

