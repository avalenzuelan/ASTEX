# This file contains classes that generates for the completeness tests
from astropy.io import fits
import pandas as pd
import numpy as np

class Grid():
    """ This class creates the grid to place the artificial stars.
    The positions are given in pixels. This class creates the positions with
    in pixels based on the specifications of each *.fits image"""

    def __init__(self, path_to_image, extention = 0):

        self.path_to_image = path_to_image
        self.extention = extention
        self._image = self._open_image(self.path_to_image,
                                      ext = self.extention)
        self._header = self._open_header(self.path_to_image,
                                   ext = self.extention)
        self.imsize = self._get_imsize()
        self._nartstars_pergrid = None
        self.stars_pos = None

    def _open_image(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].data
    
    def _open_header(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].header
    
    @property
    def nartstars_pergrid(self):
        return self._nartstars_pergrid

    def _get_imsize(self):
        return self.header['NAXIS1'], self.header['NAXIS2']
    
    def _create_stars_pos(self, x, y):
        coords = np.stack((x,y), axis = 1)
        self.stars_pos = coords 
        return coords
    

    @property
    def image(self):
        return self._image
    
    @property
    def header(self):
        return self._header
    
    def create_grid(self, psf=float, dist=2, hexshape=True):
        """ This method creates the grid. By default, the grid is hexagonal,
        but can be square if hexshape == False. This part of the code is based 
        the fortran code 'crowding.for' provided by Manuela."""

        def mask_border(psf, xmax, ymax):
            """ This subroutine creates a mask that eliminates stars that are
             too close to the border of the image """
            min_border = psf + 1
            xmax_border, ymax_border =  xmax - (psf + 1), ymax - (psf + 1)
            mask = (xv > min_border) & (xv < xmax_border) \
                & (yv > min_border) & (yv < ymax_border)
            return mask

        ymax, xmax = self.imsize
        xr, yr = np.random.uniform(), np.random.uniform() 
        xr, yr = xr*(dist*psf), yr*(dist*psf)
        nx = int((xmax - xr) // (dist*psf))
        ny = int((ymax - yr) // (dist*psf))

        #this command creates the grid given the input parameters
        xv, yv = np.meshgrid(np.arange(0,xmax, xmax/nx), 
                             np.arange(0,ymax, ymax/ny), 
                             sparse=False, indexing='xy')

        #this line is the one that moves every other row to the right, in order
        #to create hexagons
        if hexshape:
            xv[::2, :] += (xv[0][1] - xv[0][0])/2
        
        mask = mask_border(psf, xmax, ymax)
        xv, yv = xv[mask], yv[mask]

        self._nartstars_pergrid = len(xv)
        self._create_stars_pos(xv, yv)
        return xv, yv
    
    def move_random_1px(self)
        """ This method moves each position randomly within one pixel, to
        ensure that the completeness of one position is not bias because 
        a star falls in exacly the same position of the pixel each time"""
        return    
class ArtificialStars():

    def __init__(self, grid, path_to_magfile, random_mags = False):
        
        self.grid = grid
        self._magfile = self._open_magfile(path_to_magfile)

    def _open_magfile(self, path_to_magfile):
        df = pd.read_csv(path_to_magfile,
                         skiprows = [0,1,2],
                         delim_whitespace = True,
                         usecols = [3, 13],
                         names = ['H', 'Ks'])
        return df

    @property
    def magfile(self):
        return self._magfile
    
    @magfile.setter
    def magfile(self, value):
        #do something
        return self._magfile
    
    def random_magnitudes(self, min_mag, max_mag):
        #TBD - create random magnitudes with length equal to the number of stars
        #in the grid. Flat random? or with a slope?
        return

    def separate_magfile(self, grid):
        """ This method is used to split the magnitudes files in lenghts 
        equal to the lengths of the grid. This recieves as an input the grid
        object, and returns a list of N dataframes ready to be merged with the
        positions on each iteration """

        nstars_grid = grid.nartstars_pergrid
        nstars = len(self.magfile)
        n_grids = nstars // nstars_grid
        print(n_grids)
        return

    def create_pos_files(self, magfile, band):
        return