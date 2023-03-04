# This file contains classes that generates for the completeness tests
from astropy.io import fits
import numpy as np

class Grid():
    """ This class creates the grid to place the artificial stars.
    The positions are given in pixels. This class creates the *.pos file with
    the magnitude and positions in pixels, and add them to a *.fits image"""

    def __init__(self, path_to_image, extention = 0):

        self.path_to_image = path_to_image
        self.extention = extention
        self.image = self._open_image(self.path_to_image,
                                      ext = self.extention)
        self.header = self._open_header(self.path_to_image,
                                   ext = self.extention)
        self.imsize = self._get_imsize()
        self.nartstars = 0

    def _open_image(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].data
    
    def _open_header(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].header
    
    def _get_imsize(self):
        return self.header['NAXIS1'], self.header['NAXIS2']

    def get_image(self):
        return self.image
    
    def get_header(self):
        return self.header
    
    def create_grid(self, psf=float, dist=int, hexshape=True):
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

        
        xv, yv = np.meshgrid(np.arange(0,xmax, xmax/nx), 
                             np.arange(0,ymax, ymax/ny), 
                             sparse=False, indexing='xy')

        #this line is the one that moves every other row to the right, in order
        #to create hexagons
        if hexshape:
            xv[::2, :] += (xv[0][1] - xv[0][0])/2
        
        mask = mask_border(psf, xmax, ymax)
        xv, yv = xv[mask], yv[mask]

        self.nartstars = len(xv)
        return xv, yv

    