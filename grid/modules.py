# This file contains classes that generates for the completeness tests
from astropy.io import fits
import pandas as pd
import numpy as np
import re
import os

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
        self.original_stars_pos = None
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
        self.original_stars_pos = coords 
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
        the fortran code 'crowding.for' provided by Manuela.
        
        psf: float, the size of the psf fitting radius in pixels used in DAOPHOT
        photometry. This is used to avoid placing stars too close to each other.

        dist: float, the distance between the centers of the hexagons in pixels.

        hexshape: bool, if True, the grid is hexagonal, if False, the grid is
        square.
        """

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
    
    def move_random_1px(self, seed=None):
        """ This method moves each position randomly within one pixel, to
        ensure that the completeness of one position is not bias because 
        a star falls in exacly the same position of the pixel each time"""

        #TBD - this method should be improved to move the stars within a
        #radius of 1 pixel, not just 1 pixel in each direction

        #Create a copy of the original positions
        if seed is not None:
            np.random.seed(seed)
        self.stars_pos = self.original_stars_pos.copy()

        #Move the positions randomly within one pixel
        self.stars_pos[:,0] += np.random.uniform(-0.5,0.5, len(self.stars_pos))
        self.stars_pos[:,1] += np.random.uniform(-0.5,0.5, len(self.stars_pos))  
        return self.stars_pos
    
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
        #TBD - create random magnitudes with length uequal to the number of stars
        #in the grid. Flat random? or with a slope?
        return

    def separate_magfile(self, Grid, magfolder = '.', magfile = 'magfile',
                         use_leftover_stars = False,
                         remove_previous_files = False):

        """ This method is used to split the magnitudes files in lenghts 
        equal to the lengths of the grid. This recieves as an input the grid
        object, and returns a list of N files ready to be merged with the
        positions on each iteration """

        if remove_previous_files:
            os.system('rm %s/%s*' % (magfolder, magfile))   
        
        nstars_grid = Grid.nartstars_pergrid
        nstars = len(self.magfile)
        n_grids = nstars // nstars_grid
        n_leftover_stars = nstars % nstars_grid
        print(nstars, n_grids, n_leftover_stars)

        for i in range(n_grids):
            filename = magfile + '_%i.dat' % (i+1)
            path_to_save = os.path.join(magfolder, filename)
            if i == 0:
                df = self.magfile.iloc[:nstars_grid]
                df.to_csv(path_to_save, sep = ' ', index = False,
                          float_format='%.3f',)
            else:
                df = self.magfile.iloc[i*nstars_grid:(i+1)*nstars_grid]
                df.to_csv(path_to_save, sep = ' ', index = False,
                          float_format='%.3f',)
        
        #This part of the code creates a file with the leftover stars if the 
        #use_leftover_stars flag is True

        if use_leftover_stars:
            filename = magfile + '_%i.dat' % (n_grids+1)
            path_to_save = os.path.join(magfolder, filename)
            df = self.magfile.iloc[-n_leftover_stars:]
            df.to_csv(path_to_save, sep = ' ', index = False,
                                float_format='%.3f',)

        return True

    def create_addstar_files(self, Grid, band, magfolder = '.',
                             addstarfolder = '.',
                             magfile = 'magfile', addstarfile = 'addstar',
                             shuffle_mags = True):
        """ This method creates the addstar files for each grid. This method
        recieves as an input the magfile created by the separate_magfile method
        and the band to be used (H or Ks for HAWK-I fields)
        
        This method creates a file for each grid, with the positions and
        magnitudes of the stars in the grid. It reads each file in the magfolder
        and assigns a magnitude to each position.
        """

        hexgrid = Grid
        #I need to sort the files in the magfolder, because the last file is
        #the one with the leftover stars

        magfiles_list = os.listdir(magfolder)
        magfiles_list.sort(key=lambda s: int(re.search(r'\d+', s).group()))
        for file in magfiles_list:
            if file.startswith(magfile):
                path_to_open = os.path.join(magfolder, file)
                filename = file.replace(magfile, addstarfile)
                #add the band to the filename
                filename = filename.replace('.dat', '_%s.dat' % band)
                path_to_save = os.path.join(addstarfolder, filename)
                df = pd.read_csv(path_to_open, delim_whitespace = True)

                #Shuffle the magnitudes
                if shuffle_mags:
                    df = df.sample(frac=1).reset_index(drop=True)
                df['x'] = hexgrid.stars_pos[:,0]
                df['y'] = hexgrid.stars_pos[:,1]
                df.to_csv(path_to_save, sep = ' ', index = True,
                            float_format='%.3f',
                            columns = ['x', 'y', band])
        
        return True
    
    def create_addstar_fits(self, fits_file):
        """
        This method creates the fits files using the addstar files created
        by the create_addstar_files method. It uses DAOPHOT's addstar command
        to create the fits files.

        fits_file: string, path to the fits file where the stars will be added
        """

        #TBD - this method should be improved to create the fits files using
        #the addstar command from DAOPHOT
        return True