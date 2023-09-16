# This file contains classes that generates for the completeness tests
from astropy.io import fits
import pandas as pd
import numpy as np
import re
import os

class Grids():

    """ 

    This class creates the grid to place the artificial stars.
    The positions are given in pixels. This class creates the positions with
    in pixels based on the specifications of a *.fits image
    
    Then, it uses a *.mch file to transform the positions from pixels to the
    coordinate system of another image. This is done using the transformation
    matrix provided by DAOMATCH/DAOMASTER. This class also creates the grid in 
    the coordinate system of the other image.

    """

    def __init__(self, 
                 path_to_folder, 
                 master_filename,
                 list_of_filenames = [],
                 mch_filename = None,
                 extention = 0):

        self.path_to_folder = path_to_folder
        self.filename = master_filename
        self.path_to_image = os.path.join(path_to_folder, 
                                          master_filename + '.fits')
        self.extention = extention
        self._image = self._open_image(self.path_to_image,
                                      ext = self.extention)
        self._header = self._open_header(self.path_to_image,
                                   ext = self.extention)
        self.imsize = self._get_imsize()
        self.list_of_filenames = list_of_filenames
        self.list_of_paths = [os.path.join(path_to_folder, filename) \
                                for filename in self.list_of_filenames]
        self.mch_filename = mch_filename
        self.mch_path = os.path.join(path_to_folder, self.mch_filename)
        self.mch = self._read_mch()

        #output attributes
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
        self.stars_pos = {self.filename : coords}
        self.original_stars_pos = coords 
        return coords
    
    def _read_mch(self, 
                  path_to_mch = None):
        
        if path_to_mch is None:
            path_to_mch = self.mch_path

        #read the mch file, edit the index, and rename the columns
        mch = pd.read_csv(path_to_mch, delim_whitespace=True, header=None)
        mch.drop(1, axis=1, inplace=True)
        mch[0] = mch[0].apply(lambda x: x[1:-4])
        mch.columns = ['FILE', 'A', 'B', 'C', 'D', 'E', 'F', 'IDK1', 'IDK2']
        mch.set_index('FILE', inplace=True)

        return mch

    @property
    def image(self):
        return self._image
    
    @property
    def header(self):
        return self._header
    
    def create_grid(self, 
                    psf : float, 
                    dist = 2, 
                    hexshape = True):
        
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

        stars_pos = self.original_stars_pos.copy()

        #Move the positions randomly within one pixel
        stars_pos[:,0] += np.random.uniform(-0.5,0.5, len(stars_pos))
        stars_pos[:,1] += np.random.uniform(-0.5,0.5, len(stars_pos))

        self.stars_pos[self.filename] = stars_pos 
        return stars_pos
    
    def transform_to(self, 
                     to_filename : str):

        """ This method transforms the positions in pixels to positions in
        the coordinate of other image. This method recieves as an input the 
        mch file, which contains the transformation matrix from pixels to 
        privided by DAOMATCH. This method returns a file with the positions 
        in the new reference frame """

        mch = self.mch

        #name the constants and prepare the arrays
        mch_row = mch.loc[to_filename]
        a, b, c, d, e, f = mch_row['A'], mch_row['B'], mch_row['C'], \
            mch_row['D'], mch_row['E'], mch_row['F']
        xv = self.stars_pos[self.filename][:,0] 
        yv = self.stars_pos[self.filename][:,1]
        #transform the positions by solving the system of equations
        # Ax = b

        A = np.array([[c, e], [d, f]])
        b = np.array([xv - a, yv - b])
        x_new, y_new = np.linalg.solve(A, b)

        #save the new positions in the dictionary
        self.stars_pos[to_filename] = np.stack((x_new, y_new), axis = 1)
        return x_new, y_new


class Magnitudes():

    """ This class reads the photometry catalogues and creates a pandas
    dataframe with the magnitudes of the stars in the grid. This class
    recieves as an input the path to the folder where the photometry
    catalogues are located, the name of the photometry catalogue, and
    a list with the names of the alf files. The photometry catalogue
    should be a file with the magnitudes of the stars in the grid, and
    the alf files should be the files with the magnitudes of the stars
    in the image. This class returns a pandas dataframe with the
    magnitudes of the stars in the grid """

    def __init__(self,
                 path_to_folder : str,
                 phot_cat_name : str,
                 master_filename : str,
                 list_of_filenames : list) -> None:
        
        self.path_to_folder = path_to_folder
        self.phot_catalogue_name = phot_cat_name
        self.path_to_phot_cat = os.path.join(path_to_folder, phot_cat_name)
        self.master_filename = master_filename
        self.filenames = list_of_filenames
        self.filenames.insert(0, self.master_filename)
        self.alf_names = [filename + '.alf' for filename in self.filenames]
        self.path_to_alfs = [os.path.join(path_to_folder, alf_name) \
                                for alf_name in self.alf_names]
        self.phot_cat = self._open_phot_catalogue(self.path_to_phot_cat)

    def _open_phot_catalogue(self, path_to_phot_catalogue):
        names = ['ID', 'RA', 'Dec', 'H', 'eH', 'Ks', 'eKs', 'CHI', 'SHARP']
        df = pd.read_csv(path_to_phot_catalogue, 
                         delim_whitespace=True,
                         names=names,
                         skiprows=[0,1,2,3],
                         usecols=['ID'])
        df.set_index('ID', inplace=True)
        return df
    
    def get_instrumental_magnitudes(self, filename):
        
        names = ['ID', 'X', 'Y', 'MAG', 'eMAG', 'SKY', 'ITER', 'CHI', 'SHARP']  

        for alf in self.path_to_alfs:
            if filename in alf:
                alf_mag = pd.read_csv(alf, 
                                      delim_whitespace=True, 
                                      names=names,
                                      skiprows=2)
                alf_mag.set_index('ID', inplace=True)
                #self.phot_cat[filename] = alf_mag.loc[self.phot_cat.index].MAG
        
        return alf_mag
    


class ArtificialStars():

    def __init__(self, 
                 grid, 
                 path_to_folder,
                 mag_filename, 
                 random_mags = False):
        
        self.grid = grid
        self.path_to_magfile = os.path.join(path_to_folder, mag_filename)
        self._magfile = self._open_magfile(self.path_to_magfile)

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
                             magfile = 'magfile', 
                             addstarfile = 'addstar',
                             shuffle_mags = True,
                             remove_previous_files = False):
        """ This method creates the addstar files for each grid. This method
        recieves as an input the magfile created by the separate_magfile method
        and the band to be used (H or Ks for HAWK-I fields)
        
        This method creates a file for each grid, with the positions and
        magnitudes of the stars in the grid. It reads each file in the magfolder
        and assigns a magnitude to each position.
        """
        if remove_previous_files:
            os.system('rm %s/%s*' % (addstarfolder, addstarfile))

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
    
    def create_fits_artstars(self, 
                             path_to_folder, 
                             fits_file,
                             phots_per_adu = 1.0,
                             n_file = 0, 
                             addstarfolder = '.',
                             addstarfile = ''
                             ):
        """
        This method creates the fits files using the addstar files created
        by the create_addstar_files method. It uses DAOPHOT's addstar command
        to create the fits files.

        - fits_file: string, path to the fits file where the stars will be added
        - phots_per_adu: float, the number of photons per adu in the image. This
        is used to estimate the Poisson noise in the artificial stars. Can be
        calculated as GAIN/QE, where GAIN is the gain of the detector and QE
        is the quantum efficiency of the detector. If you want the artificial
        stars to not have any noise, set this to a very high number, like 99_999
        
        """
        origin_path = os.getcwd()
        dirpath = os.path.dirname(path_to_folder)

        # To run daophot, I need to change to the directory where the fits file
        # and psf file are located
        os.chdir(dirpath)
        filepath = os.path.splitext(fits_file)[0]

        path_to_addstar = os.path.join(addstarfolder, addstarfile)

        psf_filepath = filepath + '.psf'
        new_filepath = filepath + '_add' + str(n_file) + '.fits'
        cmd_addstar = 'daophot << EOF\n\nat %s\naddstar\n%s\n100\n%1.2f\n' + \
            '%s\n%s\nEOF\n' \
            % (filepath, psf_filepath, phots_per_adu, path_to_addstar, \
               new_filepath)
        print(cmd_addstar)
        os.system(cmd_addstar)

        #Return to the original directory
        os.chdir(origin_path)
        return True