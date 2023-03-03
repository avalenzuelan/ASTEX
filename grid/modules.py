# This file contains classes that generates for the completeness tests
from astropy.io import fits

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

    def _open_image(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].data
    
    def _open_header(self, path_to_image, ext = 0):
        hdulist = fits.open(path_to_image)
        return hdulist[ext].header
    
    def get_image(self):
        return self.image
    
    def get_header(self):
        return self.header
    

