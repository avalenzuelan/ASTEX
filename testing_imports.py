from grid import Grid, ArtificialStars
import matplotlib.pyplot as plt

def main():
   global hexgrid, artstars 
   hexgrid = Grid(path_to_image='Data/NPL054/blank_fits.fits')
   xv, yv =  hexgrid.create_grid(20, 2)
   print('Each grid supports %i stars' % hexgrid.nartstars_pergrid)

   fig, ax = plt.subplots(figsize=(8, 8))
   ax.scatter(xv, yv, s=1)
   ax.axis('equal')
   plt.show()

   artstars = ArtificialStars(hexgrid, 'Data/NPL054/NPL054_c1.dat.rot')
   print(artstars.magfile)
   artstars.separate_magfile(hexgrid)

if __name__ == '__main__':
   main()