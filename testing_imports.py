from grid import Grid, ArtificialStars
import matplotlib.pyplot as plt

def main():
    global hexgrid, artstars 
    hexgrid = Grid(path_to_image='Data/NPL054/blank_fits.fits')
    xv, yv =  hexgrid.create_grid(20, 2)
    hexgrid.move_random_1px()
    print('Each grid supports %i stars' % hexgrid.nartstars_pergrid)

    # fig, ax = plt.subplots(figsize=(8, 8))
    # ax.scatter(xv, yv, s=1)
    # ax.axis('equal')
    # plt.show()

    artstars = ArtificialStars(hexgrid, 'Data/NPL054/NPL054_c1.dat.rot')
    # print('Creating magfiles')
    # artstars.separate_magfile(hexgrid, magfolder='test_files/magfiles')

    print('Creating addstarfiles')
    artstars.create_addstar_files(hexgrid, 'Ks', magfolder='test_files/magfiles',
                                    addstarfolder='test_files/addstarfiles')
    artstars.create_addstar_files(hexgrid, 'H', magfolder='test_files/magfiles',
                                    addstarfolder='test_files/addstarfiles')
    
    

if __name__ == '__main__':
    main()