from grid import Grid, ArtificialStars
import matplotlib.pyplot as plt
import os
import re

def main():
    global hexgrid, artstars 
    hexgrid = Grid(path_to_image='../Data/NPL054/c1/NPL054_Ks_160s_c1.fits')
    xv, yv =  hexgrid.create_grid(20, 2)
    hexgrid.move_random_1px()
    print('Each grid supports %i stars' % hexgrid.nartstars_pergrid)

    # fig, ax = plt.subplots(figsize=(8, 8))
    # ax.scatter(xv, yv, s=1)
    # ax.axis('equal')
    # plt.show()

    artstars = ArtificialStars(hexgrid, '../Data/NPL054/c1/NPL054_c1.dat.rot')
    # print('Creating magfiles')
    # artstars.separate_magfile(hexgrid, magfolder='Data/NPL054/c1/magfiles',
    #                           remove_previous_files=True)

    # print('Creating addstarfiles')
    # artstars.create_addstar_files(hexgrid, 'Ks', 
    #                                 magfolder='../Data/NPL054/c1/magfiles',
    #                                 addstarfolder='Data/NPL054/c1/addfiles',
    #                                 remove_previous_files=True)
    
    # artstars.create_addstar_files(hexgrid, 'H', 
    #                                 magfolder='../Data/NPL054/c1/magfiles',
    #                                 addstarfolder='Data/NPL054/c1/addfiles',
    #                                 remove_previous_files=False)
    
    addfiles_list = os.listdir('../Data/NPL054/c1/addfiles')
    addfiles_list.sort(key=lambda s: int(re.search(r'\d+', s).group()))
    for addfile in os.listdir('../Data/NPL054/c1/addfiles'):
        if addfile.startswith('addstar') and \
            addfile.endswith('_Ks.dat'):
            
            n_file = int(re.search(r'\d+', addfile).group())
            print(addfile, n_file)
            artstars.create_fits_artstars('../Data/NPL054/c1/', 
                                        'NPL054_Ks_160s_c1.fits',
                                        n_file=n_file,
                                        phots_per_adu=1.88,
                                        addstarfolder='addfiles',
                                        addstarfile=addfile,)
            
            artstars.create_fits_artstars('../Data/NPL054/c1/', 
                                        'NPL054_Ks_10s_c1.fits',
                                        n_file=n_file,
                                        phots_per_adu=1.88,
                                        addstarfolder='addfiles',
                                        addstarfile=addfile,)
            
        elif addfile.startswith('addstar') and \
            addfile.endswith('_H.dat'):
            
            n_file = int(re.search(r'\d+', addfile).group())
            print(addfile, n_file)
            artstars.create_fits_artstars('../Data/NPL054/c1/', 
                                        'NPL054_H_480s_c1.fits',
                                        n_file=n_file,
                                        phots_per_adu=1.88,
                                        addstarfolder='addfiles',
                                        addstarfile=addfile,)
            
            artstars.create_fits_artstars('../Data/NPL054/c1/', 
                                        'NPL054_H_10s_c1.fits',
                                        n_file=n_file,
                                        phots_per_adu=1.88,
                                        addstarfolder='addfiles',
                                        addstarfile=addfile,)
        # break

if __name__ == '__main__':
    main()