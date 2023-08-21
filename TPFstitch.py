import numpy as np
from lightkurve import KeplerTargetPixelFile
import pandas as pd
import Construct_FITS_file as CFF
from lightkurve import search_targetpixelfile
from tqdm import tqdm
import sys
import K2_objects as K2Ob

def find_tpf_pos(tpf_ref, tpf, tpf_id, tpf_grid):
    """ Find where you tpf is relative to the reference tpf (tpf_ref). Saves the tpf id into a grid"""

    ## Check to see if tpf is within 10 pixels of tpf_ref
    if np.abs(tpf_ref.column - tpf.column) > 10 or np.abs(tpf_ref.row - tpf.row) > 10:
        print(f"TPF {tpf_id} does not boarder the reference TPF")
        return tpf_grid

    flagx = False   # share x-axis
    flagy = False   # share y-axis
    corner_flag = False


    # find if they share an x-axis or a y-axis
    if tpf_ref.column == tpf.column:
        flagx = True

    elif tpf_ref.row == tpf.row:
        flagy = True

    else:
        corner_flag = True


    # find if its up/down
    if flagx:
        if tpf_ref.row > tpf.row:
            tpf_grid[0][1] = tpf_id
        else:
            tpf_grid[2][1] = tpf_id

        return tpf_grid

    # find if its left/right
    elif flagy:
        if tpf_ref.column > tpf.column:
            tpf_grid[1][0] = tpf_id
        else:
            tpf_grid[1][2] = tpf_id

        return tpf_grid

    elif corner_flag:
        ii, jj = corners(tpf_ref, tpf)

        tpf_grid[ii][jj] = tpf_id
        return tpf_grid
    else:
        print("Tpfs don't connect")
        return tpf_grid


def corners(tpf_ref, tpf):

    flagleft = False
    flagright = False

    if tpf_ref.column > tpf.column:
        flagleft = True
    else:
        flagright = True

    if flagleft:
        if tpf_ref.row > tpf.row:
            ii, jj = 0,0
        else:
            ii, jj = 2,0

    elif flagright:
        if tpf_ref.row > tpf.row:
            ii, jj = 0,2
        else:
            ii, jj = 2,2

    return ii, jj


def load_TPF(tpf_id, K2_object, download_MAST, local_path):

    __, campaign = K2Ob.TPF_ids(K2_object)

    if download_MAST:
        print(f'Loading TPF{tpf_id} from MAST')
        pixelfile = search_targetpixelfile(f"ktwo{tpf_id}")
        tpf = pixelfile.download()

    else:
        print(f'Loading TPF{tpf_id} from local directory')
        tpf = KeplerTargetPixelFile(local_path + f'ktwo{tpf_id}-c{campaign}_lpd-targ.fits.gz')

    return tpf


def append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'FLUX'):
    superstamp = []
    shape = np.shape(tpf_grid)
    # print('TPF_grid', tpf_grid)
    # print('append_arrays', shape)
    n_rows = shape[0]
    n_cols = shape[1]

    row_big = []

    for idx in range(tpf_ref.shape[0]):

        for jj in range(n_rows):
            row_temp = tpf_store[int(tpf_grid[jj][0])].hdu[1].data[variable][idx]
            # print(f'tpf_id: {int(tpf_grid[jj][0])}, at ({jj, 0})')
            ii = 1
            while ii <= n_cols-1:
                # print(f'tpf_id: {int(tpf_grid[jj][ii])}, at ({jj, ii})')
                row_temp = np.append(row_temp, tpf_store[int(tpf_grid[jj][ii])].hdu[1].data[variable][idx], axis = 1)
                ii += 1

            if jj == 0:
                superstamp_temp = row_temp
            else:
                superstamp_temp = np.append(superstamp_temp, row_temp, axis = 0)
        del(row_temp)
        superstamp.append(superstamp_temp)

    return superstamp

def create_tpf_list(tpf_id_ref, K2_object):
    """ Create a list of the tpf ids that boarder the middle tpf. currently assumes a grid of 9 tpfs (i.e. pixel size of 15x15). Can write this to be more generalised"""

    if K2_object == 'M19':
        tpf_list = tpf_id_ref + np.array([-1, +1, -9, +9, -10, -8, +8, +10])

    elif K2_object == 'M9':
        tpf_list = K2Ob.M9BoarderTPFs(tpf_id_ref)

    ## check to see if tpfs in tpf list are associated with the object. if not the code stops
    tpf_list = checkTPFs(list(tpf_list), K2_object)

    return np.array(tpf_list)

def checkTPFs(tpf_list, K2_object):
    """ Check if the TPFs in tpf_list are associated with the object. """

    tpfs, campaign = K2Ob.TPF_ids(K2_object)

    tpf_list_copy = tpf_list.copy()

    for tpf_id in tpf_list:
        if tpf_id not in tpfs:
            tpf_list_copy.remove(tpf_id)

    return tpf_list_copy


def TPFstitch(tpf_id_ref, K2_object = None, tpf_list = None, superstamp_shape = None, file_name = 'Test.fits', download_MAST = True, local_path = '../TPFs/'):
    """ main function to run code. tpf_list is a list of tpf_ids that you want to include in the 'Superstamp'.
    If not specified, it will find a list of a 9x9 tpf grid around the centre tpf.
    Object is to check if the TPFs in tpf_list are associated with the object. """

    tpf_ref = load_TPF(tpf_id_ref, K2_object, download_MAST, local_path)

    ## define the tpf list is not input. This only works with a 3x3 tpf tpf_grid
    ## if a tpf list is input and not a superstamp shape, it checks to make sure you did want a 3x3 superstamp. Option to change shape is needed
    if tpf_list != None:
        print('Stitching the following TPFs into a new superstamp:', tpf_list)
        print('------------------------------------------------------')

        if superstamp_shape == None:
            print('You need to input a superstamp shape. Please try again.')
            sys.exit()

        while len(tpf_list)+1 != superstamp_shape[0]*superstamp_shape[1]:
            print(f'Superstamp shape of {superstamp_shape} does not equal the number of tpf ids provided')
            print('What superstamp shape did you want?')
            n = int(input('Number of rows:'))
            m = int(input('Number of columns:'))
            superstamp_shape = (n,m)
            print('\n')
        print(f'Making a superstamp with shape {superstamp_shape}')
        print('\n')


    elif K2_object != None:
        tpf_list = create_tpf_list(tpf_id_ref, K2_object)
        print('Stitching the following TPFs into a new superstamp:', tpf_list)
        print('------------------------------------------------------')

        ### is superstamp_shape necessary?


    else:
        print('You must either input a K2 object or list of tpfs')

    tpf_grid = np.zeros((3,3))

    tpf_grid[1][1] = tpf_id_ref

    tpf_store = {tpf_id_ref: tpf_ref}

    for tpf_id in tpf_list:
        # create a grid which contains the locations of each tpf
        tpf_temp = load_TPF(tpf_id, K2_object, download_MAST, local_path)

        tpf_store.update({tpf_id: tpf_temp})

        tpf_grid = find_tpf_pos(tpf_ref, tpf_temp, tpf_id, tpf_grid)


    ### if superstamp shape != (3,3) need to remove the zero-rows
    tpf_grid = tpf_grid[:,~np.all(tpf_grid == 0, axis=0)]
    tpf_grid = tpf_grid[~np.all(tpf_grid == 0, axis=1)]

    # if np.shape(tpf_grid) != superstamp_shape:   # check to make sure the tpf_grid shape == the superstamp shape
    #
    #     print(f'Something went wrong with the tpf grid shape. tpf grid shape = {np.shape(tpf_grid)}')
    #     sys.exit()

    print('\n')
    data = {}

    time = tpf_ref.hdu[1].data['TIME']         ## Append time
    data.update({'TIME': time})
    print('Time array has been appended')

    timecorr = tpf_ref.hdu[1].data['TIMECORR']         ## Append timecorr
    data.update({'TIMECORR': timecorr})
    print('Time correction array has been appended')

    cadenceno = tpf_ref.hdu[1].data['CADENCENO']         ## Append cadenceno
    data.update({'CADENCENO': cadenceno})
    print('Cadence number array has been appended')

    raw_counts = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'RAW_CNTS')         ## Append raw counts
    data.update({'RAW_CNTS': raw_counts})
    print('Raw counts array has been appended')

    flux = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'FLUX')                   ## Append flux
    data.update({'FLUX': flux})
    print('Flux array has been appended')

    flux_err = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'FLUX_ERR')           ## Append flux_err
    data.update({'FLUX_ERR': flux_err})
    print('Flux error array has been appended')

    flux_bkg = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'FLUX_BKG')           ## Append flux_bkg
    data.update({'FLUX_BKG': flux_bkg})
    print('Flux background array has been appended')

    flux_bkg_err = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'FLUX_BKG_ERR')   ## Append flux_bkg_err
    data.update({'FLUX_BKG_ERR': flux_bkg_err})
    print('Flux background error array has been appended')

    cosmic_rays = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'COSMIC_RAYS')     ## Append cosmic rays
    data.update({'COSMIC_RAYS': cosmic_rays})
    print('Cosmic rays array has been appended')

    quality = tpf_ref.hdu[1].data['QUALITY']         ## Append quality
    data.update({'QUALITY': quality})
    print('Quality array has been appended')

    poscorr1 = tpf_ref.hdu[1].data['POS_CORR1']         ## Append poscorr1
    data.update({'POS_CORR1': poscorr1})
    print('Position correction 1 array has been appended')

    poscorr2 = tpf_ref.hdu[1].data['POS_CORR2']         ## Append poscorr2
    data.update({'POS_CORR2': poscorr2})
    print('Position correction 2 array has been appended')

    rb_level = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'RB_level')           ## Append rb level
    data.update({'RB_LEVEL': rb_level})
    print('RB levels array has been appended')

    print('\n')
    print('New Superstamp has shape:', np.shape(rb_level))
    ### Construst new fits file -------------------------------------------------------------------------------------------


    tpf_shape = (np.shape(rb_level)[1],np.shape(rb_level)[2])
    bottom_tpf = tpf_store[int(tpf_grid[0][0])]
    CFF.main(data, bottom_tpf.hdu, new_tpf_name = file_name, tpf_shape = tpf_shape, verbose = False)

    return
