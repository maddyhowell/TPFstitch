import numpy as np
from lightkurve import KeplerTargetPixelFile
import pandas as pd
from .Construct_FITS_file import main as CFF
from lightkurve import search_targetpixelfile
from tqdm import tqdm
import sys
import K2_objects as K2Ob
from .null_TPFs import generate_null_tpf

def find_tpf_pos(tpf_ref, tpf, tpf_id, tpf_grid, check_flag):
    """ Find where you tpf is relative to the reference tpf (tpf_ref). Saves the tpf id into a grid"""

    ## Check to see if tpf is within 10 pixels of tpf_ref
    if check_flag:
        column_sum = np.shape(tpf_ref)[1] + np.shape(tpf)[1]
        row_sum = np.shape(tpf_ref)[2] + np.shape(tpf)[2]
        if np.abs(tpf_ref.column - tpf.column) > column_sum or np.abs(tpf_ref.row - tpf.row) > row_sum:
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


def load_TPF(tpf_id, download_MAST, local_path, campaign, verbose = True):

    if download_MAST and len(str(tpf_id).split('.')) == 1:

        pixelfile = search_targetpixelfile(f"ktwo{tpf_id}")
        tpf = pixelfile.download()
        if verbose:
            print(f'Loading TPF{tpf_id} from MAST')

    else:
        if len(str(tpf_id).split('.')) == 1:   ## if a path to a tpf is provided and its not in MAST format, do this
            tpf = KeplerTargetPixelFile(local_path + f'ktwo{int(tpf_id)}-c{campaign}_lpd-targ.fits.gz')
            if verbose:
                print(f'Loading TPF{int(tpf_id)} from local directory')

        elif 'fits' in tpf_id.split('.'):   ## if a tpf id is provided, do this
            tpf = KeplerTargetPixelFile(local_path + tpf_id)
            if verbose:
                print(f'Loading TPF with file name {tpf_id} from local directory')

        else:
            print('Something went wrong. Either not a TPF id or not an accepted file')
            ## at the moment it wont even reach here cos it will get an error on the statement above

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

def create_tpf_list(tpf_id_ref, K2_object, campaign):
    """ Create a list of the tpf ids that boarder the middle tpf. currently assumes a grid of 9 tpfs (i.e. pixel size of 15x15). Can write this to be more generalised"""

    if K2_object == 'M19' and campaign == '111':
        tpf_list = int(tpf_id_ref) + np.array([-1, +1, -9, -8, -7, +7, +8, +9])
    elif K2_object == 'M19' and campaign == '112':
        tpf_list = int(tpf_id_ref) + np.array([-1, +1, -9, +9, -10, -8, +8, +10])
    ## todo: the tpf list for M9 C111
    elif K2_object == 'M9':
        tpf_list = K2Ob.M9BoarderTPFs(tpf_id_ref, campaign)

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


def TPFstitch(tpf_id_ref, K2_object = None, tpf_list = None, file_name = 'Test.fits', download_MAST = False, local_path = '/home/user1/Documents/phd_third_year/TPFs/', campaign = '112', check_flag = True, superstamp_shape = None):
    """ main function to run code. tpf_list is a list of tpf_ids that you want to include in the 'Superstamp'.
    If not specified, it will find a list of a 9x9 tpf grid around the centre tpf.
    Object is to check if the TPFs in tpf_list are associated with the object.
    superstamp_shape used when there is tpf gaps in the resulting stitched tpf"""


    ## define the tpf list is not input. This only works with a 3x3 tpf tpf_grid
    ## if a tpf list is input and not a superstamp shape, it checks to make sure you did want a 3x3 superstamp. Option to change shape is needed

    tpf_ref = load_TPF(tpf_id_ref, download_MAST, local_path, campaign)
    tpf_id_ref_temp = tpf_id_ref

    if K2_object == None and tpf_list == None:   # no list of tpf id provided but also no K2 object provided so cant make a list of tpf ids
        print('You must either input a K2 object or list of tpfs')
        sys.exit()

    elif K2_object != None and tpf_list == None:   ## tpf id is provided but no list of tpfs ids provided so need to make one
        tpf_list = create_tpf_list(tpf_id_ref, K2_object, campaign)
        print('Stitching the following TPFs into a new superstamp:', tpf_list)
        print('------------------------------------------------------')
    elif tpf_list != None:
        tpf_id_ref = 1

    tpf_grid = np.zeros((3,3))

    tpf_grid[1][1] = tpf_id_ref

    tpf_store = {tpf_id_ref: tpf_ref}

    for idx, tpf_id in enumerate(tpf_list):
        # create a grid which contains the locations of each tpf
        tpf_temp = load_TPF(tpf_id, download_MAST, local_path, campaign)

        if type(tpf_id) != np.int64:  # If a tpf is a file path then give a random idex for tpf id to store
            tpf_id = idx + 2

        tpf_store.update({tpf_id: tpf_temp})

        tpf_grid = find_tpf_pos(tpf_ref, tpf_temp, tpf_id, tpf_grid, check_flag)

    if superstamp_shape == None:
        tpf_grid = tpf_grid[:,~np.all(tpf_grid == 0, axis=0)]
        tpf_grid = tpf_grid[~np.all(tpf_grid == 0, axis=1)]
        # print(tpf_grid)
    else:
        if superstamp_shape != np.shape(tpf_grid):
            m = 0  ### TODO: need to code this up to force a superstamp shape if not 3x3

    ## find where no tpf data exists
    no_tpfs = []
    for ii,row in enumerate(tpf_grid):
        for jj,element in enumerate(row):
            if int(element) == 0:
                no_tpfs.append((ii,jj))

    ## create null tpfs
    tpf_ref_col, tpf_ref_row = tpf_ref.column, tpf_ref.row
    for idx,thing in enumerate(no_tpfs):
        filename = generate_null_tpf(tpf_ref, thing, (tpf_ref_col, tpf_ref_row))

        tpf_temp = KeplerTargetPixelFile(filename)
        print(f'Loading null TPF with file name {filename} from local directory')

        tpf_id = 100 + idx
        tpf_store.update({tpf_id: tpf_temp})

        tpf_grid[thing[0]][thing[1]] = tpf_id

    # print(tpf_grid)
    ## to make sure HDU hasn't been rewritten in previous step
    tpf_ref = load_TPF(tpf_id_ref_temp, download_MAST, local_path, campaign, verbose = False)

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
    # print(np.shape(flux))

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

    # rb_level = append_arrays(tpf_ref, tpf_grid, tpf_store, variable = 'RB_level')           ## Append rb level
    # data.update({'RB_LEVEL': rb_level})
    # print('RB levels array has been appended')

    print('\n')
    print('New Superstamp has shape:', np.shape(flux))

    ### Construst new fits file -------------------------------------------------------------------------------------------

    tpf_shape = (np.shape(flux)[1],np.shape(flux)[2])
    bottom_tpf = tpf_store[int(tpf_grid[0][0])]

    print('Bottom_tpf:', int(tpf_grid[0][0]))
    CFF(data, bottom_tpf.hdu, new_tpf_name = file_name, tpf_shape = tpf_shape, verbose = False)

    return
