import numpy as np
from .Construct_FITS_file import main as CFF

"""Below creates a tpf full of np.nans but for the same number of cadences as the inputed existing tpf
Things to consider, the new null tpf will be the same size as the existing tpf inputted
Need to shift the pixels in a way (compared to existing tpf I guess?) , need to also consider ra and decs
!!!!! BE CAREFUL: TPF MUTABLE!!! EVERYTIME I RUN CFF.MAIN IT CHANGES TPF HDU!!!!
"""

def compass_position(idx_pos):
    if idx_pos == (0,0):
        compass_pos = 'SW'
    elif idx_pos == (0,1):
        compass_pos = 'S'
    elif idx_pos == (0,2):
        compass_pos = 'SE'
    elif idx_pos == (1,0):
        compass_pos = 'W'
    elif idx_pos == (1,2):
        compass_pos = 'E'
    elif idx_pos == (2,0):
        compass_pos = 'NW'
    elif idx_pos == (2,1):
        compass_pos = 'N'
    elif idx_pos == (2,2):
        compass_pos = 'NE'

    return compass_pos


def find_radecs(pixel_col, pixel_row, tpf):
    coords = [pixel_col - tpf.column, pixel_row - tpf.row]
    coords = np.vstack(coords).T

    radecs = tpf.wcs.all_pix2world(coords, 0)[0]

    return radecs

def generate_null_tpf(reference_tpf, position, original_tpfpos):
    data = {}

    time = reference_tpf.hdu[1].data['TIME']
    data.update({'TIME': time})

    timecorr = reference_tpf.hdu[1].data['TIMECORR']
    data.update({'TIMECORR': timecorr})

    cadenceno = reference_tpf.hdu[1].data['CADENCENO']
    data.update({'CADENCENO': cadenceno})

    raw_counts = np.empty(np.shape(reference_tpf))
    raw_counts[:] = np.nan
    data.update({'RAW_CNTS': raw_counts})

    flux = np.empty(np.shape(reference_tpf))
    flux[:] = np.nan
    data.update({'FLUX': flux})

    flux_err = np.empty(np.shape(reference_tpf))
    flux_err[:] = np.nan
    data.update({'FLUX_ERR': flux_err})

    flux_bkg = np.empty(np.shape(reference_tpf))
    flux_bkg[:] = np.nan
    data.update({'FLUX_BKG': flux_bkg})

    flux_bkg_err = np.empty(np.shape(reference_tpf))
    flux_bkg_err[:] = np.nan
    data.update({'FLUX_BKG_ERR': flux_bkg_err})

    cosmic_rays = np.empty(np.shape(reference_tpf))
    cosmic_rays[:] = np.nan
    data.update({'COSMIC_RAYS': cosmic_rays})

    quality = np.empty(np.shape(reference_tpf)[0])
    quality[:] = np.nan
    data.update({'QUALITY': quality})

    poscorr1 = np.empty(np.shape(reference_tpf)[0])
    poscorr1[:] = np.nan
    data.update({'POS_CORR1': poscorr1})

    poscorr2 = np.empty(np.shape(reference_tpf)[0])
    poscorr2[:] = np.nan
    data.update({'POS_CORR2': poscorr2})

    # rb_level = np.empty(np.shape(reference_tpf))
    # rb_level[:] = np.nan
    # data.update({'RB_LEVEL': rb_level})
    # print(np.shape(rb_level))

    tpf_shape = (np.shape(flux)[1],np.shape(flux)[2])

    idx_col = reference_tpf.get_header(1)['1CRPX4']     # [pixel] reference pixel along image axis 1
    idx_row = reference_tpf.get_header(1)['2CRPX4']     # [pixel] reference pixel along image axis 2

    ## below need to fix such that it can boarder any edge of the inputed tpf
    position_dic = {'N': (0, + tpf_shape[1]), 'S': (0, -tpf_shape[1]), 'E': (+tpf_shape[0], 0), 'W': (-tpf_shape[0], 0),
                'NE': (+tpf_shape[0], +tpf_shape[1]), 'SE': (+tpf_shape[0], -tpf_shape[1]),
                'SW': (-tpf_shape[0], -tpf_shape[1]), 'NW': (-tpf_shape[0], +tpf_shape[1])}


    compass_pos = compass_position(position)

    pixel_col = original_tpfpos[0] + position_dic[compass_pos][0]  # CCD col
    pixel_row = original_tpfpos[1] + position_dic[compass_pos][1]   # CCD row

    radecs = find_radecs(pixel_col, pixel_row, reference_tpf)
    # print(radecs)

    pixel_shifts = {'pos': position,
                    'pixel_col': pixel_col, 'pixel_row': pixel_row,
                    'idx_col': idx_col, 'idx_row': idx_row,
                    'ra': radecs[0], 'dec': radecs[1]}

    tpf_temp = reference_tpf

    filename = f'null_tpf{compass_pos}.fits'
    CFF(data, tpf_temp.hdu, new_tpf_name = filename, tpf_shape = tpf_shape, verbose = False, pixel_shifts = pixel_shifts)

    return filename
