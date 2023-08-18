import numpy as np
from astropy.io import fits


## Written so that it can take any shaped tpf (I think)!

def make_BinTableHDU(data, existing_fitsfile, tpf_shape = (15,15)):
    """ data is in the form of a list of tpfs for difference cadences. existing fits file is a fitsfile from the bottom left TPF of your stitched TPF"""

    cols = []

    temp = np.zeros(tpf_shape)

    format_shape = len(temp.flatten())

    print('Bin', tpf_shape)

    tpf_shape = (tpf_shape[1], tpf_shape[0])   # FITS ARE STUPID AND HAVE THE DIMENSIONS FLIPPED!

    # print(data['flux'])

    cols.append(fits.Column(name='TIME', format='D', disp='D14.7', array=data['TIME']))
    cols.append(fits.Column(name='TIMECORR', format='E', disp='E14.7', array=data['TIMECORR']))
    cols.append(fits.Column(name='CADENCENO', format='J', disp='I10', array=data['CADENCENO']))
    cols.append(fits.Column(name='RAW_CNTS', format=f'{format_shape}J', disp='I8', array=data['RAW_CNTS'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='FLUX', format=f'{format_shape}E', disp='E14.7', array=data['FLUX'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='FLUX_ERR', format=f'{format_shape}E', disp='E14.7', array=data['FLUX_ERR'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='FLUX_BKG', format=f'{format_shape}E', disp='E14.7', array=data['FLUX_BKG'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='FLUX_BKG_ERR', format=f'{format_shape}E', disp='E14.7', array=data['FLUX_BKG_ERR'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='COSMIC_RAYS', format=f'{format_shape}E', disp='E14.7', array=data['COSMIC_RAYS'], dim=str(tpf_shape)))
    cols.append(fits.Column(name='QUALITY', format='J', disp='E14.7', array=data['QUALITY']))
    cols.append(fits.Column(name='POS_CORR1', format='E', disp='E14.7', array=data['POS_CORR1']))
    cols.append(fits.Column(name='POS_CORR2', format='E', disp='E14.7', array=data['POS_CORR2']))
    cols.append(fits.Column(name='RB_LEVEL', format=f'{format_shape}E', disp='E14.7', array=data['RB_LEVEL'], dim=str(tpf_shape)))

    coldefs = fits.ColDefs(cols)
    hdu = fits.BinTableHDU.from_columns(coldefs)

    naxis1_orig = hdu.header['NAXIS1']
    naxis2_orig = hdu.header['NAXIS2']

    ## use the header (and wcs) from an existing MAST TPF. Use bottom left TPF
    hdu.header = existing_fitsfile[1].header

    ## change the dimensions to match the new tpf shape
    dim_index = [4,5, 6, 7, 8, 9, 13]
    for index in dim_index:
        hdu.header[f'TDIM{index}'] = str(tpf_shape)

    hdu.header['NAXIS1'] = naxis1_orig
    hdu.header['NAXIS2'] = naxis2_orig

    return hdu, coldefs

def make_PrimaryHDU(existing_fitsfile):
    """ make the primary HDU. Copied straight from an exisiting TPF """

    hdu = fits.PrimaryHDU()
    hdu.header = existing_fitsfile[0].header

    return hdu

def make_imageHDU(existing_fitsfile, tpf_shape = (15,15)):
    """ make the primary HDU. Copied straight from an exisiting TPF. Dimensions are changed to match the new stitched TPF. This is used for the aperture mask in LightKurve. Assume that the pipeline mask is the entirety of the stitched TPF """

    temp = np.zeros(tpf_shape)
    temp[:] = 3

    format_shape = len(temp.flatten())

    hdu = fits.ImageHDU(name = 'APERTURE')

    hdu.data = temp

    hdu.header = existing_fitsfile[2].header

    # print('Image', tpf_shape)

    ## Change certain keys in header to match the new dimensions
    hdu.header['NPIXSAP'] = format_shape
    hdu.header['NAXIS1'] = tpf_shape[0]
    hdu.header['NAXIS2'] = tpf_shape[1]    # this could be switched with the above line?

    return hdu


def main(data, existing_fitsfile, new_tpf_name = 'test.fits', tpf_shape = (15,15), verbose = True):

    hdu_BinTable, coldefs = make_BinTableHDU(data, existing_fitsfile, tpf_shape = tpf_shape)

    hdu_Primary = make_PrimaryHDU(existing_fitsfile)

    hdu_Image = make_imageHDU(existing_fitsfile, tpf_shape = tpf_shape)

    new_hdul = fits.HDUList()

    new_hdul.append(hdu_Primary)  #This has to come first
    new_hdul.append(hdu_BinTable)
    new_hdul.append(hdu_Image)

    ## to print the header info for hdu list
    if verbose:
        new_hdul.info()

    new_hdul.writeto(new_tpf_name, overwrite=True)   ## be careful, will overwrite tpf with the same name
    print(f'New fits file has been written called {new_tpf_name}')

    return new_hdul
