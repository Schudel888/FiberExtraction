#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
from __future__ import division #Must be first line of code in the file
import os
import sys
import numpy as np
import isolate
import config
import rht
from astropy.io import fits
import itertools

#-----------------------------------------------------------------------------------------
# Definitions
#-----------------------------------------------------------------------------------------

DIRECTORY = 'D:/' #"/Volumes/DataDavy/GALFA/lowell_fibers_project/"
ALL_CHANNELS = range(0,4)+range(10,25)+range(33,41)
GOOD_CHANNELS = [0]+range(10,25)+range(33,37)

prefix = DIRECTORY+'SC_241.66_28.675.best_'
insuffix = '_xyt_w75_s15_t70.fits'
plot_dir = DIRECTORY+'Plots/'

TABLENAME = prefix+'table.fits'

def quick_hdu_list(i):
    #READONLY
    assert isinstance(i, int)
    assert i in ALL_CHANNELS 
    assert i in GOOD_CHANNELS
    isolate.SILENT = True
    xyt_filename = prefix + str(i) + insuffix
    assert isolate.is_xyt_filename(xyt_filename)
    filaments_filename = isolate.filament_filename_from_xyt_filename(xyt_filename) 
    assert isolate.is_filament_filename(filaments_filename) 
    isolate.SILENT = False
    return config.default_open(filaments_filename, mode='readonly', memmap=True, checksum=False)

def generate_table(available_integers, output_tablename=TABLENAME):
    assert isinstance(output_tablename, str)
    if os.path.isfile(output_tablename):
        print 'A Table is already saved as', output_tablename
        if 'y' not in raw_input('Recompute anyway? ([no]/yes):  '):
            print 'Aborting tablemaker.py'
            return

    assert isinstance(available_integers, list)
    assert all([x in GOOD_CHANNELS for x in available_integers])
    available_integers.reverse()

    #Extracting a premade header
    my_header = config.default_open(DIRECTORY+'SC_241.66_28.675.best.fits', checksum=False)[0].header.copy(strip=True)
    try:
        del my_header['BUNIT'], my_header['BSCALE'], my_header['BZERO'], my_header['BLANK'], my_header['OBSFREQ'], my_header['O_BSCALE']
    except Exception:
        pass

    #Establishing which keywords will be included. GALFA is interpreted for each channel: (GALFA ==> GALFAx)
    float_keywords = ['B_MIN', 'B_MAX', 'DIAG', 'OFF_DIAG', 'GALFA_ONOFF_AVG', 'GALFA_ONOFF_MED', 'GALFA_MED', 'GALFA_AVG'] #RA, DEC
    int_keywords = ['LITPIX','MIN_X','MAX_X','MIN_Y','MAX_Y']

    output_hdulist = fits.HDUList(fits.PrimaryHDU(data=None, header=my_header))

    #Establishing the data shape of the output table #TODO Spectra
    typelist = [('channel', int), ('hdu', int), ('plotme', bool)] #, ('center_column', float), ('center_velocity', float), ('v_width', float), ('spectra', (float, 16))]
    typelist.extend(zip(float_keywords, [float]*len(float_keywords)))
    typelist.extend(zip(int_keywords, [int]*len(int_keywords)))

    #OPENS ALL HDUS BEFOREHAND TO ALLOCATE THE RIGHT SIZE TABLE
    hdu_lists = map(quick_hdu_list, available_integers) #TODO Dangerous Memory Usage
    PROBLEM_SIZE = sum(map(len, hdu_lists))-isolate.skip*len(hdu_lists)
    TABLE_arr = np.empty(PROBLEM_SIZE, dtype=typelist)
    TABLE = TABLE_arr.view(np.recarray)
    print ''

    #Loop through ALL hdu_lists and all hdus therein. Updates values.
    offset = 0
    while (offset < (PROBLEM_SIZE - 1)):
    #for i, hdu_list in itertools.izip(available_integers, hdu_lists):
        i, hdu_list = (available_integers.pop(), hdu_lists.pop())
        #hdu_list = quick_hdu_list(i)
        
        backproj_hdu = hdu_list.pop(0) #TODO: POP ALL SKIPPABLE HDUS
        backproj = backproj_hdu.data
        N = len(hdu_list) #Just the number of filaments

        if output_hdulist[0].data is None:
            my_data = np.zeros_like(backproj)
            my_data.fill(np.nan)
            output_hdulist[0].data = my_data
        assert output_hdulist[0].data is not None

        prog = float(offset)/float(PROBLEM_SIZE)
        for j, hdu in enumerate(hdu_list):
            rht.update_progress(float(prog+(j/PROBLEM_SIZE)), message=str(i)+': '+str(j))
            hdr = hdu.header

            TABLE.channel[offset+j] = i
            TABLE.hdu[offset+j] = j+isolate.skip
            TABLE.plotme[offset+j] = False
            #TABLE.spectra[offset+j] = calculate_spectra(hdu) #TODO spectra calculation

            for key in float_keywords+int_keywords:
                TABLE[key][offset+j] = hdr.get(key.replace('GALFA', 'GALFA'+str(i)), np.nan)
                
            inset = output_hdulist[0].data[hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1]
            inset[config.Cloud.nonzero_data_from_HDU(hdu, transpose=True)] = 0.0

        offset += N

        #TODO

        try:
            hdu_list.close()
            del hdu_list
        except Exception:
            pass

    rht.update_progress(1.0)

    #Put Table into HDUList
    print TABLE.shape
    output_hdulist.append(fits.BinTableHDU.from_columns(TABLE, header=fits.Header(dict(typelist)) )
    
    #Output HDUList to File
    output_hdulist.writeto(output_tablename, output_verify='silentfix', clobber=True, checksum=True)
    try:
        output_hdulist.flush()
    except Exception:
        pass
    try:
        output_hdulist.close()
    except Exception:
        pass

    #TODO
    return 

#-----------------------------------------------------------------------------------------
# Plotting, Indexing, and Viewing Library
#-----------------------------------------------------------------------------------------
class DataTable:

    def __init__(self, filename=TABLENAME):
        hdu_list = fits.open(filename, memmap=True, checksum=False, mode='update')
        self.recarray = hdu_list[1].data.view(np.recarray)
        self.dictionary = hdu_list[1].header

    def plot_key(self, key=None):
        if key is None or key not in self.dictionary:
            f = config.passive_constant(1)
        else:
            f = operator.

    def plot_cut(self, cut=config.passive_constant(True)):
        


    #TODO

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    '''
    if len(sys.argv) > 1: 
        available_integers = range(18,22) #[0]+range(10,25)+range(33,37)
    else:
    '''
    list_of_channels = range(18, 23)#range(12,25)
    filename=prefix+'table18_22.fits'
    #TODO

    print ''
    print 'Starting tablemaker.py'
    print '{'

    generate_table(list_of_channels, filename)

    table = DataTable(filename)

    print '} '
    print 'Finishing tablemaker.py '
    print ''

    #TODO
    exit()





