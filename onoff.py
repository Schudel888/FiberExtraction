#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
#import os
import numpy as np
import isolate
import config
#-----------------------------------------------------------------------------------------
# Plotting Functions
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    #all_available_integers = range(0,4)+range(10,25)+range(33,41)
    available_integers = [18,19,20,21,22]#[0]+range(10,25)+range(33,37)
    
    DIRECTORY = 'D:/'
    prefix = DIRECTORY+'SC_241.66_28.675.best_'
    insuffix = '_xyt_w75_s15_t70.fits'
    plot_dir = DIRECTORY+'Plots/'

    #outsuffix = '_xyt_w75_s15_t70_filaments.fits'
    #galfa = DIRECTORY+'SC_241.66_28.675.best.fits'

    available_cuts = [('_50', lambda h: h.header['B_MIN'] > 50.0)] #[(cut_name, _cut)]

    print ''
    print 'Starting onoff.py'
    print '{'

    for i in available_integers:
        print 'Started:', str(i)

        isolate.SILENT = True
        xyt_filename = prefix + str(i) + insuffix
        assert isolate.is_xyt_filename(xyt_filename)

        filaments_filename = isolate.isolate_all(xyt_filename, force=False) 
        #filaments_filename = isolate.filament_filename_from_xyt_filename(xyt_filename) #prefix + str(i) + outsuffix
        assert isolate.is_filament_filename(filaments_filename) 

        #Open the file once, but allow for multiple uses
        hdu_list = config.default_open(filaments_filename, mode='update')
        assert isolate.is_filament_hdu_list(hdu_list)

        _key = 'GALFA'+str(i)

        config.include(config.GALFAx(i))

        isolate.SILENT = False
        isolate.update_onoff_key(hdu_list, key=_key, force=True)

        NEWKEY = _key+'_ONOFF'
        config.include({NEWKEY: (
            None, 
            config.identity,
            np.log10,
            ['_AVG', '_MED', '_TOT']
        )})
        print 'Finished Updating', NEWKEY
        _out_name = plot_dir+filaments_filename[len(DIRECTORY):-5]+'_'+NEWKEY

        SHOW = False
        try:
            isolate.plot(hdu_list, key=NEWKEY, out_name=_out_name+'.png', show=SHOW)
            for cut_name, _cut in available_cuts:        
                isolate.plot(hdu_list, key=NEWKEY, out_name=_out_name+cut_name+'.png', show=SHOW, cut=_cut)
            print 'Finished Plotting', NEWKEY
        except Exception:
            print 'Warning: Failed to Plot', NEWKEY
            pass

        try:
            hdu_list.close()
        except Exception:
            pass

        try:
            config.exclude(NEWKEY)
            config.exclude(_key)
        except Exception:
            pass

        print 'Finished:', str(i)
        print ''
        #TODO

    print '} '
    print 'Finishing onoff.py '
    print ''

    #TODO
    exit()

