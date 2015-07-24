#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
import os
import isolate

import config
#-----------------------------------------------------------------------------------------
# Plotting Functions
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    all_available_integers = range(0,4)+range(10,25)+range(33,41)
    available_integers = [0]+range(10,25)+range(33,37)
    
    DIRECTORY = 'D:/'
    prefix = DIRECTORY+'SC_241.66_28.675.best_'
    insuffix = '_xyt_w75_s15_t70.fits'
    plot_dir = DIRECTORY+'Plots/'

    #outsuffix = '_xyt_w75_s15_t70_filaments.fits'
    #galfa = 'D:/SC_241.66_28.675.best.fits'

    #TODO DANGEROUS VARIABLES
    DO_PLOTS = False
    FORCE = False
    isolate.SILENT = True

    available_cuts = [('_50', lambda h: h.header['B_MIN'] > 50.0)] #[(cut_name, _cut)]

    def available_keys(integer):
        #return ['COLDENS']
        #return ['B']
        return ['GALFA'+str(integer)]

    #TODO


    print ''
    print 'Starting plotter.py'
    print '{'

    for i in available_integers:
        print 'Started:', str(i)

        xyt_filename = prefix + str(i) + insuffix
        assert isolate.is_xyt_filename(xyt_filename)

        filaments_filename = isolate.isolate_all(xyt_filename, force=FORCE) #isolate.filament_filename_from_xyt_filename(xyt_filename) #prefix + str(i) + outsuffix
        assert isolate.is_filament_filename(filaments_filename) 

        KEYS = available_keys(i)

        if len(KEYS) == 0 and DO_PLOTS:
            #Only open the file once in readonly mode
            isolate.plot(filaments_filename)

        elif len(KEYS) > 0:
            #Open the file once, but allow for multiple uses
            hdu_list = config.default_open(filaments_filename, mode='update')
            assert isolate.is_filament_hdu_list(hdu_list)
            if DO_PLOTS:
                isolate.plot(hdu_list)

            for _key in KEYS:

                isolate.update_key(hdu_list, _key, force=FORCE)

                if DO_PLOTS:
                    _out_name = plot_dir+filaments_filename[len(DIRECTORY):-5]+'_'+_key
                    print _out_name
                    isolate.plot(hdu_list, key=_key, out_name=_out_name+'.png')
                
                    for cut_name, _cut in available_cuts:        
                        isolate.plot(hdu_list, key=_key, out_name=_out_name+cut_name+'.png', cut=_cut)
                                            
                        #TODO

                #TODO

            try:
                hdu_list.close()
            except Exception:
                pass

        print 'Finished:', str(i)
        print ''
        #TODO

    print '} '
    print 'Finishing plotter.py '
    print ''

    #TODO
    exit()

