#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
import os
import isolate

#-----------------------------------------------------------------------------------------
# Plotting Functions
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    available_integers = range(0,4)+range(10,25)+range(33,41)
    prefix = 'D:/SC_241.66_28.675.best_'
    insuffix = '_xyt_w75_s15_t70.fits'
    #outsuffix = '_xyt_w75_s15_t70_filaments.fits'

    plot_dir = 'D:/Plots/'
    DO_PLOTS = False
    FORCE = True

    available_cuts = [('_50', lambda h: h.header['B_MIN'] > 50.0)]

    #galfa = 'D:/SC_241.66_28.675.best.fits'

    def available_keys(integer):
        '''
        keys = ['']
        keys.append('B')
        keys.append('GALFA'+str(integer))
        '''
        return ['B']

    #TODO

    print ''
    print 'Starting plotter.py'
    print '{'

    #isolate.SILENT = True

    for i in available_integers[2:]:
        print 'Started:', str(i)

        xyt_filename = prefix + str(i) + insuffix
        assert isolate.is_xyt_filename(xyt_filename)

        filaments_filename = isolate.isolate_all(xyt_filename, force=FORCE)
        #isolate.filament_filename_from_xyt_filename(xyt_filename) #prefix + str(i) + outsuffix
        assert isolate.is_filament_filename(filaments_filename) 
        if DO_PLOTS:
            isolate.plot(filaments_filename)

        #TODO

        for _key in available_keys(i):

            isolate.update_key(filaments_filename, _key, force=FORCE)
            if DO_PLOTS:
                isolate.plot(filaments_filename, key=_key, out_name=plot_dir+filaments_filename[:-5]+'_'+_key+'.png')
            
                for cut_name, _cut in available_cuts:        
                    isolate.plot(filaments_filename, key=_key, out_name=plot_dir+filaments_filename[:-5]+'_'+_key+cut_name+'.png', cut=_cut)
                                        
                    #TODO

            #TODO


        print 'Finished:', str(i)
        #TODO

    print '} '
    print 'Finishing plotter.py '
    print ''

    #TODO
    exit()

