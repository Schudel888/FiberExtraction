#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
#import os
#import config
import isolate

#-----------------------------------------------------------------------------------------
# Plotting Functions
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":


    prefix = 'D:/SC_241.66_28.675.best_'
    insuffix = '_xyt_w75_s15_t70.fits'
    outsuffix = '_xyt_w75_s15_t70_filaments.fits'
    
    #galfa = 'D:/SC_241.66_28.675.best.fits'

    for i in [17, 18, 19, 20, 21, 22]:
        for key in ['']: #'GALFA'+str(i)]:
            filename = prefix + str(i) + outsuffix
            isolate.update_key(filename, key, force=False)
            isolate.plot(filename, key=key, out_name=filename[:-5]+'_'+key+'.png')
            isolate.plot(filename, key=key, out_name=filename[:-5]+'_'+key+'_50.png', cut=(lambda h: h.header['B_MIN'] > 50.0))

            '''
            try:

            except Exception:
                
            '''
            #os.system('python isolate.py '+inprefix+str(i)+outsuffix)
            
            #os.system('python isolate.py '+inprefix+str(i)+insuffix)

