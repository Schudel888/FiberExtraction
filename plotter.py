#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
import isolate

inprefix = 'D:/SC_241.66_28.675.best_'
insuffix = '_xyt_w75_s15_t70_filaments.fits'

#outs = [('LITPIX', '_filaments_litpix.png'), 
outs= [('COLDENS','_average_coldens.png')]

ls = [18, 19, 20, 21]

for i in ls:
    for j,k in outs:
        outname = 'D:/'+str(i)+k
        filename = inprefix+str(i)+insuffix
        key = j
        isolate.plot(outname, filename, key)


