#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
from __future__ import division #Must be first line of code in the file
from astropy.io import fits
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

import scipy.ndimage.filters as filta
import math
#import os
import sys
import string
#import tempfile 
#import shutil
#import time 
#import fnmatch
import copy
import itertools
import operator

import matplotlib.pyplot as plt
import numpy as np
import rht

#-----------------------------------------------------------------------------------------
# Initialization: Program Settings
#-----------------------------------------------------------------------------------------

SUFFIX = '_filaments.'

DEBUG = True

#-----------------------------------------------------------------------------------------
# Initialization: Object Definitions
#-----------------------------------------------------------------------------------------

class Cloud:
    #See https://docs.python.org/2/tutorial/datastructures.html for List-like Sytax and Usage
    def make_mask(self):
        mask_shape = (1+max([point[0] for point in self.points])-self.min_x, 1+max([point[1] for point in self.points])-self.min_y)
        self.mask = np.zeros(mask_shape, dtype=int)
        for point in self.points:
            self.mask[point[0]-self.min_x][point[1]-self.min_y] = 1 

    def to_ImageHDU(self):
        hdr = fits.Header()
        hdr['MIN_X'] = (self.min_x, 'Lower-left x-coordinate of mask in backprojection')
        hdr['MIN_Y'] = (self.min_y, 'Lower-left y-coordinate of mask in backprojection')
        self.make_mask()
        hdr['AREA'] = (self.mask.size, 'Area covered by this mask')
        hdr['LITPIX'] = (np.count_nonzero(self.mask), 'Number of nonzero pixels in the mask')
        return fits.ImageHDU(data=self.mask, header=hdr)
    
    def __init__(self, list_of_points):
        #Expects a python list of two-integer tuples, corresponding the the x,y coordinate of the points original location in the backprojection 

        if isinstance(list_of_points, tuple):
            list_of_points = [list_of_points]
        if isinstance(list_of_points, set):
            list_of_points = list(list_of_points)
        if isinstance(list_of_points, Cloud):
            list_of_points = list_of_points.points
        
        def proper_formatting(given):
            if not isinstance(given, list):
                raise TypeError('Cloud cannot be constructed from the given datatype: '+str(type(given)))
            if len(given) is 0:
                raise ValueError('Cloud must contain at least one point')
            for point in given:
                if not isinstance(point, tuple):
                    raise TypeError('All points in a cloud must be coordinate tuples')
                if len(point) != 2:
                    raise ValueError('Points must be tuples of length 2')
                if not (isinstance(point[0], np.int64) and isinstance(point[1], np.int64)):
                    raise TypeError('Points must contain integer coordinates'+repr(point))
            return True
        #assert proper_formatting(list_of_points)

        self.points = list(set(list_of_points))
        self.min_x = min([point[0] for point in self.points])
        self.min_y = min([point[1] for point in self.points])    
        #self.make_mask()

#-----------------------------------------------------------------------------------------
# Rough Code
#-----------------------------------------------------------------------------------------

'''
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
plt.ion()
canvas = np.zeros_like(backproj)
plt.plot(canvas)
for cloud in list_of_Clouds:
    canvas[cloud.min_x:][cloud.min_y:] += cloud.mask
    plt.clf()
    plt.plot(canvas)
'''

'''
#STDOUT Progress Reporting
while not done:
    sys.stdout.write('\rSay {0} to {1}!'.format(zeroth, first))
    sys.stdout.flush()
    #DO this_tuff
    #done?
sys.stdout.write('\rDone with {0}!'.format(this_stuff))
sys.stdout.flush()
print ''
'''

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------

def show(filaments_filename):
    assert filaments_filename.endswith('.fits')
    assert '_xyt' in filaments_filename
    assert SUFFIX in filaments_filename
    print 'Accessing: '+filaments_filename+' '
    hdu_list = filter(lambda h: int(h.header['LITPIX'])>5, fits.open(filaments_filename, mode='readonly', memmap=True, save_backup=False, checksum=True)[1:]) #Allows for reading in very large files!
    
    NPlots = len(hdu_list)-1
    npages = int(math.ceil(NPlots/4)) 
    for page in range(npages):
        #plt.title('Filaments in: '+filaments_filename)
        for figure in range(4*page, min(4*(page+1), NPlots)):
            plt.subplot(2,2, (figure%4)+1)
            plt.spy(hdu_list[figure+1].data, origin='lower')
            #plt.contour(hdu_list[figure+1].data)
        plt.show()
        plt.cla()



def isolate_all(xyt_filename, BINS=6):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    header = hdu_list[0].header
    wlen = header['WLEN']
    ntheta = header['NTHETA']
    frac = header['FRAC']
    naxis1, naxis2 = hdu_list[0].header['NAXIS1'], hdu_list[0].header['NAXIS2']
    
    Hi = hdu_list[1].data['hi'] 
    Hj = hdu_list[1].data['hj'] 
    Hthets = hdu_list[1].data['hthets']
    C = np.zeros_like(Hi)
    for x in range(len(Hi)):
        C[x] = int((rht.theta_rht(Hthets[x], original=True)*BINS)//np.pi)
    #plt.plot(np.bincount(C))
    #plt.show()

    del Hthets
    
    #Set Assignment
    unprocessed = list()
    #plt.ion()
    for bin in range(BINS):
        delimiter = np.nonzero(C == bin)[0]
        raw_points = map(list, zip(Hi[delimiter],Hj[delimiter],np.zeros_like(delimiter)))
        raw_map = np.negative(np.ones((naxis1, naxis2), dtype=np.int64))
        problem_size = len(raw_points)
        message='Step '+str(bin+1)+'/'+str(BINS)+': (N='+str(problem_size)+')'

        for i, point in enumerate(raw_points):
            point[2] = i
            raw_map[point[0]][point[1]] = i

        extent = [(-1,-1), (-1, 0), (-1, 1), (0, -1), (-2, -2), (-2, -1), (-2, 0), (-2, 1), (-2, 2), (-1, -2), (-1, 2), (0,-2)] #[(-1,-1), (-1, 0), (-1, 1), (0, -1)] #
        #O(N)
        raw_points.sort(key=operator.itemgetter(0,1))
        for i, coord in enumerate(raw_points):
            rht.update_progress((i/problem_size), message=message)
            for relative_coord in extent:
                try:
                    coord[2] = raw_points[raw_map[coord[0]+relative_coord[0]][coord[1]+relative_coord[1]]][2]
                    break
                except Exception:
                    continue

        finished_map = np.zeros_like(raw_map)
        cloud = 1

        raw_points.sort(key=operator.itemgetter(2), reverse=True)
        representative = tuple(raw_points.pop())
        new_cloud = list()
        while len(raw_points)>0:
            next = raw_points.pop()
            if next[2] == representative[2]:
                new_cloud.append((next[0], next[1]))
            else:
                new_cloud.append((representative[0], representative[1]))
                if len(new_cloud) >= int(frac*wlen):
                    out_cloud = list()
                    while len(new_cloud) > 0:
                        out_point = new_cloud.pop()
                        finished_map[out_point] = cloud
                        out_cloud.append(out_point)
                    cloud += 1
                    unprocessed.append(out_cloud)
                '''
                else:
                    for out_point in new_cloud:
                        C[delimiter[raw_map[out_point]]] += 1
                '''
                representative = next
                #del new_cloud
                #new_cloud = list()
                #plt.imshow(finished_map)
                #plt.draw()

        rht.update_progress(1.0, final_message='Finished joining '+str(problem_size)+' points! Time Elapsed:')
        plt.imshow(finished_map) #.astype(np.float64)/cloud)
        #plt.draw()
        plt.show()
    plt.ioff()

    unprocessed.sort(key=len, reverse=True)
    if DEBUG:
        print map(len, unprocessed)

    #Convert lists of two-integer tuples into ImageHDUs
    list_of_Clouds = map(Cloud, unprocessed)
    list_of_HDUs = map(Cloud.to_ImageHDU, list_of_Clouds) #map(Cloud, unprocessed))
    #list_of_HDUs.sort(key=lambda hdu: hdu.header['AREA'], reverse=True)

    #Output HDUList to File
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    output_hdulist = fits.HDUList(list_of_HDUs)
    output_hdulist.insert(0, hdu_list[0].copy()) #fits.PrimaryHDU(data=backproj, header=fits.Header())) #header=header[6:-2])) #TODO Introduces Errors in Reading FITS File 
    output_hdulist.writeto(output_filename, output_verify='silentfix', clobber=True, checksum=True)

    print 'Results successfully output to '+output_filename
    return output_filename

        
#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    #Interpret Arguments
    parser = ArgumentParser(description="Run Fiber Isolation on 1+ _xyt FITS files", usage='%(prog)s [options] file(s)', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help="FITS file(s)")
    if len(sys.argv) == 1: # no arguments given, so add -h to get help msg
        sys.argv.append('-h')
    args = parser.parse_args()

    #Do Processing
    for xyt_filename in args.files: # loop over input files
        show(isolate_all(xyt_filename))

    #Cleanup and Exit
    exit()




			