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

#import scipy.ndimage
import math
#import os
import sys
import string
#import tempfile 
#import shutil
#import time 
#import fnmatch

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

def is_within_radius(clouda, cloudb, radius=6):
    for pointa in clouda:
        for pointb in cloudb:
            if (math.hypot((pointa[0]-pointb[0]), (pointa[1]-pointb[1])) <= radius):
                return True
    return False


def is_adjacent(clouda, cloudb):
    for pointa in clouda:
        for pointb in cloudb:
            if (-1<=(pointa[0]-pointb[0])<=1 and -1<=(pointa[1]-pointb[1])<=1):
                return True
    return False

class Point:

    def __init__(self, x, y, b):
        self.x = x #x coordinate in backprojection
        self.y = y #y coordinate in backprojection
        self.b = b #integrated bacprojection power

class Cloud:
    '''
    def can_join(self, other):
        return isinstance(other, Cloud) and is_adjacent(self.Points, other.points) 
    '''
    def __init__(self, list_of_points):
        self.Points = set(filter(lambda p: isinstance(p, Point), list_of_points))
        self.min_x = min([point.x for point in self.Points])
        self.min_y = min([point.y for point in self.Points])    
        mask_shape = (1+max([point.x for point in self.Points])-self.min_x, 1+max([point.y for point in self.Points])-self.min_y)
        self.mask = np.zeros(mask_shape)
        for point in self.Points:
            self.mask[point.x-self.min_x][point.y-self.min_y] = float(point.b)
        '''
        self.mask = np.zeros( 1+max([point.x for point in self.Points])-self.min_x, 1+max([point.y for point in self.Points])-self.min_y, *self.Points[0].datashape)
        for point in self.Points:
            self.mask[point.x-self.min_x][point.y-self.min_y][:] = point.data[:]
        '''


#-----------------------------------------------------------------------------------------
# Rough Code
#-----------------------------------------------------------------------------------------
'''
def isolate_some(backproj, can_join=is_adjacent, join=combine_lists):
	masks = [] #list of cloud masks of the form (x, y, 2Dimage)
	
    processed = list() #Empty list of lists
    points = np.transpose(np.nonzero(backproj)) #Iterable of points
    while (size(points) > size(processed)):
        new_cloud = list(points.pop())
        for old_cloud in processed:
            if can_join(new_cloud, old_cloud):
                join()




	from collections import deque
	U = deque() #universal set of all possible clouds, each represented by a set of points
	
	for point in np.transpose(np.nonzero(backproj)):
		cloud = [point]  
		unionable = filter(lambda x: can_union(x,cloud), U)
		for y in unionable:
			U.remove(y)
		U.append(reduce(lambda x,y: x.extend(y)))
	
	old_count = 0
	while old_count < size(U):
		new_count = 0
		a = U.pop()
		while (new_count + old_count) < size(U):
			for i in range(old_count:size(U)-new_count):
				if can_union(a, U[i]):
					

		for (i, b) in enumerate(U[count:]):
			if can_union(a, b):
				a.extend(b)

def isolate_all(xyt_filename):
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    header = hdu_list[0].header
    data = hdu_list[1].data
    Hi = data['hi'] 
    Hj = data['hj'] 
    Hthets = data['hthets']



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
    
#Command Line Set Refinement
processed = list()
while len(unprocessed) > 0:
    sys.stdout.write('\rRefining {0} sets... '.format(len(unprocessed)))
    sys.stdout.flush()
    new_cloud = unprocessed.pop()
    #no_match = True
    matches = list()
    #for old_cloud in unprocessed:
    for i, old_cloud in enumerate(unprocessed):
        if can_join(old_cloud, new_cloud):
            old_cloud.extend(new_cloud) #old_cloud grows to include all points in new_cloud
            #no_match = False
            matches.append(i)
    #if no_match:
    if len(matches) is 0:    
        processed.append(new_cloud)
    else:
        #unprocessed = map(list, map(set, unprocessed))
        for i in matches:
            unprocessed[i] = list(set(unprocessed[i]))
sys.stdout.write('\rIdentified {0} clouds!'.format(len(processed)))
sys.stdout.flush()
print ''


'''

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------

def isolate_all(xyt_filename):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    #header = hdu_list[0].header
    backproj = hdu_list[0].data
    #Hi = hdu_list[1].data['hi'] 
    #Hj = hdu_list[1].data['hj'] 
    #Hthets = hdu_list[1].data['hthets']

    #Define Adjacency Function
    can_join=is_adjacent #is_within_radius

    #Set Assignment
    coords = np.nonzero(backproj) 
    raw_points = zip(coords[0],coords[1]) #List of two-integer tuples
    del coords 
    problem_size = len(raw_points)
    rht.update_progress(0.0, message='Unioning '+str(problem_size)+' points:')
    unprocessed = list()
    while len(raw_points) > 0:
        new_cloud = [raw_points.pop()] #List of two-integer tuples
        matches = list() #List of integers
        for i, old_cloud in enumerate(unprocessed):
            if can_join(old_cloud, new_cloud):
                matches.append(i)
        while len(matches) > 0:
            new_cloud.extend(unprocessed.pop(matches.pop()))
        unprocessed.append(new_cloud) #list(set(new_cloud))) #TODO DUPLICATES
        progress = math.pow(1.0-(len(raw_points)/problem_size), 2) #O(n**2)
        if 0.0 < progress < 1.0:
            rht.update_progress(progress=progress, message='Unioning '+str(len(raw_points))+' points:') 
    unprocessed.sort(key=len, reverse=True)
    rht.update_progress(1.0, final_message='Finished unioning '+str(problem_size)+'into '+str(len(unprocessed))+' sets! Time Elapsed:')
    
    if DEBUG:
        debug_data=map(len, unprocessed)
        print debug_data

    #Convert lists of two-integer tuples into ImageHDUs
    def make_a_Point((x,y)):
        return Point(x, y, backproj[x][y])

    def make_a_Cloud(list_of_points):
        return Cloud(map(make_a_Point, list_of_points))

    def make_an_ImageHDU(cloud):
        hdr = fits.Header()
        hdr['MIN_X'] = cloud.min_x
        hdr['MIN_Y'] = cloud.min_y
        return fits.ImageHDU(data=cloud.mask, header=hdr)
    list_of_HDUs = map(make_an_ImageHDU, map(make_a_Cloud, unprocessed))
    #list_of_HDUs.sort(key=lambda hdu: hdu.data.size, reverse=True)
    
    #Output HDUList to File
    output_hdulist = fits.HDUList(list_of_HDUs)
    output_hdulist.insert(0, hdu_list[0].copy()) #fits.PrimaryHDU(data=backproj, header=header)) #TODO Introduces Errors in Reading FITS File 
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    output_hdulist.writeto(output_filename, output_verify='silentfix', clobber=True, checksum=True)
    print 'Results successfully output to '+output_filename
    #return output_filename
        
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
        isolate_all(xyt_filename)

    #Cleanup and Exit
    exit()




			