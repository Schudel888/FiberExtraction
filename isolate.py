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
#import string
#import tempfile 
#import shutil
#import time 
#import fnmatch

#import matplotlib.pyplot as plt
import numpy as np
import rht

#-----------------------------------------------------------------------------------------
# Initialization: Object Definitions
#-----------------------------------------------------------------------------------------

def combine_lists(lista, listb):
    pass

def is_adjacent(clouda, cloudb, radius=2):
    for pointa in clouda:
        for pointb in cloudb:
            if pointa is pointb:
                return True 
            if (math.hypot((pointa[0]-pointb[0]), (pointa[1]-pointb[1])) <= radius):
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
        return isinstance(other, Cloud) and is_adjacent(self.points, other.points) 
    '''
    def __init__(self, list_of_points):
        self.Points = set(filter(lambda p: isinstance(p, Point), list_of_points))
        self.min_x = min([point.x for point in self.points])
        self.min_y = min([point.y for point in self.points])    
        self.mask = np.zeros( 1+max([point.x for point in self.points])-self.min_x, 1+max([point.y for point in self.points])-self.min_y)
        for point in self.Points:
            self.mask[point.x-self.min_x][point.y-self.min_y] = point.b
        '''
        self.mask = np.zeros( 1+max([point.x for point in self.points])-self.min_x, 1+max([point.y for point in self.points])-self.min_y, *self.points[0].datashape)
        for point in self.points:
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
'''

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------

def backproj_from_xyt(xyt_filename):
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename

    print 'Accessing: '+xyt_filename+'... '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    #header = hdu_list[0].header
    backproj = hdu_list[0].data
    #Hi = hdu_list[1].data['hi'] 
    #Hj = hdu_list[1].data['hj'] 
    #Hthets = hdu_list[1].data['hthets']
    return backproj


def isolate_all(backproj, can_join=is_adjacent):
    assert backproj.ndim == 2 #2D Numpy Array 

    coords = np.nonzero(backproj) #Iterable of points
    raw_points = zip(coords[0],coords[1])
    
    problem_size = len(raw_points)
    print 'Assigning sets to '+str(problem_size)+' datapoints... '
    rht.update_progress(0.0)
    unprocessed = list()
    while len(raw_points) > 0:
        new_cloud = list()
        new_cloud.append(raw_points.pop())
        no_match = True
        for old_cloud in unprocessed:
            if can_join(old_cloud, new_cloud):
                old_cloud.extend(new_cloud) #old_cloud grows to include all points in new_cloud
                no_match = False
        if no_match:
            unprocessed.append(new_cloud)
        try:
            rht.update_progress(1.0-len(raw_points)/problem_size)
        except:
            pass
    try:
        rht.update_progress(1.0)
    except:
        pass

    print 'Refining '+str(len(unprocessed))+' sets... '
    rht.update_progress(0.0)
    processed = list()
    while len(unprocessed) > 0:
        new_cloud = unprocessed.pop()
        no_match = True
        for old_cloud in unprocessed:
            if can_join(old_cloud, new_cloud):
                old_cloud.extend(new_cloud) #old_cloud grows to include all points in new_cloud
                no_match = False
        if no_match:
            processed.append(new_cloud)
        try:
            rht.update_progress(len(processed)/(len(processed)+len(unprocessed)))
        except:
            pass
    try:
        rht.update_progress(1.0)
    except:
        pass

    print '# of Sets: '+str(len(processed))
    DEMO = True

    def make_a_Point((x,y)):
        return Point(x, y, backproj[x][y])

    def make_a_Cloud(list_of_points):
        return Cloud(map(make_a_Point, list_of_points))

    list_of_Clouds = map(make_a_Cloud, processed)

    if not DEMO:
        return list_of_Clouds

    else:
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
        isolate_all(backproj_from_xyt(xyt_filename))

    #Cleanup and Exit
    exit()




			