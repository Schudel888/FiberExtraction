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

def merge_colinear_clouds(clouda, cloudb, are_colinear, radius=3.0):
    complements_clouda = copy.deepcopy(cloudb)
    complements_cloudb = list()
    any_in_range = False
    none_perpendicular = True
    for pointa in clouda:
        is_complementary = list()
        for pointb in cloudb:
            if (math.hypot((pointa[0]-pointb[0]), (pointa[1]-pointb[1])) <= radius):
                any_in_range = True
                if pointa != pointb and are_colinear(pointa, pointb):
                    is_complementary.append(True)
                else:
                    none_perpendicular = False
                    is_complementary.append(False)
                    try:
                        complements_clouda.remove(pointb)
                    except ValueError:
                        pass 
        if any(is_complementary) and all(is_complementary):
            complements_cloudb.append(pointa)
    
    if (any_in_range and none_perpendicular):
        return True, 
    elif (not any_in_range):
        return False
    else:
        clouda.extend(complements_clouda)
        cloudb.extend(complements_cloudb)
        return False
    


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
            #plt.spy(hdu_list[figure+1].data, origin='lower')
            plt.contour(hdu_list[figure+1].data)
        plt.show()
        plt.cla()



def isolate_all(xyt_filename):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    header = hdu_list[0].header
    wlen = header['WLEN']
    #ntheta = header['NTHETA']
    backproj = hdu_list[0].data
    #backproj = hdu_list[0].data[wlen:min(header['NAXIS1'], 500), wlen:min(header['NAXIS2'], 500)]
    
    Hi = hdu_list[1].data['hi'] 
    Hj = hdu_list[1].data['hj'] 
    Hthets = hdu_list[1].data['hthets']
    theta_rhts = np.zeros_like(backproj)
    U = np.zeros(len(Hi))
    V = np.zeros(len(Hi))
    C = np.zeros(len(Hi))
    for x in range(len(Hi)):
        theta_rhts[Hi[x], Hj[x]], U[x], V[x] = rht.theta_rht(Hthets[x], original=True, uv=True)
        C[x] = theta_rhts[Hi[x], Hj[x]]
    
    if DEBUG:
        plt.quiver(Hi, Hj, U, V, C, cmap=plt.get_cmap('hsv'))
        plt.show()
    else:
        del Hthets, Hi, Hj, U, V, C

    def are_colinear(pta, ptb, threshold=0.75):
        return (threshold <= math.cos(theta_rhts[pta]-theta_rhts[ptb])**2)
    def are_nearby(pta, ptb, radius=3.0):
        return (math.hypot((pta[0]-ptb[0]), (pta[1]-ptb[1])) <= radius)
    def all_nearby_are_colinear(clouda, cloudb):
        any_in_range = False
        for pointa in clouda:
            for pointb in cloudb:
                if are_nearby(pointa, pointb):
                    any_in_range = True
                    if pointa == pointb:
                        continue
                    elif not are_colinear(pointa, pointb):
                        return False
        return any_in_range

    #Set Assignment
    coords = np.nonzero(backproj) 
    raw_points = zip(coords[0],coords[1]) #List of two-integer tuples
    del coords 
    problem_size = len(raw_points)
    rht.update_progress(0.0, message='Unioning '+str(problem_size)+' points:')
    
    unprocessed = list()
    while len(raw_points) > 0:
        new_point = raw_points.pop()
        matches = list()
        new_clouds = list()
        for i, old_cloud in enumerate(unprocessed):
            for old_point in old_cloud:
                if are_nearby(new_point, old_point) and are_colinear(new_point, old_point):
                    old_cloud.append(new_point) 
                    matches.append(i)
                    break 
        N = len(matches)
        if N == 0:
            new_clouds.append([new_point])
        elif N == 1:
            #unprocessed[matches[0]].append(new_point)
            pass
        else:
            temp_clouds = list()
            while len(matches) > 0:
                temp_clouds.append(unprocessed.pop(matches.pop()))
            del matches

            def f(temp):
                X = len(temp)
                if X > 2:
                    if all([all_nearby_are_colinear(*combo) for combo in itertools.combinations(temp, 2)]):
                        '''
                        new_cloud = [temp.pop()]
                        while len(temp) > 0:
                            new_cloud.extend(temp.pop())
                            new_cloud.remove(new_point)
                        new_clouds.append(new_cloud)
                        '''
                        new_cloud = set()
                        while len(temp) > 0:
                            new_cloud.update(temp.pop())
                        new_clouds.append(list(new_cloud))
                    else:
                        map(f, map(list, itertools.combinations(temp, X-1)))
                elif X==2:
                    a = temp.pop()
                    b = temp.pop()
                    if all_nearby_are_colinear(a,b):
                        a.extend(b)
                        a.remove(new_point)
                        new_clouds.append(a)
                    else:
                        if a not in new_clouds:
                            new_clouds.append(a)
                        if b not in new_clouds:
                            new_clouds.append(a)

            f(temp_clouds)



        '''
        new_cloud = [raw_points.pop()]
        matches = list()
        for i, old_cloud in enumerate(unprocessed):
            match = list()
            for old_point in old_cloud:
                for new_point in new_cloud:
                    if are_nearby(new_point, old_point) and are_colinear(new_point, old_point):
                        match.append(True)
                        new_cloud.append(old_point)
                    else:
                        match.append(False)
            if all(match):
                matches.append(i)

        #List of two-integer tuples

        matches = list()
        for i, other_point in enumerate(raw_points):
            if are_nearby(new_point, other_point) and are_colinear(new_point, other_point):
                matches.append(i)

        temp_points = list()
        while len(matches) > 0:
            temp_points.append(raw_points.pop(matches.pop()))

        new_cloud = 

        assert len(matches) == 0
        for i, old_cloud in enumerate(unprocessed):
            if all_nearby_are_colinear(new_cloud, old_cloud):
                matches.append(i)
        while len(matches) > 0:
            new_cloud.extend(unprocessed.pop(matches.pop()))
        unprocessed.append(new_cloud)
        '''


        unprocessed.extend(new_clouds)
        progress = math.pow(1.0-(len(raw_points)/problem_size), 2) #O(n**2)
        if 0.0 < progress < 1.0:
            rht.update_progress(progress=progress, message='Unioning '+str(len(raw_points))+' points:') 
    
    rht.update_progress(1.0, final_message='Finished unioning '+str(problem_size)+' points into '+str(len(unprocessed))+' sets! Time Elapsed:')
    unprocessed.sort(key=len, reverse=True)
    if DEBUG:
        print map(len, unprocessed)
        

    #Convert lists of two-integer tuples into ImageHDUs
    list_of_Clouds = map(Cloud, unprocessed)
    list_of_HDUs = map(Cloud.to_ImageHDU, list_of_Clouds) #map(Cloud, unprocessed))
    #list_of_HDUs.sort(key=lambda hdu: hdu.header['AREA'], reverse=True)

    #Output HDUList to File
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    if DEBUG:
        output_hdulist = fits.HDUList(list_of_HDUs)
        output_hdulist.insert(0, fits.PrimaryHDU(data=backproj, header=fits.Header())) #header=header[6:-2])) #hdu_list[0].copy()) #TODO Introduces Errors in Reading FITS File 
        output_hdulist.writeto(output_filename, output_verify='silentfix', clobber=True, checksum=True)
    else:
        shdu = fits.StreamingHDU(output_filename, header)

    print 'Results successfully output to '+output_filename
    return output_filename

    '''
    #Unfinished Methods.......................................................................
    def xyt(point):
        (x,y) = point
        #This would go a lot faster if the RHT sorted output by x coordinate
        for i in range(ntheta):
            if (Hi[i] == x) and (Hj[i] == y):
                return Hthets[i]
        return np.zeros(ntheta)

    SIGMA = int(ntheta/16.0) #TODO What fraction of pi should be covered in 1StandardDeviation?
    xyt_double_helix = rht.all_thetas(wlen, np.linspace(0.0, np.pi, ntheta), True)
    filtered_xyt_double_helix = filta.gaussian_filter1d(xyt_double_helix, SIGMA, axis=2, mode='wrap')
    def weights(displacement):
        try:
            arr = filtered_xyt_double_helix[displacement[0]+wlen//2][displacement[1]+wlen//2]
        except Exception:
            arr = np.zeros(ntheta)
            temp_theta = math.atan2(math.abs(displacement[1]), displacement[0])
            temp_index = int(math.floor((ntheta-1)*temp_theta/math.pi))
            arr[temp_index] = 1 #math.sqrt(ntheta/math.pi) #TODO Normalization???
            filta.gaussian_filter1d(arr, SIGMA, mode='wrap', output=arr)
        finally:
            temp_sum = np.sum(np.multiply(arr,arr))/2.0
            return np.divide(arr, temp_sum) #TODO Scaling Area to 1?

    def synergy_between(pointa, pointb):
        if pointa is pointb:
            return 0.0
        displacement = (pointa[0]-pointb[0], pointa[1]-pointb[1])
        r2 = (displacement[0]**2 + displacement[1]**2)/(wlen**2)
        athets = 


    print 'Beginning Fiber Refinement...'
    '''

        
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




			