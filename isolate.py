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

import math
import sys
import string
import copy
import itertools
import operator
import collections
import networkx as nx

#import matplotlib
#matplotlib.rcParams['image.origin'] = 'lower'
#matplotlib.rcParams['figure.figsize'] = 6,6
#matplotlib.rcParams['figure.dpi'] = 80

import matplotlib.pyplot as plt
import numpy as np
import rht

#-----------------------------------------------------------------------------------------
# Initialization: Program Settings
#-----------------------------------------------------------------------------------------

SUFFIX = '_filaments.'

#-----------------------------------------------------------------------------------------
# Initialization: Object Definitions
#-----------------------------------------------------------------------------------------

class Cloud:
    def make_mask(self):
        mask_shape = (1+self.max_x-self.min_x, 1+self.max_y-self.min_y)
        self.mask = np.zeros(mask_shape, dtype=int)
        for point in self.points:
            self.mask[point[0]-self.min_x][point[1]-self.min_y] = 1 

    def to_ImageHDU(self):
        hdr = fits.Header()
        hdr['MIN_X'] = (self.min_x, 'Lower-left x-coordinate of mask in backprojection')
        hdr['MAX_X'] = (self.max_x, 'Upper-right x-coordinate of mask in backprojection')
        hdr['MIN_Y'] = (self.min_y, 'Lower-left y-coordinate of mask in backprojection')
        hdr['MAX_Y'] = (self.max_y, 'Upper-right y-coordinate of mask in backprojection')
        self.make_mask()
        hdr['AREA'] = (self.mask.size, 'Area covered by this mask')
        hdr['LITPIX'] = (len(self.points), 'Number of nonzero pixels in the mask')
        return fits.ImageHDU(data=self.mask, header=hdr)
    
    def __init__(self, list_of_points):
        #Expects a python list of two-integer tuples, corresponding the the x,y coordinate of the points original location in the backprojection 

        if isinstance(list_of_points, tuple):
            if len(list_of_points) > 2:
                list_of_points = list(list_of_points)
            elif len(list_of_points) == 2:
                list_of_points = [list_of_points]
        if isinstance(list_of_points, set):
            list_of_points = list(list_of_points)
        if isinstance(list_of_points, Cloud):
            self = list_of_points
            return 

        self.points = list(set(list_of_points))

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
                if int(point[0])!=point[0] or int(point[1])!=point[1]:
                    raise TypeError('Points must contain integer coordinates'+repr(point))
            return True
        #assert proper_formatting(list_of_points)

        '''
        self.min_x = min([point[0] for point in self.points])
        self.min_y = min([point[1] for point in self.points])
        self.max_x = max([point[0] for point in self.points])
        self.max_y = max([point[1] for point in self.points])
        '''
        self.points.sort(key=operator.itemgetter(0))    
        self.min_x, self.max_x = self.points[0][0], self.points[-1][0]
        self.points.sort(key=operator.itemgetter(1))    
        self.min_y, self.max_y = self.points[0][1], self.points[-1][1]
        

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

    '''
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
    '''
    hdu_list = fits.open(filaments_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    display = np.zeros((hdu_list[0].header['NAXIS1'], hdu_list[0].header['NAXIS2']))
    
    plt.ion()
    plt.imshow(display)
    plt.draw()
    
    for i, hdu in enumerate(hdu_list):
        if i==0:
            #First HDU is the Backprojection
            continue

        hdr = hdu.header
        min_x, min_y = hdr['MIN_X'], hdr['MIN_Y']
        #max_x, max_y = hdr['MAX_X'], hdr['MAX_Y']
        
        mask = np.copy(np.nonzero(hdu.data))
        
        mask[0] += min_x 
        mask[1] += min_y
        for coord in zip(mask[0], mask[1]):
            display[coord] = i  
        
        plt.cla()
        plt.clf()
        plt.imshow(display)
        plt.draw()
    
    print 'cla'
    plt.cla()
    print 'clf'
    plt.clf()
    print 'close'
    plt.close()
    print 'ioff'
    plt.ioff()
    print 'imshow'
    
    plt.imshow(display)
    print 'show'
    plt.show()
    print 'done'
    
def isolate_all(xyt_filename, BINS=6, DEBUG = True):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    #ntheta = hdu_list[0].header['NTHETA']
    wlen = hdu_list[0].header['WLEN']
    frac = hdu_list[0].header['FRAC']
    naxis1 = hdu_list[0].header['NAXIS1']
    naxis2 = hdu_list[0].header['NAXIS2']
    Hi = hdu_list[1].data['hi'] 
    Hj = hdu_list[1].data['hj'] 

    
    #Compute TheteRHT for all pixels given, then bin by theta
    def func(Hthet):
        return int((rht.theta_rht(Hthet, original=True)*BINS)//np.pi)
    C = np.asarray(map(func, hdu_list[1].data['hthets']))
    del func

    def rel_add(*tuples):
        return tuple(map(sum, zip(*tuples)))
    '''
    if not DEBUG:
        plt.ion()
    '''
    #Set Assignment
    unprocessed = list()    
    for bin in range(BINS):
        delimiter = np.nonzero(C == bin)[0]
        raw_points = zip(Hi[delimiter],Hj[delimiter])
        del delimiter
        problem_size = len(raw_points)
        message='Step '+str(bin+1)+'/'+str(BINS)+': (N='+str(problem_size)+')'

        point_dict = dict([x[::-1] for x in enumerate(raw_points)])
        set_dict = collections.defaultdict(list)

        extent = [(-1,-1), (-1, 0), (-1, 1), (0, -1)] #[(-1, 1), (-1,-1), (-1, 0), (0, -1), (-2, -2), (-2, -1), (-2, 0), (-2, 1), (-2, 2), (-1, -2), (-1, 2), (0,-2)] # #
        for i, coord in enumerate(raw_points):
            rht.update_progress(0.15*(i/problem_size), message=message)
            for rel_coord in extent:
                try:
                    j = point_dict[rel_add(coord, rel_coord)]
                    set_dict[point_dict[coord]].append(j)
                except Exception:
                    continue
        
        del extent 
        G = nx.from_dict_of_lists(set_dict) #Undirected graph made using set_dict as an adjacency list
        del set_dict 
        
        sources = range(problem_size)
        flags = np.ones((problem_size), dtype=int)
        while len(sources) > 0: 
            source = sources.pop()
            if not flags[source]:
                continue
            else:
                rht.update_progress(0.15+0.15*(1.0-len(sources)/problem_size), message=message)
                try:
                    for member in nx.descendants(G, source):
                        #sources.remove(member)
                        flags[member] = False
                        point_dict[raw_points[member]] = source
                except nx.NetworkXError:
                    #Assume we hit an isolated pixel and move on
                    pass

        #************************************************************************************************
        finished_map = np.negative(np.ones((naxis1, naxis2), dtype=np.int64))
        #for pt in raw_points:
            #finished_map[pt] = point_dict[pt]

        histogram = np.bincount(map(point_dict.get, raw_points))
        '''
        if DEBUG:
            plt.plot(histogram)
            plt.show()
        '''
        mask = np.nonzero(histogram >= int(frac*wlen))[0]
        del histogram

        first = True
        for set_id in mask:
            prog = 0.3+0.7*(1.0-len(point_dict)/problem_size)
            if 0.0 < prog < 1.0:
                rht.update_progress(prog, message=message)
                pass
            
            out_cloud = list()
            other_dict = dict()
            
            while len(point_dict) > 0:
                temp = point_dict.popitem()
                if set_id == temp[1]:
                    out_cloud.append(temp[0])
                    finished_map[temp[0]] = set_id
                elif first and (temp[1] not in mask):
                    #del temp
                    continue
                else:
                    other_dict[temp[0]]=temp[1]
            
            first = False
            point_dict = other_dict
            unprocessed.append(out_cloud)

        rht.update_progress(1.0, final_message='Finished joining '+str(problem_size)+' points! Time Elapsed:')
        '''
        try:
            plt.imshow(finished_map)#+1)
            if not DEBUG:
                plt.cla()
                plt.clf()
                plt.draw()
            else:
                plt.show()
        except Exception:
            print 'Failed to Display Theta Slice'
        '''
    '''    
    if not DEBUG:
        plt.cla()
        plt.clf()
        plt.close()
        plt.ioff()
    '''
    unprocessed.sort(key=len, reverse=True)
    if DEBUG:
        print map(len, unprocessed)

    #Convert lists of two-integer tuples into ImageHDUs
    output_hdulist = fits.HDUList(map(Cloud.to_ImageHDU, map(Cloud, unprocessed)))
    del unprocessed

    #Output HDUList to File
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    #output_hdulist.insert(0, fits.ImageHDU(data=backprojection_with_sets, header=fits.Header()))
    output_hdulist.insert(0, hdu_list[0].copy()) #TODO Introduces Errors in Reading FITS File 
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
        #show(isolate_all(xyt_filename))
        isolate_all(xyt_filename)

    #Cleanup and Exit
    #exit()




			