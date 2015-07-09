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
import time
import string
import copy
import itertools
import operator
import collections
import networkx as nx

import numpy as np
import rht

#-----------------------------------------------------------------------------------------
# Initialization: Program Settings
#-----------------------------------------------------------------------------------------

SUFFIX = '_filaments.'

#-----------------------------------------------------------------------------------------
# Initialization: Object Definitions
#-----------------------------------------------------------------------------------------

class Progress:
    # Create progress meter that looks like: 
    # message + ' ' + '[' + '#'*p + ' '*(length-p) + ']' + time_message
    stream = sys.stdout
    all_progress = list()
    TEXTWIDTH = 79
    length = int(0.55 * TEXTWIDTH)

    def update(self, progress=None):
        return progress
        if progress is None:
            if self._incrementing:
                self._step += 1
            else:
                self._step -= 1
        elif isinstance(progress, int) and (self._problem_size >= progress >= 0):
            self._step = progress
        elif isinstance(progress, float) and (1.0 >= progress >= 0.0):
            self._step = int(self._problem_size*progress)
        else:
            return #TODO

        now = time.time()
        if (self._stop_time is None) or (self._step == 0) or (self._step == self._problem_size) or ((now - self._last_update) > 1.0):
            if self._incrementing:
                prog = self._step/self._problem_size
            else:
                prog = 1.0 - (self._step/self._problem_size)
            
            self._stop_time = self._start_time + (now - self._start_time)/prog
        
            sec_remaining = int(self._stop_time - now)
            if sec_remaining >= 3600:
                time_message = ' < ' + str(sec_remaining//3600  +1) + 'hrs'
            elif sec_remaining >= 60:
                time_message = ' < ' + str(sec_remaining//60  +1) + 'min'
            else:
                time_message = ' < ' + str(sec_remaining +1) + 'sec'
            messlen = Progress.TEXTWIDTH-(Progress.length+3)-len(time_message)
            p = int(Progress.length*prog)
            if 0 <= p < Progress.length:
                Progress.stream.write('\r{2} [{0}{1}]{3}'.format('#'*p, ' '*(Progress.length-p), string.ljust(self._message, messlen)[:messlen], time_message))
                Progress.stream.flush()
            else:
                final_offset = Progress.TEXTWIDTH-len(time_message)
                final_message = string.ljust('Finished:', final_offset)[:final_offset]
                Progress.stream.write('\r{0}{1}'.format(final_message, time_message))
                Progress.stream.flush()
                try:
                    Progress.all_progress.remove(self)
                    if len(Progress.all_progress) == 0:
                        print ''
                except Exception:
                    pass
                
    def __init__(self, problem_size, message='Progress:', incrementing=True):
        assert isinstance(problem_size, int)
        assert problem_size > 0
        assert isinstance(message, str)
        assert isinstance(incrementing, bool)

        self._start_time = time.time()
        self._stop_time = None
        self._message = message
        self._incrementing = incrementing
        self._problem_size = problem_size
        if self._incrementing:
            self._step = 0
        else:
            self._step = self._problem_size
        self._last_update = self._start_time
        Progress.all_progress.append(self)


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
        hdr['DIAG'] = (math.hypot(*self.mask.shape), 'Diagonal size of mask')
        return fits.ImageHDU(data=self.mask, header=hdr)
    
    def __init__(self, list_of_points):
        #Expects a python list of two-integer tuples, corresponding the the x,y coordinate of the points original location in the backprojection 
        if isinstance(list_of_points, Cloud):
            self = list_of_points
            return 
        
        elif isinstance(list_of_points, tuple):
            if len(list_of_points) > 2:
                self.points = list(set(list_of_points))
            elif len(list_of_points) == 2:
                self.points = [list_of_points]
        
        elif isinstance(list_of_points, set):
            self.points = list(list_of_points)

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
        assert proper_formatting(self.points)

        self.points.sort(key=operator.itemgetter(0))    
        self.min_x, self.max_x = self.points[0][0], self.points[-1][0]
        self.points.sort(key=operator.itemgetter(1))    
        self.min_y, self.max_y = self.points[0][1], self.points[-1][1]
        

#-----------------------------------------------------------------------------------------
# Rough Code
#-----------------------------------------------------------------------------------------

'''
#STDOUT Progress Reporting
#Eqivalent to map(do, things)
constant = 100/len(things)
for i, thing in things:
    sys.stdout.write('\r'+str(i*constant)+'%... ')
    sys.stdout.flush()
    do(thing)
sys.stdout.write('\rDone!')
sys.stdout.flush()
print ''
'''

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------
skip = 1 #Number of HDUs that do not correspond to filaments in output files

source_data = dict()
source_data['COLDENS'] = np.asarray(fits.open('D:/LAB_corrected_coldens.fits', mode='readonly', memmap=False, save_backup=False, checksum=True)[0].data)

def filament_average(filaments_filename, key, external_source=None):
    assert filaments_filename.endswith('.fits')
    assert '_xyt' in filaments_filename
    assert SUFFIX in filaments_filename
    assert isinstance(key, str) #TODO for now?

    print 'Accessing: '+filaments_filename+' '
    hdu_list = fits.open(filaments_filename, mode='update', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!

    if key in hdu_list[skip].header:
        print key+' keyword found in filament header...'
        if raw_input('Overwrite? [n]/y') != 'y':
            return hdu_list, key

    if external_source is not None:
        correlation_data = external_source
    elif key in source_data:
        correlation_data = source_data[key]
    else:
        raise KeyError('No source corresponding to the key: '+key)
    
    for i, hdu in enumerate(hdu_list[skip:]):
        hdr = hdu.header
        hdr[key] = np.nanmean(correlation_data[hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1][np.nonzero(hdu.data.T)])

    hdu_list.flush()
    #hdu_list.close()
    return hdu_list, key

def show(filaments_filename, key=None):
    import matplotlib
    #matplotlib.rcParams['backend'] = 'TkAgg'
    matplotlib.rcParams['image.origin'] = 'lower'
    #matplotlib.rcParams['figure.figsize'] = 6,6
    #matplotlib.rcParams['figure.dpi'] = 180
    from matplotlib import pyplot as plt

    if isinstance(filaments_filename, str):
        assert filaments_filename.endswith('.fits')
        assert '_xyt' in filaments_filename
        assert SUFFIX in filaments_filename
        print 'Accessing: '+filaments_filename+' '

        hdu_list = fits.open(filaments_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    
    elif isinstance(filaments_filename, fits.HDUList):
        hdu_list = filaments_filename

    else:
        print 'Unknown input in show'
        return 

    display = np.zeros((hdu_list[0].header['NAXIS1'], hdu_list[0].header['NAXIS2']))

    if key is None:
        for i, hdu in enumerate(hdu_list[skip:]):
            hdr = hdu.header            
            display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = i
    elif isinstance(key, str):
        for hdu in hdu_list[skip:]:
            hdr = hdu.header            
            display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = hdr[key]
    else:
        print 'Unable to show data using the given key: '+str(key)
        return

    plt.imshow(display.T)
    plt.show()

def plot(out_name, filaments_filename, key=None):
    import matplotlib
    #matplotlib.rcParams['backend'] = 'TkAgg'
    matplotlib.rcParams['image.origin'] = 'lower'
    #matplotlib.rcParams['figure.figsize'] = 6,6
    #matplotlib.rcParams['figure.dpi'] = 200
    from matplotlib import pyplot as plt

    if isinstance(filaments_filename, str):
        assert filaments_filename.endswith('.fits')
        assert '_xyt' in filaments_filename
        assert SUFFIX in filaments_filename
        print 'Accessing: '+filaments_filename+' '

        hdu_list = fits.open(filaments_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    
    elif isinstance(filaments_filename, fits.HDUList):
        hdu_list = filaments_filename

    else:
        print 'Unknown input in plot'
        return 

    display = np.zeros((hdu_list[0].header['NAXIS1'], hdu_list[0].header['NAXIS2']))

    if key is None:
        for i, hdu in enumerate(hdu_list[skip:]):
            hdr = hdu.header            
            display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = i
    elif isinstance(key, str):
        for hdu in hdu_list[skip:]:
            hdr = hdu.header            
            display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = hdr[key]
    else:
        print 'Unable to plot data using the given key: '+str(key)
        return

    plt.imshow(display.T)
    plt.colorbar(aspect=3)
    plt.savefig(out_name, dpi=500, format='png')
    plt.clf()

    if key is not None and isinstance(filaments_filename, str):
        lol = map(lambda h: h.header[key], hdu_list[skip:])
        plt.hist(lol, bins=100, histtype='stepfilled', log=(key == 'COLDENS'), label=key+' from '+filaments_filename)
        hist_name = string.rstrip(out_name, '.png')+'_hist'+'.png'
        plt.savefig(hist_name, dpi=500, format='png')
        plt.clf() 

def isolate_all(xyt_filename, BINS=6, DEBUG = True):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert '_xyt' in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = fits.open(xyt_filename, mode='readonly', memmap=True, save_backup=False, checksum=True) #Allows for reading in very large files!
    ntheta = hdu_list[0].header['NTHETA']
    wlen = hdu_list[0].header['WLEN']
    frac = hdu_list[0].header['FRAC']
    naxis1 = hdu_list[0].header['NAXIS1']
    naxis2 = hdu_list[0].header['NAXIS2']
    original = hdu_list[0].header['ORIGINAL']
    Hi = hdu_list[1].data['hi'] 
    Hj = hdu_list[1].data['hj'] 

    #Compute TheteRHT for all pixels given, then bin by theta
    C = np.multiply(np.asarray(map(rht.theta_rht,hdu_list[1].data['hthets'])), BINS/np.pi).astype(np.int_)
    
    def rel_add((a,b), (c,d)):
        return a+c,b+d
    '''
    if not DEBUG:
        plt.ion()
    '''
    #Set Assignment
    #unprocessed = list()
    list_of_HDUs = list() 
    search_pattern = [(-1,-1), (-1, 0), (-1, 1), (0, -1)] #[(-1, 1), (-1,-1), (-1, 0), (0, -1), (-2, -2), (-2, -1), (-2, 0), (-2, 1), (-2, 2), (-1, -2), (-1, 2), (0,-2)]
    for bin in range(BINS):
        delimiter = np.nonzero(C == bin)[0]
        raw_points = zip(Hi[delimiter],Hj[delimiter])
        del delimiter
        problem_size = len(raw_points)
        message='Step '+str(bin+1)+'/'+str(BINS)+': (N='+str(problem_size)+')'
        progress_bar = Progress(problem_size, message=message, incrementing=True)

        point_dict = dict([x[::-1] for x in enumerate(raw_points)])
        set_dict = collections.defaultdict(list)

        for coord in raw_points:
            #rht.update_progress(0.3*(i/problem_size), message=message)
            progress_bar.update()
            for rel_coord in search_pattern:
                try:
                    j = point_dict[rel_add(coord, rel_coord)]
                    set_dict[point_dict[coord]].append(j)
                except Exception:
                    continue
        
        G = nx.from_dict_of_lists(set_dict) #Undirected graph made using set_dict as an adjacency list
        del set_dict 
        
        progress_bar = Progress(problem_size, message=message, incrementing=False)
        sources = range(problem_size)
        flags = np.ones((problem_size), dtype=np.int_)
        while len(sources) > 0: 
            source = sources.pop()
            if not flags[source]:
                continue
            else:
                #rht.update_progress(0.3+0.3*(1.0-len(sources)/problem_size), message=message)
                progress_bar.update(len(sources))
                try:
                    for member in nx.descendants(G, source):
                        flags[member] = False
                        point_dict[raw_points[member]] = source
                        #TODO Remove members from G if that would speed up subsequent calls?
                except nx.NetworkXError:
                    #Assume we hit an isolated pixel (never made it into G) and move on
                    pass
        del sources, flags, G

        #************************************************************************************************
        #finished_map = np.negative(np.ones((naxis1, naxis2), dtype=np.int64))
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

        progress_bar = Progress(problem_size, message=message, incrementing=False)
        mask_dict = dict([x[::-1] for x in enumerate(mask)])
        out_clouds = collections.defaultdict(list)

        while len(point_dict) > 0:
            temp = point_dict.popitem()
            try:
                #Keying into mask_dict is the only operation that ought to throw an exception 
                out_clouds[mask_dict[temp[1]]].append(temp[0])
                progress_bar.update(len(point_dict))
                #rht.update_progress(0.6+0.399*(1.0-len(point_dict)/problem_size), message=message)
            except Exception:
                continue

        while len(out_clouds) > 0:
            #unprocessed.append(out_clouds.popitem()[1])
            list_of_HDUs.append(Cloud(out_clouds.popitem()[1]).to_ImageHDU())

        #rht.update_progress(1.0, final_message='Finished joining '+str(problem_size)+' points! Time Elapsed:')
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
    #unprocessed.sort(key=len, reverse=True)
    #if DEBUG:
    #print map(len, unprocessed)

    #Convert lists of two-integer tuples into ImageHDUs
    list_of_HDUs.sort(key=lambda h: h.header['DIAG'], reverse=True)
    output_hdulist = fits.HDUList(list_of_HDUs)
    #output_hdulist = fits.HDUList(map(Cloud.to_ImageHDU, map(Cloud, unprocessed)))
    #del unprocessed

    #Output HDUList to File
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    #TODO output_hdulist.insert(0, fits.ImageHDU(data=backprojection_with_sets, header=fits.Header()))
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
        if SUFFIX not in xyt_filename:
            show(isolate_all(xyt_filename))
            #isolate_all(xyt_filename)
        else:
            show(*filament_average(xyt_filename, key='COLDENS'))
            #show(xyt_filename)



    #Cleanup and Exit
    #exit()




			