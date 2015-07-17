#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
# Initialization: Imports
#-----------------------------------------------------------------------------------------
from __future__ import division #Must be first line of code in the file
from astropy.io import fits
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sys
import string
import operator
import collections
import networkx as nx
import datetime
import scipy.stats
import numpy as np
import math

import matplotlib
#matplotlib.rcParams['backend'] = 'TkAgg'
matplotlib.rcParams['image.origin'] = 'lower'
#matplotlib.rcParams['figure.figsize'] = 6,6
#matplotlib.rcParams['figure.dpi'] = 200
from matplotlib import pyplot as plt
from matplotlib import gridspec

import rht #Latest version @https://github.com/seclark/RHT
import config

#-----------------------------------------------------------------------------------------
# Initialization: Program Settings
#-----------------------------------------------------------------------------------------

SUFFIX = '_filaments.'

skip = 1 #Number of HDUs that do not correspond to filaments in output files

#-----------------------------------------------------------------------------------------
# Initialization: Class Definitions
#-----------------------------------------------------------------------------------------
class Cloud:
    functions = {
        'MIN_X': (operator.attrgetter('min_x'), 'Lower-left x-coordinate of mask in backprojection'),
        'MAX_X': (operator.attrgetter('max_x'), 'Upper-right x-coordinate of mask in backprojection'),
        'MIN_Y': (operator.attrgetter('min_y'), 'Lower-left y-coordinate of mask in backprojection'),
        'MAX_Y': (operator.attrgetter('max_y'), 'Upper-right y-coordinate of mask in backprojection'),
        'AREA': (config.chained([operator.attrgetter('mask'), operator.attrgetter('size')]), 'Area covered by this mask'),
        'LITPIX': (config.chained([operator.attrgetter('points'), len]), 'Number of nonzero pixels in the mask'),
        'DIAG': (config.chained([operator.attrgetter('mask'), config.bridge(math.hypot, operator.attrgetter('shape'))]), 'Diagonal size of filament, corner-to-corner major-axis'),
        'OFF_DIAG': (lambda self: 4*len(self.points)/(np.pi*math.hypot(*self.mask.shape)), 'Off-diagonal size of filament, computed minor-axis')
    }
 
    '''
    @staticmethod
    def nonzero_data_from_HDU(hdu):
        #assert all([key in imageHDU.header])
        if isinstance(hdu, fits.ImageHDU) or isinstance(hdu, fits.BinTableHDU):
            return np.nonzero(Cloud(hdu).mask)
        else:
            raise ValueError('Cannot Create Cloud from anything but Image and BinTable HDUs')
    '''
    def as_ImageHDU(self):
        hdr = fits.Header()
        for k,v in Cloud.functions.iteritems():
            hdr[k] = (v[0](self), v[1])
        return fits.ImageHDU(data=self.mask, header=hdr)

    def as_BinTableHDU(self):
        hdr = fits.Header()
        for k,v in Cloud.functions.iteritems():
            hdr[k] = (v[0](self), v[1])

        #column = fits.Column(name=None, format=None, unit=None, null=None, bscale=None, bzero=None, disp=None, start=None, dim=None, array=None, ascii=None)
        xy = np.nonzero(self.mask)
        xs = fits.Column(name='xs', format='1I', array=xy[0])
        ys = fits.Column(name='ys', format='1I', array=xy[1])
        #ntheta = hthets.shape[1]
        #Hthets = fits.Column(name='hthets', format=str(int(ntheta))+'E', array=hthets)
        cols = fits.ColDefs([xs, ys])

        return fits.BinTableHDU(data=cols, header=hdr)

    def as_HDU(self, sparse=False):
        '''
        hdr = fits.Header()
        for k,v in Cloud.functions.iteritems():
            hdr[k] = (v[0](self), v[1])
        '''
        if sparse:
            return self.as_BinTableHDU()
        else:
            return self.as_ImageHDU()

    def __init__(self, list_of_points):
        #Expects a python list of two-integer tuples, corresponding the the x,y coordinate of the points original location in the backprojection 
        if isinstance(list_of_points, Cloud):
            self = list_of_points
            return 
            '''
            elif isinstance(list_of_points, fits.ImageHDU):
                self = from_ImageHDU(list_of_points)
                return
            '''
        elif isinstance(list_of_points, tuple):
            if len(list_of_points) > 2:
                self.points = list(set(list_of_points))
            elif len(list_of_points) == 2:
                self.points = [list_of_points]
        
        elif isinstance(list_of_points, set):
            self.points = list(list_of_points)
        else:
            self.points = list_of_points #TODO 

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
        assert proper_formatting(self.points) #TODO

        self.points.sort(key=operator.itemgetter(0))    
        self.min_x, self.max_x = self.points[0][0], self.points[-1][0]
        self.points.sort(key=operator.itemgetter(1))    
        self.min_y, self.max_y = self.points[0][1], self.points[-1][1]
        self.mask = np.zeros((1+self.max_x-self.min_x, 1+self.max_y-self.min_y), dtype=np.int_)
        for point in self.points:
            self.mask[point[0]-self.min_x][point[1]-self.min_y] = 1 

#-----------------------------------------------------------------------------------------
# File Handling
#-----------------------------------------------------------------------------------------
def handle(fileobj, mode='readonly'):

    if isinstance(fileobj, str):
        assert fileobj.endswith('.fits')
        assert rht.xyt_suffix in fileobj
        assert SUFFIX in fileobj
        print 'Accessing: '+fileobj+' '
        hdu_list = config.default_open(fileobj, mode)
        filaments_filename = fileobj

    elif isinstance(fileobj, fits.HDUList):
        hdu_list = fileobj
        filaments_filename = hdu_list.filename()
        assert filaments_filename.endswith('.fits')
        assert rht.xyt_suffix in filaments_filename
        assert SUFFIX in filaments_filename

    else:
        raise ValueError(repr(fileobj))

    #xyt_filename = string.rstrip(filaments_filename, SUFFIX+'fits')+'.fits' TODO xyt?

    return filaments_filename, hdu_list

#-----------------------------------------------------------------------------------------
# Post-Processing Functions
#-----------------------------------------------------------------------------------------
def filament_properties(filaments_filename, key, external_source=None):
    #Computes the average value of the dataset indicated by key or external_source for each filament
    #Saves the value to the associated header entry of each filament
    assert isinstance(key, str)
    filaments_filename, hdu_list = handle(filaments_filename, mode='update')

    if key in hdu_list[skip].header:
        print key+' keyword found in filament header...'
        if 'y' not in raw_input('Overwrite? ([no]/yes):  '):
            return hdu_list, key

    if external_source is not None:
        correlation_data = external_source
    elif key in config.source_data:
        correlation_data = config.source_data[key]
    else:
        raise KeyError('No source to average over for key: '+key)

    for i, hdu in enumerate(hdu_list[skip:]):
        #for suffix, func in config.methods.iteritems():
        for suffix in config.applicable_methods[key]:        
            func = config.methods[suffix]
            hdr = hdu.header
            hdr[key+suffix] = func(correlation_data[hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1][np.nonzero(hdu.data.T)])
        hdr[key] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    hdu_list.flush()
    return hdu_list, key

def plot(filaments_filename, key=None, out_name=None, show=True, cut=(lambda x: True)):
    if not show and out_name is None:
        'Unable to plot data without showing or saving it'
        return 


    filaments_filename, hdu_list = handle(filaments_filename)
    backproj = hdu_list.pop(0)
    #TODO make sure to pop skip hdus!s
    hdu_list = filter(cut, hdu_list) #TODO

    figure_args = {'figsize':(8,7), 'facecolor':'white','dpi':250}
    
    if key is None:
        display = np.zeros((backproj.header['NAXIS1'], backproj.header['NAXIS2']))
        for i, hdu in enumerate(hdu_list):
            hdr = hdu.header            
            display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = i
        fig = plt.figure(**figure_args)
        ax1 = fig.add_axes([0.0,0.0,1.0,1.0])
        ax1.imshow(display.T)
        key = 'Filaments'

    elif isinstance(key, str) and key in hdu_list[skip].header:
        displays = dict()
        #for suffix in config.methods.iterkeys():
        for suffix in config.applicable_methods[key]:        
            displays[key+suffix] = np.zeros((backproj.header['NAXIS1'], backproj.header['NAXIS2']))
        for hdu in hdu_list:
            hdr = hdu.header
            #for suffix in config.methods.iterkeys(): 
            for suffix in config.applicable_methods[key]:           
                displays[key+suffix][hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][np.nonzero(hdu.data)] = hdr[key+suffix] #config.post_processing[key](hdr[key+suffix]) # 
        #fig, (ax1, ax2) = plt.subplots(1,2,sharey='row', **figure_args)
        fig = plt.figure(**figure_args)
        Nplots = len(config.applicable_methods[key])
        r=4
        gs = gridspec.GridSpec(Nplots*r, 5*r)
        suffixes = config.applicable_methods[key]
        for i in range(Nplots):
            title = key+suffixes[i]
            ax1 = plt.subplot(gs[1+r*i:r*(i+1),0:-4]) #plt.subplot2grid((1,10), (0,0), colspan=9)
            ax2 = plt.subplot(gs[1+r*i:r*(i+1)-1,-1:-2]) #plt.subplot2grid((1,10), (0,1), colspan=1)
            raw_data = config.post_processing[key](map(operator.itemgetter(title), map(operator.attrgetter('header'), hdu_list)))
            ax2.hist(raw_data, bins=50, orientation='horizontal', histtype='stepfilled') #range=(18, 25))
            plt.colorbar(ax1.imshow(config.post_processing[key](displays[title].T), cmap = "YlOrRd"), ax=ax2, fraction=0.10)
            plt.title(title, fontsize=8)
    else:
        print 'Unable to plot data using the given key: '+str(key)
        return
        
    
    plt.suptitle(key+' from '+filaments_filename, fontsize=12)
    if out_name is not None and isinstance(out_name, str):
        plt.savefig(out_name, dpi=500, format='png')
    if show:
        plt.show()
    plt.clf() 

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------

def isolate_all(xyt_filename, BINS=6, generateHDU=Cloud.as_ImageHDU): #Cloud.as_BinTableHDU):

    #Read in RHT Output from filename_xyt??.fits
    assert xyt_filename.endswith('.fits')
    assert rht.xyt_suffix in xyt_filename
    assert SUFFIX not in xyt_filename
    print 'Accessing: '+xyt_filename+' '
    hdu_list = config.default_open(xyt_filename)
    ntheta = hdu_list[0].header['NTHETA']
    wlen = hdu_list[0].header['WLEN']
    frac = hdu_list[0].header['FRAC']
    naxis1 = hdu_list[0].header['NAXIS1']
    naxis2 = hdu_list[0].header['NAXIS2']
    original = hdu_list[0].header['ORIGINAL']
    Hi = hdu_list[1].data['hi'] 
    Hj = hdu_list[1].data['hj'] 

    #Compute TheteRHT for all pixels given, then bin by theta
    B = map(rht.theta_rht,hdu_list[1].data['hthets']) #List of theta_rht values
    C = np.multiply(np.asarray(B), BINS/np.pi).astype(np.int_)
    
    def rel_add((a,b), (c,d)):
        return a+c,b+d
    
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
        #progress_bar = Progress(problem_size, message=message, incrementing=True)

        point_dict = dict([x[::-1] for x in enumerate(raw_points)])
        set_dict = collections.defaultdict(list)
        theta_dict = dict()

        for coord in raw_points:
            #rht.update_progress(0.3*(i/problem_size), message=message)
            #progress_bar.update()
            theta_dict[coord] = B[point_dict[coord]]
            for rel_coord in search_pattern:
                try:
                    j = point_dict[rel_add(coord, rel_coord)]
                    set_dict[point_dict[coord]].append(j)
                except Exception:
                    continue
        
        G = nx.from_dict_of_lists(set_dict) #Undirected graph made using set_dict as an adjacency list
        del set_dict 
        
        #progress_bar = Progress(problem_size, message=message, incrementing=False)
        sources = range(problem_size)
        flags = np.ones((problem_size), dtype=np.int_)
        while len(sources) > 0: 
            source = sources.pop()
            if not flags[source]:
                continue
            else:
                #rht.update_progress(0.3+0.3*(1.0-len(sources)/problem_size), message=message)
                #progress_bar.update(len(sources))
                try:
                    for member in nx.descendants(G, source):
                        flags[member] = False
                        point_dict[raw_points[member]] = source
                        G.remove_node(member) #TODO Remove members from G if that would speed up subsequent calls?
                except nx.NetworkXError:
                    #Assume we hit an isolated pixel (never made it into G) and move on
                    pass
        del sources, flags, G

        #************************************************************************************************
        #finished_map = np.negative(np.ones((naxis1, naxis2), dtype=np.int64))
        #for pt in raw_points:
            #finished_map[pt] = point_dict[pt]

        histogram = np.bincount(map(point_dict.get, raw_points))
        mask = np.nonzero(histogram >= int(frac*wlen))[0]
        del histogram

        #progress_bar = Progress(problem_size, message=message, incrementing=False)
        mask_dict = dict([x[::-1] for x in enumerate(mask)])
        out_clouds = collections.defaultdict(list)

        while len(point_dict) > 0:
            temp = point_dict.popitem()
            try:
                #Keying into mask_dict is the only operation that ought to throw an exception 
                out_clouds[mask_dict[temp[1]]].append(temp[0])
                #progress_bar.update(len(point_dict))
                #rht.update_progress(0.6+0.399*(1.0-len(point_dict)/problem_size), message=message)
            except Exception:
                continue

        while len(out_clouds) > 0:
            #cloud = out_clouds.popitem()[1]
            #unprocessed.append(out_clouds.popitem()[1])
            list_of_HDUs.append(Cloud(out_clouds.popitem()[1]).as_HDU()) #TODO Incorporate theta_dict
            #list_of_HDUs.append(generateHDU(Cloud(out_clouds.popitem()[1])))
        #rht.update_progress(1.0, final_message='Finished joining '+str(problem_size)+' points! Time Elapsed:')
        
    #Convert lists of two-integer tuples into ImageHDUs
    #unprocessed.sort(key=len, reverse=True)
    #output_hdulist = fits.HDUList(map(Cloud.as_ImageHDU, map(Cloud, unprocessed)))
    #del unprocessed
    list_of_HDUs.sort(key=config.chained([operator.attrgetter('header'), operator.itemgetter('DIAG')]), reverse=True) #lambda h: h.header['DIAG'], reverse=True)
    output_hdulist = fits.HDUList(list_of_HDUs)
    del list_of_HDUs

    #Output HDUList to File
    output_filename = string.join(string.rsplit(xyt_filename, '.', 1), SUFFIX)
    output_hdulist.insert(0, hdu_list[0].copy()) #TODO Introduces Errors in Reading FITS File 
    output_hdulist.writeto(output_filename, output_verify='silentfix', clobber=True, checksum=True)

    print 'Results successfully output to '+output_filename
    return output_filename

#-----------------------------------------------------------------------------------------
# Command Line Mode
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    #Interpret Arguments
    parser = ArgumentParser(description="Run Fiber Isolation on 1 or more RHT output (.fits) files", usage='%(prog)s [options] file(s)', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help="FITS file(s)")
    if len(sys.argv) == 1: # no arguments given, so add -h to get help msg
        sys.argv.append('-h')
    args = parser.parse_args()

    #Do Processing
    for filename in args.files: # loop over input files
        if SUFFIX not in filename:
            #This is hopefully an XYT RHT output file...
            plot(isolate_all(filename))
            #isolate_all(xyt_filename)
        else:
            #This would hopefully be a formerly run _filament.fits file...
            #plot(filename)
            #plot(*filament_properties(filename, key='COLDENS'), out_name=filename[:-5]+'_COLDENS.png')
            #plot(*filament_properties(filename, key='GALFA0'), out_name=filename[:-5]+'_GALFA0.png')
            #plot(*filament_properties(filename, key='B'))
            plot(filename, key='GALFA0', out_name=filename[:-5]+'_GALFA0_50.png', cut=lambda h: h.header['B_MIN'] > 50.0)

    #Cleanup and Exit
    #exit()




			