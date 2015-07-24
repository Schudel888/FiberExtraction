#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel
#for GALFA0 conversion http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php
#-----------------------------------------------------------------------------------------
# Initialization: Imports
#-----------------------------------------------------------------------------------------
from __future__ import division #Must be first line of code in the file
from astropy.io import fits
import operator
import numpy as np
import scipy.stats
import math
import datetime

import myfavoritefiber

print ''
print 'Starting config.py'
print '{'

#-----------------------------------------------------------------------------------------
# Initialization: Convenience Functions
#-----------------------------------------------------------------------------------------
identity = lambda x:x

def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def passive_constant(value):
    #Given passive_constant(value)(x) --> returns value #Ignores x
    return lambda x: value

def active_constant(function):
    #Given active_constant(function)(x) --> returns function() #Ignores x
    return lambda x: function()

def stubborn_constant(function, value):
    #Given stubborn_constant(function, value)(x) --> returns function(value) #Ignores x
    return lambda x: function(value)

def chained(functions):
    #Given chained([f,g,h,..])(x) --> returns ..h(g(f(x)))
    assert isinstance(functions, list)
    z = len(functions)
    if z == 0:
        return identity
    elif z == 1:
        return functions[0]
    else:
        return lambda x: functions[-1]( chained(functions[0:-1])(x) )

def bridged(funca, funcb):
    #funcb must return a tuple of length n
    #funca must expect n arguments
    return lambda x: funca(*funcb(x))

def rel_add((a,b), (c,d)):
    return a+c,b+d

def default_open(filename, mode='readonly'):
    print 'Accessing: '+filename+' '
    #try:
    return fits.open(filename, mode=mode, memmap=True, save_backup=False, checksum=True, ignore_missing_end=True)
    #except Exception:
    #return fits.open(filename, mode=mode, memmap=False, save_backup=False, checksum=True, ignore_missing_end=True)

def GALFAx(integer):
    assert isinstance(integer, int)
    assert 0<=integer<=40
    gx = {'GALFA'+str(integer): 
    (
        'D:/SC_241.66_28.675.best.fits', 
        lambda x: np.nan_to_num(default_open(x)[0].data[integer]).clip(0.0, np.inf)*1.823E18,
        np.log10,
        ['_AVG', '_MED', '_TOT']
    )}
    return gx

'''
def GALFAxN(integer):
    #arXiv:0810.1283 #http://sites.google.com/site/galfahi/
    assert isinstance(integer, int)
    return 'GALFA'+str(integer), (

        'D:/SC_241.66_28.675.best.fits', 
        lambda x: np.nan_to_num(default_open(x)[0].data[22]).clip(0.0, np.inf)*1.823E18,
        np.log10,
        ['_AVG', '_MED', '_TOT']
    )

def ONOFF(integer):
    assert isinstance(integer, int)

'''

#-----------------------------------------------------------------------------------------
# Initialization: Class Definitions
#-----------------------------------------------------------------------------------------
class Cloud:
    functions = {
        'MIN_X': (operator.attrgetter('min_x'), 'Lower-left x-coordinate of mask in backprojection'),
        'MAX_X': (operator.attrgetter('max_x'), 'Upper-right x-coordinate of mask in backprojection'),
        'MIN_Y': (operator.attrgetter('min_y'), 'Lower-left y-coordinate of mask in backprojection'),
        'MAX_Y': (operator.attrgetter('max_y'), 'Upper-right y-coordinate of mask in backprojection'),
        'AREA': (chained([operator.attrgetter('mask'), operator.attrgetter('size')]), 'Area covered by this mask'),
        'LITPIX': (chained([operator.attrgetter('points'), len]), 'Number of nonzero pixels in the mask'),
        'DIAG': (chained([operator.attrgetter('mask'), bridged(math.hypot, operator.attrgetter('shape'))]), 'Diagonal size of filament, ~major-axis'),
        'OFF_DIAG': (lambda x: 4.0*float(len(x.points))/(np.pi*math.hypot(*x.mask.shape)), 'Off-diagonal size of filament, ~minor-axis')
    }

    @staticmethod
    def on_and_off_masks_from_HDU(hdu, transpose=False):
        if isinstance(hdu, fits.BinTableHDU) or isinstance(hdu, fits.ImageHDU):
            cloud = Cloud(hdu)
            on_mask = cloud.mask
            if transpose:
                on_mask = on_mask.T
            off_mask = myfavoritefiber.off_fiber(on_mask)
            #off_mask #TODO 
            return on_mask, off_mask 
        else:
            raise ValueError('Your hdu could not be resolved to a known type in on_and_off_masks_from_HDU()')

    @staticmethod
    def nonzero_data_from_HDU(hdu, transpose=False):
        if isinstance(hdu, fits.BinTableHDU):
            if transpose:
                return (hdu.data['xs'], hdu.data['ys'])
            else:
                return (hdu.data['ys'], hdu.data['xs'])
        elif isinstance(hdu, fits.ImageHDU):
            if transpose:
                return np.nonzero(hdu.data.T)
            else:
                return np.nonzero(hdu.data)
        else:
            raise ValueError('Your hdu could not be resolved to a known type in nonzero_data_from_HDU()')

    def as_ImageHDU(self):
        hdr = fits.Header()
        for k,v in Cloud.functions.iteritems():
            hdr[k] = (v[0](self), v[1])
        hdr['SPARSE'] = False
        return fits.ImageHDU(data=self.mask, header=hdr)

    def as_BinTableHDU(self):
        hdr = fits.Header()
        for k,v in Cloud.functions.iteritems():
            hdr[k] = (v[0](self), v[1])
        hdr['SPARSE'] = True
        #column = fits.Column(name=None, format=None, unit=None, null=None, bscale=None, bzero=None, disp=None, start=None, dim=None, array=None, ascii=None)
        xy = np.nonzero(self.mask)
        xs = fits.Column(name='xs', format='1I', array=xy[0])
        ys = fits.Column(name='ys', format='1I', array=xy[1])
        #ntheta = hthets.shape[1]
        #Hthets = fits.Column(name='hthets', format=str(int(ntheta))+'E', array=hthets)
        cols = fits.ColDefs([xs, ys])

        return fits.BinTableHDU(data=cols, header=hdr)

    def as_HDU(self, sparse=False):
        if sparse:
            return self.as_BinTableHDU()
        else:
            return self.as_ImageHDU()

    def __init__(self, list_of_points):
        #Expects a python list of two-integer tuples, corresponding the the x,y coordinate of the points original location in the backprojection 
        if isinstance(list_of_points, Cloud):
            self.points = list_of_points.points
        elif isinstance(list_of_points, fits.BinTableHDU) or isinstance(list_of_points, fits.ImageHDU):
            min_coord = (list_of_points.header['MIN_X'], list_of_points.header['MIN_Y'])
            rel_coords = zip(*Cloud.nonzero_data_from_HDU(list_of_points, transpose=False))
            self.points = [rel_add(min_coord, rel_coord) for rel_coord in rel_coords]
            del min_coord, rel_coords
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
# Initialization: Useful Data
#-----------------------------------------------------------------------------------------
#Known methods that can produce a single value from a masked area of source data
methods = {
    '_AVG':np.nanmean,
    '_MED':scipy.stats.nanmedian,
    '_TOT':np.nansum,    
    '_MIN':np.nanmin,
    '_MAX':np.nanmax
}

'''
'_ONOFF':,'_VWIDTH':, '_VMAX', 
'''
#None for v[0] in source indicates all functions should take a Cloud object as their argument
#All else for v[0] assumes a 2Dimensional mask of data
sources = {
    'COLDENS': (
        'D:/LAB_corrected_coldens.fits', 
        chained([default_open, operator.itemgetter(0), operator.attrgetter('data')]), 
        np.log10,
        ['_AVG', '_MED', '_TOT']
    ),
    'GALFA': (
        'D:/SC_241.66_28.675.best.fits', 
        chained([default_open, operator.itemgetter(0), operator.attrgetter('data'), operator.itemgetter(np.s_[16:25]), lambda x: np.nansum(x, axis=0), lambda x: np.multiply(x, 1.823E18)]),
        np.log10,
        ['_AVG', '_MED', '_TOT']
    ),
    'B': (
        'D:/SC_241_2d_bs.npy',
        np.load,
        identity,
        ['_MAX', '_MIN']
    ),
    'L': (
        'D:/SC_241_2d_ls.npy',
        np.load,
        identity,
        ['_MAX', '_MIN']
    ),
    '': (
        None,
        identity,
        identity,
        list(Cloud.functions.iterkeys())
    )
}

source_data = dict()
applicable_methods = dict()
post_processing = dict() #collections.defaultdict(identity)

def clear():
    source_data.clear()
    applicable_methods.clear()
    post_processing.clear()

def exclude(key):
    if key in sources:
        del sources[key]
    if key in source_data:
        del source_data[key]
    if key in post_processing:
        del post_processing[key]
    if key in applicable_methods:
        del applicable_methods[key]
        
def include(dictionary):
    assert isinstance(dictionary, dict)
    for k,v in dictionary.iteritems():
        try:
            source_data[k] = v[1](v[0])
            post_processing[k] = v[2]
            applicable_methods[k] = v[3]
        except Exception as e:
            print 'WARNING: Unable to find source data for '+str(k)+' keyword using '+str(v[0])+' as an initial input.. (Modify config.py for better results!)'
            print repr(e)
'''
include({'': (
    None,
    identity,
    identity,
    list(Cloud.functions.iterkeys())
)})
'''
#include(sources)
#for i in [0]+range(10,25)+range(33,37):
#include(GALFAx(i))

print '} '
print 'Finishing config.py '
print ''



