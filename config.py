#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel
#for GALFA0 conversion http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php

from astropy.io import fits
import operator
import numpy as np
import scipy.stats
import collections
import math
import datetime

#-----------------------------------------------------------------------------------------
# Initialization: Convenience Functions
#-----------------------------------------------------------------------------------------
identity = lambda x:x

def passive_constant(value):
    #Given active_constant(value)(x) --> returns value #Ignores x
    return lambda x: value

def active_constant(function):
    #Given active_constant(f)(x) --> returns f() #Ignores x
    return lambda x: function()

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

def default_open(filename, mode='readonly'):
    return fits.open(filename, mode=mode, memmap=True, save_backup=False, checksum=True, ignore_missing_end=True)

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
        'DIAG': (chained([operator.attrgetter('mask'), bridged(math.hypot, operator.attrgetter('shape'))]), 'Diagonal size of filament, corner-to-corner major-axis'),
        'OFF_DIAG': (lambda self: 4*len(self.points)/(np.pi*math.hypot(*self.mask.shape)), 'Off-diagonal size of filament, computed minor-axis')
    }
 
    @staticmethod
    def nonzero_data_from_HDU(hdu, transpose=False):
        '''
        try:
            #TODO
            sparse = hdu.header['SPARSE']
        except Exception:
            sparse = None
        '''
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
        elif isinstance(list_of_points, tuple):
            if len(list_of_points) > 2:
                self.points = list(set(list_of_points))
            elif len(list_of_points) == 2:
                self.points = [list_of_points]
        
        elif isinstance(list_of_points, set):
            self.points = list(list_of_points)
        else:
            try:
                self.points = zip(*nonzero_data_from_HDU(hdu, transpose=False))
            except Exception:
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
    '_MAX':np.nanmax,
    '':chained([active_constant(datetime.datetime.now), operator.methodcaller('strftime', "%Y-%m-%d %H:%M:%S")])
}

sources = {
    'COLDENS': (
        'D:/LAB_corrected_coldens.fits', 
        chained([default_open, operator.itemgetter(0), operator.attrgetter('data')]), 
        np.log10,
        ['_AVG', '_MED', '_TOT','']
    ),
    'GALFA0': (
        'D:/SC_241.66_28.675.best.fits', 
        chained([default_open, operator.itemgetter(0), operator.attrgetter('data'), operator.itemgetter(np.s_[16:25]), lambda x: np.nansum(x, axis=0), lambda x: np.multiply(x, 1.823E18)]),
        np.log10,
        ['_AVG', '_MED', '_TOT','']
    ), 
    'B': (
        'D:/SC_241_2d_bs.npy',
        np.load,
        identity,
        ['_MAX', '_MIN']
    ),
    '': (
        None,
        identity,
        identity,
        Cloud.functions.iterkeys()
    )
}#operator.methodcaller('sum', axis=0)

source_data = dict()
applicable_methods = dict()
post_processing = dict() #collections.defaultdict(identity)
for k,v in sources.iteritems():
    try:
        source_data[k] = v[1](v[0])
        post_processing[k] = v[2]
        applicable_methods[k] = v[3]
    except Exception as e:
        print 'WARNING: Unable to find source data for '+str(k)+' keyword using '+str(v[0])+' as an initial input.. (Modify config.py for better results!)'
        print repr(e)

'''
possible_keys = list()
for key in sources.iterkeys():
    for suffix in sources[key][3]:
        possible_keys.append(key+suffix) 

def is_complete(cloud_HDU, fix=True):
    try:
        hdr = cloud_HDU.header
        for key in possible_keys:
            assert key in hdr
        return True
    except Exception:
        if not fix:
            return False
        else:
            print 'Unable to fix Cloud in config.is_complete'
            #TODO
            return False
'''




