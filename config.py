#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel
#for GALFA0 conversion http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php

from astropy.io import fits
import operator
import numpy as np
import scipy.stats
import collections

identity = lambda x:x

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

def bridge(funca, funcb):
    #funcb must return a tuple of length n
    #funca must expect n arguments
    return lambda x: funca(*funcb(x))

def default_open(filename, mode='readonly'):
    return fits.open(filename, mode=mode, memmap=True, save_backup=False, checksum=True, ignore_missing_end=True)

#Known methods that can produce a single value from a masked area of source data
methods = {
    '_AVG':np.nanmean,
    '_MED':scipy.stats.nanmedian,
    '_TOT':np.nansum,    
    '_MIN':np.nanmin,
    '_MAX':np.nanmax
}

sources = {
    'COLDENS': (
        'D:/LAB_corrected_coldens.fits', 
        chained([default_open, operator.itemgetter(0), operator.attrgetter('data')]), 
        np.log10,
        ['_AVG', '_MED', '_TOT']
    ),
    'GALFA0': (
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
    )
}#operator.methodcaller('sum', axis=0)

source_data = dict()
applicable_methods = dict()
post_processing = collections.defaultdict(identity)
for k,v in sources.iteritems():
    try:
        source_data[k] = v[1](v[0])
        post_processing[k] = v[2]
        applicable_methods[k] = v[3]
    except Exception as e:
        print 'WARNING: Unable to find source data for '+str(k)+' keyword using '+str(v[0])+' as an initial input.. (Modify config.py for better results!)'
        print repr(e)





