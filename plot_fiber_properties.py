from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pyfits
import pylab
import time
import copy
import math
import os.path
import cPickle as pickle
from matplotlib import rc
from astropy.io import fits
from astropy import wcs
import scipy.optimize
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import astropy.coordinates as coord
from astropy.io import fits
from astropy import units as u
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import random

def get_filament_data(channel, verbose = False):

    root = "/Volumes/DataDavy/GALFA/lowell_fibers_project/"
    fn = root+"SC_241.66_28.675.best_"+str(channel)+"_xyt_w75_s15_t70_filaments.fits"
    
    hdu = fits.open(fn)
    if verbose == True:
        print(hdu[0].header)
    
    return hdu
    
def get_onoff_stats(channel):

    hdu = get_filament_data(channel)
    kword = "GALFA"+str(channel)
    
    onoff_avg = np.zeros(len(hdu))
    onoff_med = np.zeros(len(hdu))
    bmax = np.zeros(len(hdu))
    bmin = np.zeros(len(hdu))
    hdu_num = np.zeros(len(hdu))
    
    for i in xrange(len(hdu)):
        try: 
             onoff_avg[i] = hdu[i].header[kword + str("_ONOFF_AVG")]
             onoff_med[i] = hdu[i].header[kword + str("_ONOFF_MED")]
             bmax[i] = hdu[i].header["B_MAX"]
             bmin[i] = hdu[i].header["B_MIN"]
             hdu_num[i] = i
             
        except KeyError:
            onoff_avg[i] = None
            onoff_med[i] = None
            bmax[i] = None
            bmin[i] = None
            hdu_num[i] = None
    
    return onoff_avg, onoff_med, bmax, bmin, hdu_num



def plot_onoff_column_avg_hist(chans):

    fig = plt.figure(facecolor = "white")
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    for _c in chans:
        onoff_avg, onoff_med, bmax, bmin, hdu_num = get_onoff_stats(_c)
        
        bmean = (bmax + bmin)/2.0
        
        #ax1.hist(onoff_avg[~np.isnan(onoff_avg)], label = r""+str(_c), bins = 200, histtype = "step")
        #ax2.hist(onoff_med[~np.isnan(onoff_med)], label = r""+str(_c), bins = 200, histtype = "step")
        
        ax1.hist(onoff_avg[(np.isnan(onoff_avg) == False) & (bmean > 50) & (onoff_avg > 1E18)], label = r""+str(_c), bins = 10, histtype = "step")
        ax2.hist(onoff_med[(np.isnan(onoff_avg) == False) & (bmean > 50) & (onoff_avg > 1E18)], label = r""+str(_c), bins = 10, histtype = "step")
    
    
    ax1.legend()
    ax2.legend()
    ax1.set_title(r"$\mathrm{Average}$ $\mathrm{on - off}$ $N_{HI}$ $b > 70$", size = 20)
    ax2.set_title(r"$\mathrm{Median}$ $\mathrm{on - off}$ $N_{HI}$ $b > 70$", size = 20)
    #ax1.semilogx()
    #ax2.semilogx()
    ax1.set_xscale('log')
    ax2.set_xscale('log')
            
            
            