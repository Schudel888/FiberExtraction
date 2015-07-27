from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pylab
import time
import copy
from matplotlib import rc
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import norm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

"""Analysis of a single fiber"""
def fibermask():
    fibermask = np.load("/Users/susanclark/Documents/Fibers/Data/GALFA/fibermask.npy")
    return fibermask

def gromask(mask):
    #Shift up, then L/R
    maskup = np.roll(mask, 1, axis=0)
    maskup1 = np.roll(maskup, 1, axis=1)
    maskup2 = np.roll(maskup, -1, axis=1)
    
    #Shift down, then L/R
    maskdown = np.roll(mask, -1, axis=0)
    maskdown1 = np.roll(maskdown, 1, axis=1)
    maskdown2 = np.roll(maskdown, -1, axis=1)
    
    #Shift right, then U/D
    maskright = np.roll(mask, 1, axis=1)
    maskright1 = np.roll(maskright, 1, axis=0)
    maskright2 = np.roll(maskright, -1, axis=0)
    
    #Shift left, then U/D
    maskleft = np.roll(mask, -1, axis=1)
    maskleft1 = np.roll(maskleft, 1, axis=0)
    maskleft2 = np.roll(maskleft, -1, axis=0)

    gromask = np.ceil((mask+maskup1+maskup2+maskdown1+maskdown2+maskright1+maskright2+maskleft1+maskleft2)/9.0)
    
    return gromask
    
def off_fiber(mask):
    masksum = np.sum(mask)
    off = copy.copy(mask)
    
    while not off.all() and np.sum(off - mask) < masksum:
        off = gromask(off)
    
    off = off - mask
        
    return off

 
"""Display fiber on/off"""
def showspectra(): 
    from astropy.io import fits
 
    #Obtain masks
    mask = fibermask()
    offmask = off_fiber(mask)
    
    """
    #Import galfa cube
    hdulist = fits.open("/Users/susanclark/Documents/Fibers/Data/GALFA/zenstrip.fits")
    galfa = hdulist[0].data
    
    #Axis 0 is velocity
    galfa1 = galfa[:, 180:310, 1900:2150]
    del galfa
    del hdulist
    print galfa1.shape   
    
    np.save("/Users/susanclark/Documents/Fibers/Data/GALFA/zenstrip_slice.npy", galfa1)
    """
    galfa = np.load("/Users/susanclark/Documents/Fibers/Data/GALFA/zenstrip_slice.npy")
    
    #maskindx = np.where(mask > 0)
    
    onfiber = np.sum(galfa[:, mask>0], axis=1)
    
    offfiber = np.sum(galfa[:, offmask>0], axis=1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #Brightness Temperature pixel scale (from zenstrip hdr)
    #bscale = 0.00662270000000
    
    #First pixel velocity (vlsr)
    crval3 = -74716.4687500 #m/s, I assume. Does not say in zenstrip hdr.
    cdelt3 = 736.117187500
    
    vels = np.linspace(crval3, crval3+cdelt3*len(onfiber), len(onfiber))
    vels = vels/1000 #To km/s
    
    """
    ax.plot(vels, onfiber*bscale/np.sum(mask), color="black")
    ax.plot(vels, offfiber*bscale/np.sum(offmask), ":", color="black")
    """
    ax.set_ylabel(r'$T_B$ $(K)$', size=20)
    ax.set_xlabel(r'$V_{LSR}$ $(km/s)$', size=20)

    ax.set_xlim(np.min(vels), np.max(vels))
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    """
    ax.legend(['On Fiber', 'Off Fiber'])
    """
    ax.plot(vels, (onfiber/np.sum(mask))-(offfiber/np.sum(offmask)), color="black")
    
    pylab.savefig("/Users/susanclark/Documents/Fibers/Paper/onoffdiff.pdf", bbox_inches=0)
    
#Figure showing fiber and on-off spectrum
def paperfig():
    fig = plt.figure()
    #ax1 = fig.add_subplot(121)
    #ax2 = fig.add_subplot(122)
    ax1 = fig.add_subplot(111)
    
    #Obtain masks
    mask = fibermask()
    offmask = off_fiber(mask)
    
    galfa = np.load("/Users/susanclark/Documents/Fibers/Data/GALFA/zenstrip_slice.npy")
    
    #maskindx = np.where(mask > 0)
    
    onfiber = np.sum(galfa[:, mask>0], axis=1)
    offfiber = np.sum(galfa[:, offmask>0], axis=1)
    
    #First pixel velocity (vlsr)
    crval3 = -74716.4687500 #m/s, I assume. Does not say in zenstrip hdr.
    cdelt3 = 736.117187500
    
    vels = np.linspace(crval3, crval3+cdelt3*len(onfiber), len(onfiber))
    vels = vels/1000 #To km/s
    
    ax1.set_ylabel(r'$T_B$ $(K)$', size=20)
    ax1.set_xlabel(r'$v_{LSR}$ $(km/s)$', size=20)

    #ax.set_xlim(np.min(vels), np.max(vels))
    ax1.set_xlim(-20, 20)
    
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
   
    plotme = (onfiber/np.sum(mask))-(offfiber/np.sum(offmask))
    ax1.plot(vels, (onfiber/np.sum(mask))-(offfiber/np.sum(offmask)), color="black")
    
    print 'column density is', (1.82*10**18)*np.sum(plotme)
    
    hdulist = fits.open('/Users/susanclark/Documents/Fibers/Data/destripe_zenith_WRONG_NAXIS3.fits')
    galfa = hdulist[0].data
    
    galfa4 = galfa[4, :, :]
    galfa4[galfa4 < 0] = None
    galfa4 = galfa4[180:310, 1900:2150]
    
    galfa5 = galfa[5, :, :]    
    galfa5[galfa5 < 0] = None
    galfa5 = galfa5[180:310, 1900:2150]
    
    minvel = -6.9936876
    maxvel = -1.10475
    
    maxk = 0.6
    mink = -0.1
    ax1.set_ylim(mink, maxk)
    
    ax1.plot([minvel, minvel], [mink, maxk], ':', lw=2, color='black')
    ax1.plot([maxvel, maxvel], [mink, maxk], ':', lw=2, color='black')
    
    mu, sigma = norm.fit(plotme)
    print mu, sigma
    
    #Fit a Gaussian to the data
    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 1., 1.]

    coeff, var_matrix = curve_fit(gauss, vels, plotme, p0=p0)
    hist_fit = gauss(vels, *coeff)
    ax1.plot(vels, hist_fit, color='green')

    print coeff[0], coeff[1], 'sigma: ', coeff[2]    
    #ax2.imshow(np.log10(galfa4+galfa5), cmap='gray', vmin=-0.8, vmax=-0.5)
    #print np.min(np.log10(galfa4+galfa5)), np.max(np.log10(galfa4+galfa5)) 

def bothchan():
    from starclass import Star, lineSeg, xy_to_radec_GALFA, radec_to_xy_GALFA
    fig = plt.figure()
    ax = fig.add_subplot(111)

    hdulist = fits.open('/Users/susanclark/Documents/Fibers/Data/destripe_zenith_WRONG_NAXIS3.fits')
    galfa = hdulist[0].data
    
    mask = fibermask()
    offmask = off_fiber(mask)
    
    galfa4 = galfa[4, :, :]
    galfa4[galfa4 < 0] = None
    galfa4 = galfa4[180:310, 1900:2150]
    
    galfa5 = galfa[5, :, :]    
    galfa5[galfa5 < 0] = None
    galfa5 = galfa5[180:310, 1900:2150]
    
    ax.imshow((galfa4+galfa5), cmap='gray')
    ax.contour(mask, levels=[0], color='blue')
    ax.contour(offmask, levels=[0], color='red')
    
    ytickps = []
    ytickls = []
    xtickps = []
    xtickls = []
        
    decs = [26, 27, 28, 30, 32] #added 27
    for i in xrange(len(decs)):
        x, y = radec_to_xy_GALFA(231.67510335, decs[i])
        print x, y
        ytickps = np.append(ytickps, y-180)
       
    ras = [209, 210, 211, 212, 213]
    for i in xrange(len(ras)):
        x, y = radec_to_xy_GALFA(ras[i], 23.82504765)
        print x, y
        xtickps = np.append(xtickps, x-1900)
    
    #old
    #ras = [197, 204, 211, 218, 225, 232]
    #ytickls = ['24$^o$', '26$^o$', '28$^o$', '30$^o$']
    #xtickls = [r'15$^h$28$^m$', r'15$^h$00$^m$', r'14$^h$32$^m$', r'14$^h$04$^m$', r'13$^h$36$^m$', r'13$^h$08$^m$']
    
    ytickls = ['24$^o$', '25$^o$', '26$^o$', '28$^o$', '30$^o$']
    #xtickls = [r'15$^h$28$^m$', r'15$^h$00$^m$', r'14$^h$32$^m$', r'14$^h$08$^m$', r'14$^h$04$^m$', r'13$^h$36$^m$', r'13$^h$08$^m$']
    xtickls = [r'13$^h$56$^m$', r'14$^h$08$^m$', r'14$^h$04$^m$', r'14$^h$08$^m$', r'14$^h$12$^m$']

    ax.set_xticks(xtickps)
    ax.set_xticklabels(xtickls, size=25)
    ax.set_yticks(ytickps)
    ax.set_yticklabels(ytickls, size=25)
    
    ax.set_ylim(0, 310-180)
    ax.set_xlim(0, 2150-1900)
            

def showfiber():
    from makefigs import popfig

    #Obtain masks
    mask = fibermask()
    offmask = off_fiber(mask)
    
    print mask.shape
    print offmask.shape

    hdulist = fits.open('/Users/susanclark/Documents/Fibers/Data/destripe_zenith_WRONG_NAXIS3.fits')
    galfa = hdulist[0].data
    galfa = galfa[5, :, :]    
    galfa[galfa < 0] = None
    
    galfa = galfa[180:310, 1900:2150]
    print galfa.shape
    
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    
    on = copy.copy(galfa)
    on[mask < 1] = 0
    ax1.imshow(np.log10(on), cmap='gray')
    
    off = copy.copy(galfa)
    off[offmask < 1] = 0
    ax2.imshow(np.log10(off), cmap='gray')

    ax3.imshow(np.log10(galfa), cmap='gray')
    
    ipoints, jpoints, hthets = loadRHT_galfa_combined_inplace(101, 11, 70, 5)
    backp1 = popfig(np.sum(hthets,axis=1), jpoints, ipoints, 600, 7800)
    backp = backp1[180:310, 1900:2150]
    
    ax4.imshow(backp, cmap='gray')
    

def single():
    from makefigs import popfig, findindex

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax2 = fig.add_subplot(122)
    #ax2.set_aspect('equal')

    """
    hdulist = fits.open('/Users/susanclark/Documents/Fibers/Data/destripe_zenith_WRONG_NAXIS3.fits')
    galfa = hdulist[0].data
    galfa = galfa[5, :, :]    
    galfa[galfa < 0] = None 
    """
    
    wlen = 101
    smr = 11
    thresh = 70
    ch = 5
    
    root = '/Users/susanclark/Documents/Fibers/Data/GALFA/'
    ipoints = np.load(root+'galfa_hi_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_masked_combined.npy')
    jpoints = np.load(root+'galfa_hj_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_masked_combined.npy')
    hthets = np.load(root+'galfa_hthets_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_masked_combined.npy')  
    thetas = np.load('/Users/susanclark/Documents/thets_w'+str(wlen)+'.npy')
    
    print 'shape of hthets', hthets.shape
    
    nthetas = len(thetas)
    th1 = 50
    th2 = 80
    
    indx = findindex(jpoints, ipoints, 600, 7800, 2164, 259)
    indx2 = findindex(jpoints, ipoints, 600, 7800, 2164, 260)
    
    print len(indx[0]), len(indx2[0])
    print indx, indx2
    indx = np.int(indx[0])
    print indx
    
    counts = np.zeros((7800, 600), np.float_)
    for i in xrange(7800):
        for j in xrange(600):
            indx = findindex(jpoints, ipoints, 600, 7800, i, j)
            #must untuple the indx:
            counts[i,j] = len(indx[0]) 
    
   
    
    #ax1.plot(hthets[indx, :], color='red')
    #ax1.plot(hthets[indx2, :], color='black')
    
    #bb1 = popfig(np.sum(hthets[:,th1:th2], axis=1), jpoints, ipoints, 600, 7800)

    """
    ax1.imshow(np.log10(galfa), cmap='gray', vmin=-2, vmax=0.001)
    ax2.imshow(bb1, cmap='gray')#, vmin = -0.25, vmax = 1.5)
    
    ax1.set_xlim(2000, 2400)
    ax1.set_ylim(50, 500)
    ax2.set_xlim(2000, 2400)
    ax2.set_ylim(50, 500)
    """
    return counts

#Load in combined galfa data THIS WORKS
def loadRHT_galfa_combined_inplace(wlen, smr, thresh, ch):

    root = '/Users/susanclark/Documents/Fibers/Data/GALFA/'
    ipoints0 = np.load(root+'galfa_hi_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_0to4500_masked.npy')
    jpoints0 = np.load(root+'galfa_hj_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_0to4500_masked.npy')
    hthets0 = np.load(root+'galfa_hthets_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_0to4500_masked.npy')
    
    ipoints1 = np.load(root+'galfa_hi_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_4300to7800_masked.npy')
    jpoints1 = np.load(root+'galfa_hj_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_4300to7800_masked.npy')
    hthets1 = np.load(root+'galfa_hthets_w'+str(wlen)+'_s'+str(smr)+'_t'+str(thresh)+'_'+str(ch)+'_4300to7800_masked.npy')
    
    #Shift all second-chunk x positions to begin at 4300
    ipoints1 = ipoints1 + 4300
    
    chunk0 = np.where(ipoints0 < 4400)
    chunk1 = np.where(ipoints1 >= 4400)
    
    print len(ipoints0), len(jpoints0)
    
    
    print len(ipoints0[chunk0[0]]), ipoints0[chunk1[0]].shape
    
    ipoints = np.concatenate((ipoints0[chunk0[0]], ipoints1[chunk1[0]]))
    jpoints = np.concatenate((jpoints0[chunk0[0]], jpoints1[chunk1[0]]))
    hthets = np.concatenate((hthets0[chunk0[0]], hthets1[chunk1[0]]))
    
    return ipoints, jpoints, hthets

    
    