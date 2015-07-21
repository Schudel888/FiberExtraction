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
import numpy as np

import matplotlib
matplotlib.rcParams['backend'] = 'TkAgg'
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['figure.dpi'] = 250
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
# File Handling
#-----------------------------------------------------------------------------------------
def is_xyt_filename(fileobj):
    try:
        assert isinstance(fileobj, str)    
        assert fileobj.endswith('.fits')
        assert rht.xyt_suffix in fileobj
        assert SUFFIX not in fileobj
        #TODO
        return True
    except Exception:
        return False

def is_filament_filename(fileobj):
    try:
        assert isinstance(fileobj, str)    
        assert fileobj.endswith('.fits')
        assert rht.xyt_suffix in fileobj
        assert SUFFIX in fileobj
        #TODO
        return True
    except Exception:
        return False

def is_filament_hdu_list(fileobj):
    try:
        assert isinstance(fileobj, fits.HDUList)
        filaments_filename = fileobj.filename()
        assert filaments_filename.endswith('.fits')
        assert rht.xyt_suffix in filaments_filename
        assert SUFFIX in filaments_filename
        #TODO
        return True
    except Exception:
        return False

#-----------------------------------------------------------------------------------------
# Post-Processing Functions
#-----------------------------------------------------------------------------------------
def update_key(filaments_filename, key, correlation_data=None, force=False):
    #Computes the value of the dataset indicated by key (or correlation_data) for each filament
    #Uses the functions indicated by key from config
    #Saves the value to the associated header entry of each filament
    assert isinstance(key, str)
    
    def do_update(hdu_list, key, correlation_data, force):
        print 'Updating key:', key, 'in', hdu_list.filename()

        if not force and key in hdu_list[skip].header:
            print key+' keyword found in filament header...'
            if 'y' not in raw_input('Overwrite? ([no]/yes):  '):
                return hdu_list, key

        #correlation_data can be None if and only if the functions you will call expect Cloud objects
        #else, it must correspond to a source_data and applicable_methods entry
        if correlation_data is None: 
            if key in config.source_data:
                correlation_data = config.source_data[key]
            else:
                raise KeyError('No source_data for key: '+key+' in config.source_data')
            
        if key not in config.applicable_methods:
            raise KeyError('No applicable_methods for key: '+key+' in config.applicable_methods')

        #for i, hdu in enumerate(hdu_list[skip:]):
        for hdu in hdu_list[skip:]:
            hdr = hdu.header
            if correlation_data is not None:
                #Assumes correlation_data can be indexed into using ndarray notation [LowerLeft to UpperRight+1]
                #Assumes config.Cloud.nonzero_data_from_HDU will return the pixel coords offset properly for the above masked region
                #Assumes func can take this weird ndarray view as input and return a scalar value
                for suffix in config.applicable_methods[key]:        
                    func = config.methods[suffix]
                    hdr[key+suffix] = func(correlation_data[hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1][config.Cloud.nonzero_data_from_HDU(hdu, transpose=True)])
                hdr[key] = config.timestamp()
            else:
                #Assumes all func require a config.Cloud object to work
                tempCloud = config.Cloud(hdu)
                for suffix in config.applicable_methods[key]:        
                    func = config.Cloud.functions[suffix][0]
                    hdr[key+suffix] = func(tempCloud)
            
        hdu_list.flush()

    if is_filament_filename(filaments_filename):
        #with config.default_open(filaments_filename, mode='update') as hdu_list:
        hdu_list = config.default_open(filaments_filename, mode='update')
        do_update(hdu_list, key, correlation_data, force)
        hdu_list.close()
        return filaments_filename

    elif is_filament_hdu_list(filaments_filename):
        hdu_list = filaments_filename
        assert hdu_list.fileinfo(0)['filemode'] == 'update'
        do_update(hdu_list, key, correlation_data, force)
        return hdu_list, key

    else:
        print 'Unknown input encountered in update_key()...'
        return None

def update_all_keys(filaments_filename, force=False):
    #Updates the properties corresponding to all known sources

    def do_update_all(hdu_list, force):
        print 'Updating all keys in', hdu_list.filename()

        if not force:
            print 'This is a long operation that involves updating ALL keys from config.sources'
            if 'y' not in raw_input('Continue? ([no]/yes):  '):
                return 'Aborted: update_all_keys()'

        exceptions = ''
        for key in config.sources.iterkeys():
            try:
                update_key(hdu_list, key, force=force) #TODO Don't need to do anything with the output?
            except Exception as e:
                print e
                if len(exceptions) == 0:
                    exceptions += 'except '+key
                else:
                    exceptions += ', '+key
        if len(exceptions) > 0 and string.count(exceptions, ',') > 0:
            exceptions = string.join(string.rsplit(exceptions, ',', 1), ' and')
        return 'All properties '+exceptions+'are now up to date in '+hdu_list.filename()
        
    if is_filament_filename(filaments_filename):
        #with config.default_open(filaments_filename, mode='update') as hdu_list:
        hdu_list = config.default_open(filaments_filename, mode='update')
        print do_update_all(hdu_list, force)
        return filaments_filename

    elif is_filament_hdu_list(filaments_filename):
        hdu_list = filaments_filename
        assert hdu_list.fileinfo(0)['filemode'] == 'update'
        print do_update_all(hdu_list, force)
        return hdu_list, key

    else:
        print 'Unknown input encountered in update_all_key()...'
        return None


def plot(filaments_filename, key=None, out_name=None, show=True, cut=config.passive_constant(True)):
    if not show and out_name is None:
        'Unable to plot data without showing or saving it in isolate.plot()'
        return None

    def do_plot(hdu_list, key, out_name, show, cut):
        try:
            filaments_filename = hdu_list.filename()
            print 'Plotting key:', key, 'in', filaments_filename
            backproj = hdu_list.pop(0)
            #TODO make sure to pop skip hdus!s
            hdu_list = filter(cut, hdu_list) #TODO TURNS AN HDULIST INTO A LIST ~_~

            if key is None:
                figure_args = {'figsize':(8,6), 'facecolor':'white','dpi':250}
                fig = plt.figure(**figure_args)
                #display = np.zeros((backproj.header['NAXIS1'], backproj.header['NAXIS2']))
                display = np.empty((backproj.header['NAXIS1'], backproj.header['NAXIS2'])).fill(np.nan)
                N = len(hdu_list)
                for i, hdu in enumerate(hdu_list):
                    hdr = hdu.header            
                    display[hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1][config.Cloud.nonzero_data_from_HDU(hdu, transpose=False)] = N-i
                ax1 = fig.add_axes([0.05,0.05,0.90,0.90])
                ax1.imshow(display.T)
                key = 'Filaments'

            elif isinstance(key, str) and key in hdu_list[0].header: #skip].header:
                displays = dict()
                datasets = dict()
                titles = [key+suffix for suffix in config.applicable_methods[key]] 

                #for suffix in suffixes:    
                for title in titles:
                    displays[title] = np.zeros_like(backproj.data)
                    displays[title].fill(np.nan)
                    datasets[title] = config.post_processing[key]([hdu.header[title] for hdu in hdu_list]) #TODO watch out for python list handling

                #for hdu in hdu_list:
                for i, hdu in enumerate(hdu_list):      
                    hdr = hdu.header
                    #for suffix in suffixes:
                    for title in titles:
                        try:
                            #mask = displays[key+suffix][hdr['MIN_X']:hdr['MAX_X']+1, hdr['MIN_Y']:hdr['MAX_Y']+1]
                            #displays[key+suffix][hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1][config.Cloud.nonzero_data_from_HDU(hdu, transpose=True)].fill(datasets[key+suffix][i])
                            displays[title][hdr['MIN_Y']:hdr['MAX_Y']+1, hdr['MIN_X']:hdr['MAX_X']+1][config.Cloud.nonzero_data_from_HDU(hdu, transpose=True)] = datasets[title][i] #.fill(datasets[title][i]) 
                        except Exception as e:
                            print e
                            print 'Failed to plot data from (1-indexed) HDU:', i+1+skip, title

                Nplots = len(titles)
                figure_args = {'figsize':(8,2.5*Nplots-0.5), 'facecolor':'white','dpi':250}
                fig = plt.figure(**figure_args)
                r=4
                gs = gridspec.GridSpec(Nplots*r, 5*r)
                for i, title in enumerate(titles):
                    ax1 = plt.subplot(gs[1+r*i:r*(i+1),0:-4]) 
                    ax2 = plt.subplot(gs[1+r*i:r*(i+1)-1,-3:-1])
                    ax2.hist(datasets[title], bins=50, orientation='horizontal', histtype='stepfilled')
                    plt.colorbar(ax1.imshow(displays[title], cmap = "YlOrRd"), ax=ax2, fraction=0.10)
                    plt.title(title, fontsize=8)
            
            else:
                print 'Unable to plot data using the given key: '+str(key)
                return
        
            plt.suptitle(key+' from '+filaments_filename, fontsize=12)
            if out_name is not None and isinstance(out_name, str):
                plt.savefig(out_name, dpi=500, format='png')
            if show:
                plt.show()

        except Exception as error:
            print error

        finally:
            plt.cla()
            plt.clf() 
            plt.close()

    if is_filament_filename(filaments_filename):
        #with config.default_open(filaments_filename, mode='readonly') as hdu_list:
        hdu_list = config.default_open(filaments_filename, mode='readonly')
        do_plot(hdu_list, key=key, out_name=out_name, show=show, cut=cut)
        hdu_list.close()
        return filaments_filename

    elif is_filament_hdu_list(filaments_filename):
        hdu_list = filaments_filename
        do_plot(hdu_list, key=key, out_name=out_name, show=show, cut=cut)
        return hdu_list, key

    else:
        print 'Unknown input encountered in plot()...'
        return None

#-----------------------------------------------------------------------------------------
# Bulk Fiber Isolation Functions
#-----------------------------------------------------------------------------------------

def isolate_all(xyt_filename, BINS=6, generateHDU=config.Cloud.as_ImageHDU): #config.Cloud.as_BinTableHDU):

    assert is_xyt_filename(xyt_filename)
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
                    j = point_dict[config.rel_add(coord, rel_coord)]
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
            list_of_HDUs.append(config.Cloud(out_clouds.popitem()[1]).as_HDU()) #TODO Incorporate theta_dict
            #list_of_HDUs.append(generateHDU(config.Cloud(out_clouds.popitem()[1])))
        #rht.update_progress(1.0, final_message='Finished joining '+str(problem_size)+' points! Time Elapsed:')
        
    #Convert lists of two-integer tuples into ImageHDUs
    #unprocessed.sort(key=len, reverse=True)
    #output_hdulist = fits.HDUList(map(config.Cloud.as_ImageHDU, map(config.Cloud, unprocessed)))
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
            #plot(*update_key(filename, key='COLDENS'), out_name=filename[:-5]+'_COLDENS.png')
            #plot(*update_key(filename, key='GALFA0'), out_name=filename[:-5]+'_GALFA0.png')
            #plot(*update_key(filename, key='B'))
            #plot(filename, key='GALFA0', out_name=filename[:-5]+'_GALFA0_50.png', cut=lambda h: h.header['B_MIN'] > 50.0)
            #plot(*update_key(filename, key='', force=True))
            update_all_keys(filename)

    #Cleanup and Exit
    #exit()




			