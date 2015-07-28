def radecs_to_xy_GALFA(ras, decs, region = "SC_241"):

    #Set wcs transformation
    w = set_wcs_galfa(region = region)
    
    #Transformation
    radec = zip(ras, decs)
    xy = w.wcs_world2pix(radec, 1)
    
    xs = xy[:,0]
    ys = xy[:,1]
    
    # Check if we're in the past-zero region of t2s4, and shift center pixel if we are.
    if region == "t2s4":
        ramin, ramax, decmin, decmax = getcorners(region = region, verbose = False)
        #if ((ra < ramin) and (ra >= 0)) or (ra > 360):
        if ((ra > ramax) and (ra <= 360)):
            #print "shifting"
            crpix = w.wcs.crpix
            crpix[0] = 360*60+crpix[0]
            w.wcs.crpix = crpix
            
            xy = w.wcs_world2pix(radec, 1)
    
            x = xy[0,0]
            y = xy[0,1]
    
    return xs, ys

def radec_to_xy_GALFA(ra, dec, region = "SC_241"):

    #Set wcs transformation
    w = set_wcs_galfa(region = region)
    
    #Transformation
    radec = [[ra, dec]]
    xy = w.wcs_world2pix(radec, 1)
    
    x = xy[0,0]
    y = xy[0,1]
    
    # Check if we're in the past-zero region of t2s4, and shift center pixel if we are.
    if region == "t2s4":
        ramin, ramax, decmin, decmax = getcorners(region = region, verbose = False)
        #if ((ra < ramin) and (ra >= 0)) or (ra > 360):
        if ((ra > ramax) and (ra <= 360)):
            #print "shifting"
            crpix = w.wcs.crpix
            crpix[0] = 360*60+crpix[0]
            w.wcs.crpix = crpix
            
            xy = w.wcs_world2pix(radec, 1)
    
            x = xy[0,0]
            y = xy[0,1]
    
    return x, y
    
def xys_to_radec_GALFA(xs, ys, region = "SC_241"):

    #Set wcs transformation
    w = set_wcs_galfa(region = region)
    
    #Transformation
    xy = zip(xs, ys)
    radec = w.wcs_pix2world(xy, 1)
    
    ras = radec[:,0]
    decs = radec[:,1]
    
    """
    # If nan's exist, shift center pixel.
    if (np.isnan(ra)) or (np.isnan(dec)):
        #print "shifting"
        crpix = w.wcs.crpix
        crpix[0] = 360*60+crpix[0]
        w.wcs.crpix = crpix
        
        #Transformation
        xy = [[x, y]]
        radec = w.wcs_pix2world(xy, 1)
    
        ra = radec[0,0]
        dec = radec[0,1]
    """
    
    return ras, decs
    

def xy_to_radec_GALFA(x, y, region = "SC_241"):

    #Set wcs transformation
    w = set_wcs_galfa(region = region)
    
    #Transformation
    xy = [[x, y]]
    radec = w.wcs_pix2world(xy, 1)
    
    ra = radec[0,0]
    dec = radec[0,1]
    
    # If nan's exist, shift center pixel.
    if (np.isnan(ra)) or (np.isnan(dec)):
        #print "shifting"
        crpix = w.wcs.crpix
        crpix[0] = 360*60+crpix[0]
        w.wcs.crpix = crpix
        
        #Transformation
        xy = [[x, y]]
        radec = w.wcs_pix2world(xy, 1)
    
        ra = radec[0,0]
        dec = radec[0,1]
    
    return ra, dec
    
def radec_to_lb(ra, dec):

    #Transformation. Conforms to astropy 0.4.3     
    obj = coord.SkyCoord(ra=np.float(ra)*u.degree, dec=np.float(dec)*u.degree, frame = "icrs")
    obj = obj.galactic
    
    l = obj.l.degree
    b = obj.b.degree
    
    return l, b 

def radecs_to_lb(ras, decs):

    #Transformation. Conforms to astropy 0.4.3     
    obj = coord.SkyCoord(ras, decs, unit = "deg", frame = "icrs")
    obj = obj.galactic
    
    ls = obj.l.degree
    bs = obj.b.degree
    
    return ls, bs 
    
def lb_to_radec(l, b):

    #Transformation. Conforms to astropy 0.4.3     
    obj = coord.SkyCoord(l=np.float(l)*u.degree, b=np.float(b)*u.degree, frame = "galactic")
    obj = obj.icrs
    
    ra = obj.ra.degree
    dec = obj.dec.degree
    
    return ra, dec  

def set_wcs_galfa(region = "SC_241"):

    if region == "t2n3":
        fn = "/Users/susanclark/ACTPol/t2n3.shallowcube.5km_bins.v6_spcor8_12_16avg_with_av_db.fits"
        hdr = fits.getheader(fn)
    if region == "t2s4":
        fn = "/Users/susanclark/ACTPol/t2s4.shallowcube.5km_bins.v6_spcor8_12_16fixed_with_av_db_x_cleanx_clean.fits"
        hdr = fits.getheader(fn)
    if region == "SC_241":    
        fn = "/Volumes/DataDavy/GALFA/SC_241/cleaned/SC_241.66_28.675.best.fits"
        hdr = fits.getheader(fn)
    
    """Set wcs transformation"""
    w = wcs.WCS(fn, naxis=2)
    #w.naxis = 2
    #w = wcs.WCS(naxis=2)
    #w.wcs.crpix = [hdr["CRPIX1"], hdr["CRPIX2"]]
    #w.wcs.cdelt = np.array([hdr["CDELT1"], hdr["CDELT2"]]) #pixel size in wcu
    #w.wcs.crval = [hdr["CRVAL1"],hdr["CRVAL2"]]
    #w.wcs.ctype = [hdr["CTYPE1"], hdr["CTYPE2"]]
    
    return w