"""
mosfire_debugging
=================

Some troubleshooting for when MOSFIRE DRP breaks
"""

import sys, os

from chun_codes import systime

from os.path import exists

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

from MOSFIRE import Options, Wavelength, IO, CSU, Fit, Detector

__version__ = 0.1

#from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength, Extract

flatops = Options.flat
waveops = Options.wavelength

def make_pixel_flat(data, results, options, outfile, inputs, lampsOff=None):
    '''
    Convert a flat image into a flat field
    '''

    def pixel_min(y): return int(np.floor(np.min(y)))
    def pixel_max(y): return int(np.ceil(np.max(y)))

    def collapse_flat_box(dat):
        '''Collapse data to the spectral axis (0)'''
        v = np.median(dat, axis=0).ravel()

        return v

    flat = np.ones(shape=Detector.npix)

    hdu = fits.PrimaryHDU((data/flat).astype(np.float32))
    hdu.header.set("version", __version__, "DRP version")
    i = 0
    for flatname in inputs:
        nm = flatname.split("/")[-1]
        hdu.header.set("infile%2.2i" % i, nm)
        i += 1

    slitno = 0
    for result in results[0:-1]:
        slitno += 1

        hdu.header.set("targ%2.2i" % slitno, result["Target_Name"])

        bf = result["bottom"]
        tf = result["top"]
        try:
            hpps = result["hpps"]
        except:
            error( "No half power points for this slit")
            hpps = [0, Detector.npix[0]]

        xs = np.arange(hpps[0], hpps[1])

        top = pixel_min(tf(xs))
        bottom = pixel_max(bf(xs))

        hdu.header.set("top%2.2i" % slitno, top)
        hdu.header.set("bottom%2.2i" % slitno, bottom)

        log.info( "%s] Bounding top/bottom: %i/%i" % (result["Target_Name"],
                                                      bottom, top))

        v = collapse_flat_box(data[bottom:top,hpps[0]:hpps[1]])

        x2048 = np.arange(Options.npix)
        v = np.poly1d(np.polyfit(xs,v,
            options['flat-field-order']))(xs).ravel()

        for i in np.arange(bottom-1, top+1):
            flat[i,hpps[0]:hpps[1]] = v

    log.info("Producing Pixel Flat...")
    for r in range(len(results)-1):
        theslit = results[r]

        try:
            bf = theslit["bottom"]
            tf = theslit["top"]
        except:
            pdb.set_trace()

        for i in range(hpps[0], hpps[1]):
            top = int(np.floor(tf(i)))
            bottom = int(np.ceil(bf(i)))

            data[top:bottom, i] = flat[top:bottom,i]

    hdu.data = (data/flat).astype(np.float32)
    bad = np.abs(hdu.data-1.0) > 0.5
    hdu.data[bad] = 1.0
    hdu.data = hdu.data.filled(1)
    if os.path.exists(outfile):
            os.remove(outfile)
    hdu.writeto(outfile)
    log.info("Done.")

def fit_edge_poly(xposs, xposs_missing, yposs, order):
    '''
    fit_edge_poly fits a polynomial to the measured slit edges.
    This polynomial is used to extract spectra.

    fit_edge_poly computes a parabola, and fills in missing data with a
    parabola

    input-
    xposs, yposs [N]: The x and y positions of the slit edge [pix]
    order: the polynomial order
    '''

    # First fit low order polynomial to fill in missing data
    fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, 2))

    xposs = np.append(xposs, xposs_missing)
    yposs = np.append(yposs, fun(xposs_missing))

    # Remove any fits that deviate wildly from the 2nd order polynomial
    ok = np.abs(yposs - fun(xposs)) < 1
    if not ok.any():
            error("Flat is not well illuminated? Cannot find edges")
            raise Exception("Flat is not well illuminated? Cannot find edges")

    # Now refit to user requested order
    fun = np.poly1d(Fit.polyfit_clip(xposs[ok], yposs[ok], order))
    res = fun(xposs[ok]) - yposs[ok]
    sd = np.std(res)
    ok = np.abs(res) < 2*sd


    # Check to see if the slit edge funciton is sane,
    # if it's not, then we fix it.
    pix = np.arange(2048)
    V = fun(pix)
    if np.abs(V.max() - V.min()) > 10:
        log.info ("Forcing a horizontal slit edge")
        fun = np.poly1d(np.median(yposs[ok]))


    return (fun, res, sd, ok)

def find_edge_pair(data, y, roi_width, edgeThreshold=450):
    '''
    find_edge_pair finds the edge of a slit pair in a flat

    data[2048x2048]: a well illuminated flat field [DN]
    y: guess of slit edge position [pix]

    Keywords:

    edgeThreshold: the pixel value below which we should ignore using
    to calculate edges.

    Moves along the edge of a slit image
            - At each location along the slit edge, determines
            the position of the demarcations between two slits

    Outputs:
    xposs []: Array of x positions along the slit edge [pix]
    yposs []: The fitted y positions of the "top" edge of the slit [pix]
    widths []: The fitted delta from the top edge of the bottom [pix]
    scatters []: The amount of light between slits


    The procedure is as follows
    1: starting from a guess spatial position (parameter y), march
        along the spectral direction in some chunk of pixels
    2: At each spectral location, construct a cross cut across the
        spatial direction; select_roi is used for this.
    3: Fit a two-sided error function Fit.residual_disjoint_pair
        on the vertical cross cut derived in step 2.
    4: If the fit fails, store it in the missing list
        - else if the top fit is good, store the top values in top vector
        - else if the bottom fit is good, store the bottom values in bottom
          vector.
    5: In the vertical cross-cut, there is a minimum value. This minimum
        value is stored as a measure of scattered light.

    Another procedure is used to fit polynomials to these fitted values.
    '''

    def select_roi(data, roi_width):
        v = data[y-roi_width:y+roi_width, xp-2:xp+2]
        v = np.median(v, axis=1) # Axis = 1 is spatial direction
        return v



    xposs_top = []
    yposs_top = []
    xposs_top_missing = []

    xposs_bot = []
    yposs_bot = []
    xposs_bot_missing = []
    yposs_bot_scatters = []

    #1
    rng = np.linspace(10, 2040, 50).astype(np.int)
    for i in rng:
        xp = i
        #2
        v = select_roi(data, roi_width)
        xs = np.arange(len(v))

        # Modified from 450 as the hard coded threshold to one that
        # can be controlled by a keyword
        if (np.median(v) < edgeThreshold):
            log.warn('median < edgeThreshold : %i' % i)
            xposs_top_missing.append(xp)
            xposs_bot_missing.append(xp)
            continue

        #3
        ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
        fit_ok = 0 < ff[4] < 4

        if fit_ok:
            (sigma, offset, mult1, mult2, add, width) = ff[0]

            xposs_top.append(xp)
            yposs_top.append(y - roi_width + offset + width)

            xposs_bot.append(xp)
            yposs_bot.append(y - roi_width + offset)

            between = offset + width/2
            if 0 < between < len(v)-1:
                start = np.max([0, between-2])
                stop = np.min([len(v),between+2])
                yposs_bot_scatters.append(np.min(v[start:stop])) # 5

                if False:
                    pl.figure(2)
                    pl.clf()
                    tmppix = np.arange(y-roi_width, y+roi_width)
                    tmpx = np.arange(len(v))
                    pl.axvline(y - roi_width + offset + width, color='red')
                    pl.axvline(y - roi_width + offset, color='red')
                    pl.scatter(tmppix, v)
                    pl.plot(tmppix, Fit.fit_disjoint_pair(ff[0], tmpx))
                    pl.axhline(yposs_bot_scatters[-1])
                    pl.draw()

            else:
                yposs_bot_scatters.append(np.nan)

        else:
            xposs_bot_missing.append(xp)
            xposs_top_missing.append(xp)
            log.info("Skipping wavelength pixel): %i" % (xp))


    return map(np.array, (xposs_bot, xposs_bot_missing, yposs_bot, xposs_top,
        xposs_top_missing, yposs_top, yposs_bot_scatters))

def find_and_fit_edges(data, header, bs, options,edgeThreshold=450):
    '''
    Given a flat field image, find_and_fit_edges determines the position
    of all slits.

    The function works by starting with a guess at the location for a slit
    edge in the spatial direction(options["first-slit-edge"]).

    Starting from the guess, find_edge_pair works out in either direction,
    measuring the position of the (e.g.) bottom of slit 1 and top of slit 2:


    ------ pixel y value = 2048

    Slit 1 data

    ------ (bottom)
    deadband
    ------ (top)

    Slit N pixel data ....

    ------- (bottom) pixel = 0

    --------------------------------> Spectral direction


    1. At the top of the flat, the slit edge is defined to be a pixel value
    2. The code guesses the position of the bottom of the slit, and runs
            find_edge_pair to measure slit edge locations.
    3. A low-order polynomial is fit to the edge locations with
            fit_edge_poly
    4. The top and bottom of the current slit, is stored into the
            result list.
    5. The top of the next slit is stored temporarily for the next
            iteration of the for loop.
    6. At the bottom of the flat, the slit edge is defined to be pixel 4.


    options:
    options["edge-order"] -- The order of the polynomial [pixels] edge.
    options["edge-fit-width"] -- The length [pixels] of the edge to
            fit over

    '''

    # TODO: move hardcoded values into Options.py
    # y is the location to start
    y = 2034
    DY = 44.25

    toc = 0
    ssl = bs.ssl

    slits = []

    top = [0., np.float(Options.npix)]

    start_slit_num = int(bs.msl[0]['Slit_Number'])-1
    if start_slit_num > 0:
        y -= DY * start_slit_num

    # Count and check that the # of objects in the SSL matches that of the MSL
    # This is purely a safety check
    numslits = np.zeros(len(ssl))
    for i in xrange(len(ssl)):
        slit = ssl[i]
        M = np.where(slit["Target_Name"] == bs.msl["Target_in_Slit"])

        numslits[i] = len(M[0])
    numslits = np.array(numslits)


    if (np.sum(numslits) != CSU.numslits) and (not bs.long_slit) and (not bs.long2pos_slit):
        error ("The number of allocated CSU slits (%i) does not match "
                " the number of possible slits (%i)." % (np.sum(numslits),
                    CSU.numslits))
        raise Exception("The number of allocated CSU slits (%i) does not match "
                " the number of possible slits (%i)." % (np.sum(numslits),
                    CSU.numslits))

    # if the mask is a long slit, the default y value will be wrong. Set instead to be the middle
    if bs.long_slit:
        y = 1104

    # now begin steps outline above
    results = []
    result = {}

    result["Target_Name"] = ssl[0]["Target_Name"]

    # 1
    result["top"] = np.poly1d([y])

    ''' Nomenclature here is confusing:

        ----- Edge  -- Top of current slit, bottom of prev slit
        . o ' Data
        ===== Data
        .;.;' Data
        ----- Edge  -- Bottom of current slit, top of next slit
    '''

    topfun = np.poly1d([y])
    xposs_top_this = np.arange(10,2000,100)
    yposs_top_this = topfun(xposs_top_this)

    initial_edges = np.array([2034])
    edge = 2034

    # build an array of values containing the lower edge of the slits

    for target in xrange(len(ssl)):
    # target is the slit number
        edge -= DY * numslits[target]
        initial_edges=np.append(initial_edges,edge)

    # collapse the 2d flat along the walenegth axis to build a spatial profile of the slits
    vertical_profile = np.mean(data, axis=1)

    # build an array containing the spatial positions of the slit centers, basically the mid point between the expected
    # top and bottom values of the slit pixels
    spatial_centers = np.array([])
    for k in np.arange(0,len(initial_edges)-1):
        spatial_centers = np.append(spatial_centers,(initial_edges[k]+initial_edges[k+1])/2)
    #slit_values=np.array([])
    #for k in np.arange(0, len(spatial_centers)):
    #    slit_values = np.append(slit_values,np.mean(vertical_profile[spatial_centers[k]-3:spatial_centers[k]+3]))

    for target in xrange(len(ssl)):

        y -= DY * numslits[target]
        y = max(y, 1)
        # select a 6 pixel wide section of the vertical profile around the slit center
        threshold_area = vertical_profile[spatial_centers[target]-3:spatial_centers[target]+3]
        log.info('threshold area : %f' % np.mean(threshold_area))
        # uses 80% of the ADU counts in the threshold area to estimate the threshold to use in defining the slits
        edgeThreshold = np.mean(threshold_area)*0.8
        #if edgeThreshold > 450:
        #    edgeThreshold = 450

        log.info("[%2.2i] Finding Slit Edges for %s ending at %4.0i. Slit "
                 "composed of %i CSU slits" % ( target,
                    ssl[target]["Target_Name"], y, numslits[target]))
        log.info("[%2.2i] Threshold used is %.1f" % (target,edgeThreshold))

        ''' First deal with the current slit '''
        hpps = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(target+1)[0], header, bs)
        ''' This might be where the bug is for J2 -- CL & SM '''
        log.info("y = %i | hpps = %.1f  %.1f " % (y, hpps[0], hpps[1]))
        if y == 1:
            xposs_bot = [1024]
            xposs_bot_missing = []
            yposs_bot = [4.25]
            botfun = np.poly1d(yposs_bot)
            ok = np.where((xposs_bot > hpps[0]) & (xposs_bot < hpps[1]))
        else:
            (xposs_top_next, xposs_top_next_missing, yposs_top_next, xposs_bot,
                xposs_bot_missing, yposs_bot, scatter_bot_this) = find_edge_pair(
                    data, y, options["edge-fit-width"],edgeThreshold=edgeThreshold)

            print xposs_top_next, xposs_top_next_missing, yposs_top_next
            print xposs_bot, xposs_bot_missing, yposs_bot, scatter_bot_this

            ok = np.where((xposs_bot > hpps[0]) & (xposs_bot < hpps[1]))
            ok2 = np.where((xposs_bot_missing > hpps[0]) & (xposs_bot_missing <
                hpps[1]))
            xposs_bot = xposs_bot[ok]
            xposs_bot_missing = xposs_bot_missing[ok2]
            yposs_bot = yposs_bot[ok]
            if len(xposs_bot) == 0:
                botfun = np.poly1d(y-DY)
            else:
                (botfun, bot_res, botsd, botok) =  fit_edge_poly(xposs_bot,
                         xposs_bot_missing, yposs_bot, options["edge-order"])


        bot = botfun.c.copy()
        top = topfun.c.copy()

        #4
        result = {}
        result["Target_Name"] = ssl[target]["Target_Name"]
        result["xposs_top"] = xposs_top_this
        result["yposs_top"] = yposs_top_this
        result["xposs_bot"] = xposs_bot
        result["yposs_bot"] = yposs_bot
        result["top"] = np.poly1d(top)
        result["bottom"] = np.poly1d(bot)
        result["hpps"] = hpps
        result["ok"] = ok
        results.append(result)

        #5
        if y == 1:
            break


        next = target + 2
        if next > len(ssl): next = len(ssl)
        hpps_next = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(next)[0],
                    header, bs)

        ok = np.where((xposs_top_next > hpps_next[0]) & (xposs_top_next <
            hpps_next[1]))
        ok2 = np.where((xposs_top_next_missing > hpps_next[0]) &
            (xposs_top_next_missing < hpps_next[1]))

        xposs_top_next = xposs_top_next[ok]
        xposs_top_next_missing = xposs_top_next_missing[ok2]
        yposs_top_next = yposs_top_next[ok]

        if len(xposs_top_next) == 0:
            topfun = np.poly1d(y)
        else:
            (topfun, topres, topsd, ok) = fit_edge_poly(xposs_top_next,
                xposs_top_next_missing, yposs_top_next, options["edge-order"])

        xposs_top_this = xposs_top_next
        xposs_top_this_missing = xposs_top_next_missing
        yposs_top_this = yposs_top_next

    results.append({"version": options["version"]})

    return results

def main(path, maskname, band, silent=False, verbose=True):

    '''
    Main function for testing

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 18 June 2018

    Modified by Chun Ly, 19 June 2018
     - Create Table summarizing results of find_and_fit_edges()
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    flat_file = path+'combflat_2d_'+band+'.fits'

    flatlist = path+'Flat.txt'
    flatlist = IO.list_file_to_strings(flatlist)

    bpos = np.ones(92) * -1

    for fname in flatlist:
        hdr, dat, bs = IO.readmosfits(fname, flatops)
        try: bs0
        except: bs0 = bs

        if np.any(bs0.pos != bs.pos):
            print "bs0: "+str(bs0.pos)+" bs: "+str(bs.pos)
            log.warn("Barset do not seem to match")

        if hdr["filter"] != band:
            log.warn("Filter name %s does not match header filter name "
                     "%s in file %s" % (band, hdr["filter"], fname))

        for i in xrange(len(bpos)):
            b = hdr["B{0:02d}POS".format(i+1)]
            if bpos[i] == -1:
                bpos[i] = b
            else:
                if bpos[i] != b:
                    log.warn("Bar positions are not all the same in "
                             "this set of flat files")

    bs = bs0

    (header, data) = IO.readfits(flat_file, use_bpm=True)
    log.info('Finding slit edges in {}'.format(path))

    results = find_and_fit_edges(data, header, bs, flatops,
                                 edgeThreshold=450.0)

    n_targets = len(results)-1
    Name      = [''] * n_targets
    yposs_top = np.zeros(n_targets)
    yposs_bot = np.zeros(n_targets)
    hpps      = np.zeros((2,n_targets))
    top       = np.zeros(n_targets)
    bottom    = np.zeros(n_targets)

    for ii in range(n_targets):
        tmp = results[ii]
        Name[ii] = tmp['Target_Name']
        if len(tmp['yposs_top']) > 0:
            yposs_top[ii] = np.max(tmp['yposs_top'])
        if len(tmp['yposs_bot']) > 0:
            yposs_bot[ii] = np.min(tmp['yposs_bot'])
        hpps[:,ii] = results[ii]['hpps']
        #top[ii]    = np.max(tmp['top'])
        #bottom[ii] = np.min(tmp['bottom'])

    arr0  = [Name, yposs_top, yposs_bot, hpps[0,:], hpps[1,:]]
    name0 = ('Name', 'yposs_top', 'yposs_bot', 'hpps0', 'hpps1')
    tab0  = Table(arr0, names=name0)
    tab0.pprint(max_lines=-1, max_width=-1)

    out = path+"pixelflat_2d_%s_CLtest.fits" % (band)
    make_pixel_flat(data, results, flatops, out, flatlist, lampsOff=True)

    if silent == False: log.info('### End main : '+systime())

    return results
#enddef

