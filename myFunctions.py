# Imports

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import desispec.io
from scipy.signal import medfilt
from astropy.table import Table
from desispec.interpolation import resample_flux
from desispec.resolution import Resolution
from desispec.spectra import Spectra
from desisim.templates import BGS, ELG, LRG
import redrock.templates

# testing github on new unix

# Convert nanomaggies to magnitudes

def nmtom(data):

    nm = 22.5 - (2.5 * (np.log10(data)))
    return nm

# Plot random model spectra using desisim

def plotspec(nspec, wave, flux):

    size = (25, nspec*3) # sets size based on number of input spectra

    if nspec == 1:

        fig, ax = plt.subplots(1,1,figsize=(12,5))
        ax.plot(wave, flux)
        ax.set_ylabel('Flux')
        ax.set_xlabel('Wavelength (A)')
        ax.set_xlim(3600, 10000)

    else:

        fig, axs = plt.subplots(int(np.ceil(nspec/3)),3,figsize=size)

        ax = axs.flatten()

        for spec in range(nspec): # Plots the spectra in a subplot

            ax[spec].plot(wave, flux[spec])
            ax[spec].ylabel('Flux')
            ax[spec].xlabel('Wavelength (A)')
            ax[spec].xlim(3600, 10000)

# Write data to a fits file based on file location

def writefits(fileloc, name, data, unit, hdrname=None, overwrite=False):

    if name[-5:] != '.fits': # Checks to see if the filename is a fits file
        name += '.fits'

    infile = os.path.join(fileloc, name)

    hdr = fits.Header() # Initiates and sets up a header
    hdr['EXTNAME'] = hdrname
    hdr['BUNIT'] = unit

    if os.path.exists(infile) == True or overwrite == False: # Checks to see if the file already exists and if so appends instead of overwriting

            fits.append(infile, data, header=hdr)

    else: # If file does not exist, writes a new file

            fits.writeto(infile, data, header=hdr, overwrite=True)

# Plot spectra using quickspectra

def plotobsspectra(outfile, specnum=0, nfilter=11, plotall=False, specrange=1, truewave=None, trueflux=None, rrtitle=False, rrfile_path=None, ax=None):

    spectra = desispec.io.read_spectra(outfile) #Read the output fits spectra fill

    if rrtitle != False and rrfile_path != None: # Set up a title as redrock redshift

        zbest = Table.read(rrfile_path, 'ZBEST')
        zb = zbest['Z']

    if plotall == True: # Plots all the spectrum in the range provided, range must be less than number of spectra read in

        size = (25, specrange*2)

        if ax is None:

            fig, axs = plt.subplots(int(np.ceil(specrange/3)),3,figsize=size)
            ax = axs.flatten()

        for spec in range(specrange):

            ax[spec].axhline(0, color='k', alpha = 0.2)
            ax[spec].plot(spectra.wave['b'], medfilt(spectra.flux['b'][spec], nfilter), 'b', alpha=0.5)
            ax[spec].plot(spectra.wave['r'], medfilt(spectra.flux['r'][spec], nfilter), 'r', alpha=0.5)
            ax[spec].plot(spectra.wave['z'], medfilt(spectra.flux['z'][spec], nfilter), 'k', alpha=0.5)

            ymin = ymax = 0.0

            for x in ['b', 'r', 'z']: # Sets domain and range
                tmpmin, tmpmax = np.percentile(spectra.flux[x][spec], [1, 99])
                ymin = min(tmpmin, ymin)
                ymax = max(tmpmax, ymax)

            ax[spec].set_ylabel('Flux [1e-17 erg/s/cm2/A]')
            ax[spec].set_xlabel('Wavelength [A]')
            ax[spec].set_ylim(ymin, ymax)

            if truewave is not None and trueflux is not None: # Plots original spectra model on top of the noise
                ax[spec].plot(truewave, trueflux[spec], 'k-', alpha = 0.7)

            if rrtitle != False: # Adds the redrock redshift as a title to each graph in provided range
                ax[spec].title.set_text('z = {}'.format(zb[spec]))

    else: # Only plots the spectra number listed in specnum

        ymin = ymax = 0.0

        for x in ['b', 'r', 'z']: # Sets domain and range of plot
            tmpmin, tmpmax = np.percentile(spectra.flux[x][specnum], [1, 99])
            ymin = min(tmpmin, ymin)
            ymax = max(tmpmax, ymax)

        if ax is None:
            fig, ax = plt.subplots(figsize=(15,6))
            ax.set_ylabel('Flux [1e-17 erg/s/cm2/A]')
            ax.set_xlabel('Wavelength [A]')
            ax.set_ylim(ymin, ymax)

        ax.axhline(0, color='k', alpha = 0.2)

        ax.plot(spectra.wave['b'], medfilt(spectra.flux['b'][specnum], nfilter), 'b', alpha=0.5)
        ax.plot(spectra.wave['r'], medfilt(spectra.flux['r'][specnum], nfilter), 'r', alpha=0.5)
        ax.plot(spectra.wave['z'], medfilt(spectra.flux['z'][specnum], nfilter), 'k', alpha=0.5)

        if truewave is not None and trueflux is not None: # Plots original model spectrum
            ax.plot(truewave, trueflux[specnum], 'k-', alpha=0.7)

        if rrtitle != False: # Adds the reshift as the title
            ax.set_title('z = {}'.format(zb[specnum]))

#Given Flux Ratio and magnitude calculate second magnitude

def maggen(mag, fratio=1):

    newmag = -2.5*np.log10(fratio)+mag

    return newmag

# Simulate and combine an ELG and BGS spectra based on given seeds and red shifts

def combspec(ELGmag, BGSmag, ELGseed, BGSseed, ELGrShifts=None, BGSrShifts=None, nrShifts=None, returnmeta=False, sepflux=False):

    BGSmags = []
    ELGmags = []
    i = 0

    if ELGrShifts is not None and BGSrShifts is not None:

        while i < len(BGSrShifts):

            BGSmags.append(BGSmag)
            ELGmags.append(ELGmag)
            i += 1

        makeBGS = LRG()
        makeELG = ELG()

        fluxBGS, waveBGS, metaBGS, objmetaBGS = makeBGS.make_templates(seed=BGSseed, nmodel=len(BGSrShifts), redshift=BGSrShifts, mag=BGSmags, nocolorcuts=True)
        fluxELG, waveELG, metaELG, objmetaELG = makeELG.make_templates(seed=ELGseed, nmodel=len(ELGrShifts), redshift=ELGrShifts, mag=ELGmags, nocolorcuts=True)

        wave = waveBGS

        flux = fluxBGS + fluxELG

        if returnmeta == True and sepflux == True:
            return wave, flux, metaBGS, metaELG, fluxBGS, fluxELG

        elif returnmeta == True:
            return wave, flux, metaBGS, metaELG

        elif sepflux == True:
            return wave, flux, fluxBGS, fluxELG

        else:
            return wave, flux

    if nrShifts is not None:

        while i < nrShifts:

            BGSmags.append(BGSmag)
            ELGmags.append(ELGmag)
            i += 1

        makeBGS = BGS()
        makeELG = ELG()

        fluxBGS, waveBGS, metaBGS, objmetaBGS = makeBGS.make_templates(seed=BGSseed, nmodel=nrShifts, mag=BGSmags, nocolorcuts=True)
        fluxELG, waveELG, metaELG, objmetaELG = makeELG.make_templates(seed=ELGseed, nmodel=nrShifts, mag=ELGmags, nocolorcuts=True)

        wave = waveBGS
        flux = fluxBGS + fluxELG

        if returnmeta == True and sepflux == True:
            return wave, flux, metaBGS, metaELG, fluxBGS, fluxELG

        elif returnmeta == True:
            return wave, flux, metaBGS, metaELG

        elif sepflux == True:
            return wave, flux, fluxBGS, fluxELG

        else:
            return wave, flux

# Takes in redrock templates and returns the data for the background ELG in a spectra class form

def rrtemp_to_spectra(infile, nrshifts=None, tempfluxes=None, tempwaves=None, zbest=None):

    '''
    returns a list of files of spectra class objects setup to run through redrock

    infile: must be a quickspectra infile with keys IVAR, MASK, RESOLUTION
    fileloc: file location
    tempfluxes: must be a list of redshift template fluxes
    tempwaves: must be a list of redshift template waves
    zbest: if tempfluxes and tempwaves are not provided the zbest data must be so the templates can be found
    '''

    fileloc = os.path.dirname(infile)

    spectra = desispec.io.read_spectra(infile)

    if tempfluxes == None and tempwaves == None:

        tempfile = redrock.templates.find_templates()[0]
        rrtemp = redrock.templates.Template(tempfile, wave=spectra.wave)

        tempfluxes = []
        tempwaves = []

        for ii in range(len(zbest)):
            ncoeff = rrtemp.flux.shape[0]
            coeff = zbest['COEFF'][ii][:ncoeff]

            tempfluxes.append(rrtemp.flux.T.dot(coeff))
            tempwaves.append(rrtemp.wave * (1+zbest[ii]['Z']))

    if nrshifts == None:

        nrshifts = len(tempfluxes)

    spectra_fibermap = Table.read(infile, 'FIBERMAP')

    ivarb = fits.getdata(infile, 'B_IVAR')
    ivarr = fits.getdata(infile, 'R_IVAR')
    ivarz = fits.getdata(infile, 'Z_IVAR')

    maskb = fits.getdata(infile, 'B_MASK')
    maskr = fits.getdata(infile, 'R_MASK')
    maskz = fits.getdata(infile, 'Z_MASK')

    resb = fits.getdata(infile, 'B_RESOLUTION')
    resr = fits.getdata(infile, 'R_RESOLUTION')
    resz = fits.getdata(infile, 'Z_RESOLUTION')

    spectra_fibermaps = []
    specfiles = []
    reszbests = []

    for row in spectra_fibermap:

        spectra_fibermaps.append(Table(row))

    specdata = []

    for ii in range(nrshifts):

        ivar = {'b': np.array([ivarb[ii]]), 'r': np.array([ivarr[ii]]), 'z': np.array([ivarz[ii]])}
        mask = {'b': np.array([maskb[ii]]), 'r': np.array([maskr[ii]]), 'z': np.array([maskz[ii]])}
        res = {'b': np.array([resb[ii]]), 'r': np.array([resr[ii]]), 'z': np.array([resz[ii]])}

        netflux = {'b': None, 'r': None, 'z': None}
        bands = []

        for band in spectra.bands:

            R = Resolution(spectra.resolution_data[band][0])
            txflux = R.dot(resample_flux(spectra.wave[band], tempwaves[ii], tempfluxes[ii]))

            netflux[band] = np.array([spectra.flux[band][ii] - txflux])

        spec = Spectra(spectra.bands, spectra.wave, netflux, ivar, resolution_data=res, mask=mask,
                       fibermap=spectra_fibermaps[ii], meta=None, single=True)

        residualoutfile = os.path.join(fileloc, 'residualdata-spectra-class-{}-{}'.format(infile[len(fileloc)+1:-5], ii))
        specfile = desispec.io.write_spectra(outfile=residualoutfile, spec=spec)
        specfiles.append(specfile)

    return specfiles

# Takes the coefficients and y-interecepts of multiple linear equations and finds their intersection point(s)

def find_intersect(coeffs, intercepts):

    '''
    returns a list of points in the form of a tuple of (x, y)

    coeffs are the coefficients of x in y=mx+b
    intercepts are the y intercepts,
    if the intercepts are the x intercepts, the function will return a tuple in the form (y, x)
    '''

    numeq = len(coeffs)
    x = []
    y = []

    for i in range(numeq):

        for k in range(numeq):

            if i != k:

                netInt = intercepts[i] - intercepts[k]
                netCoeff = coeffs[k] - coeffs[i]

                x0 = netInt / netCoeff
                x.append(x0)

                y0 = (coeffs[i] * x0) + intercepts[i]
                y.append(y0)

    sorty, yindex = np.unique(y, return_index=True)
    sortx, xindex = np.unique(x, return_index=True)

    xindex = np.sort(xindex)
    yindex = np.sort(yindex)

    newx, newy = [], []
    for i in range(len(xindex)):
        newx.append(x[i])
    for i in range(len(yindex)):
        newy.append(y[i])

    points = []
    for i in range(len(newx)):
        point = (newx[i], newy[i])
        points.append(point)

    return points
