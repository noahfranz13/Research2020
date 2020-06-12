import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import desispec.io
from scipy.signal import medfilt

# Convert nanomaggies to magnitudes

def nmtom(data):

    nm = 22.5 - (2.5 * (np.log10(data)))
    return nm

# Plot spectra using desisim

def plotspec(nspec, wave, flux):
    
    size = (25, nspec*3)
        
    plt.figure(figsize=size)
        
    for spec in range(nspec):
        
        plt.subplot(np.ceil(nspec/2),3,spec+1)
        plt.plot(wave, flux[spec])
        plt.ylabel('Flux')
        plt.xlabel('Wavelength (A)')
        
# Write data to a fits file based on file location

def writefits(fileloc, name, data, unit, hdrname=None):
    
    if name[-5:] != '.fits':
        name += '.fits'
        
    infile = os.path.join(fileloc, name)
    
    hdr = fits.Header()
    hdr['EXTNAME'] = hdrname
    hdr['BUNIT'] = unit
    
    if os.path.exists(infile) == True:
        
        fits.append(infile, data, header=hdr)
    
    else:
        
        fits.writeto(infile, data, header=hdr, overwrite=True)
        
#Plot spectra using quickspectra

def plotobsspectra(outfile, specnum=0, nfilter=11, plotall=False, specrange=1, truewave=None, trueflux=None):
    
    spectra = desispec.io.read_spectra(outfile)
    
    if plotall == True:
        
        size = (25, specrange*3)
        
        plt.figure(figsize=size)
        
        for spec in range(specrange):
            
            plt.subplot(np.ceil(specrange/2),3,spec+1)
            plt.axhline(0, color='k', alpha = 0.2)
            plt.plot(spectra.wave['b'], medfilt(spectra.flux['b'][spec], nfilter), 'b', alpha=0.5)
            plt.plot(spectra.wave['r'], medfilt(spectra.flux['r'][spec], nfilter), 'r', alpha=0.5)
            plt.plot(spectra.wave['z'], medfilt(spectra.flux['z'][spec], nfilter), 'k', alpha=0.5)
            
            ymin = ymax = 0.0
            
            for x in ['b', 'r', 'z']:
                tmpmin, tmpmax = np.percentile(spectra.flux[x][spec], [1, 99])
                ymin = min(tmpmin, ymin)
                ymax = max(tmpmax, ymax)
            
            
            plt.ylabel('Flux [1e-17 erg/s/cm2/A]')
            plt.xlabel('Wavelength [A]')
            plt.ylim(ymin, ymax)
            
            if truewave is not None and trueflux is not None:
                plt.plot(truewave, trueflux[spec], 'k-')
            
            
    else:
        
        plt.axhline(0, color='k', alpha = 0.2)
        
        plt.plot(spectra.wave['b'], medfilt(spectra.flux['b'][specnum], nfilter), 'b', alpha=0.5)
        plt.plot(spectra.wave['r'], medfilt(spectra.flux['r'][specnum], nfilter), 'r', alpha=0.5)
        plt.plot(spectra.wave['z'], medfilt(spectra.flux['z'][specnum], nfilter), 'k', alpha=0.5)
        
        ymin = ymax = 0.0
            
        for x in ['b', 'r', 'z']:
            tmpmin, tmpmax = np.percentile(spectra.flux[x][specnum], [1, 99])
            ymin = min(tmpmin, ymin)
            ymax = max(tmpmax, ymax)
        
        plt.ylabel('Flux [1e-17 erg/s/cm2/A]')
        plt.xlabel('Wavelength [A]')
        plt.ylim(ymin, ymax)
        
        if truewave is not None and trueflux is not None:
            plt.plot(truewave, trueflux[specnum], 'k-')
