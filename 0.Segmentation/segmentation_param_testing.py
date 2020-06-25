import photutils
import os
import sys
import numpy as np

import astropy
from astropy.visualization import SqrtStretch
from photutils import detect_threshold
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.functional_models import Gaussian2D

from pyspark.sql.types import Row
from pyspark import SparkContext
from pyspark.sql import SparkSession

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
from photutils import deblend_sources

import time
import matplotlib.pyplot as plt

import shutil
from math import pi, log

import warnings
warnings.filterwarnings("ignore")

# Define data paths
DR2_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2/"
mosaic_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/"

hdul = fits.open(mosaic_path + "P004+38_mosaic-blanked.fits")
data = hdul[0].data


def getImage(gsigma, nsigma, npixels):
    print('Creating: n%s_g%s_np%s.png'%(nsigma, gsigma, npixels))
    # Define noise threshold
    threshold = detect_threshold(data, nsigma=nsigma)

    # Sigma?
    sigma = gsigma * gaussian_fwhm_to_sigma  # FWHM = 3.
    # Use kernel (3x3) to find borders of sources
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    # sources = object with label of sources for each pixel in the fits file
    sources = detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)

    segm_deblend = deblend_sources(data, sources, npixels=npixels,
                               filter_kernel=kernel, nlevels=32, # kernel should be same as that of thresholding
                               contrast=0.001)  

    
    # Create nice figure
    norm = ImageNormalize(stretch=SqrtStretch())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
    ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Data')
    cmap = sources.make_cmap(random_state=12345)
    ax2.imshow(sources, origin='lower', cmap=cmap)
    ax2.set_title('Segmentation Image')
    
    #print(hdul[0]._header['NAXIS1'])
    #print(hdul[0]._header['NAXIS2'])
    x = hdul[0]._header['NAXIS1']
    y = hdul[0]._header['NAXIS2']

    plt.figure(figsize=(x/1000, y/1000), dpi=100)
    plt.imshow(sources, origin='lower', cmap=cmap)

    plt.savefig('./images/n%s_g%s_np%s.png'%(nsigma, gsigma, npixels), dpi=1000)#, bbox_inches='tight', pad_inches=0)
    plt.clf()
    
    #---distribution of sources size---#
    sources_array = np.asarray(segm_deblend)
    num_sources = np.max(sources_array)

    # return the size of each label
    _, size = np.unique(segm_deblend, return_counts=True)

    # count occurances of each size
    size, counts = np.unique(size, return_counts=True)

    fig, ax = plt.subplots(1,2,figsize=(10,6))
    ax[0].scatter(size[:-1], np.log2(counts[:-1]))
    ax[1].scatter(size[:-1], np.log2(counts[:-1]))
    ax[1].set_xlim([0,100])
    plt.savefig('./distri/n%s_g%s_np%s.png'%(nsigma, gsigma, npixels))
    plt.close()



for gsigma in range(10, 11): # orignal = (2, 11); current (4, 5); later (5, 11)
    #nsigma defines threshold, larger the value higher the threshold
    for nsigma in range(2, 11): # original = (2, 11); current (6, 11); later (2, 11)
        # npixels defines the number of pixels a source at least contains
        for npixels in range(12, 17):
            #temp.append(getImage.remote(gsigma, nsigma, npixels))
            getImage(gsigma, nsigma, npixels)
