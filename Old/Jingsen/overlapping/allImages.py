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

# from pyspark.sql.types import Row
# from pyspark import SparkContext
# from pyspark.sql import SparkSession

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
from photutils import deblend_sources

import time
import matplotlib.pyplot as plt

import shutil
from math import pi, log
from tqdm import tqdm
from PIL import Image

import warnings
warnings.filterwarnings("ignore")

# Define data paths
DR2_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2/"
mosaic_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/"
#write_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_writable/dr2_images"
write_path = './images_upper'

# hdul = fits.open(mosaic_path + "P004+38_mosaic-blanked.fits")
# data = hdul[0].data




def getImage(gsigma, nsigma, npixels, data, name):
    print('Creating: %s.png'%(name))
    # Define noise threshold
    threshold = detect_threshold(data, nsigma=nsigma)

    # Sigma?
    sigma = gsigma * gaussian_fwhm_to_sigma  # FWHM = 3.
    # Use kernel (3x3) to find borders of sources
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    sources = detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)

    cmap = sources.make_cmap(random_state=12345)


    x = hdul[0]._header['NAXIS1']
    y = hdul[0]._header['NAXIS2']

    fig = plt.figure(frameon=False)
    fig.set_size_inches(x/1000,y/1000)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(sources, origin='upper', cmap=cmap)
    fig.savefig('%s/%s.png'%(write_path,name), dpi=1000)

    # plt.figure(figsize=(x/1000, y/1000), dpi=1000)
    # plt.imshow(sources, origin='upper', cmap=cmap)
    # plt.axis('off')
    # plt.savefig('%s/%s.png'%(write_path,name), dpi=1000)#, bbox_inches='tight', pad_inches=0)

    plt.clf()


files = [i for i in os.listdir(mosaic_path) if 'fits' in i]
for f in tqdm(files):
    hdul = fits.open(mosaic_path + f)
    data = hdul[0].data
    name = f.split('_')[0]
    getImage(3,3,16,data,name)

