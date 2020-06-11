import photutils
import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
#from pyspark import SparkContext
#from pyspark.sql import SQLContext
from photutils import detect_threshold

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import os
import csv
import numpy as np
import plotly
from io import StringIO

import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.express as px

import matplotlib.pyplot as plt
import pandas as pd
import requests
from tqdm import tqdm

DR2_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2/"
mosaic_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/"

hdul = fits.open(mosaic_path + "P004+38_mosaic-blanked.fits")
data = hdul[0].data

# gsimga defines a Gaussuan kernel, which is used to smooth(convolve) the noise
# this parameter has the smallest influence
for gsigma in range(1,10):
    #nsigma defines threshold, larger the value higher the threshold
    for nsigma in range(3,10):
        # npixels defines the number of pixels a source at least contains
        for npixels in range(5,15):

            threshold = detect_threshold(data, nsigma=nsigma)

            # define a Gaussian kernal
            sigma = gsigma * gaussian_fwhm_to_sigma  # FWHM = 3.
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            kernel.normalize()

            sources = detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)

            

            # plot
            fig, ax1 = plt.subplots(1, 1, figsize=(10,10), frameon=False)
            ax1.set_axis_off()
            cmap = sources.make_cmap(random_state=12345)
            plt.imshow(sources, origin='lower', cmap=cmap)
            plt.savefig('./images/n%s_g%s_np%s.png'%(nsigma, gsigma, npixels), dpi=96 * 10, bbox_inches='tight', pad_inches=0)


            # calcualte the size distribution
            '''
            sources_array = np.asarray(sources)
            num_sources = np.max(sources_array)

            # return the size of each label
            _, size = np.unique(sources, return_counts=True)

            # count occurances of each size
            size, counts = np.unique(size, return_counts=True)

            fig, ax = plt.subplots(1,2,figsize=(10,6))
            ax[0].scatter(size[:-1], np.log2(counts[:-1]))
            ax[1].scatter(size[:-1], np.log2(counts[:-1]))
            ax[1].set_xlim([0,100])
            plt.savefig('./distri/n%s_g%s_np%s.png'%(nsigma, gsigma, npixels))
            '''
