import pandas as pd
import numpy as np
import os
from PIL import Image, ImageDraw
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from photutils import detect_threshold
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel

import os, fnmatch

IMAGE_PATH = '/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/' # './images_upper'
WRITE_PATH = './temp'

# given two rows of component information, generate a plot contains two sub-plots represented by the two rows
def getTwoImgs(r1, r2, name):
    num_pxl = max(r1['total_pixels'], r2['total_pixels'])
    r1_crop, rect1, cmap1 = crop(r1, num_pxl)
    r2_crop, rect2, cmap2 = crop(r2, num_pxl)

    fig, ax = plt.subplots(1,2)
    ax[0].add_artist(rect1)
    ax[0].imshow(r1_crop, cmap=cmap1)
    ax[0].scatter(r1_crop.shape[0]// 2, r1_crop.shape[1] // 2, s=4, c='red')
    ax[0].set_title("%s;#pxl:%s"%(r1['label'], r1['total_pixels']))
    
    ax[1].add_artist(rect2)
    ax[1].imshow(r2_crop, cmap=cmap2)
    ax[1].scatter(r1_crop.shape[0]// 2, r1_crop.shape[1] // 2, s=4, c='red')
    ax[1].set_title("%s;#pxl:%s"%(r2['label'], r2['total_pixels']))

    plt.savefig(os.path.join(WRITE_PATH, '%s.png'%(name)))


def crop(r, num_pxl):
    file = ''.join(fnmatch.filter(os.listdir(IMAGE_PATH), r['mosaic'] + '*.fits'))
    print(file)
    hdul = fits.open(os.path.join(IMAGE_PATH, file))
    
    data = hdul[0].data
    threshold = detect_threshold(data, nsigma=3.)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    # Use kernel (3x3) to find borders of sources
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    # sources = object with label of sources for each pixel in the fits file
    sources = detect_sources(data, threshold, npixels=16, filter_kernel=kernel)
    
    r_img = sources._data
    # print(r_img)
    
    w, h = r_img.shape

    x = r['brightest_pixel_x']
    y = r['brightest_pixel_y']
    # x = r['center_of_mass_x']
    # y = r['center_of_mass_y']

    # draw a 5x5 rectangle around the brightest pixel
    # draw = ImageDraw.Draw(r_img)
    # draw.rectangle(((x-5, y-5), (x+5, y+5)), outline='red')
    rect = plt.Rectangle((x, y), 5, 5, color='red', fill=False)

    # crop cout a sub-image with a size of num_pxl by num_pxl
    # the briestest pixel should be at the center of the sub-image
    left = x - num_pxl if x - num_pxl > 0 else 0
    top = y - num_pxl if y - num_pxl > 0 else 0
    right = x + num_pxl if x + num_pxl < w else w
    bottom = y + num_pxl if y + num_pxl < h else h

    print(left, top, right, bottom)
    print('Ra:', r['brightest_pixel_RA'])
    print('Dec:', r['brightest_pixel_DEC'])
    print('Noise Threshold:', threshold)
    r_crop = r_img[top:bottom, left:right]
    cmap = sources.make_cmap(random_state=12345)
    return r_crop, rect, cmap



df = pd.read_csv('./my_csv_0016.csv')
for i in range(0, len(df), 2):
    r1 = df.iloc[i]
    r2 = df.iloc[i+1]
    getTwoImgs(r1, r2, i)


