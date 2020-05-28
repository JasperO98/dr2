import photutils
import os
import sys
import numpy as np

from photutils import detect_threshold
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.functional_models import Gaussian2D

from pyspark.sql.types import Row
from pyspark import SparkContext
from pyspark.sql import SparkSession

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
import time

import shutil
from math import pi, log, sqrt

import warnings
warnings.filterwarnings("ignore")



# The function to create the catalogue
def create_sub_matrices(hdul, sources):
    data = hdul.data
    min_max_coords = np.full((sources.max_label, 2, 2), ((np.inf, np.inf), (-np.inf, -np.inf)))
    for i in range(np.shape(data)[0]):
        for j in range(np.shape(data)[1]):
            lab = sources.data[i, j] - 1  # Label of pixel i,j
            if lab + 1 != 0:
                # Define min max coords
                min_coords = min_max_coords[lab][0]
                max_coords = min_max_coords[lab][1]
                # Check if smallest x coordinate
                if i < min_coords[0]:
                    min_coords[0] = i
                # Check if biggest x coordinate
                elif i > max_coords[0]:
                    max_coords[0] = i
                # Check if smallest y coordinate
                if j < min_coords[1]:
                    min_coords[1] = j
                # Check if biggest y coordinate
                elif j > max_coords[1]:
                    max_coords[1] = j
    
    min_max_coords = min_max_coords.astype(int)
    source_intensities = []
    for i in range(len(min_max_coords)):
        # Define min and max coords
        min_coords = min_max_coords[i][0]
        max_coords = min_max_coords[i][1] + 1
        # Mask of labels
        mask = sources.data[min_coords[0]:max_coords[0],
                            min_coords[1]:max_coords[1]]
        # Only take labels of current source and factor it against the intensities
        source_intensities.append((mask == i+1) * data[min_coords[0]:max_coords[0],
                                                       min_coords[1]:max_coords[1]])
    return np.column_stack([
        [hdul.header['OBJECT'] + '_' + str(l) for l in range(len(min_max_coords))],
        min_max_coords[:, 0, :],
        min_max_coords[:, 1, :],
        np.array(source_intensities)
    ])


def create_catalogue_submatrix(hdul, label, min_coords, max_coords, mat, deconv):
    if mat.size == 0:
        return (label, ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN',
                        'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])

    conv = pi * pow((hdul.header['BMAJ'] * 3600), 2) / (4 * log(2)) 
    w = WCS(hdul.header)
    
    # matrix filled with index values so first row for x = [0, 1, 2, ...] and y = [0, 0, 0, ..]
    # Used to fit gaussian and to calculate norm center of mass
    y, x = np.mgrid[:mat.shape[0], :mat.shape[1]]
    
    # Get num pixels
    total_pixels = np.sum(mat > 0)
    
    # Get sum intensity
    sum_mat = np.sum(mat)
    integrated_intensity = conv * sum_mat
    
    # Brightest pixel (account for local submatrix coordinates)
    brightest_pixel_x, brightest_pixel_y = np.unravel_index(np.argmax(mat, axis=None), mat.shape)
    brightest_pixel = mat[(brightest_pixel_x, brightest_pixel_y)]
    brightest_pixel_x += min_coords[0] # Add the min x coord of source
    brightest_pixel_y += min_coords[1] # Add the min y coord of source
    brightest_pixel_RA, brightest_pixel_DEC = w.all_pix2world(brightest_pixel_x,
                                                              brightest_pixel_y, 0, ra_dec_order=True)
    
    # center of mass = sum(( coords(x, y) * mat(inten) )) / integrated_intensity
    center_of_mass_x, center_of_mass_y = np.sum(np.sum(np.array(( (x + min_coords[0]),
                                                                  (y + min_coords[1]) )) * mat, axis=1), axis=1) / sum_mat 
   
    # center_of_mass = np.around(center_of_mass, 0).astype(int)

    center_of_mass_RA, center_of_mass_DEC = w.all_pix2world(center_of_mass_x,
                                                            center_of_mass_y, 0, ra_dec_order=True)
    
    # Define parameters of gaussian
    x_m = center_of_mass_x - min_coords[0]
    y_m = center_of_mass_y - min_coords[1]
    x_s = max_coords[0] - min_coords[0]
    y_s = max_coords[1] - min_coords[1]
    
    # Define model
    mod = Gaussian2D(
        amplitude=brightest_pixel, # Max intensity of source
        x_mean=x_m,    # X center of mass
        y_mean=y_m,    # Y center of mass
        x_stddev=x_s,  # X of submatrix
        y_stddev=y_s   # Y of submatrix
    )
    
    # Setting resutrictions for model
    mod.x_mean.min=0   # x_mean Min = 0
    mod.y_mean.min=0   # y_mean Min = 0
    mod.x_mean.max=x_s # x_mean Max = X of submatrix
    mod.y_mean.max=y_s # y_mean Max = Y of submatrix

    mod.x_stddev.min=0   # x_stdv Min = 0
    mod.y_stddev.min=0   # y_stdv Min = 0
    mod.x_stddev.max=x_s # x_stdv Max = X of submatrix
    mod.y_stddev.max=y_s # y_stdv Max = Y of submatrix
    
    fitter = LevMarLSQFitter() # Use least square fitter
    
    try:
        best_fit_gauss = fitter(mod, x, y, mat) # Fit model
    except:
        return (label, [total_pixels, x_s, y_s,
                        integrated_intensity, brightest_pixel, brightest_pixel_x, brightest_pixel_y,
                        float(brightest_pixel_RA), float(brightest_pixel_DEC), center_of_mass_x, center_of_mass_y,
                        center_of_mass_RA, center_of_mass_DEC, 'NaN', 'NaN', 'NaN', 'NaN',
                        'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']
       )
    
    # best_fit_gauss = fitter(mod, x, y, mat) # Fit model

    # Define centers of gaussian fit
    center_of_gaus_fit_x = best_fit_gauss.x_mean.value + min_coords[0]
    center_of_gaus_fit_y = best_fit_gauss.y_mean.value + min_coords[1]
    center_of_gaus_fit_RA, center_of_gaus_fit_DEC = w.all_pix2world(center_of_gaus_fit_x,
                                                                    center_of_gaus_fit_y, 0, ra_dec_order=True)
    
    # Define axis and theta of fit
    fit_x_axis = best_fit_gauss.x_stddev.value
    fit_y_axis = best_fit_gauss.y_stddev.value
    fit_theta  = best_fit_gauss.theta.value
    
    deconv_x = deconv(fit_x_axis)
    deconv_y = deconv(fit_y_axis)
    
    # Integrated intensity of fit
    integrated_intensity_fit = conv * np.sum(best_fit_gauss(x, y))
    
    # Residual / Source (sum intensities)
    ratio_residual = np.sum((mat - best_fit_gauss(x, y)).clip(min=0)) / sum_mat
    
    return (label, [total_pixels, x_s, y_s,
                    integrated_intensity, brightest_pixel, brightest_pixel_x, brightest_pixel_y,
                    float(brightest_pixel_RA), float(brightest_pixel_DEC), center_of_mass_x, center_of_mass_y,
                    float(center_of_mass_RA), float(center_of_mass_DEC), center_of_gaus_fit_x, center_of_gaus_fit_y,
                    float(center_of_gaus_fit_RA), float(center_of_gaus_fit_DEC), fit_x_axis, fit_y_axis, fit_theta,
                    deconv_x, deconv_y, integrated_intensity_fit, ratio_residual]
           )


# Transform lines to csv string
def toCSVLine(obj, data):
    # return obj + [np.array2string(t, precision=10, separator=',', max_line_width=np.inf)[1:-1] for t in data]
    return obj + ',' + str(data).replace(' ', '').strip('[]')


def main():
    # Ignore the warnings
    warnings.filterwarnings("ignore")
    
    # Create spark context 132.229.44.220:7077
    sc = SparkContext(appName="dr2", master='spark://132.229.44.220:7077')
    
    # Define data paths
    DR2_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2/"
    mosaic_path = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/"
    out = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_writable/"
    
    # Get paths of fits files
    fits_files = [mosaic_path + f for f in os.listdir(mosaic_path)]
    
    print('Analyzing:')
    for f in range(len(fits_files)):
        print(fits_files[f])
    
    # Put the paths in an RDD and determine number of partitions
    # More partitions == more cpu and faster (can crash when partitions are too high)
    file_paths = sc.parallelize(fits_files, 8)
    
    # Store paths of fits files in rdd
    fits_content = file_paths.map(lambda file: fits.open(file)[0])
    
    # Calculate noise threshold
    fits_thresh = fits_content.map(lambda content: (content, detect_threshold(content.data , nsigma=3.)))


    # Sigma?
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.

    # Use kernel (3x3) to find borders of sources
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    # sources = object with label of sources for each pixel in the fits file
    fits_sources= fits_thresh.map(lambda ft: (ft[0],
                                                      detect_sources(ft[0].data, ft[1], npixels=16, filter_kernel=kernel)
                                                     )
                                         )
    
    # RDD with fits file as key and value = list ( list( label, min_x, min_y, max_x, max_y, intensity of source ), ...)
    sub_matrix = fits_sources.map(lambda fts: (fts[0], list(create_sub_matrices(fts[0], fts[1]))))
    sub_matrix.getNumPartitions()
    
    deconv = lambda axis: sqrt( (pow( axis * 2 * sqrt(2 * log(2)), 2) - 16).clip(min=0) )
    # RDD with fits file as key and value = list( label, min_x, min_y, max_x, max_y, intensity of source ) for each source
    catalogue_data = sub_matrix.flatMap(lambda x : map(lambda y: create_catalogue_submatrix(x[0], y[0],
                                                                                        (y[1], y[2]),
                                                                                        (y[3], y[4]), y[5], deconv),
                                                   x[1]
                                                  )
                                   )
    
    # Header for csv
    header = ["label", "total_pixels", "x_pixels", "y_pixels", 
              "integrated_intensity", "brightest_pixel", "brightest_pixel_x", "brightest_pixel_y",
              "brightest_pixel_RA", "brightest_pixel_DEC", "center_of_mass_x", "center_of_mass_y",
              "center_of_mass_RA", "center_of_mass_DEC", "center_of_gaus_fit_x", "center_of_gaus_fit_y",
              "center_of_gaus_fit_RA", "center_of_gaus_fit_DEC", "fit_x_axis", "fit_y_axis", "fit_theta",
              "deconv_x", "deconv_y", "integrated_intensity_fit", "ratio_residual"
             ]
                          
    # Save the RDDs as one (coaslesce=1) csv 
    if os.path.isdir(out + 'catalogue_v3'):
        shutil.rmtree(out + 'catalogue_v3')


    start = time.time()
    catalogue_data.map(lambda c: toCSVLine(c[0], c[1])).coalesce(1, shuffle = True).saveAsTextFile(out + 'catalogue_v3')
    print('Execution time:', ((time.time()-start) / 60 ), 'minutes')
    
    # Exit and stop the sparkcontext
    sc.stop()

main()



