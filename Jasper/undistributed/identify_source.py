import sys
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from photutils import detect_sources
import warnings


def create_catalogue_data(data, sources):
	# Calculate integrated intensity, center of mass and brightest pixel per source
	integrated_intensity = np.zeros(sources.max_label)
	center_of_mass = np.zeros((sources.max_label, 2))
	brightest_pixel = np.zeros((sources.max_label, 3))

	for i in range(np.shape(data)[0]):
		for j in range(np.shape(data)[1]):
			lab = sources.data[i, j] - 1  # Label of pixel i,j
			val = data[i, j]  # Intensity value of pixel i,j
			if lab + 1 != 0:
				integrated_intensity[lab] += val  # Sum for each label the intensity of each pixel
				center_of_mass[lab, :] += val * np.array([i, j])  # Sum the intensity times pixel position (i, j)
				if val > brightest_pixel[lab][2]:
					brightest_pixel[lab] = np.array([i, j, val])  # Keep the highest value per source
	norm_center_of_mass = center_of_mass / np.array([integrated_intensity, integrated_intensity]).T
	return np.array([
        np.arange(1, sources.max_label+1), # Source label
        integrated_intensity, # Integrated intensity
        brightest_pixel[:,0], # x pixel
        brightest_pixel[:,1], # y pixel
        brightest_pixel[:,2], # intensity
        norm_center_of_mass[:,0], # norm center of mass x
        norm_center_of_mass[:,1] # norm center of mass x
    ]).T

def write_catalogue_data(data, out, file):
	np.save(out + '.'.join(file.split('/')[-1].split('.')[:-1]) + '.tmp', data)


def main():
	# ignore warnings
	warnings.filterwarnings("ignore")
	
	# Get params
	file = sys.argv[1]
	out = sys.argv[2] + '/'
	
	# Open fits file and retrieve the intensity values
	hdul = fits.open(file)
	data = hdul[0].data

	# Noise threshold
	threshold = photutils.detect_threshold(data, nsigma=2.)

	# Sigma
	sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.

	# Use kernel (3x3) to find borders of sources
	kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
	kernel.normalize()

	# Detect sources and save a label for each pixel
	sources = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
	# Get catologue data
	catalogue_data = create_catalogue_data(data, sources)

	# Write catalogue data
	write_catalogue_data(catalogue_data, out, file)


main()
