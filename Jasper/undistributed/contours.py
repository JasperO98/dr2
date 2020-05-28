# --------------------------------------------------------------------
# Program to make VLBA images for 2019 A+A Paper
# Created by Michael Keim (michaelkeim2468@gmail.com)
# --------------------------------------------------------------------


# General usage
import os
import sys
import time

# For calculations
import numpy as np
from math import log10, floor
from PIL import Image

# For plotting
import matplotlib.pyplot as mpl
import aplpy
from astropy.io import fits
import montage_wrapper as montage
from matplotlib import rcParams
df = 1.5
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 20*df
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker


# --------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------


def magnitude(x):
	return int(log10(abs(x)))

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# --------------------------------------------------------------------
# Read in data
# --------------------------------------------------------------------


data_name          = 'specifications.csv'
path               = np.genfromtxt(data_name, delimiter=', ', usecols=0 , dtype=str)
name               = np.genfromtxt(data_name, delimiter=', ', usecols=1 , dtype=str)
extension          = np.genfromtxt(data_name, delimiter=', ', usecols=2 , dtype=str)
standard_deviation = np.genfromtxt(data_name, delimiter=', ', usecols=3 , dtype=np.float64) # Units?
first_contour      = np.genfromtxt(data_name, delimiter=', ', usecols=4 , dtype=np.float64)
x_center           = np.genfromtxt(data_name, delimiter=', ', usecols=5 , dtype=np.float64) # Pixels
y_center           = np.genfromtxt(data_name, delimiter=', ', usecols=6 , dtype=np.float64) # Pixels
x_width            = np.genfromtxt(data_name, delimiter=', ', usecols=7 , dtype=np.float64) # Milliarcseconds
y_width            = np.genfromtxt(data_name, delimiter=', ', usecols=8 , dtype=np.float64) # Milliarcseconds
comp_a_x_center    = np.genfromtxt(data_name, delimiter=', ', usecols=9 , dtype=np.float64) # Pixels
comp_a_y_center    = np.genfromtxt(data_name, delimiter=', ', usecols=10, dtype=np.float64) # Pixels
comp_a_x_degrees   = np.genfromtxt(data_name, delimiter=', ', usecols=11, dtype=np.float64) # Degrees
comp_a_y_degrees   = np.genfromtxt(data_name, delimiter=', ', usecols=12, dtype=np.float64) # Degrees
comp_a_x_ra        = np.genfromtxt(data_name, delimiter=', ', usecols=13, dtype=str)        # hms
comp_a_y_dec       = np.genfromtxt(data_name, delimiter=', ', usecols=14, dtype=str)        # dms
comp_b_x_center    = np.genfromtxt(data_name, delimiter=', ', usecols=15, dtype=np.float64) # Pixels
comp_b_y_center    = np.genfromtxt(data_name, delimiter=', ', usecols=16, dtype=np.float64) # Pixels
comp_b_x_degrees   = np.genfromtxt(data_name, delimiter=', ', usecols=17, dtype=np.float64) # Degrees
comp_b_y_degrees   = np.genfromtxt(data_name, delimiter=', ', usecols=18, dtype=np.float64) # Degrees
comp_b_x_ra        = np.genfromtxt(data_name, delimiter=', ', usecols=19, dtype=str)        # hms
comp_b_y_dec       = np.genfromtxt(data_name, delimiter=', ', usecols=20, dtype=str)        # dms
stddev_nrnd        = np.genfromtxt(data_name, delimiter=', ', usecols=21, dtype=np.float64) # Units?
pmax               = np.genfromtxt(data_name, delimiter=', ', usecols=22, dtype=np.float64)
offset             = np.genfromtxt(data_name, delimiter=', ', usecols=23, dtype=np.int64)


# --------------------------------------------------------------------
# Plot images
# --------------------------------------------------------------------


offset_cal = ''
for i in range(len(name)):

	# Make all same size
	x_width[i], y_width[i] = 100., 100.

	# Except A few C Band
	if i == 10 or i == 8 or i == 6 or i == 4:
	#if name[i][10] == 'C':
		x_width[i], y_width[i] = 50., 50.
	if i > 7:
		if name[i][10] == 'C':
			x_width[i], y_width[i] = 40., 40.
		else:
			x_width[i], y_width[i] = 80., 80.			
	if i == 1 or i == 2:
		x_width[i], y_width[i] = 80., 80.
	if i == 4:
		x_width[i], y_width[i] = 40., 40.
	if i == 5:
		x_width[i], y_width[i] = 80., 80.
	if i == 0:
		x_width[i], y_width[i] = 80., 80.
	if i == 1:
		x_width[i], y_width[i] = 60., 60.
	if i == 2:
		x_width[i], y_width[i] = 60., 60.
	if i == 10:
		x_width[i], y_width[i] = 30., 30.
	if i == 9:
		x_width[i], y_width[i] = 60., 60.
	if i == 8:
		x_width[i], y_width[i] = 30., 30.
	if i == 7:
		x_width[i], y_width[i] = 80., 80.

	# Get fits file
	file = path[i] + '/' + name[i] + '.' + extension[i]
	fig = mpl.figure(figsize=(10, 10))

	

	# Add axes and title
	# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	ax = fig.add_axes([0.1, 0.1118, 0.7765, 0.776])
	title = name[i].replace('_L',' at 1.55 GHz')
	title = title.replace('_C',' at 4.96 GHz')
	title = title[0:5] + '+' + title[5:len(title)] if comp_a_y_degrees[i] > 0. else title[0:5] + '-' + title[5:len(title)]
	ax.set_title(title, fontsize=22.5*df)
	ax.set_xlabel('Relative R.A. (mas)', fontsize=22.5*df)
	ax.set_ylabel('Relative Dec. (mas)', fontsize=22.5*df)
	ax.set_xlim(-x_width[i]/2., x_width[i]/2.)
	ax.set_ylim(-y_width[i]/2., y_width[i]/2.)
	if i > 7 and i != 8 and i != 9:
		if name[i][10] == 'C':
			ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
			ax.set_yticks([-15, -10, -5, 0, 5, 10, 15])
		else:
			ax.set_xticks([-30, -15, 0, 15, 30])	
			ax.set_yticks([-30, -15, 0, 15, 30])
	if i == 2 or i == 9 or i == 1:
		ax.set_xticks([-20, -10, 0, 10, 20])
		ax.set_yticks([-20, -10, 0, 10, 20])	
	if i == 4:
		ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
		ax.set_yticks([-15, -10, -5, 0, 5, 10, 15])
	if i == 5:
		ax.set_xticks([-30, -15, 0, 15, 30])
		ax.set_yticks([-30, -15, 0, 15, 30])
	if i == 0:
		ax.set_xticks([-30, -15, 0, 15, 30])
		ax.set_yticks([-30, -15, 0, 15, 30])
	if i == 7:
		ax.set_xticks([-30, -15, 0, 15, 30])
		ax.set_yticks([-30, -15, 0, 15, 30])
	if i == 8 or i == 10:
		ax.set_xticks([-10, -5, 0, 5, 10])
		ax.set_yticks([-10, -5, 0, 5, 10])

	# Add image
	f = aplpy.FITSFigure(file, figure=fig, subplot=[0.1, 0.1, 0.8, 0.8])
	hdul = fits.open(file)
	data = hdul[0].data
	order = np.max([abs(magnitude(np.max(data[0][0]))), abs(magnitude(np.min(data[0][0])))])

	# Recenter L Band, which are all odd indicies
	if (i % 2) != 0:
		x_p0, y_p0  = x_center[i], y_center[i]
		x_d0, y_d0 = f.pixel2world(x_center[i], y_center[i])
		x_center[i], y_center[i] = f.world2pixel(x_center[i-1], y_center[i-1])
		x_pf, y_pf  = x_center[i], y_center[i]
		x_df, y_df = f.pixel2world(x_center[i], y_center[i])
		if comp_b_x_center[i] > 0.:
			comp_a_x_center[i]  += x_pf - x_p0
			comp_a_y_center[i]  += y_pf - y_p0
			comp_b_x_center[i]  += x_pf - x_p0
			comp_b_y_center[i]  += y_pf - y_p0
		# if i == 0:
		# 	comp_a_x_degrees[i] += x_df - x_d0
		# 	comp_a_y_degrees[i] += y_df - y_d0
		# 	comp_b_x_degrees[i] += x_df - x_d0
		# 	comp_b_y_degrees[i] += y_df - y_d0


	# Add colorbar
	# f.show_colorscale(pmin=0., pmax=100.)
	pmax[0] = 99.85
	f.show_colorscale(pmin=0., pmax=pmax[i], cmap='gist_heat')
	f.add_colorbar()
	# f.colorbar.set_axis_label_text('Flux Density (Jy/beam)')
	# f.colorbar.set_axis_label_font(fontproperties= font_manager.FontProperties(family = 'serif', size = 22.5*df))
	# f.colorbar.set_axis_label_pad(((np.max(offset)-offset[i])/1.386))
	# fig.colorbar(, ax=ax, format=ticker.FuncFormatter(fmt))
	f.colorbar.set_axis_label_text('Intensity (Jy/beam)')
	f.colorbar.set_axis_label_font(fontproperties= font_manager.FontProperties(family = 'serif', size = 22.5*df))
	f.colorbar.set_axis_label_pad(20)
	# f.colorbar.set_ticks([10])
	# f.colorbar._ticklabel_fontproperties
	# f.colorbar._colorbar_axes.figure.show(labels=True)

	# Add scalebar
	f.add_scalebar(10.*0.00000027777778)
	f.scalebar.set_color('white')
	f.scalebar.set_corner('top right')
	f.scalebar.set_label('10 mas')

	# Recenter
	x_center[i], y_center[i] = f.pixel2world(x_center[i], y_center[i])
	f.recenter(x_center[i], y_center[i], width=x_width[i]*0.00000027777778, height=y_width[i]*0.00000027777778)

	# Add beam
	f.add_beam(color='0.5')
	f.beam.set_frame(True)

	# Add contours
	# levels = np.array([-4,4,5,6,7,10,20,50,100,200,400,800,1600],dtype='float64')
	# levels *= stddev_nrnd[i]
	# levels = np.concatenate(([-(2.)**(j) for j in range(11, -1, -1)], [(2.)**(j) for j in range(0, 12, 1)]))
	# levels *= standard_deviation[i] * first_contour[i]
	# levels *= stddev_nrnd[i] * 3.
	# f.show_contour(file, colors='black', levels=levels)
	# levels = np.array([-4, 4, 5, 6, 7, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 ],dtype='float64')
	levels = np.concatenate(([-(np.sqrt(2.))**(j) for j in range(17, -1, -1)], [(np.sqrt(2.))**(j) for j in range(0, 18, 1)]))
	print(levels)
	levels *= 3.# 4.
	print(levels)
	levels *= stddev_nrnd[i]
	f.show_contour(file, colors='lightgrey', levels=levels)


	# Remove ticks
	f.axis_labels.hide_x()
	f.axis_labels.hide_y()
	f.tick_labels.hide()
	f.tick_labels.hide_x()
	f.tick_labels.hide_y()

	# Add component labels
	if comp_b_x_center[i] > 0.:
		# f.show_markers(comp_a_x_degrees[i], comp_a_y_degrees[i], marker='+', color='black', s=200, facecolor='darkred', edgecolor='black', lw=3.5, zorder=1, alpha=0.9)
		# f.show_markers(comp_b_x_degrees[i], comp_b_y_degrees[i], marker='+', color='black', s=200, facecolor='darkred', edgecolor='black', lw=3.5, zorder=1, alpha=0.9)
		comp_a_x_center[i], comp_a_y_center[i] = f.pixel2world(comp_a_x_center[i], comp_a_y_center[i])
		comp_b_x_center[i], comp_b_y_center[i] = f.pixel2world(comp_b_x_center[i], comp_b_y_center[i])
		# f.add_label(comp_a_x_center[i], comp_a_y_center[i], 'W', color='darkred', size= 30)
		# f.add_label(comp_b_x_center[i], comp_b_y_center[i], 'E', color='darkred', size= 30)
		f.add_label(comp_a_x_center[i], comp_a_y_center[i], 'W', color='white', size= 30)
		f.add_label(comp_b_x_center[i], comp_b_y_center[i], 'E', color='white', size= 30)
		factorax = 1.
		if i == 2 or i == 3 or i == 6 or i == 7:
			factorax = 2.
		comp_a_x_center[i] = comp_a_x_center[i] + factorax*2.0*0.00000027777778*x_width[i]/100. if comp_a_x_center[i] < comp_a_x_degrees[i] else comp_a_x_center[i] - factorax*2.0*0.00000027777778*x_width[i]/100.
		comp_a_y_center[i] = comp_a_y_center[i] + 2.0*0.00000027777778*x_width[i]/100. if comp_a_y_center[i] < comp_a_y_degrees[i] else comp_a_y_center[i] - 2.0*0.00000027777778*x_width[i]/100.
		comp_b_y_center[i] = comp_b_y_center[i] + 2.3*0.00000027777778*x_width[i]/100. if comp_b_y_center[i] < comp_b_y_degrees[i] else comp_b_y_center[i] - 2.3*0.00000027777778*x_width[i]/100.
		comp_b_x_center[i] = comp_b_x_center[i] + 2.3*0.00000027777778*x_width[i]/100. if comp_b_x_center[i] < comp_b_x_degrees[i] else comp_b_x_center[i] - 2.3*0.00000027777778*x_width[i]/100.
		dx = comp_b_x_degrees[i]-comp_a_x_degrees[i]
		dy = comp_b_y_degrees[i]-comp_a_y_degrees[i]
		ls = np.sqrt(dx*dx + dy*dy)
		ls /= 0.00000027777778
		print("Linear Size = " + str(ls))
		f.show_lines([np.array( [ [comp_a_x_degrees[i], comp_a_x_center[i]], [comp_a_y_degrees[i], comp_a_y_center[i]] ] )], color='lightgrey')
		f.show_lines([np.array( [ [comp_b_x_degrees[i], comp_b_x_center[i]], [comp_b_y_degrees[i], comp_b_y_center[i]] ] )], color='lightgrey')

	# Save figure
	fig.savefig(name[i]+'.png', bbox_inches='tight')
	f.close()
	# if 1 == 0: #if comp_b_x_center[i] > 0.:
	# 	os.system('open ' + name[i] + '.png')
	# 	break
	# if comp_b_x_center[i] != -1:
	# 	continue

	im = Image.open(name[i] + '.png')
	xim, yim = im.size
	# print(np.max(offset), offset[i])
	# if np.max(offset) != offset[i]:
	background = Image.new('RGBA', (np.max(offset), yim), (255, 255, 255, 255))
	off = (int(round(((np.max(offset) - xim) / 2), 0)), int(round(((yim - yim) / 2),0)))
	background.paste(im) #, (0, 0, np.max(offset), yim))
	background.save(name[i] + '.png')
	x2im, y2im = background.size
	print(0.93/(y2im/x2im+2))
	print(0.93-2.*0.93/(y2im/x2im+2))

	offset_cal = offset_cal + str(x_center[i]) +', ' +str(y_center[i])+'\n'
	#print(0.93/(yim/xim+2))
	#print(0.93-2.*0.93/(yim/xim+2))
print(offset_cal)



