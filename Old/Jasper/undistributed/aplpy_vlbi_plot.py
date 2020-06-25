# !/bin/python
# Plotting the VLBI images in aplpy

import matplotlib
import os
import re
import aplpy
import numpy as np
import glob
from astropy.io import fits
import time
import sys
import montage_wrapper as montage

start = time.time()

# Input
file_name = sys.argv[1]
stddev = float(sys.argv[2])
multi = float(sys.argv[3])
x = int(sys.argv[4])
y = int(sys.argv[5])
dx = float(sys.argv[6])
dy = float(sys.argv[7])
a, b = 1, 1
try:
    ax = int(sys.argv[8])
    ay = int(sys.argv[9])
    amx = float(sys.argv[10])
    amy = float(sys.argv[11])
    aname = str(sys.argv[12])
except:
    a = 0
try:
    bx = int(sys.argv[13])
    by = int(sys.argv[14])
    bmx = float(sys.argv[15])
    bmy = float(sys.argv[16])
    bname = str(sys.argv[17])
except:
    b = 0

# Names
name = str(file_name)
name = name.replace('.fits','')
name = name.replace('.FITS','')
j = 0
while True:
    if name.find('/', j+1, len(name)) == -1:
        if j != 0:
            directory = name[0:j+1]
        name = name[j+1:len(name)]
        break
    else:
        j = name.find('/', j+1, len(name))

# Create Figure
fig = aplpy.FITSFigure(file_name)
# fig.show_grayscale()

# Add Contours
levs = np.array([-2048, -1024, -512, -256, -128, -64, -32, -16, -8, -4, -2, -1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048], np.float64)
levs *= multi*stddev
fig.show_contour(file_name,colors='black',levels=levs)

# fig.add_colorbar()
# fig.colorbar.set_location('right')
# fig.colorbar.set_pad(0.07)
# fig.colorbar.set_axis_label_text('Flux (Jy/beam)')
fig.tick_labels.set_xformat('hh:mm:ss.sss')
fig.tick_labels.set_yformat('dd:mm:ss.sss')
fig.axis_labels.set_xtext('Right Ascension (J2000)')
fig.axis_labels.set_ytext('Declination (J2000)')

# Resize
x, y = fig.pixel2world(x, y)
dx *= 0.00000027777778
dy *= 0.00000027777778
fig.recenter(x, y, width=dx, height=dy)
fig.tick_labels.set_font(size='x-large',stretch='ultra-condensed')
fig.axis_labels.set_font(size='x-large')
# fig.add_grid()
# fig.grid.set_color('black')
# fig.grid.set_alpha(0.5)
#fig.axis_labels.hide()
if a == 1:
    ax, ay = fig.pixel2world(ax, ay)
    aname = "A ("+aname+")" 
    fig.add_label(ax, ay, aname, color='black', size= 15)
    fig.show_markers(amx, amy, marker='+', color='black', s=200, facecolor='darkred', edgecolor='black', lw=1.5, zorder=1, alpha=0.9)
if b == 1:
    bx, by = fig.pixel2world(bx, by)
    bname = "B ("+bname+")"
    fig.add_label(bx, by, bname, color='black', size= 15)
    fig.show_markers(bmx, bmy, marker='+', color='black', s=200, facecolor='darkred', edgecolor='black', lw=1.5, zorder=1, alpha=0.9)

# Add Beam
try:
    fig.add_beam(color='0.5')
    fig.beam.set_frame(True)
except(KeyError):
    print('No Beam info in the header.')

# Save
# fig.set_title(name)
fig.save(directory+name+'.png')
fig.close()

end = time.time()
print('Completed in: ',end-start)
