#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits
from photutils import detect_threshold,detect_sources
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

DR2_path = "../../../../data/mostertrij/data/LoTSS_DR2"

mosaic_path = "../../../../data/mostertrij/data/LoTSS_DR2_mosaic/"

listOfFiles = []
for (directory_path, directory_names, filenames) in os.walk(DR2_path):
    listOfFiles += [os.path.join(directory_path, file) for file in filenames if file.endswith('mosaic-blanked'+'.' + 'fits')]
    
listOfFiles
data = fits.open(listOfFiles[2])[0].data
sigma_1, sigma_2 = 7,10
# gsimga defines a Gaussuan kernel, which is used to smooth(convolve) the noise
# this parameter has the smallest influence
for gsigma in range(sigma_1,sigma_2):
    #nsigma defines threshold, larger the value higher the threshold
    for nsigma in range(sigma_1,sigma_2):
        # npixels defines the number of pixels a source at least contains
        for npixels in range(12,13):

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
            #plt.savefig('./../n%s_g%s_np%s.png'%(nsigma, gsigma, npixels), dpi=96 * 10, bbox_inches='tight', pad_inches=0)


# ## Deblending 

# In[2]:


from photutils import deblend_sources
segm_deblend = deblend_sources(data, sources, npixels=12,
                               filter_kernel=kernel, nlevels=32, # kernel should be same as that of thresholding
                               contrast=0.001)  
fig, ax1 = plt.subplots(1,figsize=(10,10), frameon=False)
ax1.set_axis_off()
cmapdbS = segm_deblend.make_cmap(random_state=12345)
plt.imshow(segm_deblend, origin='lower', cmap=cmapdbS)


# In[3]:


segm_deblend


# In[4]:


sources


# In[5]:


from photutils import deblend_sources
segm_deblend = deblend_sources(data, sources, npixels=12,
                               filter_kernel=kernel, nlevels=32, # kernel should be same as that of thresholding
                               contrast=0.001)  
fig, ax1 = plt.subplots(1,figsize=(10,10), frameon=False)
ax1.set_axis_off()
cmapdbS = segm_deblend.make_cmap(random_state=12345)
plt.imshow(segm_deblend, origin='lower', cmap=cmapdbS)


# ## Validating Deblending 

# In[6]:


segm_deblend


# In[7]:


print(sources)


# In[8]:


sources[55]


# In[9]:


from photutils import source_properties
cat_sources = source_properties(data, sources)
cat_segmented = source_properties(data, segm_deblend)
cat_segmented.to_table()


# In[10]:


import seaborn as sns
catList_sources = list(cat_sources.area)
catList_segmented =list(cat_segmented.area)
Areas_sources =[]
Areas_segmented = []
for i in catList_sources:
    Areas_sources.append(float(str(i)[:-5]))
for i in catList_segmented:
    Areas_segmented.append(float(str(i)[:-5]))

f, ax = plt.subplots(figsize=(20, 8))
sns.distplot(Areas_sources, label="Sources size distribution")
sns.distplot(Areas_segmented, label = "Sources size distribution after segmentation")
plt.legend()


# In[11]:


val = 1439  
x_min = int(str(cat_sources[val].bbox_xmin)[:-6])
x_max = int(str(cat_sources[val].bbox_xmax)[:-6])
y_min = int(str(cat_sources[val].bbox_ymin)[:-6])
y_max = int(str(cat_sources[val].bbox_ymax)[:-6])
plt.imshow(data[y_min:y_max,x_min:x_max], origin='lower', cmap ='Greys_r')


# In[12]:


plt.imshow(sources[1439])


# In[13]:


plt.imshow()


# In[ ]:


plt.imshow(segm_deblend, origin='lower', cmap=cmapdbS)
plt.xlim([x_min,x_max])
plt.ylim([y_min,y_max])


# In[ ]:


for i in range(3000):
    if (x_min +3) >int(str(cat_segmented[i].bbox_xmin)[:-6]) > (x_min -3 ):
        print(i)


# In[ ]:


int(str(cat_segmented[1439].bbox_xmin)[:-6])


# In[ ]:




