# Introduction
This repo is made for the course seminar distributed data mining.
In this course the potential of distributed computing is explored.
In this repo distribution is done using PySpark on the fs.dslc.liacs.nl server to find radio wave sources in space.
This is done by processing FITS files which contain "snapshots" of a particular place in the sky.
The end goal for this project was to produce a catalogue of the radio sources found in the  FITS files.
This catalogue contains their properties and location in the sky and is created in the first step "1.Catalogue_Creation".

# 1. Creating the catalogue
This directory contains 3 scripts, namely example.ipynb, create_catalogue_noSpark.py and create_catalogue.py.
The script example.ipynb is a jupyter notebook contains the steps taken to identify sources in a single mosaic.
This script also shows most ouputs of each step in order to give an idea about what data we are dealing with.
To execute the processes in a distributed way the file create_catalogue.py should be executed.
This script processes all our fits files in a distributed fashion.
It first identifies all sources in all the FITS files and each source is then redistributed using PySpark to process it more efficiently. The script create_catalogue_noSpark.py does the same process as the previous script but without PySpark to evaluate the speed up achieved with PySpark. 

# 2. Overlap
In this directory the files overlap_mosaics.csv and overlap_pairs.csv are created by executing filter_overlaps.ipynb.
This notebook also uses PySpark to find overlapping sources in the created catalogue in step 1.
These overlapping sources are unwanted and should be removed/filtered.
This directory also contains filter_overlaps_noSpark.ipynb which is again the implementation of previous script but without PySpark.
You should read the report if you want to know more about the algorithm to detect the overlapping sources.

# 3. Artifacts.
In this directory the catalogue from the previous step is taken and the artifacts are filtered sources.
These are artificial sources detected by the algorithms in step 1, caused bu highly intense sources.
These should also be removed/filtered, since these are not actual objects in the sky.

# 4. Component_Analyses
In this step the filtered catalogue is analysed and results are produced for our report. All the results are discussed there.

# 5. Lookup_Source
This directory contains the scripts explore_catalogue.ipynb and create_mosaic_images.ipynb.
In explore_catalogue.ipynb it is possible to find sources of interest in the final catalogue.
However, this first requires you to execute create_mosaic_images.ipynb to create the actual high resolution png file that contains the mosaics where the sources are located. This allows you to see the actual detected source.

# 0. Extra scripts
In the directories Deblending and Segmentation scripts are found which test parameters for segmentation of the sources and deblending sources. At the start of the project we wanted to explore the use and effect of deblending however we found that it was not usefull to apply it directly on all sources, therefore we did not implement it in our catalogue. However, in some cases deblending is quite usefull and should be used and it can also create interesting sources to look at.
