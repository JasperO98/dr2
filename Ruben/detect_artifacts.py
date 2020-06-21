import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# INIT
CHUNKSIZE = 10**5
CATALOGUE_FILE = "./catalogue_v5.csv.tar"
HEADER_FILE = "./catalogue_v5_header.csv"
OVERLAP_FILE = "./overlap_dict.csv"
WRITE_TO = "./artifact_detection.csv"

INTENSITY_THRESHOLD = 0.8
RADIUS_THRESHOLD = 200 # arcseconds
ECCENTRICITY_THRESHOLD = 3

header = list(pd.read_csv(HEADER_FILE, sep=",", quotechar='"'))
header_dict = {x : header.index(x) for x in header}

# INTENSITY HISTOGRAM
disable = True
if not disable:
    print("Make intensities histogram")
    nbins = 100
    for ci, chunk in enumerate(pd.read_csv(CATALOGUE_FILE, compression='infer', header=0, sep=',', quotechar='"', chunksize=CHUNKSIZE, error_bad_lines=False)):
        print("Chunk", ci+1)
        br_pixels = chunk.iloc[:, header_dict['brightest_pixel']]
        br_pixels = br_pixels[np.isfinite(br_pixels)]
        
        if ci == 0:
            edges = np.histogram_bin_edges(br_pixels, bins=nbins)
            hist = np.histogram(br_pixels, bins=edges)[0].astype("float64")
        else:
            hist = hist + np.histogram(br_pixels, bins=edges)[0]

    hist = hist / np.sum(hist) # normalise
    hist *= 10
    hist += np.min(hist[hist > 0]) # scale-up, avoid zeroes
    hist = np.log10(hist) # logarithm

    fig, ax = plt.subplots()
    ax.set_xlabel("Intensity threshold")
    ax.set_ylabel("Log10-scaled frequency")
    ax.set_title("Floored log scale frequency of brightest pixel intensities")
    ax.plot(edges[0:-1], hist)
    plt.show()

# HIGH INTENSITY SOURCES
print("Find intense sources")
intense_sources = []
for ci, chunk in enumerate(pd.read_csv(CATALOGUE_FILE, compression='infer', header=0, sep=',', quotechar='"', chunksize=CHUNKSIZE, error_bad_lines=False)):
    print("Chunk", ci+1)
    intense = np.full(len(chunk.index), INTENSITY_THRESHOLD)
    br_pixels = chunk.iloc[:, header_dict['brightest_pixel']]
    br_indices = br_pixels >= intense
    br_chunk = chunk.iloc[np.where(br_indices == True)[0], :]
    for ri, row in br_chunk.iterrows():
        intense_sources.append(row)
    
print("Intense sources:", len(intense_sources))
intense_sources_dict = {} # convert to dictionary indexed by fits files
for src in intense_sources:
    fits = src[header_dict['label']].split("_")[0]
    try:
        intense_sources_dict[fits].append(src)
    except KeyError:
        intense_sources_dict[fits] = [src]

# FIND ARTIFACTS
print("Generate overlap dict")
overlap_file = pd.read_csv(OVERLAP_FILE, sep=",", quotechar='"') # find overlap mosaics
overlap_dict = {} # make a better dictionary
for _, mosaic in overlap_file.iterrows():
    overlaps = mosaic["overlaps"]
    overlaps = [overlaps, ""][pd.isna(overlaps)]
    overlaps = overlaps.split(";")
    for i, ov in enumerate(overlaps):
        overlaps[i] = ov.split("_")[0]
    overlap_dict[mosaic["fitsfile"].split("_")[0]] = overlaps

print("Start artifact finding")
artifact_indexed_list = [] # index-wise list of sources marked as artifacts
RADIUS_THRESHOLD_DEGREES = RADIUS_THRESHOLD / 3600 # convert from arcseconds to deg
for ci, chunk in enumerate(pd.read_csv(CATALOGUE_FILE, compression='infer', header=0, sep=',', quotechar='"', chunksize=CHUNKSIZE, error_bad_lines=False)):
    print("Chunk", ci+1) # pass over all chunks
    eccentricity = np.exp(np.abs(np.log(chunk.iloc[:, header_dict['fit_x_axis']] / chunk.iloc[:, header_dict['fit_y_axis']])))
    passes_eccentricity = eccentricity >= ECCENTRICITY_THRESHOLD # calculate sources above eccentricity threshold

    for ri, row in chunk.iterrows():
        if ri % 1000 == 0: print("Row", ri+1)
        mark_artifact = False
        if passes_eccentricity[ri]:
            fits = row[header_dict['label']].split("_")[0]
            mosaics_to_check = overlap_dict[fits] # high-intensity sources to check in intense_sources_dict
            mosaics_to_check.append(fits) # append fits file itself
            compile_near_sources = []
            for mosaic in mosaics_to_check:
                try:
                    for src in intense_sources_dict[mosaic]:
                        compile_near_sources.append([src[header_dict['center_of_mass_RA']], src[header_dict['center_of_mass_DEC']]])
                except KeyError:
                    pass
            compile_near_sources = np.array(compile_near_sources)
            vec = np.array([row[header_dict['center_of_mass_RA']], row[header_dict['center_of_mass_DEC']]]) - compile_near_sources
            distance = np.sqrt(np.sum(np.power(vec, 2))) # distance (euclidean)
            if np.sum(distance <= RADIUS_THRESHOLD_DEGREES) > 0:
                mark_artifact = True
            ##for isrc in compile_near_sources:
            ##    vec = np.array([row[header_dict['center_of_mass_RA']], row[header_dict['center_of_mass_DEC']]]) - np.array(isrc)
            ##    distance = np.sqrt(np.sum(np.power(vec, 2))) # distance (euclidean)
            ##    if distance <= RADIUS_THRESHOLD:
            ##        mark_artifact = True
            ##        break
        artifact_indexed_list.append(mark_artifact)

# WRITE OUT
print("Writing...")
wfile = open(WRITE_TO, 'w')
wfile.write("label,is_artifact") # make csv
for ci, chunk in enumerate(pd.read_csv(CATALOGUE_FILE, compression='infer', header=0, sep=',', quotechar='"', chunksize=CHUNKSIZE, error_bad_lines=False)):
    for ri, row in chunk.iterrows():
        index = ci * CHUNKSIZE + ri
        wfile.write("\n" + str(row[header_dict['label']]) + "," + str(artifact_indexed_list[index]))
wfile.close()
