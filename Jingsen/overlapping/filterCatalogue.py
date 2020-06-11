import pandas as pd
from overlap import getOverlapDict


header = "label,total_pixels,x_pixels,y_pixels,integrated_intensity,brightest_pixel,brightest_pixel_x,brightest_pixel_y,brightest_pixel_RA,brightest_pixel_DEC,center_of_mass_x,center_of_mass_y,center_of_mass_RA,center_of_mass_DEC,center_of_gaus_fit_x,center_of_gaus_fit_y,center_of_gaus_fit_RA,center_of_gaus_fit_DEC,fit_x_axis,fit_y_axis,fit_theta,deconv_x,deconv_y,integrated_intensity_fit,ratio_residual"
header = header.split(',')
df = pd.read_csv("./cata")
df.columns = header

# drive two columns: 1) mosaic file, 2) component ID from "label"
mosaic, comp = [], []
label = list(df['label'])
for item in label:
    m, c = item.split('_')
    mosaic.append(m)
    comp.append(c)

df['mosaic'] = mosaic
df['comp'] = comp

# a dictionary; a key is one mosiac file, a value is a list of mosaic files which are overlapped with the key.
overlaps = getOverlapDict()

cnt = 0
for i, (key, value) in enumerate(overlaps.items()):
    # create a dataframe which stores all components in the mosaic file "key"
    df_key = df[df['mosaic'] == key]
    # create a dataframe which stores components from all mosaic files in the list "value"
    df_value = df[df['mosaic'].isin(value)]

    print("searching %s/%s"%(i, len(overlaps)))

    eps = 0.0016

    # iterate over all components in the "value" dataframe to match components in the "key" value
    for i1, r1 in df_value.iterrows():
        ra = r1['brightest_pixel_RA']
        dec = r1['brightest_pixel_DEC']
        a = df_key[
                df_key['brightest_pixel_RA'].between(ra-eps, ra+eps) &
                df_key['brightest_pixel_DEC'].between(dec-eps, dec+eps)
            ]

        if not a.empty:
            a.to_csv('my_csv_0016.csv', mode='a', header=False, index=False)
            r1.to_frame().transpose().to_csv('my_csv_0016.csv', mode='a', header=False, index=False)
            cnt += 1

    print("accumulated number of matches: %s"%cnt)
