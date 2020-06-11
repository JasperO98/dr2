import pandas as pd
import numpy as np
import os
from PIL import Image, ImageDraw
from matplotlib import pyplot as plt

IMAGE_PATH = './images_upper'
WRITE_PATH = './temp'

# given two rows of component information, generate a plot contains two sub-plots represented by the two rows
def getTwoImgs(r1, r2, name):
    num_pxl = max(r1['total_pixels'], r2['total_pixels'])
    r1_crop = crop(r1, num_pxl)
    r2_crop = crop(r2, num_pxl)

    fig, ax = plt.subplots(1,2)
    ax[0].imshow(r1_crop)
    ax[0].set_title("%s;#pxl:%s"%(r1['label'], r1['total_pixels']))
    ax[1].imshow(r2_crop)
    ax[1].set_title("%s;#pxl:%s"%(r2['label'], r2['total_pixels']))

    plt.savefig(os.path.join(WRITE_PATH, '%s.png'%(name)))


def crop(r, num_pxl):
    r_img = Image.open(os.path.join(IMAGE_PATH, r['mosaic']+'.png'))
    w, h = r_img.size

    x = r['brightest_pixel_x']
    y = r['brightest_pixel_y']
    # x = r['center_of_mass_x']
    # y = r['center_of_mass_y']

    # draw a 5x5 rectangle around the brightest pixel
    draw = ImageDraw.Draw(r_img)
    draw.rectangle(((x-5, y-5), (x+5, y+5)), outline='red')

    # crop cout a sub-image with a size of num_pxl by num_pxl
    # the briestest pixel should be at the center of the sub-image
    left = x - num_pxl if x - num_pxl > 0 else 0
    top = y - num_pxl if y - num_pxl > 0 else 0
    right = x + num_pxl if x + num_pxl < w else w
    bottom = y + num_pxl if y + num_pxl < h else h

    print(left, top, right, bottom)
    r_crop = r_img.crop((left, top, right, bottom))

    return r_crop



df = pd.read_csv('./my_csv_0016.csv')
for i in range(0, len(df), 2):
    r1 = df.iloc[i]
    r2 = df.iloc[i+1]
    getTwoImgs(r1, r2, i)


