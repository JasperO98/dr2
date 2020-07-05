import shutil
import os


out = "/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_writable/"
shutil.move(out + "catalogue_v6.0/part-00000", 'catalogue_v6.0.csv')
shutil.rmtree(out + 'catalogue_v6.0')
