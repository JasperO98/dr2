import pandas as pd

# generate a dictionary; 
# a key is one mosiac file, a value is a list of mosaic files which are overlapped with the key.
def getOverlapDict():
    df = pd.read_csv('./overlap_dict.csv')
    df = df[['fitsfile', 'overlaps']]

    overlaps = {}
    for index, row in df.iterrows():
        key = row['fitsfile'].split('_')[0]
        # if the file does not overlap with any other files
        if row.isnull().any().any():
            continue
        values = [i.split('_')[0] for i in row['overlaps'].split(';')]
        overlaps[key] = values
    return overlaps
