# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:16:30 2023

@author: jpv88
"""

import pandas as pd

path = r'C:\\Users\\jpv88\\Documents\\GitHub\\Image-Analysis\\annotation_gui\\'
file = 'test.csv'

barcodes_full = pd.read_csv(path + file)
barcodes = pd.DataFrame()

for (columnName, columnData) in barcodes_full.items():
    if ('checked' in columnName) and all(columnData.to_numpy()):
        gene_name = '_'.join(columnName.split('_')[:2])
        barcodes[gene_name] = barcodes_full[gene_name].to_numpy()
        
