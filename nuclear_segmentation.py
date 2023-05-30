# -*- coding: utf-8 -*-
"""
Created on Mon May 15 19:48:42 2023

@author: jpv88
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import skimage

from scipy import ndimage as ndi
from skimage.metrics import structural_similarity as ssim
from tqdm import tqdm
from os import listdir
from os.path import isfile, join
from matplotlib.widgets import Button

import pandas as pd
from ipywidgets import widgets, interactive, fixed, interact

from matplotlib.widgets import Button
import napari

path = r'E:\\05-23_phox2b_retro_tracing\\image_processing\\slice_1\\individual genes\\'
fname = 'GFP_r1.tif'

image = skimage.io.imread(path+fname)

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
Snap25_image_fns = [fname for fname in onlyfiles if 'Snap25' in fname]
Reln_image_fns = [fname for fname in onlyfiles if 'Reln' in fname]
Snap25_images = []
for fn in Snap25_image_fns:
    Snap25_images.append(skimage.io.imread(path+fn))
    
for fn in Reln_image_fns:
    Snap25_images.append(skimage.io.imread(path+fn))
    
    
    # %%

def extract_plane_window(coords, image, window_sz):
    
    # row is dimension, columns are min and max
    to_extract = np.zeros((3, 2))
    
    image_dims = np.shape(image)
    
    for i in range(2):
        i += 1
        coord = coords[i]
        to_extract[i,:] = [coord - window_sz/2, coord + window_sz/2]
        if (coord + window_sz/2) > image_dims[i]:
            to_extract[i,1] = image_dims[i]
        if (coord - window_sz/2) < 0:
            to_extract[i,0] = 0
    
    to_extract[0,:] = [coords[0], coords[0]]
    to_extract = to_extract.astype(int)
    
    z = to_extract[0,0]
    xmin = to_extract[1,0]
    xmax = to_extract[1,1]
    ymin = to_extract[2,0]
    ymax = to_extract[2,1]
    
    return image[z, xmin:xmax, ymin:ymax]

# assumes that there is only one file with the given token in the folder
def load_image(token, path):
    
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    image_fns = [fname for fname in onlyfiles if token in fname]
    
    image_fn = image_fns[0]
    
    image = skimage.io.imread(path + image_fn)
    
    return image

# loads the Snap25 image from the same round as the token, assumes only one 
# file with the given token in the folder
def load_Snap25_image(token, path):
    
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    image_fns = [fname for fname in onlyfiles if token in fname]
    image_fn = image_fns[0]
    
    round_str = image_fn.split('_')[1][:2]
    
    Snap25_token = 'Snap25_' + round_str
    
    image_fns = [fname for fname in onlyfiles if Snap25_token in fname]
    image_fn = image_fns[0]
    
    image = skimage.io.imread(path + image_fn)
    
    return image

def find_gene_names(path):
    
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    genes = []
    for f in onlyfiles:
        round_str = f.split('.')[0]
        genes.append(round_str)
        
    genes = [gene.split('_')[0] for gene in genes]
    
    return genes

# normalize pixel intensity values in im2 to match im1
def normalize_grayscale(im1, im2):
    
    newMax = np.max(im1)
    newMin = np.min(im1)
    
    Max = np.max(im2)
    Min = np.min(im2)
    
    im_norm = (im2 - Min)*(newMax-newMin)/(Max-Min) + newMin
    
    return im_norm


# %%

# minimum cell volume to be kept, in voxels
min_cell_size = 50

n_planes = np.shape(image)[0]

thresh_image = np.zeros(np.shape(image)) 

for i in range(n_planes):
    thresh = skimage.filters.threshold_otsu(image[i,:,:])
    thresh_image[i,:,:] = (image[i,:,:] > thresh)
    
thresh_image = scipy.ndimage.binary_fill_holes(thresh_image)
    
connect_image = skimage.measure.label(thresh_image)

nucleus_radius = 15

distance = ndi.distance_transform_edt(thresh_image)
coords = skimage.feature.peak_local_max(distance, 
                                        footprint=np.ones((nucleus_radius, nucleus_radius, nucleus_radius)), 
                                        labels=thresh_image)
mask = np.zeros(distance.shape, dtype=bool)
mask[tuple(coords.T)] = True
markers = skimage.measure.label(mask, connectivity=3)
markers = skimage.segmentation.expand_labels(markers, distance=1)
markers = skimage.measure.label(markers>0, connectivity=3)
labels = skimage.segmentation.watershed(-distance, markers, mask=thresh_image)

# test = napari.utils.colormaps.low_discrepancy_image(labels)
# viewer = napari.view_image(image)
# labels_layer = viewer.add_labels(markers, name='segmentation')
# labels_layer = viewer.add_labels(labels, name='segmentation')
# labels_layer = viewer.add_labels(mask, name='segmentation')
# print(np.max(labels))


n_cells = np.max(labels)
cell_sizes = np.zeros((n_cells, 2))
cell_sizes[:,0] = np.array(range(n_cells)) + 1

for i in tqdm(range(n_cells), desc='Determining cell sizes...'):
    cell_sizes[i,1] = np.sum(labels == cell_sizes[i,0])

for i in tqdm(range(n_cells), desc='Removing small cells...'):
    if cell_sizes[i,1] < min_cell_size:
        labels[labels == cell_sizes[i,0]] = 0
        
labels = skimage.measure.label(labels)
n_cells = np.max(labels)
cell_centroids = np.zeros((n_cells, 3))

for i in tqdm(range(n_cells), desc='Calculating cell centroids'):
    coords = np.argwhere(labels == i+1)
    cell_centroids[i,:] = np.round(np.mean(coords, axis=0))
    
cell_centroids = cell_centroids.astype(int)

NMIs = np.zeros(len(cell_centroids))
for i in range(len(cell_centroids)):
    window1 = extract_plane_window(cell_centroids[i,:], Snap25_images[0], 100)
    window2 = extract_plane_window(cell_centroids[i,:], Snap25_images[1], 100)
    NMIs[i] = skimage.metrics.normalized_mutual_information(window1, window2)

thresh = skimage.filters.threshold_otsu(NMIs)
cell_centroids = cell_centroids[NMIs > thresh]


df = pd.DataFrame(cell_centroids)
df = df.rename(columns={0: "z", 1: "y", 2: "x"})


df.to_csv('cell_centroids')

# %%

k = 54

window1 = extract_plane_window(cell_centroids[k,:], Snap25_images[0], 100)
window2 = extract_plane_window(cell_centroids[k,:], Snap25_images[1], 100)

fig, ax = plt.subplots(1,2)
ax[0].imshow(window1)
ax[1].imshow(window2)

# %%



gene = 'Dlk1'
centroid_idx=1

path = r'E:\\05-23_phox2b_retro_tracing\\image_processing\\slice_1\\individual genes\\'

GFP_image = load_image('GFP', path)
Snap25_GFP_image = load_Snap25_image('GFP', path)
gene_image = load_image(gene, path)
Snap25_gene_image = load_Snap25_image(gene, path)

fig = plt.figure(figsize=[8, 8])
ax1 = plt.subplot(2,2,1)
plt.axis('off')
ax2 = plt.subplot(2,2,2, sharex=ax1, sharey=ax1)
plt.axis('off')
ax3 = plt.subplot(2,2,3, sharex=ax1, sharey=ax1)
plt.axis('off')
ax4 = plt.subplot(2,2,4, sharex=ax1, sharey=ax1)
plt.axis('off')

fig.subplots_adjust(wspace=0, hspace=0)

window_sz = 100

z = df['z'][0]

ax1.imshow(skimage.exposure.adjust_log(gene_image[z,:,:]), cmap='magma')
ax2.imshow(skimage.exposure.adjust_log(Snap25_gene_image[z,:,:]), cmap='magma')
ax3.imshow(skimage.exposure.adjust_log(Snap25_GFP_image[z,:,:]), cmap='magma')
ax4.imshow(skimage.exposure.adjust_log(GFP_image[z,:,:]), cmap='magma')

ax1.set_title(str(cell_centroids[centroid_idx, :]))

x = df['x'][0]
y = df['y'][0]

axes = [ax1, ax2, ax3, ax4]

for ax in axes:
    ax.set_xlim(x - window_sz/2, x + window_sz/2)
    ax.set_ylim(y - window_sz/2, y + window_sz/2)
    
class Index:
    ind = 0

    def yes(self, event):
        self.ind += 1
        plt.draw()

    def no(self, event):
        self.ind += 1
        ax1.imshow(skimage.exposure.adjust_log(gene_image[self.ind,:,:]), cmap='magma')
        ax2.imshow(skimage.exposure.adjust_log(Snap25_gene_image[self.ind,:,:]), cmap='magma')
        ax3.imshow(skimage.exposure.adjust_log(Snap25_GFP_image[self.ind,:,:]), cmap='magma')
        ax4.imshow(skimage.exposure.adjust_log(GFP_image[self.ind,:,:]), cmap='magma')
    
ax_yes = fig.add_axes([0.7, 0.05, 0.1, 0.075])
ax_no = fig.add_axes([0.81, 0.05, 0.1, 0.075])

b_yes = Button(ax_yes,'Yes')
b_no = Button(ax_no, 'No')

callback = Index()
b_yes.on_clicked(callback.yes)
b_no.on_clicked(callback.no)

        






# %%

path = r'E:\\05-23_phox2b_retro_tracing\\image_processing\\slice_1\\individual genes\\'

GFP_image = load_image('GFP', path)
Snap25_GFP_image = load_Snap25_image('GFP', path)

gene = 'Dlk1'
gene_image = load_image(gene, path)
Snap25_gene_image = load_Snap25_image(gene, path)

fig = plt.figure(figsize=[8, 8])
ax1 = plt.subplot(2,2,1)
plt.axis('off')
ax2 = plt.subplot(2,2,2, sharex=ax1, sharey=ax1)
plt.axis('off')
ax3 = plt.subplot(2,2,3, sharex=ax1, sharey=ax1)
plt.axis('off')
ax4 = plt.subplot(2,2,4, sharex=ax1, sharey=ax1)
plt.axis('off')

num_ROIs = df.shape[0]

axes = [ax1, ax2, ax3, ax4]


def plot_ROI(centroid_idx, fig, axes):

    window_sz = 100
    
    z = df['z'][centroid_idx]
    x = df['x'][centroid_idx]
    y = df['y'][centroid_idx]
    
    [ax1, ax2, ax3, ax4] = axes

    ax1.imshow(skimage.exposure.adjust_log(gene_image[z,:,:]), cmap='magma')
    ax2.imshow(skimage.exposure.adjust_log(Snap25_gene_image[z,:,:]), cmap='magma')
    ax3.imshow(skimage.exposure.adjust_log(Snap25_GFP_image[z,:,:]), cmap='magma')
    ax4.imshow(skimage.exposure.adjust_log(GFP_image[z,:,:]), cmap='magma')

    ax1.set_title(str(cell_centroids[centroid_idx, :]))

    for ax in axes:
        ax.set_xlim(x - window_sz/2, x + window_sz/2)
        ax.set_ylim(y - window_sz/2, y + window_sz/2)

global centroid_idx 
centroid_idx = 0
    
from tkinter import *

root = Tk()  # create parent window

# use Button and Label widgets to create a simple TV remote
def markYes(df, gene):
    '''callback method used for turn_on button'''
    # use a Toplevel widget to display an image in a new window
    global centroid_idx
    plot_ROI(centroid_idx, fig, axes)
    centroid_idx += 1
    

    
# use Button and Label widgets to create a simple TV remote
def markNo(df, gene):
    '''callback method used for turn_on button'''
    # use a Toplevel widget to display an image in a new window
    global centroid_idx
    plot_ROI(centroid_idx, fig, axes)
    centroid_idx += 1
    

turn_on = Button(root, text="YES", command=markYes(df, gene))
turn_on.pack()

turn_off = Button(root, text="NO", command=markNo(df, gene))
turn_off.pack()


root.mainloop()


# %%

# 42 is bad
centroid_idx = 163

test1 = extract_plane_window(cell_centroids[centroid_idx,:], Snap25_images[0], 100)
plt.subplots()
plt.imshow(test1)
plt.show()

test2 = extract_plane_window(cell_centroids[centroid_idx,:], Snap25_images[1], 100)
plt.subplots()
plt.imshow(test2)

print(skimage.metrics.normalized_root_mse(test1, test2))

NRMSEs = np.zeros(len(cell_centroids))
for i in range(len(cell_centroids)):
    test1 = extract_plane_window(cell_centroids[i,:], Snap25_images[0], 100)
    test2 = extract_plane_window(cell_centroids[i,:], Snap25_images[1], 100)
    NRMSEs[i] = skimage.metrics.normalized_root_mse(test1, test2)
    

# %%

import imagej

fiji_path = r'C:\\Users\\jpv88\\Fiji.app'
IJ = imagej.init(fiji_path)

imp = IJ.openImage("E:/05-23_phox2b_retro_tracing/image_processing/slice_1/individual genes/GFP.tif");




    
    
    
        




