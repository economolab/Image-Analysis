# -*- coding: utf-8 -*-
"""
Created on Sat May 20 18:25:03 2023

@author: jpv88
"""

# %%

path = r'E:\\05-23_phox2b_retro_tracing\\image_processing\\slice_1\\individual genes\\'

GFP_image = load_image('GFP', path)
Snap25_GFP_image = load_Snap25_image('GFP', path)

gene = 'Dlk1'
gene_image = load_image(gene, path)
Snap25_gene_image = load_Snap25_image(gene, path)

fig, ax = plt.subplots(2,2)

centroid_idx = 10
window_sz = 100

gene_plane = extract_plane_window(cell_centroids[centroid_idx,:], gene_image, window_sz)
Snap25_gene_plane = extract_plane_window(cell_centroids[centroid_idx,:], Snap25_gene_image, window_sz)
Snap25_GFP_plane = extract_plane_window(cell_centroids[centroid_idx,:], Snap25_GFP_image, window_sz)
GFP_plane = extract_plane_window(cell_centroids[centroid_idx,:], GFP_image, window_sz)

ax[0,0].imshow(skimage.exposure.adjust_log(gene_plane), cmap='hot')
ax[0,1].imshow(skimage.exposure.adjust_log(Snap25_gene_plane), cmap='hot')
ax[1,0].imshow(skimage.exposure.adjust_log(Snap25_GFP_plane), cmap='hot')
ax[1,1].imshow(skimage.exposure.adjust_log(GFP_plane), cmap='hot')

ax[0,0].set_title(str(cell_centroids[centroid_idx, :]))


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




    
    
    
        
