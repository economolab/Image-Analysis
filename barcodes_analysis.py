# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:16:30 2023

@author: jpv88
"""

from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import JV_utils

path = r'C:\\Users\\jpv88\\Documents\\GitHub\\Image-Analysis\\annotation_gui\\'
file = 'annotations_s02.csv'

barcodes_full = pd.read_csv(path + file)
barcodes = pd.DataFrame()

file2 = 's01_annotations.csv'
barcodes_full2 = pd.read_csv(path+file2)
barcodes2 = pd.DataFrame()


for (columnName, columnData) in barcodes_full.items():
    if ('checked' in columnName) and all(columnData.to_numpy()):
        gene_name1 = '_'.join(columnName.split('_')[:2])
        gene_name2 = gene_name1.split('_')[0]
        barcodes[gene_name2] = barcodes_full[gene_name1].to_numpy()
        
for (columnName, columnData) in barcodes_full2.items():
    if ('checked' in columnName) and all(columnData.to_numpy()):
        gene_name1 = '_'.join(columnName.split('_')[:2])
        gene_name2 = gene_name1.split('_')[0]
        barcodes2[gene_name2] = barcodes_full2[gene_name1].to_numpy()
        
barcodes_Phox2b = barcodes[barcodes['Phox2b'] == 1]

barcodes_Phox2b = barcodes_Phox2b[['Dlk1', 'Slc32a1']]
        
# %%

import sklearn

import numpy as np

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(metric='euclidean', linkage='ward', 
                                n_clusters=4, compute_distances=True)

model = model.fit(barcodes_Phox2b.to_numpy())

# %%

barcodes_Phox2b_array = barcodes_Phox2b.to_numpy()
labels = model.labels_

barcodes_sets = []
for val in np.unique(labels):
    barcodes_set = barcodes_Phox2b_array[labels == val]
    barcodes_sets.append(barcodes_set)
    
heat_data = np.vstack(barcodes_sets)

fig, ax = plt.subplots()
plt.pcolormesh(np.rot90(heat_data), cmap='RdYlBu')

mpl.rcParams['image.composite_image'] = False
plt.rcParams['svg.fonttype'] = 'none'

# %%

import numpy as np
from sklearn.manifold import TSNE


model_TSNE = TSNE(n_components=2, perplexity=1, learning_rate='auto',
                init='random', n_iter=10000)

X_embedded = model_TSNE.fit_transform(barcodes_Phox2b.to_numpy())

fig, ax = plt.subplots()
for i, color in enumerate(['r', 'g', 'b']):
    x, y = X_embedded[:,0][model.labels_ == i], X_embedded[:,1][model.labels_ == i]
    ax.scatter(x, y, c=color, s=30, label=color)
    
# plt.scatter(X_embedded[:,0], X_embedded[:,1], c=model.labels_, s=30)
ax.legend()

test2 = ['Dlk1-/slc-', 'Dlk1-/Slc+', 'Dlk1+/Slc-']

mpl.rcParams['image.composite_image'] = False
plt.rcParams['svg.fonttype'] = 'none'

# %%

import sklearn

distance = sklearn.metrics.pairwise_distances(barcodes.to_numpy(),
                                              metric='manhattan')
      
clustering = sklearn.cluster.AffinityPropagation(affinity='precomputed')

clustering.fit(distance)
        
num_per_label = []
for label in np.unique(clustering.labels_):
    num_per_label.append(sum(clustering.labels_ == label))
        
        
# %% total expression level per gene

genes = barcodes.columns.tolist()

cmap = mpl.cm.viridis

reads = np.sum(barcodes, axis=0)
reads = reads/len(barcodes)
reads = reads.tolist()

genes_sorted = JV_utils.sort_list_by_list(reads, genes)
reads = sorted(reads)

genes_sorted = list(reversed(genes_sorted))
reads = list(reversed(reads))

fig, ax = plt.subplots()

colors = np.array(reads)*100/100
bars = plt.bar(genes_sorted, np.array(reads)*100, zorder=2, edgecolor='k', 
               linewidth=2, color=cmap(colors))
plt.setp(ax.get_xticklabels(), rotation=40, horizontalalignment='right', 
         rotation_mode="anchor")

x = []
for bar in bars:
    x.append(bar.xy[0])

plt.yticks(fontsize=14)
plt.xticks(fontsize=11)
plt.ylabel('% Expression', fontsize=16)
plt.title('Expression Levels in GFP Tagged Neurons', fontsize=18)
plt.grid(axis='y', which='both', ls='--', alpha=1, lw=0.2, zorder=1)

plt.ylim(0, 100)

# Create offset transform by 5 points in x direction
dx = 5/72.; dy = 0/72. 
offset = mpl.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
for label in ax.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)
    
    
plt.tight_layout()

# %%
import numpy as np

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    R = dendrogram(linkage_matrix, **kwargs)
    return R, linkage_matrix

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(metric='manhattan', linkage='average', 
                                n_clusters=5, compute_distances=True)

model = model.fit(barcodes.to_numpy())
plt.title("Hierarchical Clustering Dendrogram")
# plot the top three levels of the dendrogram
R, linkage_matrix = plot_dendrogram(model, truncate_mode="level")
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()

# %%

test = (model.labels_ == 4)
test2 = barcodes[test]


# %%

from sklearn.cluster import DBSCAN
import numpy as np

clustering = DBSCAN(eps=0.9, metric='manhattan').fit(barcodes.to_numpy())
clustering.labels_

# %%

from sklearn.ensemble import RandomForestClassifier
from sklearn.cluster import KMeans

importances = []
for _ in tqdm(range(10)):
    n_clusters = np.arange(2, 10, 1)
    for n in n_clusters:
        kmeans = KMeans(n_clusters=n).fit(barcodes)
        labels = kmeans.labels_
        forest = RandomForestClassifier()
        forest.fit(barcodes, labels)
        importances.append(forest.feature_importances_)
    
importances = np.stack(importances)
imp_means = np.mean(importances, axis=0)

cmap = mpl.cm.viridis
colors = imp_means/max(imp_means)

fig, ax = plt.subplots()
plt.bar(genes, imp_means, color=cmap(colors))
plt.errorbar(genes, imp_means, yerr=np.std(importances, axis=0), fmt='none', 
             color="k", capsize=2)
plt.ylabel('Mean Decrease in Impurity (MDI)', fontsize=13)
plt.title('Relative Importance of Gene in Random Forest Classification',  
          fontsize=13)

plt.setp(ax.get_xticklabels(), rotation=40, horizontalalignment='right', 
         rotation_mode="anchor")

plt.yticks(fontsize=14)
plt.xticks(fontsize=11)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

plt.ylim(0)