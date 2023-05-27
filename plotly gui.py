# -*- coding: utf-8 -*-
"""
Created on Sun May 21 15:18:09 2023

@author: jpv88
"""

import skimage
import plotly

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

import plotly.io as io
io.renderers.default='browser'

path = r'E:\\05-23_phox2b_retro_tracing\\image_processing\\slice_1\\individual genes\\'
fname = 'GFP_r1.tif'

image = skimage.io.imread(path+fname)

img = image[4,:,:]

img_rgb = np.array([[[255, 0, 0], [0, 255, 0], [0, 0, 255]],
                    [[0, 255, 0], [0, 0, 255], [255, 0, 0]]
                   ], dtype=np.uint8)

config = dict({'scrollZoom':True})

fig = make_subplots(rows=1, 
                    cols=2,
                    shared_xaxes='all',
                    shared_yaxes='all')

fig.add_trace(go.Heatmap(z=img), 1, 1)
fig.add_trace(go.Heatmap(z=img), 1, 2)
fig.update_traces(showscale=False)

fig.update_xaxes(range=[1100, 1200], row=1, col=1)
fig.update_yaxes(range=[2100, 2200], row=1, col=1)

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"title": "Slider switched to step: " + str(i)}],  # layout attribute
    )
    steps.append(step)

sliders = [dict(
    active=10,
    currentvalue={"prefix": "Frequency: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(sliders=sliders)

fig.show(config=config)






