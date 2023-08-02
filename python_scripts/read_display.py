import os
import numpy as np
from scipy.io import loadmat, savemat
from datetime import datetime
import matplotlib.pyplot as plt
import random




def plot_2d(x, y, z, path, l = 10): 
    batch_size = z.shape[0]

    randomlist = random.sample(range(batch_size), 5)

    w = int(batch_size/l)

    fig = plt.figure(figsize=((5+1)*l, 5*w))

    for i in randomlist:# range(batch_size):
        print("Plotting figure %s of %s"%(i+1, batch_size))
        plt.subplot(l,w,i+1)
        colorbar = 'jet'
        plt.scatter(x, y, c = z[i], cmap=colorbar, marker='.')
        plt.axis('square')
        plt.axis('off')
        plt.colorbar()

    fig.tight_layout()
    print('Saving figure\n')
    plt.savefig(path+'/conductivities.png')


path = '/pvfs2/Derick/EIT/Mine/dataset'

for i in range(1,14):
    directory = path + '/part_%s'%(i)
    print('\nWorking on directory %s'%(directory))
    data_dom = loadmat(directory +'/dataset_domain.mat')
    cond = data_dom['inputConductivity']
    x    = data_dom['x1'               ]
    y    = data_dom['x2'               ]


    plot_2d(x, y, cond, directory)



# directory = '/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-12'
# data_dom = loadmat(directory +'/dataset_domain.mat')
# cond = data_dom['inputConductivity']
# x    = data_dom['x1'               ]
# y    = data_dom['x2'               ]

# plot_2d(x, y, cond, directory)

# plt.show()