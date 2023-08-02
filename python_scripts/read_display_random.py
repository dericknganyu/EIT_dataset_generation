import os
import numpy as np
from scipy.io import loadmat, savemat
from datetime import datetime
import matplotlib.pyplot as plt
import random




def plot_2d(x, y, z, path, w = 10, nsamp = 100, name = 'conductivities', circle = False): 
    batch_size = z.shape[0]

    randomlist = random.sample(range(batch_size), nsamp)
    batch_size = nsamp

    h = int(batch_size/w)

    fig = plt.figure(figsize=((5+1)*w, 5*h))

    for samp, i in zip(randomlist, range(batch_size)):# range(batch_size):
        print("Plotting figure %s of %s : Sample %s"%(i+1, batch_size, samp))
        plt.subplot(h, w, i+1)
        if circle:
            plot_circle()
        colorbar = 'jet'
        plt.scatter(x, y, c = z[i], cmap=colorbar, marker='.')
        plt.axis('square')
        plt.axis('off')
        plt.colorbar()

    fig.tight_layout()
    print('Saving figure\n')
    plt.savefig(path+'/%s_samples_%s.png'%(nsamp, name))

def plot_circle(R = [1, 0.90]):
    th = np.arange(0, 2 * np.pi, np.pi / 100)
    for r in R: # Add more values here if needed
        xunit = r * np.cos(th)
        yunit = r * np.sin(th)
        plt.plot(xunit, yunit, linestyle='dashed')
    return



path = '/pvfs2/Derick/EIT/Mine/dataset'

for i in range(1,14):
    directory = path + '/part_%s'%(i)
    print('\nWorking on directory %s'%(directory))
    data_dom = loadmat(directory +'/dataset_domain.mat')
    cond = data_dom['inputConductivity']
    x    = data_dom['x1'               ]
    y    = data_dom['x2'               ]


    plot_2d(x, y, cond, directory, w = 20, nsamp = 400)#, circle = True, name = 'conductivities1')



# directory = '/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-12'
# data_dom = loadmat(directory +'/dataset_domain.mat')
# cond = data_dom['inputConductivity']
# x    = data_dom['x1'               ]
# y    = data_dom['x2'               ]

# plot_2d(x, y, cond, directory)

# plt.show()