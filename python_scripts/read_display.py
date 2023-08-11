import os
import numpy as np
from scipy.io import loadmat, savemat
from datetime import datetime
import matplotlib.pyplot as plt
import random

def plot_2d(x, y, z, path, w = 10, nsamp = 100, random= False, circle = False): 
    
    batch_size = z.shape[0]

    if random:
        random.seed(0)
        randomlist = random.sample(range(batch_size), nsamp)
        name = 'random_samples'
        #print(randomlist)
    else:
        randomlist = range(0, batch_size) 
        name = 'samples'   

    figname = path+'/%s_conductivities_%s.png'%(nsamp, name)

    if os.path.isfile(figname): #if file exists, no need to proceed
        print('%s already exists'%(figname))
        return 0
    #batch_size = nsamp

    h = int(nsamp/w)

    fig = plt.figure(figsize=((5+1)*w, 5*h))

    for samp, i in zip(randomlist, range(nsamp)):
        if (i+1)%25 ==0:
            print("Plotting figure %s of %s : Sample %s"%(i+1, nsamp, samp))

        plt.subplot(h, w, i+1)
        if circle:
            plot_circle()
        colorbar = 'jet'
        plt.scatter(x, y, c = z[samp], cmap=colorbar, marker='.') # !!mistake you are plotting first 400 since z[i]. Use z[samp] for random samples
        plt.axis('square')
        plt.axis('off')
        plt.title(str(samp))
        plt.colorbar()

    fig.tight_layout()
    print('Saving figure\n')
    plt.savefig(figname)

    return 1

def plot_circle(R = [1, 0.90]):
    th = np.arange(0, 2 * np.pi, np.pi / 100)
    for r in R: # Add more values here if needed
        xunit = r * np.cos(th)
        yunit = r * np.sin(th)
        plt.plot(xunit, yunit, linestyle='dashed')
    return 


ending = "domain.mat"
path = '/pvfs2/Derick/EIT/Mine/data_texture'
runs = 0
for root, _, files in os.walk(path):
    for file in files:
        if ending in file: #file.endswith(ending):
            file_path = os.path.join(root, file)
            data_dom = loadmat(file_path)
            print('\nWorking on directory %s'%(root))
            cond = data_dom['inputConductivity']
            x    = data_dom['x1'               ]
            y    = data_dom['x2'               ]


            runs += plot_2d(x, y, cond, root)

print('\nTotal of %d runs completed'%(runs))
