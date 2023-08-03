import os
import numpy as np
from scipy.io import loadmat, savemat
from datetime import datetime
import matplotlib.pyplot as plt
import random


def plot_2d(x, y, z, path, w = 10, nsamp = 100, random= False, circle = False): 
    
    batch_size = z.shape[0]
    mini, maxi = 0.2, 2

    if random:
        random.seed(0)
        randomlist = random.sample(range(batch_size), nsamp)
        name = 'random_samples'
        #print(randomlist)
    else:
        randomlist = range(0, batch_size) 
        name = 'samples'   

    figname = path+'/%s_00conductivities_%s.png'%(nsamp, name)

    if os.path.isfile(figname): #if file exists, no need to proceed
        print('%s already exists'%(figname))
        return 0
    #batch_size = nsamp

    h = int(nsamp/w)

    
    fig, axes = plt.subplots(nrows=h, ncols=w, figsize=((5+1)*w, 5*h), constrained_layout=True)
    # fig = plt.figure(figsize=((5+1)*w, 5*h))

    for samp, i, ax in zip(randomlist, range(nsamp), axes.flat):
        if nsamp%(i+1) ==0:
            print("Plotting figure %s of %s : Sample %s"%(i+1, nsamp, samp))

        plt.subplot(h, w, i+1)
        if circle:
            plot_circle()
        colorbar = 'jet'
        # plt.scatter(x, y, c = z[samp], cmap=colorbar, marker='.') # !!mistake you are plotting first 400 since z[i]. Use z[samp] for random samples
        im = ax.scatter(x, y, c = z[samp], cmap=colorbar, marker='.', vmin=mini, vmax=maxi)#imshow(np.random.random((10,10)), vmin=0, vmax=1)
        
        ax.axis('square')
        ax.axis('off')
        plt.title(str(samp))
        #im.tight_layout()
        
    fig.colorbar(im, ax=axes.ravel().tolist())
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
path = '/pvfs2/Derick/EIT/Mine/data'
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
            break
        #break
    #break    
print('\nTotal of %d runs completed'%(runs))



import numpy as np
import matplotlib.pyplot as plt

