import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

colourMap = plt.cm.magma 

def add_texture(x, y, kx, ky, cond=2, angle=0, centre = (10,0)):
    TIMESTAMP = datetime.utcnow().strftime('%Y%m%d-%H%M%S-%f')
    angle = angle*np.pi/180

    x_rot = centre[0] + x*np.cos(angle) - y*np.sin(angle)
    y_rot = centre[1] + x*np.sin(angle) + y*np.cos(angle)

    z = cond*(0.5*(2+ np.sin(kx*x_rot) + np.sin(ky*y_rot)))

    print(np.max(z), np.min(z))
    plt.imshow(z, extent = [0,1,0,1])
    plt.colorbar()
    plt.savefig(TIMESTAMP+'png')
    return


x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)

X, Y = np.meshgrid(x, y)

angle = 45 

k = np.linspace(20, 0, num=5, endpoint=False)
print(k)




for kx, ky in zip(k, k):
    #print(kx)
    #ky = ky + 20
    fig = plt.figure(figsize=((5+1)*2, 5))
    plt.subplot(1, 2, 1)
    add_texture(X, Y, kx, ky, angle=0)
    plt.subplot(1, 2, 2)
    add_texture(X, Y, kx, ky, angle=angle)
    fig.tight_layout()
    plt.show