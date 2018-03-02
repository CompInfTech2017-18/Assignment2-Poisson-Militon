import imageio
import numpy as np
from matplotlib import lines as line
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from matplotlib.lines import Line2D

im = imageio.imread('potential.png')
k=np.zeros((100,100))
black=[0,0,0]
blue=[11,210,239]
red=[255,0,0]
green=[11,148,1]
i=0
j=0
for y in im:
    for x in y:
        if np.array_equal(x, blue):
            k[i,j]=1
        elif np.array_equal(x,red):
            k[i,j]=100
        elif np.array_equal(x,green):
            k[i,j]=-100
        print(x)
        print(j)
        j+=1
    i += 1
    j = 0

print(im)
print(im.shape[0])
print(k)

#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.scatter(im, marker='.')