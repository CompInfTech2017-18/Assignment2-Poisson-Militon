# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 01:38:24 2018

@author: kazy
"""
import pylab
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

P=np.pi


a=100
V=100
interval=0.1
inf=100


def Uk(x,y,a,V,inf):
    U=0
    k=0
    for i in range(inf):
        U+=((4*V)/P*(np.sinh((2*k+1)*(a-x)*P/a)*np.sin((2*k+1)*P*y/a))/((2*k+1)*np.sinh((2*k+1)*P)))
        k=k+1    
    return U


X=np.arange(0,a,interval)
Y=np.arange(0,a,interval)
xgrid, ygrid = np.meshgrid(X,Y)

U=Uk(xgrid,ygrid,a,V,inf)

eps=0.001

N=100



class Solve:
    def __init__(self,N,M,):
        self.N=N
        self.M=M
        self.Unew=np.zeros((N,M))
        self.Uold=np.ones((N,M))
        
    def Jacoby(self):
        while (np.absolute(np.trace(self.Unew)-np.trace(self.Uold))>eps):
            self.Uold=self.Unew
            self.Unew = (                 self.Uold[0:-2,1:-1] +
                         self.Uold[1:-1,0:-2]                  + self.Uold[1:-1,2:] +
                                     self.Uold[2:  ,1:-1])/4
            self.Unew = np.insert(self.Unew,0,0,axis=0)
            self.Unew = np.insert(self.Unew,N-1,0,axis=0)
            self.Unew = np.insert(self.Unew,0,100,axis=1)
            self.Unew = np.insert(self.Unew,N-1,0,axis=1)
        return self.Unew

    def Gauss(self):
        U1 = 0
        U2 = 1
        while (np.absolute(U1 - U2) > eps):
            U1 = np.trace(self.Unew)
            self.Unew = (                       self.Unew[0:-2, 1:-1] +
                         self.Unew[1:-1, 0:-2]                       + self.Unew[1:-1, 2:] +
                                           self.Unew[2:, 1:-1]) / 4
            self.Unew = np.insert(self.Unew, 0, 0, axis=0)
            self.Unew = np.insert(self.Unew, N - 1, 0, axis=0)
            self.Unew = np.insert(self.Unew, 0, 100, axis=1)
            self.Unew = np.insert(self.Unew, N - 1, 0, axis=1)
            U2 = np.trace(self.Unew)
        return self.Unew

    def Relaxation(self,w):
        U1 = 0
        U2 = 1
        # for i in range(100):
        while (np.absolute(U1 - U2) > eps):
            U1 = np.trace(self.Unew)
            Uold1 = self.Unew
            Uold1 = np.delete(Uold1, 0, 0)
            Uold1 = np.delete(Uold1, N - 2, 0)
            Uold1 = np.delete(Uold1, 0, 1)
            Uold1 = np.delete(Uold1, N - 2, 1)
            R = (self.Unew[0:-2, 1:-1] +
                 self.Unew[1:-1, 0:-2] + self.Unew[1:-1, 2:] +
                 self.Unew[2:, 1:-1]) / 4
            self.Unew = Uold1 + w * (R - Uold1)
            self.Unew = np.insert(self.Unew, 0, 0, axis=0)
            self.Unew = np.insert(self.Unew, N - 1, 0, axis=0)
            self.Unew = np.insert(self.Unew, 0, 100, axis=1)
            self.Unew = np.insert(self.Unew, N - 1, 0, axis=1)
            R = np.insert(R, 0, 0, axis=0)
            R = np.insert(R, N - 1, 0, axis=0)
            R = np.insert(R, 0, 100, axis=1)
            R = np.insert(R, N - 1, 0, axis=1)
            U2 = np.trace(self.Unew)
            #print(Unew)
        return self.Unew

    def plotter(self):
        X1=np.linspace(0,,100)
        Y1=np.linspace(0,100,100)
        xgrid1, ygrid1 = np.meshgrid(X1,Y1)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
CS = plt.contour(xgrid, ygrid, U,colors=('indigo','purple','b','m','violet','aqua'), linewidths=0.8)
ax.plot_wireframe(xgrid, ygrid, U,color='black',linewidth=0.3)

pylab.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
k=Solve(100,100)
j=k.Gauss()
CS1 = plt.contour(xgrid1, ygrid1, j,colors=('indigo','purple','b','m','violet','aqua'), linewidths=0.8)
ax1.plot_wireframe(xgrid1, ygrid1, j,color='black',linewidth=0.3)

pylab.show()



    
    
    
 



        

            
        

    