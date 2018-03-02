# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 01:38:24 2018

@author: kazy
"""
import imageio
import pylab
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


class Solve:
    def __init__(self,N,M):
        self.eps = 0.001
        self.P = np.pi
        self.N=N
        self.M=M
        self.Unew=np.zeros((M,N))
        self.Uold=np.ones((M,N))

    def Fourie(self):
        x=self.xgrid
        y=self.ygrid
        a = 100
        V = 100
        inf = 100
        U = 0
        k = 0
        for i in range(inf):
            U += ((4 * V) / self.P * (np.sinh((2 * k + 1) * (a - x) * self.P / a) * np.sin((2 * k + 1) * self.P * y / a)) / (
                        (2 * k + 1) * np.sinh((2 * k + 1) * self.P)))
            k = k + 1
        return U

    def Jacoby(self):
        for i in range(100000):
        #while (np.absolute(np.trace(self.Unew)-np.trace(self.Uold))>self.eps):
            self.Uold=self.Unew
            self.Unew = (                 self.Uold[0:-2,1:-1] +
                         self.Uold[1:-1,0:-2]                  + self.Uold[1:-1,2:] +
                                     self.Uold[2:  ,1:-1])/4
            self.Unew = np.insert(self.Unew,0,0,axis=0)
            self.Unew = np.insert(self.Unew,self.M-1,0,axis=0)
            self.Unew = np.insert(self.Unew,0,100,axis=1)
            self.Unew = np.insert(self.Unew,self.N-1,-100,axis=1)
        return self.Unew

    def Gauss(self):
        U1 = 0
        U2 = 1
        while (np.absolute(U1 - U2) > self.eps):
            U1 = np.trace(self.Unew)
            self.Unew = (                       self.Unew[0:-2, 1:-1] +
                         self.Unew[1:-1, 0:-2]                       + self.Unew[1:-1, 2:] +
                                           self.Unew[2:, 1:-1]) / 4
            self.Unew = np.insert(self.Unew, 0, 0, axis=0)
            self.Unew = np.insert(self.Unew, self.M - 1, 0, axis=0)
            self.Unew = np.insert(self.Unew, 0, 100, axis=1)
            self.Unew = np.insert(self.Unew, self.N - 1, 0, axis=1)
            U2 = np.trace(self.Unew)
        return self.Unew

    def Relaxation(self,w):
        U1 = 0
        U2 = 1
        # for i in range(100):
        while (np.absolute(U1 - U2) > self.eps):
            U1 = np.trace(self.Unew)
            Uold1 = self.Unew
            Uold1 = np.delete(Uold1, 0, 0)
            Uold1 = np.delete(Uold1, self.M - 2, 0)
            Uold1 = np.delete(Uold1, 0, 1)
            Uold1 = np.delete(Uold1, self.N - 2, 1)
            R = (self.Unew[0:-2, 1:-1] +
                 self.Unew[1:-1, 0:-2] + self.Unew[1:-1, 2:] +
                 self.Unew[2:, 1:-1]) / 4
            self.Unew = Uold1 + w * (R - Uold1)
            self.Unew = np.insert(self.Unew, 0, 0, axis=0)
            self.Unew = np.insert(self.Unew, self.M - 1, 0, axis=0)
            self.Unew = np.insert(self.Unew, 0, 100, axis=1)
            self.Unew = np.insert(self.Unew, self.N - 1, 0, axis=1)
            R = np.insert(R, 0, 0, axis=0)
            R = np.insert(R, self.M - 1, 0, axis=0)
            R = np.insert(R, 0, 100, axis=1)
            R = np.insert(R, self.N - 1, 0, axis=1)
            U2 = np.trace(self.Unew)
            #print(Unew)
        return self.Unew

    def plotter(self, U):
        X=np.linspace(0,self.N,self.N)
        Y=np.linspace(0,self.M,self.M)
        self.xgrid, self.ygrid = np.meshgrid(X,Y)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        CS = plt.contour(self.xgrid, self.ygrid, U, colors=('indigo', 'purple', 'b', 'm', 'violet', 'aqua'), linewidths=0.8)
        ax.plot_wireframe(self.xgrid, self.ygrid, U, color='black', linewidth=0.3)
        for angle in range(0, 360):
            ax.view_init(30, angle)
            plt.draw()
            plt.pause(.001)
        pylab.show()

    def run(self,method):
        if method=='jacoby':
            solve=self.Jacoby()
            self.plotter(solve)
        elif method=='gauss':
            solve=self.Gauss()
            self.plotter(solve)
        elif method=='relax':
            solve=self.Relaxation(1)
            self.plotter(solve)
        elif method=='fourie':
            solve=self.Fourie()
            self.plotter(solve)
        elif method=='conductor':
            solve=self.conductor()
            self.plotter(solve)


    def cond_read(self, file_name):
        self.im = imageio.imread(file_name)
        self.N=self.im.shape[1]
        self.M=self.im.shape[0]
        self.Ucond1 = np.zeros((self.M ,self.N ))
        black = [0, 0, 0]
        blue = [11, 210, 239]
        self.pot_blue=0
        red = [255, 0, 0]
        self.pot_red=100
        green = [11, 148, 1]
        self.pot_green = -100
        i = 0
        j = 0
        for y in self.im:
            for x in y:
                if np.array_equal(x, blue):
                    self.Ucond1[i, j] = self.pot_blue
                elif np.array_equal(x, red):
                    self.Ucond1[i, j] = self.pot_red
                elif np.array_equal(x, green):
                    self.Ucond1[i, j] = self.pot_green
                #print(x)
                #print(j)
                j += 1
            i += 1
            j = 0
        return self.Ucond1


    def conductor(self):
        self.Ucond1 = self.cond_read('potential.png')
        self.Ucond = self.Ucond1
        self.plotter(self.Ucond)
        U1 = 0
        U2 = 1
        for t in range(1000):
            print(t)
        #while (np.absolute(U1 - U2) > self.eps):
            U1 = np.trace(self.Ucond)
            self.Ucond = (self.Ucond[0:-2, 1:-1] +
                         self.Ucond[1:-1, 0:-2] + self.Ucond[1:-1, 2:] +
                         self.Ucond[2:, 1:-1]) / 4
            self.Ucond = np.insert(self.Ucond, 0, 0, axis=0)
            self.Ucond = np.insert(self.Ucond, self.M - 1, 0, axis=0)
            self.Ucond = np.insert(self.Ucond, 0, 0, axis=1)
            self.Ucond = np.insert(self.Ucond, self.N - 1, 0, axis=1)
            #if t%100==0:
                #self.plotter(self.Ucond)
            i = 0
            j = 0
            for y in self.Ucond1:
                for x in y:
                    if (x==100):
                        self.Ucond[i, j] = 100
                    elif (x==-100):
                        self.Ucond[i, j] = -100
                    j += 1
                i += 1
                j = 0

            U2 = np.trace(self.Ucond)

        return self.Ucond


#Resh=Solve(20,1000)
#Resh.run('jacoby')
conductor=Solve(100,100)
conductor.run('conductor')
#Resh.run('gauss')
#Resh.run('fourie')
#Resh.run('relax')















