# -*- coding: utf-8 -*-
"""
@author: Sally Anne McCarthy
@This code is the simple ising model. It creates and then orders a random matrix 
@Matrix is NxN big, and this model does N**3 flips
"""


import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm
import timeit


#Parameters
start = timeit.default_timer()
N = 50         #matrix size N rows and N columns
T = 1        #Temperature of the system. Will be adapting this
kb = 1
h = 0 
J = 1
#A function that produces the initial matrix of N x N elements
#where each element is randomly 1 or -1
def initialise(N):
    seq = [-1, 1]
    lis = []
    for i in range (0, N):
        b = random.choices(seq, k = N)
        lis.append(b)
    lis_matrix = np.asmatrix(lis)    
    return lis_matrix

#Need to randomly spin one of them
#This function returns the same matrix with one random value flipped
def randspin(lis, x, y):
    if lis[x, y] == 1:
        lis[x, y] = -1
    elif lis[x, y] == -1:
        lis[x, y] = 1
    else:
        print('Something has gone wrong. Impossible value at ', x, y)
    return lis    


#This function calculates the energy change associated with flipping the value at (i,j)
def energychange(lis, i, j): 
    
    if i==0:
        lefts = lis[N-1, j]
    else:
        lefts = lis[i-1, j]
    if i == N-1:
        rights = lis[0, j]
    else:
        rights = lis[i+1, j]
    if j == 0:
        bottoms = lis[i, N-1]
    else:
        bottoms = lis[i, j-1]
    if j == N-1:
        tops = lis[i, 0]
    else:
        tops = lis[i, j+1] 
    return 2 * lis[i, j] * (lefts + rights + tops + bottoms) * J


#Average magnetism per element
def avgmagnet(lis, N):
    a = np.matrix.sum(lis)
    b = a / (N**2)
    return b

#This is the Ising function. It takes a random point in the matrix
#It calculates the energy change associated with that flip
#If that energy change is less than zero, it makes the flip
#if pflip is greater than a predefined constant, it also makes the flip
#it does this process N**3 times.     

mag_sus=[]
def Ising(N, lis, T):
    nstep = 0
    
    while nstep <= N**3.5:
        flip_constant = random.random()
        random_x = random.choice(np.arange(0, N))
        random_y = random.choice(np.arange(0, N))
        dU = energychange(lis, random_x, random_y)
        if dU <= 0:
            lis = randspin(lis, random_x, random_y)
        else:
            pflip = np.exp(-dU/(T*kb))
            if pflip >= flip_constant:
                    lis = randspin(lis, random_x, random_y)
                    
        #if nstep == 0:
        #    plt.matshow(lis, cmap=cm.coolwarm)          
        #if nstep == (0.2)*N**3:
         #   plt.matshow(lis, cmap=cm.coolwarm)
        #if nstep == (0.4)*N**3:
         #   plt.matshow(lis, cmap=cm.coolwarm)  
        #if nstep == (0.6)*N**3:
         #   plt.matshow(lis, cmap=cm.coolwarm)  
        #if nstep == (0.8)*N**3:
         #   plt.matshow(lis, cmap=cm.coolwarm)          
        nstep+=1
        #print(nstep)
        a = np.abs(np.matrix.sum(lis))
        if a == N**2:
            break
    plt.matshow(lis, cmap=cm.coolwarm)    
    
    return lis


nstep = 0
testmatr = initialise(N)
magno1 = np.matrix.sum(testmatr)
productmatr = Ising(N, testmatr, 0.1)
magno2 = np.matrix.sum(productmatr)
magsus = (magno2 - magno1) / kb * T
mag = avgmagnet(productmatr, N)
#print(np.sum(mag_sus))
#print(len(mag_sus)) 
plt.show()  
stop = timeit.default_timer()
print(stop - start)    