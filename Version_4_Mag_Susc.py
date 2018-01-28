# -*- coding: utf-8 -*-
"""
@author: Sally Anne McCarthy
This one does magnetic susceptability


"""

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm


#Parameters
N = 50       #matrix size N rows and N columns
T = 1           #Temperature of the system. Will be adapting this
kb = 1
J=1
iter_no = N**3  

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

#Energy of the configuration
def total_energy(lis, N):
    tot_energy = 0
    for i in range(0, N-1):
        for j in range(0, N-1):            
            if i==0:
                left = lis[N-1, j]
            else:
                left = lis[i-1, j]
            if i == N-1:
                right = lis[0, j]
            else:
                right = lis[i+1, j]
            if j == 0:
                bottom = lis[i, N-1]
            else:
                bottom = lis[i, j-1]
            if j == N-1:
                top = lis[i, 0]
            else:
                top = lis[i, j+1] 
            a = lis[i, j]
            b = left + right + top + bottom
            tot_energy += -a*b
    return tot_energy * 0.25

def susceptability(m1, m2, T):
    a = (m1 - m2)**2
    return a / (T * kb) 


    
    


#This is the Ising function. It takes a random point in the matrix
#It calculates the energy change associated with that flip
#If that energy change is less than zero, it makes the flip
#if pflip is greater than a predefined constant, it also makes the flip
#it does this process N**3 times.
    


def Ising(N, lis, T):
    nstep = 0
    while nstep <= N**3:
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
         #   plt.matshow(lis, cmap=cm.coolwarm)          
        #if nstep == (0.2)*N**3:
        #    plt.matshow(lis, cmap=cm.coolwarm)
        #if nstep == (0.4)*N**3:
        #    plt.matshow(lis, cmap=cm.coolwarm)  
        #if nstep == (0.6)*N**3:
        #    plt.matshow(lis, cmap=cm.coolwarm)  
        #if nstep == (0.8)*N**3:
        #    plt.matshow(lis, cmap=cm.coolwarm)          
        nstep+=1
    #plt.matshow(lis, cmap=cm.coolwarm)    
    return lis
        
progress = 0
temps = np.linspace(0.1, 5, 20)
suscep = []

for elem in temps:
    print(progress)
    matr = initialise(N)
    magno1 = np.matrix.sum(matr)
    matrp = Ising(N, matr, elem)
    magno2 = np.matrix.sum(matrp)
    mag_sus = susceptability(magno1, magno2, elem)
    suscep.append(mag_sus)
    
#    if progress == 0:
#        plt.matshow(matr)
#    if progress == 2:
#        plt.matshow(matr)

#    if progress == 4:
#        plt.matshow(matr)

#    if progress == 6:
#        plt.matshow(matr)

#    if progress == 8:
#        plt.matshow(matr)
        
#    if progress == 10:
#        plt.matshow(matr)

    progress += 1
        
        


plt.scatter(temps, suscep)
plt.xlabel('Temperature')
plt.ylabel('Magnetic Susceptability after N**3 iterations')

plt.show()