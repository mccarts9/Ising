# -*- coding: utf-8 -*-
"""
@author: Sally Anne McCarthy
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
        nstep+=1
    return lis
        
progress = 0
temps = np.linspace(0.1, 5.1, 10)
magnets = []

for elem in temps:
    print(progress)
    matr = initialise(N)
    matr = Ising(N, matr, elem)
    magspin = avgmagnet(matr, N)
    magspin = np.abs(magspin)
    magnets.append(magspin)
    if progress == 0:
        plt.matshow(matr, cmap=cm.coolwarm)
    if progress == 2:
        plt.matshow(matr, cmap=cm.coolwarm)

    if progress == 4:
        plt.matshow(matr, cmap=cm.coolwarm)

    if progress == 6:
        plt.matshow(matr, cmap=cm.coolwarm)

    if progress == 8:
        plt.matshow(matr, cmap=cm.coolwarm)
        
    if progress == 10:
        plt.matshow(matr, cmap=cm.coolwarm)

    progress += 1
        
        
        
print(len(temps))
print(len(magnets))

plt.figure()
plt.scatter(temps, magnets)
plt.xlabel('Temperature')
plt.ylabel('Average Magnetism per spin after N**3 iterations')

