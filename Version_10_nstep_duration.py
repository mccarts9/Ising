# -*- coding: utf-8 -*-
"""
Sample Size vs nsteps to duration
Laptop too unreliable to just do run times

@author: Sally Anne McCarthy
"""
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm
import timeit


#Parameters
amount = np.linspace(2, 31, 30)    #matrix size N rows and N columns
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
def energychange(lis, i, j, N):  
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

def spec_heat(e1, e2, T):
    a = (e1 - e2)**2
    return a / (T * T *  kb)   
    



#This is the Ising function. It takes a random point in the matrix
#It calculates the energy change associated with that flip
#If that energy change is less than zero, it makes the flip
#if pflip is greater than a predefined constant, it also makes the flip
#it does this process N**3.5 times.
#it stops if the matrix is homogenous 
totals = []

def Ising(N, lis, T):
    nstep = 0
    while nstep <= N**8:
        flip_constant = random.random()
        random_x = random.choice(np.arange(0, N))
        random_y = random.choice(np.arange(0, N))
        dU = energychange(lis, random_x, random_y, N)
        if dU <= 0:
            lis = randspin(lis, random_x, random_y)
        else:
            pflip = np.exp(-dU/(T*kb))
            if pflip >= flip_constant:
                    lis = randspin(lis, random_x, random_y)        
        nstep+=1
        a = np.abs(np.matrix.sum(lis))
        if a == N**2:
            break
    totals.append(nstep)
    #plt.matshow(lis, cmap=cm.coolwarm)    
    return lis
        
progress = 0
temps = np.linspace(0.1, 12.6, 100)




times = []
prog = 0    
for i in amount:
    i = np.int(i)
    matr = initialise(i)
    matrp = Ising(i, matr, T)
    prog +=1
    print('I have ', prog, 'system sizes completed out of ', len(amount))


plt.figure(1)
plt.scatter(amount, totals)
plt.xlabel('System Size -Length of one row')
plt.ylabel('NSteps for Ising Model to have equal spin states.')
plt.title('System Size and Duration of Ising Model iterations at T = 1')



log_sizes = np.log10(amount)
log_total = np.log10(totals)

coef = np.lib.polyfit(log_sizes, log_total, 1) 
line = np.lib.polyval(coef, log_sizes) 
slope = str(coef[0])
intercept = str(coef[1])

log_label = 'Y = ' + slope + 'X' +  intercept
log_label = str(log_label)


plt.figure(2)
plt.plot(log_sizes, log_total, 'ro', label = 'Log of iterations of Ising Model')
plt.plot(log_sizes, line, 'b--', label = log_label)
plt.legend(loc=2)
plt.xlabel('Log of Length of one row')
plt.ylabel('Log of Nstep duration of Ising Model')
plt.title('Log Log Plot of System Size and Duration of Ising Model at T = 1')

plt.show()