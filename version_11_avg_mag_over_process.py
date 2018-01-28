# -*- coding: utf-8 -*-
"""


@author: Sally Anne McCarthy
This captures the process of the average magnetism thing
"""




import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm


#Parameters
N = 30      #matrix size N rows and N columns
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
#it does this process N**3 times.     
    
nsteps_list = []
mags_list = []
energies_list = []
mag_sus_1_list = []
def Ising(N, lis, T):
    nstep = 0
    while nstep <= N**5:
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
       
        #sys_mag = avgmagnet(lis, N)
        #mags_list.append(sys_mag)
        sys_energy = total_energy(lis, N)
        energies_list.append(sys_energy)
        a=np.matrix.sum(lis)
        #mag_sus_1_list.append(a)
        nstep+=1
        nsteps_list.append(nstep)
        if np.abs(a) == N**2:
            break
    return lis
        
progress = 0
matr1 = initialise(N)
m0 = np.matrix.sum(matr1)
e0 = total_energy(matr1, N)
matr2 = Ising(N, matr1, T)

specific_heat_array = []
for element in energies_list:
    specific_heat = spec_heat(e0, element, T)
    specific_heat_array.append(specific_heat)
    

susceptibility_list=[]
for elem in mag_sus_1_list:
    y = susceptability(m0, elem, T)
    susceptibility_list.append(y)
   
plt.figure()
plt.scatter(nsteps_list, specific_heat_array)   
plt.xlabel('Iterations of Ising Model Loop')
plt.ylabel('Specific Heat')
plt.title('Evolution of Specific Heat over Ising Model')
plt.show() 
    
    
    
"""
plt.figure()
plt.scatter(nsteps_list, susceptibility_list)
plt.xlabel('Iterations of Ising Model Loop')
plt.ylabel('Magnetic Susceptibility')
plt.title('Evolution of Magnetic Susceptibility over Ising Model')
plt.show()


plt.figure()
plt.scatter(nsteps_list, mags_list)
plt.xlabel('Iterations of Ising Model Loop')
plt.ylabel('Average Magnetism per element')
plt.title('Evolution of Average Magnetism per particle over Ising Model')

plt.figure()
plt.scatter(nsteps_list, energies_list)
plt.xlabel('Iterations of Ising Model Loop')
plt.ylabel('System Energy')
plt.title('Evolution of System Energy over Ising Model')

"""
