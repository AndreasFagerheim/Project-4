# -*- coding: utf-8 -*-

    
import itertools
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

def Read_three_Column_File(file_name):
    with open(file_name, 'r') as data:
            MC_cycles = []
            E = []
            A = []
            M = []
            x = []
            for line in itertools.islice(data, 4, None):
                p = line.split()
                MC_cycles.append(float(p[0]))
                E.append(float(p[1]))
                M.append(float(p[2]))
                A.append(float(p[3]))
                
                #x.append(float(p[4]))
    return MC_cycles, E,M,A#cv, x

                                
plt.figure()
MC_cycles, E, M,A= Read_three_Column_File('C:\Data\configsR.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,A)


plt.grid('on')
plt.title("\n".join(wrap('Accpted configurations per Monte Carlo cycle, T = 1', 60)))

plt.ylabel('Accepted conifurations per cycle(random)')
plt.xlabel('MC cycles')
MC_cycles, E, M,A= Read_three_Column_File('C:\Data\configsU.txt')
plt.subplot(2, 1, 2)
plt.plot(MC_cycles, A)

plt.ylabel('Accepted conf /cycle (uniform)')
plt.xlabel('MC cycles')
plt.grid('on')

#plt.figure().tight_layout()

plt.savefig('acceptedConfigsT1.png')
plt.show()

plt.figure()
MC_cycles, E, M,A= Read_three_Column_File('C:\Data\configsR2.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,A)


plt.grid('on')
plt.title("\n".join(wrap('Accpted configurations per Monte Carlo cycle, T = 2.4', 60)))

plt.ylabel('Accepted conifurations per cycle(random)')
plt.xlabel('MC cycles')
MC_cycles, E, M,A= Read_three_Column_File('C:\Data\configsU2.txt')
plt.subplot(2, 1, 2)
plt.plot(MC_cycles, A)

plt.ylabel('Accepted conf /cycle (uniform)')
plt.xlabel('MC cycles')
plt.grid('on')

#plt.figure().tight_layout()

plt.savefig('acceptedConfigsT2.png')
plt.show()

