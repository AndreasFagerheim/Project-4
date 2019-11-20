# -*- coding: utf-8 -*-

    
import itertools
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

def Read_three_Column_File(file_name):
    with open(file_name, 'r') as data:
            MC_cycles = []
            E = []
            cv = []
            M = []
            x = []
            for line in itertools.islice(data, 4, None):
                p = line.split()
                MC_cycles.append(float(p[0]))
                E.append(float(p[1]))
                M.append(float(p[2]))
                #cv.append(float(p[2]))
                
                #x.append(float(p[4]))
    return MC_cycles, E,M#cv, x

                                
plt.figure()
MC_cycles, E, M = Read_three_Column_File('C:\Data\data_BT1random1.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,np.abs(E))
plt.plot([0,len(E)],[1.9959820,1.9959820],"g--")

plt.grid('on')
plt.title("\n".join(wrap('Calculated Energy and Magnetization for a  2 x 2 lattice developing with increasing MC cycles.', 60)))

plt.ylabel('E/N')
plt.xlabel('MC cycles')

plt.subplot(2, 1, 2)
plt.plot(MC_cycles, M)
plt.plot([0,len(M)],[0.9986,0.9986],"g--")
plt.ylabel('|M|/N')
plt.xlabel('MC cycles')
plt.grid('on')

#plt.figure().tight_layout()

plt.savefig('N=2T=1.png')
plt.show()
#---------------------------------------------------------------------20x20, T = 1 random--------------------------------
plt.figure()
MC_cycles, E, M = Read_three_Column_File('C:\Data\C_N20T1random.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,np.abs(E))
plt.grid('on')
plt.title("\n".join(wrap('Calculated Energy and Magnetization for a  20 x 20 lattice. Inital random spin orientation, T = 1.', 60)))
plt.ylabel('E/N')
plt.xlabel('MC cycles')

plt.subplot(2, 1, 2)
plt.plot(MC_cycles, M)
plt.ylabel('|M|/N')
plt.xlabel('MC cycles')
plt.grid('on')

plt.savefig('N=20T=1R.png')
plt.show()


#---------------------------------------------------------------------20x20, T = 1 uniform--------------------------------
plt.figure()
MC_cycles, E, M = Read_three_Column_File('C:\Data\C_N20T1Uniform.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,np.abs(E))


plt.grid('on')
plt.title("\n".join(wrap('Calculated Energy and Magnetization for a  20 x 20. Inital uniform spin orientation, T = 1.', 60)))
plt.ylabel('E/N')
plt.xlabel('MC cycles')

plt.subplot(2, 1, 2)
plt.plot(MC_cycles, M)
plt.ylabel('|M|/N')
plt.xlabel('MC cycles')
plt.grid('on')

plt.savefig('N=20T=1U.png')
plt.show()
print(len(MC_cycles))

#---------------------------------------------------------------------20x20, T = 1 random--------------------------------
plt.figure()
MC_cycles, E, M = Read_three_Column_File('C:\Data\C_N20T2random.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,np.abs(E))
plt.grid('on')
plt.title("\n".join(wrap('Calculated Energy and Magnetization for a  20 x 20 lattice. Inital random spin orientation, T = 2.4.', 60)))
plt.ylabel('E/N')
plt.xlabel('MC cycles')

plt.subplot(2, 1, 2)
plt.plot(MC_cycles, M)
plt.ylabel('|M|/N')
plt.xlabel('MC cycles')
plt.grid('on')

plt.savefig('N=20T=2.4R.png')
plt.show()


#---------------------------------------------------------------------20x20, T = 1 uniform--------------------------------
plt.figure()
MC_cycles, E, M = Read_three_Column_File('C:\Data\C_N20T2uniform.txt')
plt.subplot(2, 1, 1)
plt.plot(MC_cycles,np.abs(E))


plt.grid('on')
plt.title("\n".join(wrap('Calculated Energy and Magnetization for a  20 x 20. Inital uniform spin orientation, T = 2.4.', 60)))
plt.ylabel('E/N')
plt.xlabel('MC cycles')

plt.subplot(2, 1, 2)
plt.plot(MC_cycles, M)
plt.ylabel('|M|/N')
plt.xlabel('MC cycles')
plt.grid('on')

plt.savefig('N=20T=2.4U.png')
plt.show()
print(len(MC_cycles))