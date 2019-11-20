import itertools
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

def Read_eight_Column_File(file_name):
    with open(file_name, 'r') as data:
            MC_cycles = []
            
            A = []
           
            P = []
            for line in itertools.islice(data, 4, None):
                p = line.split()
                MC_cycles.append(float(p[0]))
                A.append(float(p[6]))
                P.append(float(p[7]))
                
                
                #x.append(float(p[4]))
    return MC_cycles, A,P#cv, x

plt.figure()
MC_cycles,A,P= Read_eight_Column_File('C:\Data\configsU2111.txt')
t2 = np.divide(P,400)

plt.hist(t2, bins= 200)

plt.xlim(-2,-1.6)
plt.title("\n".join(wrap('Probability distribution, T = 1', 60)))
plt.ylabel('"Probabilty"')
plt.xlabel('Total Energy E')

#plt.savefig('probdist1.png')
plt.show()
plt.figure()
MC_cycles,A,P2= Read_eight_Column_File('C:\Data\prob24.txt')
#MC_cycles, E, M,A= Read_eight_Column_File('C:\Data\configsU.txt')

t = np.divide(P2,400)
plt.hist(t,bins= 35)
plt.xlim(-2,-0.5)
plt.title("\n".join(wrap('Probability distribution, T = 2.4', 60)))
plt.ylabel('"Probabilety"')
plt.xlabel('Total Energy E')
#plt.plot(MC_cycles, A)
#plt.savefig('probdist2.png')

plt.show()

plt.figure()
MC_cycles,A3,P3= Read_eight_Column_File('C:\Data\prob24.txt')
plt.plot(MC_cycles, A3)
plt.ylim(0,0.5)
plt.show()



