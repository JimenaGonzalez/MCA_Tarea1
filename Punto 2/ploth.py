# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import sys
import pylab
import numpy as np
import matplotlib.pyplot as plt
import os
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																														

import random
plt.isinteractive()
True
plt.ishold
True

lst = []
t =[]
E =[]
f = open('datos.dat', 'r')
for line in f:
    lst += [line.split()]
t = [x[2] for x in lst]
E = [x[3] for x in lst]
 
  
 
plt.scatter(t,E, s=0.5)
 
plt.ylabel('E')
plt.ylabel('t')
 
plt.xlim(0,1500)
plt.ylim(0,-1 )  
plt.savefig(sys.argv[1])




   

# <codecell>


