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
a=[]
b=[]
c=[]
d=[]


############################################
f = open(sys.argv[1], 'r')
for line in f:
    lst += [line.split()]

a= [x[0] for x in lst]
b= [x[1] for x in lst]
c= [x[2] for x in lst]
d= [x[3] for x in lst]


for i in range(0,len(a)):
 
 os.system("./a.out "+str(a[i])+" "+str(b[i])+" "+str(c[i])+" "+str(d[i])+" > datos.dat ")
 
 lst = []
 q3_R =[]
 p3_R =[]
 f = open('datos.dat', 'r')
 for line in f:
    lst += [line.split()]
 q3 = [x[0] for x in lst]
 p3 = [x[1] for x in lst]
 
  
 
 plt.scatter(q3, p3, s=0.5)
 print "./a.out "+str(a[i])+" "+str(b[i])+" "+str(c[i])+" "+str(d[i])+" > datos.txt ",i

plt.ylabel('P3')
plt.xlabel('Q3') 
plt.xlim(float(sys.argv[2]),float (sys.argv[3]))
plt.ylim(float(sys.argv[4]), float (sys.argv[5]))  
plt.savefig(sys.argv[6])




   

# <codecell>


