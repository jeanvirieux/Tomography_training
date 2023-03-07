##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# 
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import ray_slowness2_gradient

#############################################
#  main program: ray tracing
#############################################
################################### call from the command line ... 
#if len(sys.argv) == 14 :
#  xsrc=float(sys.argv[1])
#  zsrc=float(sys.argv[2])
#  xmin=float(sys.argv[3])
#  xmax=float(sys.argv[4])
#  zmin=float(sys.argv[5])
#  zmax=float(sys.argv[6])
#  theta=float(sys.argv[7])
#  dtheta=float(sys.argv[8])
#  ntheta=int(sys.argv[9])
#  dtau=float(sys.argv[10])
#  nt_max=int(sys.argv[11])
#  ioption=int(sys.argv[12])
#  idebug=int(sys.argv[13])
#else:
#  sys.exit()
#  print(len(sys.argv))
#  print(' we need 13 arguments ')

infile = open('input.inp','r')
line=infile.readline()
coords = line.split()
xsrc=float(coords[0])   # source position
zsrc=float(coords[1])
  
line=infile.readline()
coords = line.split()
xmin=float(coords[0])   # box definition
xmax=float(coords[1])
zmin=float(coords[2])
zmax=float(coords[3])
  
line=infile.readline()
coords = line.split()
theta=float(coords[0])  # shooting angles
dtheta=float(coords[1])
ntheta=int(coords[2])

line=infile.readline()
coords = line.split()
dtau=float(coords[0])   # tau integration
nt_max=int(coords[1])

line=infile.readline()
coords = line.split()
ioption=int(coords[0])  # option for output and debug
idebug=int(coords[1])
infile.close()
  
print('source position: xsrc,zsrc')
print('box: xmin,xmax,zmin,zmax')
print('initial angle, angle sampling, number of angles: theta,dtheta,ntheta')
print('integration step (m): dtau')
print('max number of points along a ray: nt_max')
print('option for output: ioption')
print('option for debugging: idebug')
 

theta=3.14159/180.*theta   # radian
dtheta=3.14159/180.*dtheta

############################################################
########### model definition u2=u2_0+gammax.x+gammaz.z
############################################################

infile = open('model.inp','r')
line=infile.readline()
coords = line.split()
czero=float(coords[0])
gammax=float(coords[1])
gammaz=float(coords[2])
infile.close()

############################################################
########### saving file
############################################################
print('WARNING: do not forget to erase ray.dat')

  
############################################
# starting ray tracing
############################################
theta=theta-dtheta
for ix in range(0,ntheta):
   theta=theta+dtheta
   
   print(' shooting angle:',theta*57.3)
   
   ray_rk2=ray_slowness2_gradient.rk2_ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
                   czero,gammax,gammaz,\
                   theta,dtau,nt_max,\
                   ioption,idebug)

#   xmin=min(ray_rk2[0])
#   xmax=max(ray_rk2[0])
#   zmin=min(ray_rk2[1])
#   zmax=max(ray_rk2[1])
#   print('box',xmin,xmax,zmin,zmax)

   plt.xlim(xmin,xmax)
   plt.ylim(zmax,zmin)

   plt.xlabel('Horizontal distance - km ')
   plt.ylabel('vertical distance - km ')
   plt.title('Ray Tracing')
   
   plt.plot(ray_rk2[2],ray_rk2[3],'b')
   plt.plot(ray_rk2[0],ray_rk2[1],'r--',linewidth=4)

plt.savefig("ray_slowness_2D.pdf",dpi=300,format='pdf')  
plt.show()

  

 
