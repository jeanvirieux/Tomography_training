##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# 
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import struct

import ray_rk

#############################################
#  main program: ray tracing
#############################################
###### call from the command line is avoided ...

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
########### model.xzp grid definition xpos,zpos (binary)
########### model.vel velocity (binary)
############################################################

infile = open('model.hed','r')
line=infile.readline()
coords = line.split()
nx=int(coords[0])
nz=int(coords[1])

model_grd='model.xzp'
xzp = np.fromfile(str(model_grd), 'float32')
xpos=np.split(xzp,2)   # xpos[0] along x and xpos[1] along z

model_vel='model.vel'
velocity = np.fromfile(str(model_vel), 'float32')
vmax=np.amax(velocity)
vmin=np.amin(velocity)

u2=1./velocity
u2=u2.reshape(nz,nx)

print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
print('xgrid min,max',xpos[0].min(),xpos[0].max())
print('zgrid min,max',xpos[1].min(),xpos[1].max())
print('velocity min',vmin)
print('velocity max',vmax)
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

velocity=velocity.reshape(nz,nx)

############################################################
########### saving file
############################################################
print('WARNING: do not forget to erase ray.dat')

#################### velocity structure   
plt.imshow(velocity,origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[xmin,xmax,zmax,zmin],vmin=vmin,vmax=vmax,cmap=plt.cm.get_cmap('bwr',24),clip_on=True)
   
############################################
# starting ray tracing
############################################
theta=theta-dtheta
for ix in range(0,ntheta):
   theta=theta+dtheta
   
   print(' shooting angle:',theta*57.3)
   
   ray_rk2=ray_rk.ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
                   u2,xpos[0],xpos[1],nx,nz,\
                   theta,dtau,nt_max,\
                   ioption,idebug)

   plt.xlim(xmin,xmax)
   plt.ylim(zmax,zmin)

   plt.xlabel('Horizontal distance - km ')
   plt.ylabel('vertical distance - km ')
   plt.title('Ray Tracing')

#################### rays   
   plt.plot(ray_rk2[0],ray_rk2[1],'black',linewidth=2)
   plt.savefig("ray_2D.pdf",dpi=300,format='pdf')
   
plt.show()

  

 
