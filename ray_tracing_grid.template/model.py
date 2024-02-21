##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# 
##############################################################################
import sys
import numpy as np
import math

#############################################
#  building the mesh and the velocity structure
#############################################
###### call from the command line is avoided ...

infile = open('model.inp','r')
line=infile.readline()
coords = line.split()
xorg=float(coords[0])     # source position
zorg=float(coords[1])
  
line=infile.readline()
coords = line.split()
nx=int(coords[0])         # grid definition
nz=int(coords[1])
dx=float(coords[2])
dz=float(coords[3])
  
line=infile.readline()
coords = line.split()
czero=float(coords[0])    # velocity
gammax=float(coords[1])   # constant gradient x and z
gammaz=float(coords[2])   # first z and then x

line=infile.readline()
coords = line.split()
x_ball=float(coords[0])   # center of the ball
z_ball=float(coords[1])

line=infile.readline()    # exponential ball
coords = line.split()
xlength=float(coords[0])
zlength=float(coords[1])
dpert_v=float(coords[2])
infile.close()

############################################################
########### model definition grid and values
########### model.inp and model.vel (binary)
############################################################

xpos=np.empty(nx)
zpos=np.empty(nz)
for ix in range(0,nx):
  xpos[ix]=xorg+float(ix)*dx
for iz in range(0,nz):
  zpos[iz]=zorg+float(iz)*dz

xtot=xpos.max()-xpos.min()
ztot=zpos.max()-zpos.min()

print('min,max xgrid,xlength',xpos.min(),xpos.max(),xtot)
print('min,max zgrid,zlength',zpos.min(),zpos.max(),ztot)

outfile = open('model.hed','w')
outfile.write(str(nx)+' '+str(nz)+'\n')
outfile.write(str(dx)+' '+str(dz)+'\n')
outfile.close()
  
outfile = open('model.xzp','wb')
a = np.array(xpos,'float32')
a.tofile(outfile)
a = np.array(zpos,'float32')
a.tofile(outfile)
outfile.close()

velocity=np.empty([nz,nx])

for iz in range(0,nz):
   z=zpos[iz]-zorg
   zr=(zpos[iz]-z_ball)/zlength
   for ix in range(0,nx):
      x=xpos[ix]-xorg
      xr=(xpos[ix]-x_ball)/xlength
      arg=xr*xr+zr*zr
      velocity[iz][ix]=czero+gammax*x+gammaz*z+dpert_v*math.exp(-arg)

print('min velocity',velocity.min())
print('max velocity',velocity.max())

a = np.array(velocity,'float32')

output_file = open('model.vel', 'wb')
a.tofile(output_file)
output_file.close()

  

 
