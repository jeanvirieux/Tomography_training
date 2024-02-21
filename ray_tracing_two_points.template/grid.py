##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# 
##############################################################################
import sys
import numpy as np

import bspline

##############################################################################
#
# numerical evaluation of the slowness square and its derivatives
# inside a grid using bspline ...
#
##############################################################################

def u2(x,z,u2tab,xpos,zpos,nx,nz,iopt,ier):

#  xpos=struct.unpack("f"*nx,xp)
#  zpos=struct.unpack("f"*nz,zp)
#  u2tup=struct.unpack("f"*nx*nz,u2grid)
#  u2tab=np.asarray(u2tup)
#  u2tab=u2tab.reshape(nx,nz)

#  print('grid u2 min',u2tab.min())
#  print('grid u2 max',u2tab.max())
  
  pp=np.empty([4,4])
  ntx=0
  ntz=0

  for ix in range(0,nx):
    if x < xpos[ix]:
      ntx=ix-1
      break
    ntx=nx
  for iz in range(0,nz):
    if z < zpos[iz]:
      ntz=iz-1
      break
    ntz=nz
    
  if ntx < 1:
    ntx=1
  if ntx > nx-2:
    ntx=nx-2
  if ntz < 1:
    ntz=1
  if ntz > nz-2:
    ntz=nz-2
    
######################### value u2 always
  xu=xpos[ntx]
  zw=zpos[ntz]
  dxpas=xpos[ntx+1]-xpos[ntx]
  dzpas=zpos[ntz+1]-zpos[ntz]
  ux=(x-xu)/dxpas
  wz=(z-zw)/dzpas

#  print('ux,wz',ux,wz)
#  print('ntx,ntz',ntx,ntz)
  
  for k in range(0,4):
     for i in range(0,4):
#        print('i,k',i,k)
#        print('ix,iz',ntx-1+i,ntz-1+k,u2tab[ntx-1+i][ntz-1+k])
        pp[i][k]=u2tab[ntz-1+k][ntx-1+i]
#  print('pp',pp)
  
######################### value
  u2val=bspline.bsp2(ux,wz,pp)

######################### first derivative
  if iopt == 1:
    u2qx=bspline.bsp2s1(ux,wz,pp)/dxpas
    u2qz=bspline.bsp2t1(ux,wz,pp)/dzpas
    return u2val,u2qx,u2qz
######################### first and second derivatives  
  if iopt == 2:
    u2qx=bspline.bsp2s1(ux,wz,pp)/dxpas
    u2qz=bspline.bsp2t1(ux,wz,pp)/dzpas
    u2qxqx=bspline.bsp2s2(ux,wz,pp)/dxpas/dxpas
    u2qxqz=bspline.bsp2ts(ux,wz,pp)/dxpas/dzpas
    u2qzqz=bspline.bsp2t2(ux,wz,pp)/dzpas/dzpas
    return u2val,u2qx,u2qz,u2qxqx,u2qzqz,u2qxqz
#########################
  if iopt == 0:
    return u2val,u2val
