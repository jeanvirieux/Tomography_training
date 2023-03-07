##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# computing the ray using Runge-Kutta of 2nr order
# until xtag,ztag for drawing ... essentially
##############################################################################
import sys
import numpy as np

import grid

###############################################################
#  simple Runge-Kutta integration: order 2 ... with possible hamiltonian checking
###############################################################

def ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
            u2tab,xpos,zpos,nx,nz,\
            theta,dtau,nt_max,\
            xtag,ztag,eps,\
            ioption,idebug):
  
  import math

  ier=-999
  iopt=0  
  est_u2=grid.u2(xsrc,zsrc,u2tab,xpos,zpos,nx,nz,iopt,ier)
  u_source=math.sqrt(est_u2[0])
  
# initial conditions in the 4D space
  pxi=u_source*math.sin(theta) # slowness vector at the source
  pzi=u_source*math.cos(theta) # ...
  qxi=xsrc
  qzi=zsrc

  dtau0=dtau    # to be converted dl in m ....
  dtau1=dtau0/u_source # we jump into the tau sampling with correct dimensions km**2/s
  dtau2=0.5*dtau1  # half-step for RK2

  # ' source ',qxi,qzi,u_source*u_source,dtau1
  # ' pxi, pzi', pxi,pzi,theta*180./3.14159

  # initial ray coordinates in the 4D phase space
  qx0=qxi
  qz0=qzi
  px0=pxi
  pz0=pzi
  # initial analytical solution
  qx=qxi
  qz=qzi
  px=pxi
  pz=pzi
  
  X_pos=[]
  Z_pos=[]
  X_pos.append(qxi)
  Z_pos.append(qzi)

  # initial travel-time  (compute triangular estimation)
  #           should be improved by the Simpson rule
  temps=0.
  
  # loop for tracing ray
  nt=0
  iopt=1

### open file with appending option
  if ioption >= 1:
    outfile = open('ray.dat','a+')
  
  while nt < nt_max:
    nt=nt+1
    est_u2=grid.u2(qx0,qz0,u2tab,xpos,zpos,nx,nz,iopt,ier)
    u2=est_u2[0]
    u2qx=est_u2[1]
    u2qz=est_u2[2]
    
    if idebug == 3:
      ham=0.5*(px0**2+pz0**2-u2)/u2    # normalized hamiltonian
      print(' hamiltonian drift : error estimation ',ham, qx0,qz0,px0,pz0)
      
    if ioption/10 >= 1:
       q=(qx0,-qz0)  
       outfile.write(str(q)+"\n")
       
    X_pos.append(qx0)
    Z_pos.append(qz0)
    
    ##################################################
    ##################################################
    ##################################################
    qx1=qx0+dtau2*px0          # first step of RK2
    qz1=qz0+dtau2*pz0
    px1=px0+dtau2*0.5*u2qx
    pz1=pz0+dtau2*0.5*u2qz

    est_u2=grid.u2(qx1,qz1,u2tab,xpos,zpos,nx,nz,iopt,ier)
    
    u2_time=est_u2[0] # save this middle estimation for trapezoidal rule of time integration
    u2qx=est_u2[1]
    u2qz=est_u2[2]

    ##################################################
    ##################################################
    ##################################################
    qx1=qx0+dtau1*px1          # second step of RK2
    qz1=qz0+dtau1*pz1
    px1=px0+dtau1*0.5*u2qx
    pz1=pz0+dtau1*0.5*u2qz

  # we must update for further computation along the ray
    qx0=qx1
    qz0=qz1
    px0=px1
    pz0=pz1

    dist=(xtag-qx0)**2+(ztag-qz0)**2
    if dist < eps*eps:
      break
    
    temps=temps+u2_time*dtau1  #trapeze rule

    if idebug == 2:
      print('qx0,qz0,qx1,qz1 ',qx0,qz0,qx1,qz1,nt,nt_max)
      print(' xmin,xmax,zmin,zmax ',xmin,xmax,zmin,zmax)
      
    if qz0 <= zmin or qz0 >= zmax or qx0 <= xmin or qx0 >= xmax:
      break
    
############################################# end of the while    
  if ioption/10 >=1:
    outfile.write("     \n")   # empty line for going to the next ray
    outfile.close()               # close the file ... for this ray
  
  return X_pos,Z_pos

