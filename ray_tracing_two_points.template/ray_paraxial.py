##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# computing ray and paraxial ray using Runge-Kutta of 2nr order
# with a given shooting angle theta.
# return values are
#
# ray quantities qx_pos,qz_pos,px_pos,pz_pos,t_pos
# paraxial quantities dqx_pos,dqz_pos for expansion and 2-pts tracing
#
##############################################################################

import sys
import numpy as np
import math
import grid

###############################################################
#  simple Runge-Kutta integration: order 2 ... with possible hamiltonian checking
###############################################################


def ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
            u2tab,xpos,zpos,nx,nz,\
            theta,dtau,nt_max,\
            ioption,idebug):
  

  ier=-999
  iopt=0  
  est_u2=grid.u2(xsrc,zsrc,u2tab,xpos,zpos,nx,nz,iopt,ier)
  u_source=math.sqrt(est_u2[0])

  dtet=1.e20 # initial angle perturbation

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

  # initial paraxial trajectories (not yet paraxial rays)
  # q paraxial used for Gaussian Beam summation
  dqx1=1.; dqz1=0.; dpx1=0.; dpz1=0.
  dqx2=0.; dqz2=1.; dpx2=0.; dpz2=0.
  # p paraxial   used for ray tube estimation
  dqx3=0.; dqz3=0.; dpx3=1.; dpz3=0.
  dqx4=0.; dqz4=0.; dpx4=0.; dpz4=1.
  
  #%%%%%%%%%%%%%%%%%%%%%%%% few other quantities
  # scale the jacobian and set the kmah index
  fact=1.; kmah_current=0
  
  # point paraxial ray for shooting angle estimation through (dqx,dqz)
  dqx=0.; dqz=0.
  dqx_pos=[]
  dqx_pos.append(dqx)
  dqz_pos=[]
  dqz_pos.append(dqz)
  
  qx_pos=[]
  qz_pos=[]
  px_pos=[]
  pz_pos=[]
  qx_pos.append(qx0)
  qz_pos.append(qz0)
  px_pos.append(px0)
  pz_pos.append(pz0)
  
  # initial travel-time  (compute triangular estimation)
  #           should be improved by the Simpson rule
  temps=0.
  t_pos=[]
  t_pos.append(temps)
    
  # loop for tracing ray
  #
  nt=0
  iopt=1

### open file with appending option
  if ioption/10 >= 1:
    outfile = open('centray.dat','a+')
    parfile = open('pararay.dat','a+')
  
  while nt < nt_max:
    nt=nt+1
    
    iopt=1
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
       if ioption/10 == 2:
         dq=(qx0+dqx,-qz0-dqz) # the full paraxial ray (no scaling)
         parfile.write(str(dq)+"\n")
       if ioption/10 == 3:     # only the paraxial perturbation
         dq=(dqx,-dqz)
         parfile.write(str(dq)+"\n")
         
    ##################################################
    ##################################################
    ##################################################
    qx1=qx0+dtau2*px0          # first step of RK2
    qz1=qz0+dtau2*pz0
    px1=px0+dtau2*0.5*u2qx
    pz1=pz0+dtau2*0.5*u2qz

    iopt=2
    est_u2=grid.u2(qx1,qz1,u2tab,xpos,zpos,nx,nz,iopt,ier)
    
    u2_time=est_u2[0] # save this middle estimation for trapezoidal rule of time integration
    u2qx=est_u2[1]
    u2qz=est_u2[2]
    u2qxqx=est_u2[3]
    u2qzqz=est_u2[4]
    u2qxqz=est_u2[5]

    ##################################################
    ##################################################
    ##################################################
    qx1=qx0+dtau1*px1          # second step of RK2
    qz1=qz0+dtau1*pz1
    px1=px0+dtau1*0.5*u2qx
    pz1=pz0+dtau1*0.5*u2qz

        
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # paraxial integration  euler integration ... 
    # we use the estimation at the midpoint of integration
    # reason for computing second derivatives
    # but we can use Euler integration only
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dqx1=dqx1+dtau1*dpx1
    dqz1=dqz1+dtau1*dpz1
    dpx1=dpx1+dtau1*0.5*(dqx1*u2qxqx+dqz1*u2qxqz)
    dpz1=dpz1+dtau1*0.5*(dqx1*u2qxqz+dqz1*u2qzqz)

    dqx2=dqx2+dtau1*dpx2
    dqz2=dqz2+dtau1*dpz2
    dpx2=dpx2+dtau1*0.5*(dqx2*u2qxqx+dqz2*u2qxqz)
    dpz2=dpz2+dtau1*0.5*(dqx2*u2qxqz+dqz2*u2qzqz)

    dqx3=dqx3+dtau1*dpx3
    dqz3=dqz3+dtau1*dpz3
    dpx3=dpx3+dtau1*0.5*(dqx3*u2qxqx+dqz3*u2qxqz)
    dpz3=dpz3+dtau1*0.5*(dqx3*u2qxqz+dqz3*u2qzqz)

    dqx4=dqx4+dtau1*dpx4
    dqz4=dqz4+dtau1*dpz4
    dpx4=dpx4+dtau1*0.5*(dqx4*u2qxqx+dqz4*u2qxqz)
    dpz4=dpz4+dtau1*0.5*(dqx4*u2qxqz+dqz4*u2qzqz)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # paraxial ray estimation through initial conditions
    # point excitation
    # dqxi=0 dqzi=0   dpxi=pzi   dpzi=-pxi  point paraxial
    #                 dpxi=0.    dpzi=0.    plane paraxial 
    # This gives the derivative with respect to the initial angle
    # of the horizontal displacement.
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dqx=(dqx3*pzi-dqx4*pxi)  # at constant tau
    dqz=(dqz3*pzi-dqz4*pxi)  # paraxial ray related to
    dpx=(dpx3*pzi-dpx4*pxi)  # ray tube (solutions 3 and 4)
    dpz=(dpz3*pzi-dpz4*pxi)  # GBS and plane wave (solutions 1 and 2)

    dqx_pos.append(dqx)
    dqz_pos.append(dqz)
    
  #
  #  caustic detection   indice kmah
  #
    jacobian=dqx*pz1-dqz*px1    # do we cross the central ray by the paraxial ray
    jacobian=jacobian*fact      # fact is equal to one: differential quantity!
    if math.copysign(1.,jacobian) < 0. and kmah_current%2 == 0:
      kmah_current=kmah_current+1
    if math.copysign(1.,jacobian) > 0. and kmah_current%2 == 1:
      kmah_current=kmah_current+1
    
  # we must update for further computation along the ray
    qx0=qx1
    qz0=qz1
    px0=px1
    pz0=pz1
    
    qx_pos.append(qx0)
    qz_pos.append(qz0)
    px_pos.append(px0)
    pz_pos.append(pz0)
    
    temps=temps+u2_time*dtau1  #trapeze rule
    t_pos.append(temps)
    
    if idebug == 2:
      print('qx0,qz0,qx1,qz1 ',qx0,qz0,qx1,qz1,nt,nt_max)
      print(' xmin,xmax,zmin,zmax ',xmin,xmax,zmin,zmax)
      
    if qz0 <= zmin or qz0 >= zmax or qx0 <= xmin or qx0 >= xmax:
      break
        
############################################# end of the while    
  if ioption/10 >=1:
     outfile.write("     \n")   # empty line for going to the next ray
     if ioption/10 == 2:
       dq=(qx0+dqx,-qz0-dqz) # the paraxial ray
       parfile.write("     \n")
     if ioption/10 == 3:     # only the paraxial values
       dq=(dqx,-dqz)
       parfile.write("     \n")
     outfile.close()               # close the file ... for this ray
     parfile.close()               # close the file ... for the paraxial ray
  
  return qx_pos,qz_pos,px_pos,pz_pos,t_pos,dqx_pos,dqz_pos,nt,kmah_current

