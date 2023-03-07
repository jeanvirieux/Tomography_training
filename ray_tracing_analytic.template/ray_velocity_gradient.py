##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# computing the ray using RK2 and the analytical solution
# for constant vertical gradient of velocity: horizontal gradient should be 0
##############################################################################
import sys
import numpy as np

##############################################################################
#
# analytical model of the velocity with a vertical constant gradient
#
##############################################################################

def u2_est(qx,qz,czero,gammax,gammaz,iopt):
  
  if iopt == 0:
    u= 1./(czero+gammax*qx+gammaz*qz)
    u2=u*u
    return u2,u

  if iopt == 1:
    u=1./(czero+gammax*qx+gammaz*qz)
    u2=u*u
    u3=u2*u
    u2qx=-2*gammax*u3
    u2qz=-2.*gammaz*u3
    return u2,u2qx,u2qz
  
  if iopt == 2:
    u=1./(czero+gammax*qx+gammaz*qz)
    u2=u*u
    u3=u2*u
    u4=u3*u
    u2qx=-2*gammax*u3
    u2qz=-2.*gammaz*u3
    u2qxqx=6.*gammax*gammax*u4
    u2qxqz=6.*gammax*gammaz*u4
    u2qzqz=6.*gammaz*gammaz*u4
    return u2,u2qx,u2qz,u2qxqx,u2qzqz,u2qxqz

###############################################################
#  simple Runge-Kutta integration: order 2 ... with possible hamiltonian checking
###############################################################

def rk2_ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
            czero,gammax,gammaz,\
            theta,dtau,nt_max,\
            ioption,idebug):
  import math
  

  iopt=0  
  est_u2=u2_est(xsrc,zsrc,czero,gammax,gammaz,iopt)
  u_source=math.sqrt(est_u2[0])
  print(' source position:',xsrc,zsrc,' km,km; source velocity:',1./u_source,' km/s')
  

# initial conditions in the 4D space
  pxi=u_source*math.sin(theta) # slowness vector at the source
  pzi=u_source*math.cos(theta) # ...
  qxi=xsrc
  qzi=zsrc
  sin_theta0=pxi/u_source
  cos_theta0=math.sqrt(1.-sin_theta0**2)

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

  X_ana=[]
  Z_ana=[]
  X_ana.append(qxi)
  Z_ana.append(qzi)

  # initial travel-time  (compute triangular estimation)
  #           should be improved by the Simpson rule
  temps=0.
  
  # loop for tracing ray
  nt=0
  qtau=0.
  iopt=1

### open file with appending option
  outfile = open('ray.dat','a+')
  
  while nt < nt_max:
    nt=nt+1
    est_u2=u2_est(qx0,qz0,czero,gammax,gammaz,iopt)
    u2=est_u2[0]
    u2qx=est_u2[1]
    u2qz=est_u2[2]
    
    if idebug == 3:
      ham=0.5*(px0**2+pz0**2-u2)/u2    # normalized hamiltonian
      print(' hamiltonian drift : error estimation ',ham, qx0,qz0,px0,pz0)
      
    if ioption/10 >= 1:
       q=(qx0,-qz0)    # analytical solution 
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
    
    est_u2=u2_est(qx1,qz1,czero,gammax,gammaz,iopt)
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

    temps=temps+u2_time*dtau1  #trapeze rule

    ##################################################
    ##################################################
    ####  analytical solution for the velocity with a constant vertical gradient
    ##################################################
    ##################################################
    sin_thetar=pxi/math.sqrt(u2_time)   # find the local angle
    if sin_thetar > 0.999999:
      sin_thetar=1.
    
    cos_thetar=math.sqrt(1.-sin_thetar**2)   # restric to positive value
    qx=qxi+(cos_theta0-cos_thetar)/gammaz/pxi
    qz=qzi+(sin_thetar-sin_theta0)/gammaz/pxi

    X_ana.append(qx)
    Z_ana.append(qz)

    if idebug == 2:
      print('qx0,qz0,qx1,qz1 ',qx0,qz0,qx1,qz1,nt,nt_max)
      print(' xmin,xmax,zmin,zmax ',xmin,xmax,zmin,zmax)
      
    if qz0 <= zmin or qz0 >= zmax or qx0 <= xmin or qx0 >= xmax:
      break
    
############################################# end of the while    
  if ioption/10 >=1:
     outfile.write("     \n")   # empty line for going to the next ray
     
  outfile.close()               # close the file ... for this ray
  
  print(' reached increment:',nt,nt_max)
  return X_pos,Z_pos,X_ana,Z_ana   

