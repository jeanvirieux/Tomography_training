##############################################################################
# need python 3.*  (I guess)
# maybe conda activate base
# point source xsrc,zsrc
# xmin,xmax,zmin,zmax box
# grid description u2,xpos,zpos,nx,nz
# receiver xobs,zobs,nobs
# distance convergence eps,distance tolerance eps1, max iterations ishoot_max
# initial angle and ray evolution parameters theta,dtau,nt_max
# option for storing and printing  ioption,idebug
############################################# when calling ray_sht=ray_paraxial()
#   ray_sht[0]    # qx_pos[] 
#   ray_sht[1]    # qz_pos[]
#   ray_sht[2]    # px_pos[]
#   ray_sht[3]    # pz_pos[]
#   ray_sht[4]    # t_pos[]
#   ray_sht[5]=   # dqx_pos[]
#   ray_sht[6]=   # dqz_pos[]
#   ray_sht[7]=   # nt
#   ray_sht[8]=   # kmah
##############################################
##############################################################################
import sys
import numpy as np

import ray_paraxial

def ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
                   u2,xpos,zpos,nx,nz,\
                   xobs,zobs,nobs,eps,eps1,ishoot_max,\
                   theta,dtau,nt_max,\
                   ioption,idebug):
   
############################################
# starting two-points ray tracing
############################################
   print(' searching for shooting angles for the receiver array')
   print('######################################################')
   print(' initial shooting angle:',theta*57.3)
   tet_obs=[]
   for iobs in range(0,nobs):
      print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
      print('&&&&&&&&& receiver number',iobs+1,'receiver position',xobs[iobs],zobs[iobs])
      print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
      dist_save=1.e20
      for ishoot in range(0,ishoot_max):
   
         ray_sht=ray_paraxial.ray(xsrc,zsrc,xmin,xmax,zmin,zmax,\
                            u2,xpos,zpos,nx,nz,\
                            theta,dtau,nt_max,\
                            ioption,idebug)
         nt_retour=ray_sht[7]
         qx0=ray_sht[0]
         qz0=ray_sht[1]
         distold=(xobs[iobs]-qx0[0])**2+(zobs[iobs]-qz0[0])**2

         for it in range(1,nt_retour):  # loop over the ray
            dist=(xobs[iobs]-qx0[it])**2+(zobs[iobs]-qz0[it])**2
            if dist < distold:
              it_save=it
              distold=dist  

         print("nearest point",it_save,qx0[it_save],qz0[it_save])
         if distold < dist_save:
           dist_save=distold      # save the smallest distance

         if dist_save < eps*eps:  # we converge
           break
     
# improving the shooting
# estimation of the correction from the paraxial

         px0=ray_sht[2]
         pz0=ray_sht[3]
         dqx=ray_sht[5]
         dqz=ray_sht[6]
# ok DL= derivative. dtet
# DL=pz0*(xobs-qx0)-px0*(zobs-qz0)   derivative=pz0*dqx-px0*dqz .... projection of vector (dqx,dqz)
#                                    orthogonal to the ray direction
         dtet=(pz0[it_save]*(xobs[iobs]-qx0[it_save])-px0[it_save]*(zobs[iobs]-qz0[it_save]))/ \
          (pz0[it_save]*dqx[it_save]-px0[it_save]*dqz[it_save])  # correction paraxiale
      
         theta=theta+dtet
#         print('next theta',theta*57.27,'ishoot',ishoot)
      
      if ishoot < ishoot_max:
        ier=0
        print('selected theta',theta*57.27)
        tet_obs.append(theta)
  
      if ishoot == ishoot_max:
        if dist_save < eps1*eps1:          # tolerance convergence (often 10*eps)
          print('tolerated theta',theta*57.27)
          tet_obs.append(theta)
        else:
          print(' failure for reaching this receiver *******')
          tet_obs.append(-999.)                       # failure

   return tet_obs
