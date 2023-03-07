#########################################################################################
#                                                                                       #
#  example using a more complex model      (including rays in this case)                #
#                                                                                       #
# python HFM_all.py file name nz nx  dcarre zorg xorg zsrc xsrc  (z first and then x)   #
# python HFM_all.py nankai_double.bin 67 365 468.75 -468.75 -468.75 30000. 75000.       #
#                                                                                       #
#########################################################################################

import sys
import os
import numpy as np
import struct
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import def_HFM_all


verbose=1.         # could be zero for less output  (useless when doing it through this python)
order=2.           # minimal stencil order
tsrc=0.            # time at seed
radius=-1.         # automatic factorization
export=[]          # define the list of quantities to be evaluated
export.append(1.)  # output time field
export.append(0.)  # output flow (to be multiplied by the slowness square to get the slowness vector
export.append(0.)  # output combined sensitivities
export.append(1.)  # output rays

xrec=[];zrec=[];trec=[]
rec=8.
xrec.append(25000.)
zrec.append(20000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(40000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(60000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(80000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(100000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(120000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(140000.)
trec.append(1.)
xrec.append(25000.)
zrec.append(160000.)
trec.append(1.)

if int(rec) == 0:
  export[2]=0.
  export[3]=0.

fac=0.5

print('nbre of arguments',len(sys.argv))
##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 9:
  model_vel = sys.argv[1]
  nx=float(sys.argv[2])
  nz=float(sys.argv[3])
  hcarre=float(sys.argv[4])
  xorg=float(sys.argv[5])
  zorg=float(sys.argv[6])
  xsrc=float(sys.argv[7])
  zsrc=float(sys.argv[8])
else:
  print(len(sys.argv),' 8 arguments needed ')
  print(' need the name of the velocity file ')
  print(' nx,nz of the grid ')
  print(' hcarre of the square grid ')
  print(' xorg,zorg of the grid ')
  print(' xsrc,zsrc ')
  sys.exit()
  
nxc=int(nx)
nzc=int(nz)

print('nxc,nzx ',nxc,nzc)

xmin=xorg;xmax=xorg+nx*hcarre
zmin=zorg;zmax=zorg+nz*hcarre


print('box',xmin,xmax,zmin,zmax)

if xsrc < xmin or xsrc > xmax:
    print(' xsrc outside the box ',xsrc,'box',xmin,xmax)
    sys.exit()

if zsrc < zmin or zsrc > zmax:
    print(' zsrc outside the box ',zsrc,'box',xmin,zmax)
    sys.exit()


scale=3.
fig_x=6. * scale
fig_z=fig_x*(zmax-zmin)/(xmax-xmin) * scale

plt.rcParams["figure.figsize"] = (fig_x, fig_z)

plt.yticks(fontsize=20)
plt.xticks(fontsize=20)

velocity = np.fromfile(str(model_vel), "float64")
invvel = 0.001/velocity    # s/m

velocity = velocity.reshape(nxc,nzc)

############# pack into a buffer
slowness=struct.pack('d'*nxc*nzc,*invvel)    # flatten the list for packing 

############# calling the HFM function depending on the flag export
if int(export[0]) == 1 and int(export[1]) == 0 and int(export[2]) == 0 and int(export[3]) == 0:
  field=def_HFM_all.hfm_py(verbose,order,xorg,zorg,hcarre,nx,nz,slowness,\
                         xsrc,zsrc,tsrc,radius,\
                         rec,xrec,zrec,trec,fac,export)
  
############# calling the HFM function depending on the flag export
if int(export[0]) == 1 and int(export[1]) == 0 and int(export[2]) == 0 and int(export[3]) == 1:
  hfm_all=def_HFM_all.hfm_py(verbose,order,xorg,zorg,hcarre,nx,nz,slowness,\
                         xsrc,zsrc,tsrc,radius,\
                         rec,xrec,zrec,trec,fac,export)
  field=hfm_all[0]    # time
  nray=hfm_all[1]     # number of rays
  npts_ray=hfm_all[2] # npts_ray[nray]
  npts_geo=hfm_all[3] # total number of points
  pts_ray=hfm_all[4]  # ray points [x,z]
  
############# get back the array from the buffer   always computed
time=np.asarray(field)

time=time.reshape(nxc,nzc)

################## velocity
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.imshow(velocity,origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[zmin,zmax,xmax,xmin],cmap=cm.seismic,clip_on=True)

cb=plt.colorbar(orientation='horizontal')
cb.set_label(label='(km/s)',size=22) #'large')
cb.ax.tick_params(labelsize=20)

plt.xlabel('Range (km)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.title('P wave velocity',fontsize=28)

################## should we include rays?
if int(export[0]) == 1 and int(export[1]) == 0 and int(export[2]) == 0 and int(export[3]) == 1:
  ioffset=0
  for iray in range(0,nray):
    xpos=np.empty(int(npts_ray[iray]))
    zpos=np.empty(int(npts_ray[iray]))
    ipoint=0
    for jray in range(0,int(npts_ray[iray])):
       xpos[jray]=pts_ray[ioffset+2*jray]
       zpos[jray]=pts_ray[ioffset+2*jray+1]
       ipoint=ipoint+2
    plt.plot(zpos,xpos,'black',linewidth=2)
    ioffset=ioffset+ipoint
      
plt.show()

##################  time
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.imshow(time,origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[zmin,zmax,xmax,xmin],cmap=cm.seismic,clip_on=True)

cb=plt.colorbar(orientation='horizontal')
cb.set_label(label='(sec)',size=22) #'large')
cb.ax.tick_params(labelsize=20)

plt.xlabel('Range (km)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.title('Travel time',fontsize=28)

################## should we include rays?
if int(export[0]) == 1 and int(export[1]) == 0 and int(export[2]) == 0 and int(export[3]) == 1:
  ioffset=0
  for iray in range(0,nray):
    xpos=np.empty(int(npts_ray[iray]))
    zpos=np.empty(int(npts_ray[iray]))
    ipoint=0
    for jray in range(0,int(npts_ray[iray])):
       xpos[jray]=pts_ray[ioffset+2*jray]
       zpos[jray]=pts_ray[ioffset+2*jray+1]
       ipoint=ipoint+2
    plt.plot(zpos,xpos,'black',linewidth=2)
    ioffset=ioffset+ipoint
    
plt.show()
