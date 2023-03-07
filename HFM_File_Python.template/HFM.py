#####################################################################################
#                                                                                   #
#  example using a more complex model                                               #
#                                                                                   #
# python HFM.py file name nz nx  dcarre zorg xorg zsrc xsrc  (z first and then x)   #
# python HFM.py nankai_double.bin 67 365 468.75 -468.75 -468.75 30000. 75000.       #
#                                                                                   #
#####################################################################################

import sys
import os
import numpy as np
import struct
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import def_HFM


verbose=1.         # could be zero for less output  (useless when doing it through this python)
order=2.           # minimal stencil order
tsrc=0.            # time at seed
radius=-1.         # automatic factorization
export_time=1.     # output time field

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
plt.imshow(velocity,origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[zmin,zmax,xmax,xmin],cmap=cm.seismic,clip_on=True)

cb=plt.colorbar(orientation='horizontal')
cb.set_label(label='(km/s)',size=22) #'large')
cb.ax.tick_params(labelsize=20)

plt.xlabel('Range (km)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.title('P wave velocity',fontsize=28)

plt.show()

############# pack into a buffer
slowness=struct.pack('d'*nxc*nzc,*invvel)    # flatten the list for packing 

############# calling the HFM function
field=def_HFM.hfm_py(verbose,order,xorg,zorg,hcarre,nx,nz,slowness,xsrc,zsrc,tsrc,radius,export_time)

############# get back the array from the buffer
time=np.asarray(field)
time=time.reshape(nxc,nzc)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)

plt.imshow(time,origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[zmin,zmax,xmax,xmin],cmap=cm.seismic,clip_on=True)

cb=plt.colorbar(orientation='horizontal')
cb.set_label(label='(sec)',size=22) #'large')
cb.ax.tick_params(labelsize=20)

plt.xlabel('Range (km)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.title('Travel time',fontsize=28)

plt.show()
