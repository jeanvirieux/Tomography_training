############################################################################
# Interface in Python for the C++ code of Jean-Marie Mirebeau and co-workers
#
# simple exercise for using HFM library (without going into the adg module
# using direct call to C++ codes with multiplexing files
#
# need python and HFM executables (to be installed from
# https://github.com/Mirebeau/HamiltonFastMarching to be downloaded and compiled
# C++ compiler (> 2017 convention 11.3.0 works) + cmake (3.22.1 works) +
# make (4.3 works)
#
############################################################################
#
#  verbose: option for more or less outputs   0 or 1
#  order: stencil order 2 or 4
#  xorg,zorg: grid origin
#  hcarre: square grid step
#  nx,nz: number of intervals (number of nodes nx+1 [0,nx])
#  slowness: slowness value constant in this case
#  xsrc,zsrc: source position
#  tsrc: time value at the source
#  radius: circle where the factorization is applied (-1 ... automatic)
#  export_time: 1 for time
#
############################################################################
def hfm_py(verbose,order,xorg,zorg,hcarre,nx,nz,slowness,xsrc,zsrc,tsrc,radius,export_time):
  
  import sys
  import os
  import numpy as np
  import struct
  import time

#############################################
#  option = 1 verbose
#############################################
  Fmt = open('input_Format.txt','w')
  Fmt.write('verbosity\n')
  Fmt.write('0\n')
  Fmt.write(' \n')

  Data = open('input_Data.dat','wb')
  myByte=struct.pack("d",verbose)
  Data.write(myByte)

#############################################
#  select forward engine Isotropic2?
#############################################

  Fmt.write('model\n')
  Fmt.write('-1\n')
  Fmt.write('Isotropic2\n')
  Fmt.write(' \n')

#############################################
#  Stencil ordering
#############################################

  Fmt.write('order\n')
  Fmt.write('0\n')   # a scalar
  Fmt.write(' \n')

  myByte=struct.pack("d",order)
  Data.write(myByte)

#############################################
#  array ordering   RowMajor ... in this case
#############################################

  Fmt.write('arrayOrdering\n')
  Fmt.write('-1\n')   # 
  Fmt.write('RowMajor\n')
  Fmt.write(' \n')

#############################################
#  define origin 
#############################################

  Fmt.write('origin\n')
  Fmt.write('1\n')   # a 2D vector
  Fmt.write('2\n')   # two coordinates
  Fmt.write(' \n')

  myByte=struct.pack("d",xorg)
  Data.write(myByte)
  myByte=struct.pack("d",zorg)
  Data.write(myByte)

#############################################
#  define gridScale
#############################################

  Fmt.write('gridScale\n')
  Fmt.write('0\n')   # a scalar
  Fmt.write(' \n')

  myByte=struct.pack("d",hcarre)
  Data.write(myByte)

#############################################
#  define dims
#############################################

  Fmt.write('dims\n')
  Fmt.write('1\n')   # a 2D vector
  Fmt.write('2\n')   # two coordinates
  Fmt.write(' \n')

  myByte=struct.pack("d",nx)
  Data.write(myByte)
  myByte=struct.pack("d",nz)
  Data.write(myByte)

#############################################
#  define cost (slowness)
#############################################

  Fmt.write('cost\n')
  Fmt.write('2\n')   # two inputs
  Fmt.write(str(nx)+'\n')   # first value
  Fmt.write(str(nz)+'\n')   # second value
  Fmt.write(' \n')

  Data.write(slowness)      # already packed
############################ could be done faster
#  for ix in range(0,100):
#     for iz in range(0,100):
#        myByte=struct.pack("d",slowness)
#        Data.write(myByte)
        
#############################################
#  define seed (source)
#############################################

  Fmt.write('seeds\n')
  Fmt.write('2\n')   # two inputs
  Fmt.write('1\n')   # one point for source
  Fmt.write('2\n')   # two coordinates of this point
  Fmt.write(' \n')

  myByte=struct.pack("d",xsrc)
  Data.write(myByte)
  myByte=struct.pack("d",zsrc)
  Data.write(myByte)

#############################################
#  define seedvalue (=0)  (at the source)
#############################################

  Fmt.write('seedValues\n')
  Fmt.write('1\n')   # one input
  Fmt.write('1\n')   # one component
  Fmt.write(' \n')

  myByte=struct.pack("d",tsrc)
  Data.write(myByte)

  Fmt.write('factoringRadius\n')
  Fmt.write('0\n')   # a scalar
  Fmt.write(' \n')

  myByte=struct.pack("d",radius)
  Data.write(myByte)

#############################################
#  define export values
#############################################

  Fmt.write('exportValues\n')
  Fmt.write('0\n')   # a scalar
  Fmt.write(' \n')

  myByte=struct.pack("d",export_time)
  Data.write(myByte)

#############################################
#  close and run using the multiplex input file
#############################################

  Fmt.close()
  Data.close()


  print('using the multiplex input file')
  print('index file input_Format.txt')
  print('binary file input_Data.dat')

  print('runnning C++ code')
  os.system('./FileHFM_Isotropic2.exe')
  print('runnning C++ code OK')

  print('decoding the multiplex output file')
  print('index file output_Format.txt')
  print('binary file output_Data.dat')


#############################################
# direct access of the output data
# construction of indexes ... for extraction
#############################################

  clef=[]
  direct=np.empty([50])
  lec_dim=np.empty([3,50])

  iflag=0
  idirect=-1
  irec_hfm=1

  time.sleep(2.0)
  
  Fmt = open('output_Format.txt','r')
  while iflag == 0:
    
    line=Fmt.readline()
    coords=line.split()
    string=str(coords[0])
    line=Fmt.readline()
    coords=line.split()
    iclef=int(coords[0])
  
    if string == 'visitedUnset':    # freeze everything
       iflag=1
       iclef=-2

#    print('clef',string)
#    print('type clef',iclef)
      
    if iclef == -1:   # just an information
       line=Fmt.readline()
       coords=line.split()
       line=Fmt.readline()
       coords=line.split()   # str(coords[0]) should be zero
#       print('end',str(coords[0]))

    if iclef == 0:   # just a scalar
       idirect=idirect+1
       clef.append(string)
       direct[idirect]=irec_hfm
       irec_hfm=irec_hfm+1
       lec_dim[:,idirect]=-1
       line=Fmt.readline()
       coords=line.split()   # str(coords[0]) should be zero
#       print('end',str(coords[0]))

    if iclef == 1:  # a   
       idirect=idirect+1
       clef.append(string)
       direct[idirect]=irec_hfm
       line=Fmt.readline()
       coords=line.split()
       ix=int(coords[0])
       irec_hfm=irec_hfm+ix
       lec_dim[0,idirect]=ix 
       lec_dim[1:2,idirect]=-1
       line=Fmt.readline()
       coords=line.split()   # str(coords[0]) should be zero
#       print('end',str(coords[0]))

    if iclef == 2:
       idirect=idirect+1
       clef.append(string)
       direct[idirect]=irec_hfm
       line=Fmt.readline()
       coords=line.split()
       ix=int(coords[0])
       line=Fmt.readline()
       coords=line.split()
       iz=int(coords[0])
       irec_hfm=irec_hfm+ix*iz
       lec_dim[0,idirect]=ix   
       lec_dim[1,idirect]=iz 
       lec_dim[2,idirect]=-1 
       line=Fmt.readline()
       coords=line.split()   # str(coords[0]) should be zero
#       print('end',str(coords[0]))
       
    if iclef == 3:
       idirect=idirect+1
       clef.append(string)
       direct[idirect]=irec_hfm
       line=Fmt.readline()
       coords=line.split()
       ix=int(coords[0])
       line=Fmt.readline()
       coords=line.split()
       iz=int(coords[0])
       line=Fmt.readline()
       coords=line.split()
       it=int(coords[0])
       irec_hfm=irec_hfm+ix*iz*it
       lec_dim[0,idirect]=ix  
       lec_dim[1,idirect]=iz 
       lec_dim[2,idirect]=it 
       line=Fmt.readline()
       coords=line.split()   # str(coords[0]) should be zero
#       print('end',str(coords[0]))
      
  Fmt.close()
  idirect_save=idirect+1
#  print('save direct',idirect_save)

#######################################################
# direct-access through keys of the output files of HFM
#######################################################
  Data = open('output_Data.dat','rb')

#######################################################
# extraction of times   (values)
#######################################################

  for idirect in range(0,idirect_save):
#       print('clef news',idirect,clef[idirect]) 
       if clef[idirect] == 'values':
         nxc=int(lec_dim[0,idirect])
         nzc=int(lec_dim[1,idirect])
         print('array nx,nz',nxc,nzc)
         irec_hfm=direct[idirect]
         Data.seek(int(irec_hfm-1)*8,0)
############################################# testing the seek OK       
#       for ix in range(0,int(irec_hfm)-1):
#             myByte=Data.read(8)
#             print('myByte',struct.unpack('d',myByte))
#       myByte=Data.read(8)
#       print('myByte',struct.unpack('d',myByte))
#############################################
#       myByte=bytes(8*nxc*nzc)   #seems to be done automatically
         myByte=Data.read(8*nxc*nzc)
         field=struct.unpack('d'*nxc*nzc,myByte)
         del(myByte)
#         time=np.asarray(field)
#         time=time.reshape(nxc,nzc)

#######################################################
# extraction of other quantities if needed
#####################################################
       
#######################################################
# direct-access through keys of the output files of HFM
#######################################################       
  Data.close()
  return field
