!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing synthetic travel times for an homogeneous medium
! using configuration files of acquisition of tomo2D.
!
! Output file: fcal.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program temps_homo
  implicit none

  include 'model_mesh.h'
  include 'model_type.h'
  include 'acqui_type.h'
  include 'data_type.h'

  type(acqui_type) :: acqui

  type(data_type) :: data

  type(model_mesh) :: mesh
  type(model_type) :: model

  real(kind=4) :: temps,speed
  integer(kind=4) :: isrc,irec,init,idata

  character*255 :: receiver_file,source_file,model_file

  !=========================================== model specification
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !======================== done here for making easier analysis of other files
  write(*,*) 'enter origin (xmin,zmin) and sampling (dx,dz)'
  read(*,*) mesh%xmin,mesh%zmin,mesh%dx,mesh%dz
  write(*,*) 'enter nx,nz'
  read(*,*) mesh%nx,mesh%nz
  mesh%xmax=mesh%xmin+float(mesh%nx-1)*mesh%dx
  mesh%zmax=mesh%zmin+float(mesh%nz-1)*mesh%dz
  mesh%ndof=mesh%nx*mesh%nz
  !======================================= allocate memory for velocity array
  allocate(model%vel(mesh%ndof))

  write(*,*) ' model box xmin,xmax,zmin,zmax: ',mesh%xmin,mesh%xmax,mesh%zmin,mesh%zmax
  ! true box where should be sources and receivers
  mesh%xleft=2*mesh%dx
  mesh%xright=(mesh%nx-2)*mesh%dx ! we have two layers around the box
  mesh%zleft=2*mesh%dz
  mesh%zright=(mesh%nz-2)*mesh%dz
  write(*,*) ' Attention: sources and receivers should be between'
  write(*,*) ' along x',mesh%xleft,mesh%xright
  write(*,*) ' along z',mesh%zleft,mesh%zright
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !================================================ acquisition
  write(*,*) 'enter receiver network file name'
  read(*,'(a)') receiver_file
  write(*,*) 'enter source network file name'
  read(*,'(a)') source_file
  !============================= read receiver/source file and check if OK
  init=0
  call read_src_rec(source_file,receiver_file,acqui,mesh,model,init)
  allocate(acqui%xsrc(acqui%nsrc))
  allocate(acqui%zsrc(acqui%nsrc))
  allocate(acqui%id_src(acqui%nsrc))
  allocate(acqui%xrec(acqui%nrec))
  allocate(acqui%zrec(acqui%nrec))
  allocate(acqui%id_rec(acqui%nrec))
  init=1
  call read_src_rec(source_file,receiver_file,acqui,mesh,model,init)
  !================================================ model
  write(*,*) 'enter model file name'
  read(*,'(a)') model_file
  !============================= read model
  call read_model(model_file,mesh,model)

  data%nobs=acqui%nsrc*acqui%nrec
  !===============================================
  ! homogeneous medium
  !===============================================
  speed=model%vel(1)

  write(*,*) ' speed ',speed

  idata=0
  open(7,file='fcal.dat',status='unknown')
  write(7,*) data%nobs
  do isrc=1,acqui%nsrc
     do irec=1,acqui%nrec
        temps=dsqrt(dble((acqui%xsrc(isrc)-acqui%xrec(irec))**2+(acqui%zsrc(isrc)-acqui%zrec(irec))**2))/speed
        idata=idata+1                       ! creation of an ID for the data
        write(7,*) isrc,irec,temps,idata
     enddo
  enddo
  close(7)
  deallocate(acqui%xsrc)
  deallocate(acqui%zsrc)
  deallocate(acqui%id_src)
  deallocate(acqui%xrec)
  deallocate(acqui%zrec)
  deallocate(acqui%id_rec)

  deallocate(model%vel)
  stop
end program temps_homo
