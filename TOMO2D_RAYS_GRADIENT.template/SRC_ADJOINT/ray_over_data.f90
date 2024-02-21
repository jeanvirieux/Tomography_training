!#####################################################
! loop over data
!     for estimating travel times data%time_syn
!     for estimating the gradient grad_cur through adjoint formulation 
!     for estimating the misfit function  fcost
!
!     ioption   =0 only the forward modeling (time)
!     ioption   =1 both time and Frechet
!#####################################################
subroutine ray_over_data(fcost,model_cur,acqui,data,mesh,grad_cur,ioption)

  implicit none

  include 'acqui_type.h'
  include 'data_type.h'
  include 'model_mesh.h'

  type(model_mesh) :: mesh
  type(acqui_type) :: acqui
  type(data_type) :: data

  real(kind=4) :: fcost

  real(kind=4) :: model_cur(mesh%ndof)
  real(kind=4) :: grad_cur(mesh%ndof)

  real(kind=4),dimension(:),allocatable :: frechet_single
  real(kind=4) :: time_syn

  integer(kind=4) :: idata,index,ioption
  integer(kind=4) :: iflag_src,iflag_rec,isrc,irec

  integer(kind=4),save :: irecord=0,ix,iz
 


  fcost=0.

  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§    FORWARD
  !§§§§§§§§§§§§§§§§§§§§§§                   FORWARD
  !§§§§§§§§§                                FORWARD
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

  if(ioption == 0) then
     data%time_syn(:)=0.
     do idata=1,data%nobs
        !================================================ check if the data is related to meaningful 'Source'
        iflag_src=-999
        do isrc=1,acqui%nsrc
           if(data%idt_src(idata) == acqui%id_src(isrc)) then
              iflag_src=isrc
              goto 101
           endif
        enddo
        !============================ source not found
        write(*,*) ' Problem with the source ID at the data ', idata
        stop 'ray_over_data'
101     continue
        !================================================ check if the data is related to a meaningful 'Receiver'
        iflag_rec=-999
        do irec=1,acqui%nrec
           if(data%idt_rec(idata) == acqui%id_rec(irec)) then
              iflag_rec=irec
              goto 201
           endif
        enddo
        !============================ receiver not found
        write(*,*) ' Problem with the receiver ID at the data ',idata
        stop 'ray_over_data'
201     continue
        !======================================== compute synthetic time
        call ray_sngl_data_src_rec(fcost,idata,model_cur,acqui,iflag_src,iflag_rec,time_syn,mesh,frechet_single,ioption)
        data%time_syn(idata)=time_syn
        data%idt_dat(idata)=idata               ! create an ID for this data
     enddo
     return
  endif

  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§    INVERSION
  !§§§§§§§§§§§§§§§§§§§§§§                   INVERSION
  !§§§§§§§§§                                INVERSION
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

  if(ioption == 1) then
     allocate(frechet_single(mesh%ndof)) ! full allocation - never in real life
     grad_cur(:)=0.
     do idata=1,data%nobs
        ! set synthetic time and Frechet derivative to zero and start the recursive investigation for time and Frechet derivatives
        frechet_single(:)=0.
        !======================================== compute synthetic time, gradient and Frechet_single
        !======================================== and cumul partial derivative for gradient
        !================================================ check if the data is related to meaningful 'Source'
        iflag_src=-999
        do isrc=1,acqui%nsrc
           if(data%idt_src(idata) == acqui%id_src(isrc)) then
              iflag_src=isrc
              goto 100
           endif
        enddo
        !============================ source not found
        write(*,*) ' Problem with the source ID at the data ', idata
        stop 'ray_over_data'
100     continue
        !================================================ check if the data is related to a meaningful 'Receiver'
        iflag_rec=-999
        do irec=1,acqui%nrec
           if(data%idt_rec(idata) == acqui%id_rec(irec)) then
              iflag_rec=irec
              goto 200
           endif
        enddo
        !============================ receiver not found
        write(*,*) ' Problem with the receiver ID at the data ',idata
        stop 'ray_over_data'
        !============================ ok validation of this data
200     continue
        call ray_sngl_data_src_rec(fcost,idata,model_cur,acqui,iflag_src,iflag_rec,time_syn,mesh,frechet_single,ioption)
        !==================================== from Frechet to gradient (we should use the sparsity by storing involved index)
        !==================================== here simplified version using global loop ...   grad = - Jt * DT
        do index=1,mesh%ndof
           grad_cur(index)=grad_cur(index)-frechet_single(index)*(data%time_obs(idata)-time_syn)
        enddo
!        write(99,*) ' data frechet ',maxval(frechet_single),data%time_obs(idata)-time_syn
        !======================================== estimate the misfit function L2
        data%time_syn(idata)=time_syn
        fcost=fcost+(data%time_obs(idata)-time_syn)**2
     enddo
     grad_cur(:)=grad_cur(:)/float(data%nobs)
     write(*,*) ' gradient ',maxval(abs(grad_cur))


!!JEAN patch for adjoint 

!     open(1,file='gradient.bin',recl=4*mesh%ndof,access='direct')
     open(1,file='gradient.bin',recl=4,access='direct')
     do ix=1,101
     do iz=1,201
     index=ix+101*(iz-1)
     irecord=irecord+1
     write(1,rec=irecord) grad_cur(index)
!     write(1,rec=irecord) (grad_cur(index),index=1,mesh%ndof)
     enddo
     enddo
     close(1)




     fcost=fcost/float(data%nobs)
     deallocate(frechet_single)
     return
  endif


  write(*,*) ' ray_over_data - ioption should be 0 or 1 ',ioption
  stop
end subroutine ray_over_data

