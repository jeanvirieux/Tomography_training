!#####################################################
! loop over data
!     for estimating the gradient through Frechet derivatives  NOT DONE HERE  because CG method
!     for estimating the Frechet derivative (Big Matrix)
!     for estimating the misfit function  fcost
!
!     ioption   =0 only the forward modeling (time)
!     ioption   =1 both time and Frechet
!#####################################################
subroutine ray_over_data(fcost,model_cur,acqui,data,mesh,model,ioption)

  use lsmrdatamodule

  implicit none

  include 'acqui_type.h'
  include 'data_type.h'
  include 'model_mesh.h'
  include 'model_type.h'

  type(model_mesh) :: mesh
  type(model_type) :: model
  type(acqui_type) :: acqui
  type(data_type) :: data

  real(kind=4) :: fcost

  real(kind=dp) :: model_cur(mesh%ndof)
  real(kind=4),dimension(:),allocatable :: frechet_single
  real(kind=4) :: time_syn

  integer(kind=4) :: idata,index,ioption
  integer(kind=4) :: iflag_src,iflag_rec,isrc,irec

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
        call ray_sngl_data_src_rec(fcost,idata,model_cur,mesh,acqui,iflag_src,iflag_rec,time_syn,frechet_single,ioption)
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
     allocate(frechet_single(mesh%ndof))
     model%frechet(:,:)=0.
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
        call ray_sngl_data_src_rec(fcost,idata,model_cur,mesh,acqui,iflag_src,iflag_rec,time_syn,frechet_single,ioption)
        !======================================== fill-in the complete Frechet derivative
        do index=1,mesh%ndof
           model%frechet(idata,index)=frechet_single(index)    ! full storage
        enddo
        !======================================== estimate the misfit function L2
        data%time_syn(idata)=time_syn
        fcost=fcost+(data%time_obs(idata)-time_syn)**2
     enddo
     fcost=fcost/float(data%nobs)
     deallocate(frechet_single)
     return
  endif


  write(*,*) ' ray_over_data - ioption should be 0 or 1 ',ioption
  stop
end subroutine ray_over_data

