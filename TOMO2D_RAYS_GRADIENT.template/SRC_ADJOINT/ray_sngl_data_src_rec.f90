!#####################################################
! loop over data
!     for estimating the gradient through Frechet derivatives
!                     done and reset for each data
!                     we do not store Frechet derivatives
!     for estimating synthetic data
!     for estimating the cost function
!     
!  fcost: cost function L2 norm
!  idata: the current data under investigation
!  model_cur: current slowness model
!
!  acqui: the acquisition structure
!  data: the data structure
!  model: the model structure
!   
!#####################################################
subroutine ray_sngl_data_src_rec(fcost,idata,model_cur,acqui,isrc,irec,time_syn,mesh,frechet_single,ioption)

  implicit none

  include 'acqui_type.h'
  include 'model_mesh.h'

  type(model_mesh) :: mesh
  type(acqui_type) :: acqui

  real(kind=4) :: model_cur(mesh%ndof)
  real(kind=4) :: frechet_single(mesh%ndof)
  real(kind=4) :: time_syn

  real(kind=4) :: fcost,xleng_total,xleng_cur
  real(kind=4) :: xbeg,zbeg
  real(kind=4) :: xend,zend

  integer(kind=4) :: idata   ! current data
  integer(kind=4) :: isrc,irec
  integer(kind=4) :: ioption   ! ioption=0 only time; ioption=1 time and Frechet
  integer(kind=4) :: irecur  ! counting for recursive loops
  integer(kind=4) :: debug=1

  if(debug == 0) write(*,*) ' contribution to the gradient of the data ',idata

  !=================================== dichotomie compute Frechet and time
  time_syn=0.      ! set synthetic time to zero and start the recursive investigation for time and Frechet derivatives
  xbeg=acqui%xsrc(isrc);zbeg=acqui%zsrc(isrc)
  xend=acqui%xrec(irec);zend=acqui%zrec(irec)
  xleng_total=sngl(dsqrt(dble((xbeg-xend)**2+(zbeg-zend)**2)))
  xleng_cur=0.
  irecur=0
  call ray_segment(idata,model_cur,xbeg,zbeg,xend,zend,xleng_cur,time_syn,mesh,frechet_single,ioption,irecur)
  if(debug == 0) write(*,*) ' total length and current length ',xleng_total,xleng_cur,abs(xleng_total-xleng_cur)
  return
end subroutine ray_sngl_data_src_rec

