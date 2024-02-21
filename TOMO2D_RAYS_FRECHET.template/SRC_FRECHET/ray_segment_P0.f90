!######################################################
! moving along the ray between two points
!  - check if they belong to the same cell
!  - if yes, compute sensitivity kernel and travel time
!  - if not, divide the segment in two segments 
!           and proceed recursively
!
! we assume a straight line between the two points
! this depends on the ray sampling.
!
!============== trigger the recursive status of this subroutine
!
!
!============== compute travel time (ioption=0) 
!============== compute travel time and Frechet derivative (ioption=1)
!
!######################################################
recursive subroutine ray_segment(idata_cur,model_cur,xbeg,zbeg,xend,zend,xleng_total,time_syn,mesh,frechet_single,ioption,irecur)

  use lsmrdatamodule

  implicit none

  include 'model_mesh.h'

  TYPE (model_mesh) :: mesh

  real(kind=dp) :: model_cur(mesh%ndof)
  real(kind=dp) :: xlenteur_beg,xlenteur_end

  integer(kind=4) :: idata_cur,index
  integer(kind=4) :: ixbeg,izbeg,ixend,izend
  real(kind=4) :: xleng_segment,xleng_total

  real(kind=4) :: xbeg,zbeg
  real(kind=4) :: xend,zend
  real(kind=4) :: xmid,zmid

  integer(kind=4) :: ix_err,iz_err

  integer(kind=4) :: ioption   ! ioption=0 only time ioption =1 time and frechet
  integer(kind=4),parameter :: irecur_max=20
  integer(kind=4) :: irecur,irecur1,irecur2

  real(kind=4),dimension(mesh%ndof) :: frechet_single
  real(kind=4) :: time_syn

  !========================================= check if same cell
  !================= index of two ends of the segment
  ixbeg=(xbeg-mesh%xmin)/mesh%dx + 1
  izbeg=(zbeg-mesh%zmin)/mesh%dz + 1
  ixend=(xend-mesh%xmin)/mesh%dx + 1
  izend=(zend-mesh%zmin)/mesh%dz + 1

  call check_bout(mesh%nx,mesh%nz,ixbeg,izbeg,ix_err,iz_err)

  !================= length of the segment
  xleng_segment=sngl(dsqrt(dble((xbeg-xend)**2+(zbeg-zend)**2)))
  if(ixbeg == ixend .and. izbeg == izend) then   ! same cell for both ends of the ray segment
     xleng_total=xleng_total+xleng_segment
     !================================== compute the frechet derivative for this data
     !================================== left bottom corner
     index=mesh%nx*(izbeg-1)+ixbeg
     if(ioption == 1) then
        !================================== simple slowness description
        frechet_single(index)=frechet_single(index)+xleng_segment
     endif
     !================================== compute the contribution to the synthetic time 
     call slowness(xbeg,zbeg,ixbeg,izbeg,model_cur,mesh,xlenteur_beg)
     call slowness(xend,zend,ixbeg,izbeg,model_cur,mesh,xlenteur_end)
     !================================== trapeze rule for getting the integral of travel time
     time_syn=time_syn+xleng_segment*(xlenteur_beg+xlenteur_end)*0.5
  else
     !======================================= we may end if segment too small
     if(xleng_segment < mesh%scale_length) then
       return
     endif
     if(irecur > irecur_max) return           ! limit the recursive depth for dynamic memory allocation
     irecur=irecur+1
     !======================================= divide by two and recursive call
     xmid=0.5*(xbeg+xend)
     zmid=0.5*(zbeg+zend)
     irecur1=irecur
     call ray_segment(idata_cur,model_cur,xbeg,zbeg,xmid,zmid,xleng_total,time_syn,mesh,frechet_single,ioption,irecur1)
     irecur2=irecur
     call ray_segment(idata_cur,model_cur,xmid,zmid,xend,zend,xleng_total,time_syn,mesh,frechet_single,ioption,irecur2)
  endif
  return
end subroutine ray_segment
!============================================
!
! bilinear interpolation of slowness inside the same element of the mesh
! return the value xlenteur at the position (x,z)
!============================================
subroutine slowness(x,z,ixbeg,izbeg,model_cur,mesh,xlenteur)

  use lsmrdatamodule

  implicit none

  include 'model_mesh.h'
  TYPE (model_mesh) :: mesh

  real(kind=dp) :: model_cur(mesh%ndof)

  integer(kind=4) :: ixbeg,izbeg,index
  real(kind=4) :: x,z

  real(kind=4) :: xorg_local,zorg_local
  real(kind=dp) :: xx,zz
  real(kind=dp) :: xvalue(4,4),xlenteur

  index=mesh%nx*(izbeg-1)+ixbeg
  xvalue(1,1)=model_cur(index)
  index=mesh%nx*(izbeg-1)+ixbeg+1
  xvalue(2,1)=model_cur(index)
  index=mesh%nx*(izbeg)+ixbeg
  xvalue(1,2)=model_cur(index)
  index=mesh%nx*(izbeg)+ixbeg+1
  xvalue(2,2)=model_cur(index)

  xorg_local=float(ixbeg-1)*mesh%dx+mesh%xmin
  zorg_local=float(izbeg-1)*mesh%dz+mesh%zmin
  xx=(x-xorg_local)/mesh%dx
  zz=(z-zorg_local)/mesh%dz
  xlenteur=xvalue(1,1)*(1.0d0-xx)*(1.0d0-zz)  &
          +xvalue(2,1)*    xx *(1.0d0-zz)  &
          +xvalue(1,2)*(1.0d0-xx)*    zz  &
          +xvalue(2,2)*    xx*     zz
return
end subroutine slowness
!==========================================
!
! bound checks for linear interpolation
!
!==========================================
subroutine check_bout(nx,nz,ix,iz,ix_err,iz_err)
integer(kind=4) :: nx,nz,ix,iz,ix_err,iz_err
ix_err=0      ! OK
iz_err=0      ! OK
if(ix < 2) then
  ix_err=-1
  ix=2
endif
if(ix > nx-1) then
  ix_err=1
  ix=nx-1
endif
if(iz < 2) then
  iz_err=-1
  iz=2
endif
if(iz > nz-1) then
  iz_err=1
  iz=nz-1
endif
return
end subroutine check_bout
