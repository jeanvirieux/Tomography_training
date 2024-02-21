program model_homo

implicit none

integer(kind=4) :: nx,nz
real(kind=4) :: value

real(kind=4), dimension(:,:), allocatable :: vel

integer(kind=4) :: ix,iz

write(*,*) ' enter the number of grid points in x and in z '
read(*,*) nx,nz
write(*,*) ' enter the homogeneous value of velocity '
read(*,*) value

allocate(vel(nx,nz))
open(8,file='model_homo.bin',access='direct',recl=4*nx*nz)

open(7,file='model_homo.dat')

do iz=1,nz
do ix=1,nx
vel(ix,iz)=value
write(7,*) value
enddo
enddo
close(7)

write(8,rec=1) vel
close(8)

deallocate(vel)

stop
end program model_homo
