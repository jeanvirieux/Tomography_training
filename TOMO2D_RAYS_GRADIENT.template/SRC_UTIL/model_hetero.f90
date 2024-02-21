program model_hetero

implicit none

integer(kind=4) :: nx,nz
real(kind=4) :: dx,dz,valeur,dvalue,speed
real(kind=4) :: xorg,zorg,xleng,zleng,x,z
real(kind=4) :: dist_norm

real(kind=4), dimension(:,:), allocatable :: vel

integer(kind=4) :: ix,iz

write(*,*) ' enter the number of grid points in x and in z '
read(*,*) nx,nz
write(*,*) ' enter the step size dx,dz in meters '
read(*,*) dx,dz

allocate(vel(nx,nz))
open(8,file='model_hetero_ref.bin',access='direct',recl=4*nx*nz)

open(7,file='model_hetero_ref.dat')

write(*,*) ' enter the homoegenous value of velocity '
read(*,*) valeur

write(*,*) ' enter velocity anomaly m/s'
read(*,*) dvalue

write(*,*) ' enter position of anomaly (xorg,zorg) in meters'
read(*,*) xorg,zorg
write(*,*) ' enter extension (xleng,zleng) in meters'
read(*,*) xleng,zleng

do iz=1,nz
z=(iz-1)*dz
do ix=1,nx
x=(ix-1)*dx
dist_norm=((x-xorg)/xleng)**2+((z-zorg)/zleng)**2
speed=valeur+dvalue*exp(-dsqrt(dble(dist_norm)))
vel(ix,iz)=speed
write(7,*) speed
enddo
enddo

close(7)

write(8,rec=1) vel
close(8)

deallocate(vel)

stop
end program model_hetero
