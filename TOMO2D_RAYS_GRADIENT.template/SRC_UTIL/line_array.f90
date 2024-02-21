program line_array

write(*,*) ' enter the number of positions '
read(*,*) npos

write(*,*) ' enter the initial position of the line array '
read(*,*) xorg,zorg

write(*,*) ' enter the vertical stepping '
read(*,*) dz
write(*,*) ' enter the horizontal stepping '
read(*,*) dx

write(7,*) npos

do ipos=1,npos
write(7,*) xorg+(ipos-1)*dx,zorg+(ipos)*dz,ipos
enddo

end program line_array
