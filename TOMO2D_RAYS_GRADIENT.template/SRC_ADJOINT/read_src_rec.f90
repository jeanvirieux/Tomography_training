!##############################################################
! reading source and receiver network on unit 7
!
! check if they are inside the active box
!
!##############################################################

subroutine read_src_rec(src_name_file,rec_name_file,acqui,mesh,model,init)

  implicit none

  include 'acqui_type.h'
  include 'model_mesh.h'
  include 'model_type.h'

  TYPE(acqui_type) :: acqui
  TYPE(model_mesh) :: mesh
  TYPE(model_type) :: model

  integer(kind=4) :: isrc,irec,init

  character*255 src_name_file,rec_name_file

  if(init == 0) then
     !=========================================== reading the source file
     open(7,file=src_name_file,status='old')
     read(7,*) acqui%nsrc
     !=========================================== reading the source file
     open(77,file=rec_name_file,status='old')
     read(77,*) acqui%nrec
     return
  else   ! init=1
     do isrc=1,acqui%nsrc
        read(7,*,end=100) acqui%xsrc(isrc),acqui%zsrc(isrc),acqui%id_src(isrc)
        !=========================================== check if source inside active box
        if(acqui%xsrc(isrc) < mesh%xleft .or. acqui%xsrc(isrc) > mesh%xright .or. &
             acqui%zsrc(isrc) < mesh%zleft .or. acqui%zsrc(isrc) > mesh%zright) then
           write(*,*) 'STOP: source outside the active box',isrc,acqui%xsrc(isrc), acqui%zsrc(isrc)
           stop
        endif
     enddo
     do irec=1,acqui%nrec
        read(77,*,end=200) acqui%xrec(irec),acqui%zrec(irec),acqui%id_rec(irec)
        !=========================================== check if source inside active box
        if(acqui%xrec(irec) < mesh%xleft .or. acqui%xrec(irec) > mesh%xright .or.  &
             acqui%zrec(irec) < mesh%zleft .or. acqui%zrec(irec) > mesh%zright) then
           write(*,*) 'STOP: receiver outside the active box',irec,acqui%xrec(irec),acqui%zrec(irec)
           stop
        endif
     enddo
     close(7)
     close(77)
     return
  endif
100 continue
  write(*,*) ' error when reading the source file - read_src_rec '
  stop
200 continue
  write(*,*) ' error when reading the receiver file - read_src_rec '
  stop

end subroutine read_src_rec
