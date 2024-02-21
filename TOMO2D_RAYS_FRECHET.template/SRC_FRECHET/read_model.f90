!##############################################################
! reading model
!
!         some quantities have been red from standard input
!         we consider the velocity which will be transformed in slowness
!         we reserve space for Frechet derivative (sensitivity kernel only)
!
!##############################################################

subroutine read_model(model_name_file,mesh,model)

  implicit none

  include 'model_mesh.h'
  include 'model_type.h'

  type(model_mesh) :: mesh
  type(model_type) :: model

  character*255 model_name_file

  integer(kind=4) :: imodel

  !=========================================== reading the source file
  open(7,file=model_name_file,status='old')

  do imodel=1,mesh%ndof                        ! read values of velocity
     read(7,*) model%vel(imodel)          ! we are interested in slowness
  enddo

  close(7)

end subroutine read_model
