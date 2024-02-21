
TYPE model_mesh

SEQUENCE
                 integer(kind=4) :: nx,nz,ndof
                 real(kind=4)    :: dx,dz
		 real(kind=4)    :: xmin,xmax,zmin,zmax
		 real(kind=4)    :: xleft,xright,zleft,zright
                 real(kind=4)    :: scale_length    ! minimum length

END TYPE model_mesh
