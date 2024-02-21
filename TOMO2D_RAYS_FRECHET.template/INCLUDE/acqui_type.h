TYPE acqui_type

SEQUENCE

integer(kind=4) :: nsrc
		 integer(kind=4),allocatable,dimension(:) :: id_src
		 real(kind=4),allocatable,dimension(:)    :: xsrc,zsrc

integer(kind=4) :: nrec
		 integer(kind=4),allocatable,dimension(:) :: id_rec
		 real(kind=4),allocatable,dimension(:)    :: xrec,zrec

END TYPE acqui_type
