
TYPE data_type

SEQUENCE

integer(kind=4) :: nobs
		 integer(kind=4),allocatable,dimension(:)  :: idt_src,idt_rec,idt_dat
		 real(kind=4),allocatable,dimension(:)    :: time_obs,time_syn,rhs


END TYPE data_type
