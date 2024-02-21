
TYPE model_type

SEQUENCE

                 real(kind=4),allocatable,dimension(:)    :: vel
                 real(kind=4),allocatable,dimension(:)    :: frechet_single
                 real(kind=4),allocatable,dimension(:,:)  :: frechet

END TYPE model_type
