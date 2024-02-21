!##############################################################
! reading data on unit 7
!
!##############################################################

subroutine read_data(data_name_file,data,acqui,init)

  implicit none

  include 'data_type.h'
  include 'acqui_type.h'

  TYPE(acqui_type) :: acqui
  TYPE(data_type) :: data

  character*255 data_name_file

  integer(kind=4) :: idata,isrc,irec,iflag_src,iflag_rec,init

  if(init == 0) then
     !=========================================== reading the source file
     open(7,file=data_name_file,status='old')
     read(7,*) data%nobs
     return
  else    ! init=1 usually ...
     do idata=1,data%nobs
        read(7,*,end=300) data%idt_src(idata),data%idt_rec(idata),data%time_obs(idata),data%idt_dat(idata)

        !================================================ check if the data is meaningful 'Source'
        iflag_src=-999
        do isrc=1,acqui%nsrc
           if(data%idt_src(idata) == acqui%id_src(isrc)) then
              iflag_src=isrc
              goto 100
           endif
        enddo
        !============================ source not found
        write(*,*) ' Specified source ID in data file is not found in the source file',idata,data%idt_src(idata)
        write(*,*) ' ID sources ',acqui%id_src(:)
        stop

100     continue
        !================================================ check if the data is meaningful 'Receiver'
        iflag_rec=-999
        do irec=1,acqui%nrec
           if(data%idt_rec(idata) == acqui%id_rec(irec)) then
              iflag_rec=irec
              goto 200
           endif
        enddo
        !============================ receiver not found
        write(*,*) ' Specified receivers ID in data file is not found in the receivers file',idata,data%idt_rec(idata)
        write(*,*) ' ID receivers ',acqui%id_rec(:)
        stop

200     continue
     enddo
     close(7)
     return
  endif

300 continue
  write(*,*) ' error when reading the data file - read_data '
  stop
end subroutine read_data
