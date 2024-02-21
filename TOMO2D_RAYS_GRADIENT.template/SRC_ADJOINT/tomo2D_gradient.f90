!########################################################
! program of delayed travel-time tomography
!
! ================ CROSS-HOLE TOMOGRAPHY (straight lines)
! 
! - Rays are assumed to be straight lines between sources
!   and receivers (cross-hole tomography)
!
! - Linear problem ... as rays stay straight lines ...
!                      so a very simple ray tracing ...
!
!########################################################
program tomo2D_simple_adjoint

  implicit none

  include 'model_mesh.h'
  include 'model_type.h'
  include 'acqui_type.h'
  include 'data_type.h'

  type(acqui_type) :: acqui
  type(data_type) :: data
  type(model_mesh) :: mesh
  type(model_type) :: model

  character*255 :: receiver_file,source_file,data_file,model_file

  include 'optim_type.h'
  type(optim_type) :: optim                   ! data structure for the optimizer
  character*4 :: FLAG    

!!!================================== model%vel full model in velocity
  real(kind=4),dimension(:),allocatable :: model_cur     ! full model in slowness
  real(kind=4),dimension(:),allocatable :: grad_cur      ! gradient of the misfit function
  real(kind=4),dimension(:),allocatable :: grad_precon   ! precond. gradient



  real(kind=4) :: fcost                                    ! cost function value
  integer(kind=4) :: ioption,idata,iter_max,init,imodel,iter

  real(kind=4) :: damp

  !===================== file variables
  character(len=250)  :: config_file
  integer(kind=4)  :: unitf
  integer(kind=4) :: ierr
  !==========================================
  ! one source for this simulation
  !==========================================
  config_file = 'inputs'
  unitf=8
  open(unit = unitf, file = config_file, status = 'old', action='read', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) " error : could not open file : ", trim(adjustl(config_file))
     stop
  endif

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! forward modeling : ioption=0
  ! inverse problem  : ioption=1
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  read(unitf,*) ioption
  if(ioption == 0) then
     write(*,*) ' FORWARD MODELING: output data will be in the file fcal.dat'
  elseif(ioption == 1) then
     write(*,*) ' INVERSE PROBLEM: output model will be in the file model.dat'
  else
     write(*,*) ' unknown option: please check it'
     stop 'unrecognized option: 0 or 1 only'
  endif

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !=========================================== model specification
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !======================== done here for making easier analysis of other files
  read(unitf,*) mesh%xmin,mesh%zmin,mesh%dx,mesh%dz
  write(*,*) 'enter nx,nz'
  read(unitf,*) mesh%nx,mesh%nz
  mesh%xmax=mesh%xmin+float(mesh%nx-1)*mesh%dx
  mesh%zmax=mesh%zmin+float(mesh%nz-1)*mesh%dz
  mesh%ndof=mesh%nx*mesh%nz                       ! number of degrees of freedom
  write(*,*) ' model box xmin,xmax,zmin,zmax: ',mesh%xmin,mesh%xmax,mesh%zmin,mesh%zmax
  ! true box where should be sources and receivers
  mesh%xleft=2*mesh%dx                             ! working box # defined box
  mesh%xright=(mesh%nx-2)*mesh%dx ! we have two layers around the box
  mesh%zleft=2*mesh%dz
  mesh%zright=(mesh%nz-2)*mesh%dz
  write(*,*) ' Attention: sources and receivers should be between'
  write(*,*) ' along x',mesh%xleft,mesh%xright
  write(*,*) ' along z',mesh%zleft,mesh%zright
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !================================
  ! setup the initial model
  !================================
  !----------------------------------------------------!
  ! intial model                                       !
  !----------------------------------------------------!
  read(unitf,'(a)') model_file
  !----------------------------------------------------!
  allocate(model_cur(mesh%ndof))               ! allocate space for slowness
  allocate(model%vel(mesh%ndof))               ! allocate space for initial velocity
  allocate(grad_cur(mesh%ndof))
  allocate(grad_precon(mesh%ndof))
  !============================= initial velocity unit 7 is used for reading but closed after it
  call read_model(model_file,mesh,model) 
  write(*,*) ' OK for reading model file ',mesh%ndof
  !===================================== we consider the slowness as the parameter
  model_cur(:)= 1.d0/dble(model%vel(:))       ! put the initial model in the running model

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !================================================ acquisition
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  read(unitf,'(a)') receiver_file
  read(unitf,'(a)') source_file
  !============================= read rec/src files and check if OK
  !                              read rec/src units 7 and 77 are used for reading
  !============================= init=0 for getting values for allocation (done only in the main)
  init=0
  call read_src_rec(source_file,receiver_file,acqui,mesh,model,init)
  allocate(acqui%xsrc(acqui%nsrc))
  allocate(acqui%zsrc(acqui%nsrc))
  allocate(acqui%id_src(acqui%nsrc))
  allocate(acqui%xrec(acqui%nrec))
  allocate(acqui%zrec(acqui%nrec))
  allocate(acqui%id_rec(acqui%nrec))
  init=1
  call read_src_rec(source_file,receiver_file,acqui,mesh,model,init)
  write(*,*) ' OK for reading source and receiver files '

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !=============================================================== data
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  read(unitf,'(a)') data_file
  !============================= read data unit 7 is used for reading
  init=0
  call read_data(data_file,data,acqui,init)
  allocate(data%idt_src(data%nobs))
  allocate(data%idt_rec(data%nobs))
  allocate(data%idt_dat(data%nobs))
  allocate(data%time_obs(data%nobs))
  allocate(data%time_syn(data%nobs))
  init=1
  call read_data(data_file,data,acqui,init)
  write(*,*) ' OK for reading data file ',data%nobs

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !==================== damping and iteration for inversion (not used when doing fwd)
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  read(unitf,*) iter_max,damp
  close(unitf)                        ! end of the reading here

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !   mesh%scale_length to be set to the standard precision
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  mesh%scale_length=0.00005   ! valeur absolue ... to be better defined
  !                         defined at 0.1 msec for the box [0,100m] x [0,200m]
  !                         used for Frechet and travel times

  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§    FORWARD
  !§§§§§§§§§§§§§§§§§§§§§§                   FORWARD
  !§§§§§§§§§                                FORWARD
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  if(ioption == 0) then
     !----------------------------------------------------!
     ! computation of the synthetic time for acquisition  !
     !----------------------------------------------------!
     write(*,*) ' Forward problem store the new data in file fcal.dat '
     call ray_over_data(fcost,model_cur,acqui,data,mesh,grad_cur,ioption)
     write(*,*) ' End of travel time computation '
     unitf=8
     config_file='fcal.dat'
     open(unit = unitf, file = config_file, status = 'unknown', action='write', iostat=ierr)
     if (ierr /= 0) then
        write(*,*) " error : could not open file : ", trim(adjustl(config_file))
        stop
     endif
     write(unitf,*) data%nobs
     do idata=1,data%nobs
        write(unitf,*) data%idt_src(idata),data%idt_rec(idata),data%time_syn(idata),data%idt_dat(idata)
     enddo
     close(unitf)
     stop ' ending the forward problem '
  endif

  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§    INVERSION
  !§§§§§§§§§§§§§§§§§§§§§§                   INVERSION using adjoint formulation while considering rays
  !§§§§§§§§§                                INVERSION
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  if(ioption == 1 ) then

     !-------------------------------------------------------------------------!
     ! computation of the synthetic data  and gradient                         !
     !-------------------------------------------------------------------------!
     write(*,*) ' true gradient here '
     call ray_over_data(fcost,model_cur,acqui,data,mesh,grad_cur,ioption)
     write(*,*) ' initial non-normalized data misfit ',fcost

     write(*,*) ' initial normalized data misfit ',fcost

     write(*,*) ' travel time and gradient computation have been performed '
     unitf=8
     config_file='fcal_init.dat'
     open(unit = unitf, file = config_file, status = 'unknown', action='write', iostat=ierr)
     if (ierr /= 0) then
        write(*,*) " error : could not open file : ", trim(adjustl(config_file))
        stop
     endif
     write(unitf,*) data%nobs
     do idata=1,data%nobs
        write(unitf,*) data%idt_src(idata),data%idt_rec(idata),data%time_syn(idata),data%idt_dat(idata)
     enddo
     close(unitf)

     !----------------------------------------------------!
     ! parameter initialization                           !
     !----------------------------------------------------!
     FLAG='INIT'             ! first flag
     optim%niter_max=iter_max! maximum iteration number 
     optim%conv=1e-6         ! tolerance for the stopping criterion
     optim%print_flag=1      ! print info in output files 
     optim%debug=.false.     ! level of details for output files
     optim%l=15              ! max nber of stored models for lBFGS

     !----------------------------------------------------!
     ! optimization loop: while convergence not reached or!
     ! linesearch not failed, iterate                     !
     !----------------------------------------------------!
     iter=0
     grad_precon(:)=1.0
     write(*,*) 'model init', MAXVAL(model_cur)
     do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
!      call PSTD(mesh%ndof,model_cur,fcost,grad_cur,grad_precon,optim,FLAG)
       call LBFGS(mesh%ndof,model_cur,fcost,grad_cur,optim,FLAG)
        write(*,*) ' Reverse communication FLAG option (GRAD,NSTE,FAIL,CONV): ',FLAG
        if(FLAG.eq.'GRAD') then        
           iter=iter+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PLEASE NOTE IT IS DONE HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INSIDE THE (lBFGS) OPTIMIZATION LOOP 
           call ray_over_data(fcost,model_cur,acqui,data,mesh,grad_cur,ioption)
           write(*,*) ' iteration, fcost_data; fcost_data, model amp.',iter,fcost
           write(*,'(E20.10,2F20.10)') fcost,1./MAXVAL(model_cur),1./MINVAL(model_cur)
           if(iter > 15) goto 2000 ! ending ...
        endif
        if(FLAG == 'NSTE') then
          write(*,*) 'step length found: update the full model in slowness'
          ! here for writing intermediate model
        endif
     enddo
     if(FLAG .eq. 'FAIL') then
        ! Warning message
        write(*,*) 'convergence has been reached with the specified iteration number ',optim%niter_max
        write(*,*) 'See the convergence history in iterate_LB.dat for identifying the problem'
     endif
     if(FLAG .eq. 'CONV') then
        !Helpful console writings
        write(*,*) 'END OF TEST'
        write(*,*) 'See the convergence history in iterate_LB.dat if you want'
     endif
2000 continue

     !===================================================================
     !=============@@@@@@@@@@@@@@@@@@@@@@@ check the fcost and compute final synthetic data
     !===================================================================
     write(*,*) ' final fmisfit ',fcost
     call ray_over_data(fcost,model_cur,acqui,data,mesh,grad_cur,0)
     unitf=8
     config_file='fcal.dat'
     open(unit = unitf, file = config_file, status = 'unknown', action='write', iostat=ierr)
     if (ierr /= 0) then
        write(*,*) " error : could not open file : ", trim(adjustl(config_file))
        stop
     endif
     write(unitf,*) data%nobs
     do idata=1,data%nobs
        write(unitf,*) data%idt_src(idata),data%idt_rec(idata),data%time_syn(idata),data%idt_dat(idata)
     enddo
     close(unitf)

     !===================================================================
     !=============@@@@@@@@@@@@@@@@@@@@@@@ output the updated model
     !===================================================================
     open(7,file='model.dat',status='unknown')
     do imodel=1,mesh%ndof            
        !================================ build final velocity model
        model%vel(imodel)=1.d0/model_cur(imodel)     ! back to velocity
        write(7,*) model%vel(imodel)                       ! we are interested in slowness
     enddo
     close(7)
     !================================ write it under binary format
     open(7,file='model.bin',access='direct',recl=4*mesh%ndof)
     write(7,rec=1) model%vel(:)
     close(7)
     stop ' ending the inverse problem '
  endif ! end of the inversion step

  write(*,*) ' we should not arrive here: option should be 0 or 1 ',ioption
  stop
end program tomo2D_simple_adjoint
