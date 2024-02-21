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
!========================================================
!
! The outer loop over models in order to compute new Fréchet
! derivative is missing because of the linear relation
! of the tomography we assume: Fréchet derivatives are independent
! of the model (only lengths of straight ray segment)
!
!########################################################
program tomo2D_simple

  USE lsmrdatamodule
  USE lsmrreverse

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
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !                       interface with LSMR    double precision ...
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  type ( lsmr_keep_type ) :: keep 
  ! A variable of type lsmr_keep_type that is used to
  ! preserve internal data between reverse-communication
  ! calls. The components are private.

  type ( lsmr_options_type ) :: options
  ! A variable of type lsmr_options_type that is used to control the options.
  ! It should not be altered by the user.

  type ( lsmr_inform_type ) :: inform
  ! A variable of type lsmr_info_type that is used to hold information.
  ! It should not be altered by the user.
  real(kind=dp),dimension(:),allocatable :: data_res      ! data projection
  real(kind=dp),dimension(:),allocatable :: data_cur      ! data projection
  real(kind=dp),dimension(:),allocatable :: model_cur     ! current point of the model
  real(kind=dp),dimension(:),allocatable :: model_upd     ! updated model
  real(kind=dp) :: damp
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !                       interface with LSMR    end ...
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  real(kind=4) :: fcost                                    ! cost function value
  integer(kind=4) :: ioption,idata,iter_max,init,action,index,imodel

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
  if(ioption == 1) then
     allocate(model%frechet_single(mesh%ndof))    ! allocate space for Frechet derivative (one data)
  endif
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

  if(ioption == 1) then
     !===================@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ should be sparse not yet
     allocate(model%frechet(data%nobs,mesh%ndof))  ! allocate space for total Frechet derivative
     !===================@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  endif

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
     call ray_over_data(fcost,model_cur,acqui,data,mesh,model,ioption)
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
  !§§§§§§§§§§§§§§§§§§§§§§                   INVERSION
  !§§§§§§§§§                                INVERSION
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  if(ioption == 1 ) then

     allocate(data_cur(data%nobs))   ! allocate the vector for storing A'* data_perturbation
     allocate(data_res(data%nobs))   ! allocate the vector for rhs
     allocate(model_upd(mesh%ndof)) !

!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!
!!! THIS WILL BE THE POINT WHERE THE OUTER MAIN LOOK SHOULD BE WHEN NON-LINEAR TOMO
!!!
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


     !-------------------------------------------------------------------------!
     ! computation of the synthetic data  and Frechet derivatives              !
     !-------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PLEASE NOTE IT IS DONE HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTSIDE THE (LSRM) OPTIMIZATION LOOP
     call ray_over_data(fcost,model_cur,acqui,data,mesh,model,ioption)
     write(*,*) ' initial data misfit ',fcost
     write(*,*) ' Frechet derivatives have been computed '

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
        !====================== @@@@@@@@@@@@@@@@@@@@@@@ put residues into the data_res for inversion
        data_res(idata)=data%time_obs(idata)-data%time_syn(idata)
     enddo
     close(unitf)
     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     !=======================================================
     ! start here the loop for solving A du = dt 
     !                     dt == data_res         (used by LSMR)
     !                     A  == model%frechet_single (reverse communication will let this to the user)
     !                     du == model_upd        slowness perturbation (used by LSMR)
     !              data_cur and model_cur are used for internal storage by LSMR
     !=======================================================
     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     model_cur(:)=0.d0;data_cur(:)=0.0d0
     options%ctest= 3 ! we leave the stopping reason to Fong & Saunders criterion
     options%atol=1.e-3
     options%btol=1.e-3
     options%conlim=5000
     !     options%nout=11  ! for output log
     !     options%localSize=6*mesh%ndof ! no reorthogonalisation of the model vector yet
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !  CONJUGATE GRADIENT using lsmr with the Dongarra's 
     !  reverse communication architecture: action is controlling what it is asked to the user (four levels)
     !    action=0  for starting (user setup) and for ending (lsmr setup)
     !    action=1  from data to model
     !    action=2  from model to data
     !    action=10 cleaning
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !===============@@@@@@@@@@@@@@@@@@@@@@@@ action is set to zero for initialisation
     action=0
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     ! back here after reverse communication for doing what it is asked by LSMR_reverse
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
1000 continue
     call LSMR_reverse(data%nobs,mesh%ndof,action,data_cur,model_cur,data_res,damp,model_upd,keep,options,inform)
     !==============@@@@@@@@@@@@@@@@@@@@@@@@ return value of action will be zero when ending the CG method
     if(inform%itn >= iter_max) then
        write(*,*) ' ending with maximum iterations without convergence',iter_max,inform%itn
        action=0    ! we have reached the maximum of iterations set action to 0 for getting out
     endif
     if(action == 0) then
        write(*,*) ' reason for termination see doc of LSMR ',inform%istop
        write(*,*) ' number of iterations ',inform%itn
        write(*,*) ' norm of rhs ',inform%normb
        write(*,*) ' norm of A ',inform%norma
        write(*,*) ' condition number of A ',inform%conda
        write(*,*) ' norm of residus ',inform%normr
        write(*,*) ' norm of normal residus ',inform%normAr
        write(*,*) ' norm of solution ',inform%normx
        goto 2000
     endif
     !==============@@@@@@@@@@@@@@@@@@@@@@@@ requesting the following product  data -> model
     if(action == 1) then    ! model_cur=model_cur + model%frechet_trans * data_cur    v= v + A' u
        do index=1,mesh%ndof
           do idata=1,data%nobs
              model_cur(index)=model_cur(index)+model%frechet(idata,index)*data_cur(idata)
           enddo
        enddo
     !==============@@@@@@@@@@@@@@@@@@@@@@@@ requesting the following product model -> data
     endif
     if(action == 2) then    ! data_cur= data_cur + model%frechet * model_cur    v= v + A' u
        do idata=1,data%nobs
           do index=1,mesh%ndof
              data_cur(idata)=data_cur(idata)+model%frechet(idata,index)*model_cur(index)
           enddo
        enddo

     endif
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     ! reverse communication back to LSMR_reverse
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     goto 1000
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     ! end of the linear inversion ...
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
2000 continue


!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!
!!! THIS WILL BE THE POINT WHERE THE OUTER MAIN LOOK SHOULD BE WHEN NON-LINEAR TOMO
!!!
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


     !===================================================================
     !=============@@@@@@@@@@@@@@@@@@@@@@@ output the updated model
     !===================================================================
     open(7,file='model.dat',status='unknown')
     do imodel=1,mesh%ndof            
        !================================ build final slowness model
        model_cur(imodel)=model_upd(imodel)+1.d0/model%vel(imodel)
        !================================ deduce final velocity model
        model_upd(imodel)=1.0d0/model_cur(imodel)          ! go back to velocity
        write(7,*) model_upd(imodel)                       ! we are interested in slowness
     enddo
     close(7)
     !================================ write it under binary format
     open(7,file='model.bin',access='direct',recl=4*mesh%ndof)
     write(7,rec=1) sngl(model_upd(:))
     close(7)

     !===================================================================
     !=============@@@@@@@@@@@@@@@@@@@@@@@ check the fcost and compute final synthetic data
     !===================================================================

     call ray_over_data(fcost,model_cur,acqui,data,mesh,model,ioption)
     write(*,*) ' fmisfit ',fcost
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

     !=============================== deallocation of structure "keep" for cleaning
     action=10
     write(*,*) ' cleaning LSMR '
     call LSMR_reverse(data%nobs,mesh%ndof,action,data_cur,model_cur,data_res,damp,model_upd,keep,options,inform)
     stop ' ending the inverse problem '
  endif ! end of the inversion step

  write(*,*) ' we should not arrive here: option should be 0 or 1 ',ioption
  stop
end program tomo2D_simple
