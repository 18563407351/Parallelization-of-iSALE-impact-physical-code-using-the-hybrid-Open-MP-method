!> @mainpage iSALE3D
!! ===================================================================
!!                     iSALE-3D
!!       Lagrangian/Eulerian Hydrocode for impact
!!            modelling in three dimensions
!!
!!                by Dirk Elbeshausen
!!            Museum fuer Naturkunde Berlin
!!            dirk.elbeshausen@mfn-berlin.de
!!
!!               http://www.isale-code.de
!!
!! This code uses a solver as described in Hirt et al., 1974.
!! EOS and material models used in this code have been developed by
!! Kai Wünnemann, Gareth S. Collins, Boris Ivanov, and H. Jay Melosh.
!!
!! @authors Dirk Elbeshausen, MfN
!<----------------------------------------------------- 
program isale3d
  use mod_isale
  use mod_time
  use ptool_interface
  use mod_parallel
  use mod_identify_mat
  use mod_mpi_exchange
  use mod_regrid
  implicit none

  out   = STDOUT ! for the very first beginning...
  ioerr = STDOUT
  ! first start with initializing MPI (if required)...
  call initialize_mpi()
  ! analyze flags set by the user...
  call ANALYZE_FLAGS()
  ! just some salutation (at the moment)...
  call SETUP_CREDITS()
  ! read material and parameter files...
  call read_param()
  call setup_derive_param()
  ! acquire specific user information
  call SETUP_USER_INFO()
  ! starting space decomposition
  call SPACE_DECOMP()
  ! allocate memory, init fields, perform some pre-mainloop-calculations
  call SETUP_INITIALIZE()
  ! initialize new data structures
  ! (intermediate step until implementation of 
  ! object oriented approach is finished)
  call init_cells()
  print*,1
  call CELSET()
  if (dumpname(1:4) == "NONE") then
     ! setup grid, projectile and layer
     ! calculate pressure gradient (if required)...
     ! and place tracers...
  else
     ! restart simulation from dump-file
     if (regridmesh .eq. 1) then
       call READ_DUMP_regrid(DUMPNAME)
     else
       call READ_DUMP(DUMPNAME) 
     endif
  endif
  call ptool_text(io_all,"finished reading dump!","")
  ! here we want to be sure that all processes are synchronous...
  call BARRIER() 

  ! now write header and grid information
  ! into jpeg-file (only done by primary thread)
  !if (pid == PRIMARY) call WRITE_JPEG_INFO()
  ! write the first timestep into jpeg file (except when reading from dump)
  ! and only, if running in debug-mode. The first timestep is usually written
  ! directly after entering the main-loop. Writing here helps to find bugs
  ! introduced by the setup-routines...
#ifdef DEBUG_MODE
  if (dumpname(1:4) == "NONE" .and. writeall .ne. 1) call WRITE_JPEG_STEP(1)
#endif
  ! calculating the volume of each cell once before starting
  ! - for eulerian runs (ALE_MODE==EULERIAN) this will be the only time...
  call UPDATE_VOLUME()
  ! how many materials are in cell??
  IF (ALE_MODE /= LAGRANGIAN) call IDENTIFY_MAT(eps_min) 
  ! once before starting, exchange all ghostpoints
  call EXCHANGE_ALL() 
  ! performing an MPI_TEST (only if environment MPI_TEST is set)
  call MPITEST_MAIN() 
  ! ------------------------------------------------------------------
  ! here we are entering the main loop now...
  ! all action is taking place in subroutines called from here...
  ! note: be careful when adding new routines 
  ! -> check if you need to exchange some ghostpoints...
  ! (let's try following method: exchange ghostpoints before entering
  ! the routine...)
  ! ------------------------------------------------------------------

  ! If just performing a set-up check end here
  if (pid == PRIMARY) call WRITE_JPEG_INFO()
  call info_stop_after_setup()

  ! aquire and store start time information
  ! TMD: gettime uses the MPI-safe MPI_Wtime function if compiled with MPI
  time_start = gettime()
  time_old = time_start    

#ifdef DEBUG_MODE
  call ptool_text(io_oin,"A T T E N T I O N","RUNNING IN DEBUG MODE!")
#endif

  DO
     ! calculate the new time increment...
     call UPDATE_TIMESTEP()
     ! initialize the new timestep (something already done at end of advect)...
     call UPDATE_STATE()
     ! now calculate the new tracer positions and values (if required...)
     ! NOTE: for efficiency enhancements, we exchange tracers in timstp ...
     if (tracer_motion .eq. TR_VEL) call MOVE_TRACER()
     ! start the new iteration (write some data etc.)
     call UPDATE_CYCLE()
     ! compute the deviatoric stress tensor (if required)...
     if (calc_stress == 1) call UPDATE_VELOCITY_STRESS()      
     ! calculate the new velocity, calculate artificial viscosity
     ! consider artificial viscosity and stress part...
     call UPDATE_VELOCITY_PRESSURE()
     ! calculate specific internal energy...
     call UPDATE_ENERGY()
     ! perform 'advection':
     ! if ALE_MODE==EULERIAN (eulerian step) perform full donor cell advection...
     if(ALE_MODE == EULERIAN) then    
        call ADVECT()
     else ! completely lagrangian or mixed step...
        call MYSTOP("ONLY EULERIAN CASE AT THE MOMENT!")
     endif

     ! add ncyc, dt, time into statistics file...
     ! if statistics > 1 calculate further statistics (monitor mass conservation)
     ! finally derive some runtime statistics...
     time_act = gettime()
     time_delta = time_act - time_old
     time_old = time_act
     time_inter = time_inter + time_delta
  ENDDO

end program isale3d
