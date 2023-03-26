module mod_setup_params
  use mod_ndim
  use mod_io, only : MAXPROJ
  implicit none
  
  integer :: obj_num,proj_num,par_num
 
  type type_obj

     character(len=10) :: name         ! material name
     character(len=12) :: type         ! projectile type (shape)
     character(len=12) :: tprof        ! Thermal profile flag
     character(len=12) :: pprof        ! Porosity profile
     character(len=12) :: par_dist     ! Type of particle distribution
     character(len=12) :: bitmap       ! Name of bitmap file for meso-scale object
     
     real*8            :: vel          ! velocity
     real*8            :: angle(2)     ! incidence angles (X/Z-dir, Y/Z-dir)
     real*8            :: vol          ! volume
     real*8            :: rad(1:ndim)  ! radius
     real*8            :: mid(1:ndim)  ! coordinates of center
     real*8            :: sie          ! specific internal energy
     real*8            :: temp         ! temperature
     real*8            :: density      ! density
     real*8            :: pressure     ! pressure
     real*8            :: damage       ! initial damage of object
     real*8            :: vvel         ! initial vibrational velocity
     real*8            :: par_frac     ! fraction of space occupied by particles
     real*8            :: par_range    ! Particle size range factor about mean
     real*8            :: par_margin_h ! Particle margin on sides relative to particle size
     real*8            :: alpha_bot    ! Distension at the bottom of layer (infinite depth) for linear (exponential) model
     real*8            :: alpha_top    ! Distension at the surface of layer for linear/exp model
     real*8            :: alpha_d      ! Porosity e-folding depth for exponential porosity profile        

     integer           :: mat          ! material ID
     integer           :: off(1:ndim)  ! offset from center
     integer           :: cppr(0:ndim) ! resolution (CPPR), 0=max.
     integer           :: damprof
     integer           :: par_host_obj ! host object for particles
     integer           :: par_rgb_id   ! RGB particle identifier for bitmap reading
     
     !DErelic: pos is only used for 'layers'. Merge with offset later...
     integer           :: pos
     integer           :: angle_offset_l ! offset of start of layer slope in X direction (default is end of left extension zone)
     integer           :: angle_offset_r ! offset of end of layer slope in X direction (default is start of right extension zone)
     
     !DErelic: type (of layers) is still integer. Merge with character(len=12) :: type
     !         of standard objects later...
     integer           :: laytype

  end type type_obj

  type(type_obj), target, allocatable :: obj(:)

  integer :: lay_num
  real*8  :: v_imp_max
  integer :: col_site
  integer :: symmetry   ! set eq 1 if no halfspace in y
  real*8  :: T_surf, P_surf

#ifdef ISALE2D 
  real*8  :: pvel(maxproj,1:2)                   ! projectile velocity
#endif

#ifdef ISALE3D
  real*8, allocatable :: cmc_input(:,:,:,:)  ! for use with setup_geometry_fromfile
#endif

  ! Setup types
  character(len=*), parameter :: STYPE_DEFAULT = "DEFAULT"
  character(len=*), parameter :: STYPE_PLANET = "PLANET"
  character(len=*), parameter :: STYPE_BOUNCE = "BOUNCE"
  character(len=*), parameter :: STYPE_ROTATE = "ROTATE"
  character(len=*), parameter :: STYPE_SHEARFLOW = "SHEARFLOW"
  character(len=*), parameter :: STYPE_LANDSLIDE = "LANDSLIDE"
  character(len=*), parameter :: STYPE_MESO_PORE = "MESO_PORE"
  character(len=*), parameter :: STYPE_MESO_PART = "MESO_PART"
  character(len=*), parameter :: STYPE_MESO_BITMAP = "MESO_BMP"
  character(len=*), parameter :: STYPE_DEFORM = "DEFORM"
  character(len=*), parameter :: STYPE_IMPORT_GEOM = "IMPRT_GEOM"

  ! Object types
  character(len=*), parameter :: OBJTYPE_SPHEROID = "SPHEROID"
  character(len=*), parameter :: OBJTYPE_CUBOID = "CUBOID"
  character(len=*), parameter :: OBJTYPE_BITMAP = "BITMAP"
  character(len=*), parameter :: OBJTYPE_PLATE = "PLATE"
  character(len=*), parameter :: OBJTYPE_CYLINDER = "CYLINDER"

  real, allocatable :: LAY_POS_REL(:) ! POSITION OF EACH LAYER (relative, in % of high res heigh

  ! Layer types
  integer, parameter :: LAYTYPE_STANDARD = 0
  integer, parameter :: LAYTYPE_LANDSLIDE = 10
  integer, parameter :: LAYTYPE_MESOSCALE = 20
  integer, parameter :: LAYTYPE_MESHFILE = 30
  character(len=*), parameter :: OBJTYPE_HEMISPHERE="HEMISPHERE"
  character(len=*), parameter :: OBJTYPE_CONE    ="CONE"
  character(len=*), parameter :: OBJTYPE_TRUNCONE="TRUNCONE"

  ! Temperature profiles
  character(len=*), parameter :: TPROF_CONST = "CONST"
  character(len=*), parameter :: TPROF_COND = "COND"
  character(len=*), parameter :: TPROF_CONDCONV = "CONDCONV"
  character(len=*), parameter :: TPROF_CONDCONVCAP = "CONDCONVCAP"
  character(len=*), parameter :: TPROF_USER = "USER"

  ! Porosity and maybe other profiles
  character(len=*), parameter :: PROF_NONE   = "NONE"
  character(len=*), parameter :: PROF_LINEAR = "LINEAR"
  character(len=*), parameter :: PROF_EXP    = "EXP" 

  ! Solver types
  character(len=*), parameter :: EULERIAN = "EULER"
  character(len=*), parameter :: LAGRANGIAN = "LAGRANGE"
  character(len=*), parameter :: ALE = "ALE"

  ! Tracer motion
  character(len=*), parameter :: TR_VEL = "VELOCITY" ! Interpolate tracer velocity from node velocities
  character(len=*), parameter :: TR_MAT = "MATERIAL" ! Use material fluxes to calculate tracer velocity

  ! Tracer setup types
  character(len=*), parameter :: TR_DEFAULT = "DEFAULT" ! Default tracer setup in cartesian grid
  character(len=*), parameter :: TR_SPHERICAL = "SPHERICAL" ! Tracers arranged in spherical coords


contains
! --------------------------------------------------------------------------------
!> Subroutine to read all input parameters for projectile (object) no. 'i' and
!! store that information in the corresponding 'obj' data structure.
!!
!! Called from read_asteroid
!!
!! @authors Dirk Elbeshausen, MfN
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (layed out from read_asteroid)....DE   2014-06-02
!!
!<--------------------------------------------------------------------------------
  subroutine param_read_object(i,oo)
    use ptool_interface
    implicit none
    integer,        intent(in)    :: i  !< projectile number to read
    type(type_obj), intent(inout) :: oo !< allocated (!!) object type
    ! local variables
    integer :: ierr

    if (i > obj_num) &
         call mystop("in param_read_proj: object no. > max. number of objects!")

    !---- first, initialize object with appropriate values
    oo%vel  = 0.D0
    oo%cppr = 0

    call ptool_get_int("OBJRESH",oo%cppr(1),icol=i,serr=1) ! must be contained...
#ifdef ISALE2D
    call ptool_get_int("OBJRESV",oo%cppr(2),icol=i,defval=oo%cppr(1))
#else
    call ptool_get_int("OBJRESV",oo%cppr(3),icol=i,defval=oo%cppr(1))
    call ptool_get_int("OBJRESD",oo%cppr(2),icol=i,defval=oo%cppr(1))
#endif
    call ptool_get_dble("OBJVEL",oo%vel,icol=i,serr=1)
    
    call ptool_get_dble("ANGLE", oo%angle(1),icol=i,defval=90.D0)
#ifdef ISALE3D
    call ptool_get_dble("ANGLE2",oo%angle(2),icol=i,defval=0.D0)
#endif
    
    call ptool_get_str("OBJMAT",oo%name,ierr,icol=i,serr=1)
    call ptool_get_str("OBJTYPE",oo%type,ierr,icol=i,defstr=OBJTYPE_SPHEROID)
    call transform_objtype(oo%type)

    ! If object is to be read as a bitmap, call the bitmap reader here:
    if (oo%type == OBJTYPE_BITMAP) then
       call ptool_get_str("OBJBITMAP",oo%bitmap,ierr,icol=i,serr=1)
    else
       oo%bitmap = "NONE"
    end if
    
    ! Temperature profile flag (needed in PLANET mode or if LAY_NUM=0)
    call ptool_get_str("OBJTPROF",oo%tprof,ierr,icol=i,defstr=TPROF_CONST)

    ! Porosity profile flag (needed in OBJ_GRAD \= 0). 
    call ptool_get_str("OBJPPROF", oo%pprof, ierr, icol=i, defstr=PROF_NONE)
    oo%alpha_bot = 1.D0 
    oo%alpha_top = 1.D0 
    call ptool_get_dble("OBJALP_B", oo%alpha_bot, icol=i, defval=1.D0)
    call ptool_get_dble("OBJALP_T", oo%alpha_top, icol=i, defval=1.D0)   
    !           -- This is currently supported for 2D flat layer case.
    oo%alpha_d   = 1.D0
    call ptool_get_dble("OBJALP_D", oo%alpha_d, icol=i, defval=1.D0)   

    ! Optional: set projectile internal energy (sie), if not set, sie is 
    !           set to 0. - and will later be set depending on OBJTEMP
    call ptool_get_dble("OBJENER",oo%sie,icol=i,defval=0.D0)
    if (ndim==3 .and. oo%sie.ne. 0.D0) call ptool_text(io_all,&
         "   OBJENER input parameter not yet implemented in iSALE3D","using 0.D0")
    
    ! Optional: set initial damage of the projectile
    !           -- currently supported by iSALE-2D only
     call ptool_get_dble("OBJDAM",oo%damage,icol=i,defval=0.D0)

     call ptool_get_int("OBJDAMPROF", oo%damprof, icol=i, defval=0)

     ! Optional : set initial vibrational velocity of the projectile
     !           -- used for testing of Melsoh acfl.
     call ptool_get_dble("OBJVVEL",oo%vvel,icol=i,defval=0.D0)

     ! Optional: set projectile temperature, if not set, T_surf is used
     call ptool_get_dble("OBJTEMP",oo%temp,icol=i,defval=T_surf)

     ! Optional: set projectile density, default is zero.
     call ptool_get_dble("OBJDENS",oo%density,icol=i,defval=0.D0)
     if (ndim==3 .and. oo%density .ne. 0.D0) call ptool_text(io_all,&
          "   OBJDENS input parameter not yet implemented in iSALE3D","using default value")

     ! Optional: set projectile pressure, if not set, P_surf is used
     call ptool_get_dble("OBJPRES",oo%pressure,icol=i,defval=P_surf)

     call ptool_get_int("OBJOFF_H",oo%off(X),icol=i,defval=0)
#ifdef ISALE2D
     call ptool_get_int("OBJOFF_V",oo%off(Y),icol=i,defval=0)
#else
     call ptool_get_int("OBJOFF_V",oo%off(Z),icol=i,defval=0)
     call ptool_get_int("OBJOFF_D",oo%off(Y),icol=i,defval=0)
#endif
     
   end subroutine param_read_object
  
! --------------------------------------------------------------------------------
!> Subroutine to read all input parameters for layer no. 'i' and
!! store that information in the corresponding 'obj' data structure.
!!
!! Called from read_asteroid
!!
!! @authors Dirk Elbeshausen, MfN
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (layed out from read_asteroid)....DE   2014-06-02
!!
!<--------------------------------------------------------------------------------
  subroutine param_read_lay(i,oo)
    use ptool_interface
    implicit none
    integer,        intent(in)    :: i  !< projectile number to read
    type(type_obj), intent(inout) :: oo !< allocated (!!) object type
    ! local variables
    integer :: ierr

    !---- first, initialize object with appropriate values
    oo%vel  = 0.D0
    oo%cppr = 0

    call ptool_get_int("LAYTYPE", oo%laytype,icol=i,defval=0)
    call ptool_get_int("LAYPOS",  oo%pos,icol=i,serr=1)
    call ptool_get_str("LAYMAT",  oo%name,ierr,icol=i,serr=1)
    ! optional: initial damage for layer
    !           -- not supported for iSALE-3D
    call ptool_get_dble("LAYDAM", oo%damage,icol=i,defval=0.D0)
    call ptool_get_int("LAYDAMPROF", oo%damprof, icol=i, defval=0)
    call ptool_get_str("LAYTPROF",oo%tprof,ierr,icol=i,defstr=TPROF_CONST)

    call ptool_get_dble("LAYANG", oo%angle(1),icol=i,defval=0.D0)
    call ptool_get_int("LAYOFF_L", oo%angle_offset_l, icol=i, defval=0)
    call ptool_get_int("LAYOFF_R", oo%angle_offset_r, icol=i, defval=0)

    ! Define the profile of porosity
    call ptool_get_str("LAYPPROF", oo%pprof, ierr, icol=i, defstr=PROF_NONE)
    oo%alpha_bot = 1.D0 
    oo%alpha_bot = 1.D0 
    call ptool_get_dble("LAYALP_B", oo%alpha_bot, icol=i, defval=1.D0)
    call ptool_get_dble("LAYALP_T", oo%alpha_top, icol=i, defval=1.D0)   
    !           -- This is currently supported for 2D flat layer case.
    oo%alpha_d   = 1.D0
    call ptool_get_dble("LAYALP_D", oo%alpha_d, icol=i, defval=1.D0)

   ! Define object properties as appropriate for layers
    oo%density = 0.D0; oo%sie = 0.D0; oo%pressure = P_surf; oo%temp = T_surf

  end subroutine param_read_lay



  subroutine transform_objtype(thistype)
    use ptool_interface
    implicit none
    character(len=*), intent(inout) :: thistype

    select case (trim(thistype))
    case (OBJTYPE_SPHEROID, OBJTYPE_CUBOID, OBJTYPE_BITMAP, OBJTYPE_PLATE, & 
          OBJTYPE_HEMISPHERE, OBJTYPE_CONE, OBJTYPE_TRUNCONE)
       ! nothing to do here
    case ("ELLIPSE", "SPHERE")
       thistype = OBJTYPE_SPHEROID
    case ("BOX", "RECTANGLE")
       thistype = OBJTYPE_CUBOID
    case ("CYLINDER")
       !DE in 3D, we keep this as a special case (CYLINDER)
#ifdef ISALE2D
       thistype=OBJTYPE_CUBOID
#endif
    case default
       call ptool_text(io_all, "in transform_objtype: unknown string for objtype",thistype)
       call mystop("in transform_objtype: unknown string for objtype: "//trim(thistype))
    end select

  end subroutine transform_objtype

end module mod_setup_params

! --------------------------------------------------------------------
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! -------- setup modules ---------------------------------------------
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! --------------------------------------------------------------------

module mod_param_mesh
  use mod_ndim
  implicit none
  integer :: next(1:ndim,0:2) ! number of cells per dir
  real*8  :: extfac           ! extension factor
#ifdef ISALE2D
  real*8  :: dx,dy
#else
  real*8  :: dx,dy,dz
#endif
  real*8  :: d(1:ndim)        ! grid spacing (highres-zone)
  real*8  :: d_max(0:ndim)     ! max. grid spacing (ext. grid)

contains
  
  ! --------------------------------------------------------------------------------
  !> Subroutine to read all mesh parameters from the input file and
  !! store that information in the 'next' array.
  !!
  !! Called from read_asteroid
  !!
  !! @authors Dirk Elbeshausen, MfN
  !!
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (layed out from read_asteroid)....DE   2014-06-02
  !!
  !<--------------------------------------------------------------------------------
  subroutine param_read_mesh()
    use mod_setup_params, only : symmetry
    use ptool_interface
    use mod_shared, only : NOVALI
    implicit none
    integer, parameter      :: nowarn=-1
    call ptool_get_int("GRIDH",next(X,FIRST), icol=1)
    call ptool_get_int("GRIDH",next(X,CENTER),icol=2)
    call ptool_get_int("GRIDH",next(X,LAST),  icol=3)
#ifdef ISALE2D
    call ptool_get_int("GRIDV",next(Y,FIRST), icol=1)
    call ptool_get_int("GRIDV",next(Y,CENTER),icol=2)
    call ptool_get_int("GRIDV",next(Y,LAST),  icol=3)
#else
    call ptool_get_int("GRIDV",next(Z,FIRST), icol=1)
    call ptool_get_int("GRIDV",next(Z,CENTER),icol=2)
    call ptool_get_int("GRIDV",next(Z,LAST),  icol=3)
    call ptool_get_int("SYMMETRY",symmetry, serr=nowarn, defval=2)
    if (next(Y,FIRST)==NOVALI) then ! use same as in x-dir
       next(Y,:) = next(X,:)
       if (symmetry==2 .or. symmetry==4) next(Y,FIRST) = next(Y,FIRST)/2
    else
       call ptool_get_int("GRIDD",next(Y,FIRST),  icol=1)
       call ptool_get_int("GRIDD",next(Y,CENTER), icol=2)
       call ptool_get_int("GRIDD",next(Y,LAST)  , icol=3)
    endif
#endif
    call ptool_get_dble("GRIDEXT",extfac,defval=1.04D0)
    call ptool_get_dble("GRIDSPC",dx,icol=1,serr=1)
    call ptool_get_dble("GRIDSPC",dy,icol=2,defval=dx)
#ifdef ISALE3D
    call ptool_get_dble("GRIDSPC",dz,icol=3,defval=dx)
#endif
#ifdef ISALE2D
    call ptool_get_dble("CYL",cyl,defval=1.D0)
#endif
    call ptool_get_dble("GRIDSPCM",d_max(1),icol=1,defval=-20.D0)
    call ptool_get_dble("GRIDSPCM",d_max(2),icol=2,defval=d_max(1))
#ifdef ISALE3D
    call ptool_get_dble("GRIDSPCM",d_max(3),icol=3,defval=d_max(1))
#endif
  end subroutine param_read_mesh

end module mod_param_mesh

! --------------------------------------------------------------------

module mod_regrid
  implicit none
  integer :: regridmesh, regrid_num
  real*8  :: dxnew,dynew
#ifdef ISALE2D  
  integer :: nchl,nchr,nchb,ncht            ! Cells to add/minus from new grid
  integer :: extaddl,extaddr,extaddb,extaddt! Cells for new mesh extension
  integer :: nextaddl,nextaddr               ! Cells for new mesh extension (dynamic regrid)
  integer :: nextaddb,nextaddt               ! Cells for new mesh extension (dynamic regrid)
  integer :: nnchl(0:4), nnchr(0:4)         ! Cells to add/minus from new grid dynamically  
  integer :: nnchb(0:4), nncht(0:4)         ! Cells to add/minus from new grid dynamically   
  
  integer :: dynamic_regrid
  real*8, allocatable :: Cp(:,:,:)        ! COORDS IN X,Y
  real*8, allocatable :: Vp(:,:,:)        ! VELOCITY COMPONENTS
  real*8, allocatable :: ROp(:,:,:)       ! DENSITY
  real*8, allocatable :: MCp(:,:)         ! MASS IN CELL
  real*8, allocatable :: CMCp(:,:,:)      ! MAT. CONCENTRATION IN CELL
  real*8, allocatable :: VOIDCONCp(:,:)   ! CONCENTRATION OF VOID
  real*8, allocatable :: SIEp(:,:,:)      ! SPECIFIC INTERNAL ENERGY
  real*8, allocatable :: VOLp(:,:)        ! CELL VOLUME
  real*8, allocatable :: ALPHAp(:,:,:)    ! DISTENSION 
  real*8, allocatable :: STRESDEVp(:,:,:) ! COMP. OF STRESS DEV. 
  real*8, allocatable :: EPSTRAINp(:,:)
  real*8, allocatable :: VOLSTRAINp(:,:)
  real*8, allocatable :: PLWp(:,:)
  real*8, allocatable :: VELOp(:,:)
  real*8, allocatable :: DAMAGEp(:,:)
  real*8, allocatable :: DUMMYFp(:,:,:)   ! DUMMY FIELD #1-9
  real*8, allocatable :: INITCOORDp(:,:,:)! DUMMY FIELD #2
  real*8, allocatable :: GRAVITYp(:,:,:)  ! Gravity
  real*8              :: regrid_thresh
  integer             :: dynamic_buffer_v, ii, counter, jj 
  integer             :: dynamic_position, dynamic_coarsen
  character(len=20)   :: F_NAME, regrid_trig
  integer :: jerr
  integer :: dynamic_time(0:4), dynamic_dtsave    ! *savetime for dynamic time
#else 
  real*8, allocatable :: C_tmp(:,:,:,:)        ! COORDS IN X,Y
  real*8, allocatable :: V_tmp(:,:,:,:)        ! VELOCITY COMPONENTS
  real*8, allocatable :: RO_tmp(:,:,:,:)       ! DENSITY
  real*8, allocatable :: MC_tmp(:,:,:)         ! MASS IN CELL
  real*8, allocatable :: CMC_tmp(:,:,:,:)      ! MAT. CONCENTRATION IN CELL
  real*8, allocatable :: VOIDCONC_tmp(:,:,:)   ! CONCENTRATION OF VOID
  real*8, allocatable :: SIE_tmp(:,:,:,:)      ! SPECIFIC INTERNAL ENERGY
  real*8, allocatable :: VOL_tmp(:,:,:)        ! CELL VOLUME
  real*8, allocatable :: ALPHA_tmp(:,:,:,:)    ! DISTENSION 
  real*8, allocatable :: STRESDEV_tmp(:,:,:,:) ! COMP. OF STRESS DEV. 
  real*8, allocatable :: EPSTRAIN_tmp(:,:,:)
  real*8, allocatable :: VOLSTRAIN_tmp(:,:,:)
  real*8, allocatable :: PLW_tmp(:,:,:)
  real*8, allocatable :: VELO_tmp(:,:,:)
  real*8, allocatable :: DAMAGE_tmp(:,:,:)
  real*8, allocatable :: DUMMYF_tmp(:,:,:,:)   ! DUMMY FIELD #1-9
  real*8, allocatable :: INITCOORD_tmp(:,:,:,:)! DUMMY FIELD #2
  real*8, allocatable :: GRAVITY_tmp(:,:,:,:)  ! Gravity
  
  real*8, allocatable :: MV_tmp(:,:,:)
  real*8, allocatable :: RMV_tmp(:,:,:)
  real*8, allocatable :: ROL_tmp(:,:,:,:)
  real*8, allocatable :: P_tmp(:,:,:)
  real*8, allocatable :: T_tmp(:,:,:)
  real*8, allocatable :: Q_tmp(:,:,:)
  real*8, allocatable :: TMELT_tmp(:,:,:)
  real*8, allocatable :: CSOUND_tmp(:,:,:)
  real*8, allocatable :: ALPHAP_tmp(:,:,:,:)
  real*8, allocatable :: VOLSTRP_tmp(:,:,:)
  real*8, allocatable :: YAC_tmp(:,:,:)
  real*8, allocatable :: PVIBR_tmp(:,:,:)
  real*8, allocatable :: YIELD_VALS_tmp(:,:,:)
  real*8, allocatable :: e_rate_tmp(:,:,:)
  real*8, allocatable :: mu_tmp(:,:,:)
  integer*1, allocatable :: failure_tmp(:,:,:)
  real*8, allocatable :: grav_magnitude_tmp(:,:,:)
  real*8, allocatable :: MAXP_tmp(:,:,:)
#endif
contains
  
  subroutine read_param_regrid()
    use ptool_interface
    use mod_ndim, only : X,Y,FIRST,LAST
    use mod_param_mesh, only : next
    implicit none

    call ptool_text(io_oin,"   obtaining regridding parameters","",2)
    call ptool_get_int("REGRID",regridmesh,defval=0)
    call ptool_text_int(io_ini,"      doing regrid",regridmesh,3)
    
#ifdef ISALE2D
    ! SR: check if dynamically regridding  
    call ptool_get_int("RGRD_DYNM",dynamic_regrid,defval=0) 

    ! If regridding grab parameters defining cells to add
    ! or remove from high-resolution zone
    ! NB: extension cells will be taken from next
    if ((regridmesh .eq. 1) .and. (dynamic_regrid .eq. 0)) then
       call ptool_get_int("NEWCELLH",nchl,icol=1,defval=0)
       call ptool_get_int("NEWCELLH",nchr,icol=2,defval=0)
       call ptool_get_int("NEWCELLV",nchb,icol=1,defval=0)
       call ptool_get_int("NEWCELLV",ncht,icol=2,defval=0)
       ! Define new extension-zone cells for mesh
       extaddl = next(X,FIRST)
       extaddr = next(X,LAST)
       extaddb = next(Y,FIRST)
       extaddt = next(Y,LAST)
    end if
#endif

  end subroutine read_param_regrid

#ifdef ISALE2D  
  subroutine read_param_dynamic_regrid()
    use mod_ndim, only : X,Y,FIRST,LAST
    use mod_param_mesh, only : next
    use ptool_interface
    implicit none

    call ptool_text(io_oin,"   obtaining dynamic regridding parameters","",2)
    call ptool_get_int("RGRD_DYNM",dynamic_regrid,defval=0)

    ! SR: If dynamically regridding (dynamic_regrid = 1 or 2) grab
    ! parameters definig the number of regrids requested, parameters 
    ! defining how many cells above bottom to buffer from and the 
    ! value for the dynamic coarsening
    if ((dynamic_regrid .eq. 1) .or. (dynamic_regrid .eq. 2)) then
      call ptool_get_int("RGRD_NUM",regrid_num,defval=1)
      call ptool_get_int("RGRD_CRS",dynamic_coarsen,defval=0)
      call ptool_get_int("RGRD_BUFF",dynamic_buffer_v,defval=20)
      call ptool_get_int("RGRD_DYNT",dynamic_dtsave,defval=1)
      
      ! SR: If low number of regrids selected by the user (<5) 
      ! read parameters for every regrid 
      do ii =0, 4
        call ptool_get_int("NEWCELLH",nnchl(ii),icol=2*ii+1,defval=0)
        call ptool_get_int("NEWCELLH",nnchr(ii),icol=2*ii+2,defval=0)
        call ptool_get_int("NEWCELLV",nnchb(ii),icol=2*ii+1,defval=0)
        call ptool_get_int("NEWCELLV",nncht(ii),icol=2*ii+2,defval=0)
      end do
      
      ! Define new extension-zone cells for mesh
      call ptool_get_int("EXTCELLH",nextaddl,icol=1,defval=0)
      call ptool_get_int("EXTCELLH",nextaddr,icol=2,defval=0)
      call ptool_get_int("EXTCELLV",nextaddb,icol=1,defval=0)
      call ptool_get_int("EXTCELLV",nextaddt,icol=2,defval=0)
      
    end if 
    
    
    ! SR: If dinamically regridding using spatial threshold, grab
    ! parameters defining: how many cells above the bottom to buffer 
    ! from, use density or pressure as threshold and the threshold value
    if (dynamic_regrid.eq.1) then
       call ptool_get_int("RGRD_POS",dynamic_position,defval=-55555)
       call ptool_get_str("RGRD_TRIG",regrid_trig,jerr, defstr='DENSITY') 
       call ptool_get_dble("RGRD_THR",regrid_thresh,defval=1.10D0)
       F_NAME = 'init'
       
    ! SR: If dynamically regridding using temporal threshold, 
    ! grab parameters defining the regrid time
    elseif (dynamic_regrid.eq.2) then
      do ii =0, 4
        call ptool_get_int("RGRD_TIME",dynamic_time(ii), icol=ii+1,defval=20)
      end do
      
       F_NAME = 'init'
       
    end if

  end subroutine read_param_dynamic_regrid
#endif

end module mod_regrid

! --------------------------------------------------------------------

#ifdef ISALE2D
module mod_selfgrav
  use mod_shared, only : PI
  implicit none
  !for case to select correct brute force sum
  integer, parameter    :: SAME_X=1
  integer, parameter    :: SAME_Y=2
  integer, parameter    :: SAME_XANDY=3
  integer, parameter    :: MASS_ON_Y_AXIS=10

  ! params for determining when to calculate self gravity
  !  integer, parameter    :: SG_1ST=5        ! first cycle to calculate self-grav
  !  integer, parameter    :: SG_NEXT=50      ! calc sg every this many cycles thereafter
  !  integer, parameter    :: SG_0CMV=150     ! calc sg for all points in the mesh, not just those
  !                                           !   with mass.  Must be a multiple of 

  ! Object type for nested meshes used in self gravity calculation
  type self_grav_mesh

     real*8, allocatable :: mass(:,:)        ! Mass of each vertex in the mesh
     real*8, allocatable :: position(:,:,:)  ! Position of vertex in the mesh
     integer, allocatable :: count(:,:)      ! Count of finest-mesh vertices in mesh vertex
     integer :: resolution                   ! Size of mesh cell in finest-mesh cells
     integer :: numnodes(2)                  ! Number of nodes in mesh in each direction

  end type self_grav_mesh

  type(self_grav_mesh), target, allocatable :: SelfGravMesh(:)

  integer :: nummeshes
  integer :: selfgrav_cycle    ! User option: number of cycles per update
  integer :: finest_mesh       ! User option: mesh level for finest mesh (default = 1)

  real*8 :: barhut_param ! User option for choosing the accuracy of the self-grav alg.

  ! brute force sum parameters
  real*8, parameter     :: theta_sum_limit=PI  ! extent of azimuthal sum
  integer, parameter    :: n_theta=36 ! number of slices to divide the 2pi azimuthal ring into
  real*8, parameter     :: n_theta_div2=DBLE(n_theta)/2.0D0
  real*8, parameter     :: increment_theta=2.0D0*PI/DBLE(n_theta)

  integer, parameter    :: theta_loop_limit=n_theta/2
  real*8, save          :: cosine_arr(theta_loop_limit)
end module mod_selfgrav

#endif

! --------------------------------------------------------------------

module mod_param_time
implicit none
  real*8  :: dt          ! current time increment (calculated in timstp)
  real*8  :: time        ! current model time
  real*8  :: timescale   ! scaling factor to derive nondimensional time
  integer :: ncyc,ncycold
end module mod_param_time

module mod_particle

  use mod_ndim
  implicit none
  type type_particle
     integer :: iobj
     integer :: class
     real*8  :: size
     real*8  :: radius(1:ndim)
     real*8  :: midpoint(1:ndim)
  end type type_particle
  type(type_particle), target, allocatable :: particle(:)

end module mod_particle

module mod_probe
  implicit none
  integer, parameter :: nprobemax = 100
  integer :: nprobe
  character(len=20) :: probe_type
  integer :: probe(1:2,1:nprobemax)
  real*8 :: probe_position(1:2,1:nprobemax)

  ! Probe arrays
  integer :: nparray,maxnprobes
  integer, parameter :: nparraymax = 10
  integer :: parray(1:4,1:nparraymax)
  integer, allocatable :: parrayc(:,:,:)

contains 

#ifdef ISALE2D  
  integer function find_tracer_id(pos)

    use mod_tracer, only : tracer_num, tracer, TR_X, TR_Y
    use mod_ndim, only : X,Y
    implicit none
    real*8, intent(in) :: pos(2)
    real*8 :: distance,distance_min
    integer :: it, ii

    ! Find the tracer index of the tracer closest to pos
    distance_min = 1.D10
    do ii = 1,tracer_num
       distance = (tracer(TR_X,ii)-pos(X))**2 + (tracer(TR_Y,ii)-pos(Y))**2
       if (distance .lt. distance_min) then
          distance_min = distance
          find_tracer_id = ii
       end if
    end do

  end function find_tracer_id

  function find_cell_ids(pos,located)
    ! ===============================================================
    ! Function to find the i-j cell location of a physical position
    ! (pos). If the position lies outside the grid (left or bottom),
    ! the function returns the bottom-right cell coordinates. . .
    ! ===============================================================
    use mod_tracer, only : tracer_num, tracer, TR_X, TR_Y
    use mod_ndim, only   : X,Y,n
    use mod_fields, only : C
    implicit none
    real*8, intent(in)   :: pos(2)
    logical, intent(out) :: located
    real*8 :: distance,distance_min
    integer :: i,j
    integer, dimension(2) :: find_cell_ids

    ! If position is outside grid, return bottom right cell id
    if (pos(X) > C(X,n(X)+1,1) .or. pos(Y) < C(Y,1,1)) then
       find_cell_ids(X) = n(X)
       find_cell_ids(Y) = 1
       located = .false.
       return
    end if

    ! Find the tracer index of the tracer closest to pos
    located = .true.
    distance_min = 1.D10
    find_cell_ids(X) = n(X); find_cell_ids(Y) = n(Y)
    do i = 1,n(X)
       distance = (0.5D0*(C(X,i,1)+C(X,i+1,1))-pos(X))**2
       if (distance .lt. distance_min) then
          distance_min = distance
          find_cell_ids(X) = i
       end if
    end do
    distance_min = 1.D10
    do j = 1,n(Y)
       distance = (0.5D0*(C(Y,1,j)+C(Y,1,j+1))-pos(Y))**2
       if (distance .lt. distance_min) then
          distance_min = distance
          find_cell_ids(Y) = j
       end if
    end do

  end function find_cell_ids

#endif
  
end module mod_probe



