module CropType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing variables needed for the crop model
  !
  ! TODO(wjs, 2014-08-05) Move more crop-specific variables into here - many are
  ! currently in CNVegStateType
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use spmdMod             , only : masterproc
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use clm_varcon          , only : spval
  use clm_varctl          , only : iulog, use_crop
  use pftconMod           , only : mxnp
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC DATA TYPES:
  !
  ! Crop state variables structure
  type, public :: crop_type

     ! Note that cropplant and harvdate could be 2D to facilitate rotation
     integer , pointer :: nyrs_crop_active_patch  (:)   ! number of years this crop patch has been active (0 for non-crop patches)
     logical , pointer :: croplive_patch          (:)   ! patch Flag, true if planted, not harvested
     logical , pointer :: cropplant_patch         (:)   ! patch Flag, true if planted
     integer , pointer :: harvdate_patch          (:)   ! patch harvest date
     real(r8), pointer :: fertnitro_patch         (:)   ! patch fertilizer nitrogen
     real(r8), pointer :: gddplant_patch          (:)   ! patch accum gdd past planting date for crop       (ddays)
     real(r8), pointer :: gddtsoi_patch           (:)   ! patch growing degree-days from planting (top two soil layers) (ddays)
     real(r8), pointer :: vf_patch                (:)   ! patch vernalization factor for cereal
     real(r8), pointer :: cphase_patch            (:)   ! phenology phase
     real(r8), pointer :: latbaset_patch          (:)   ! Latitude vary baset for gddplant (degree C)
     character(len=20) :: baset_mapping
     real(r8) :: baset_latvary_intercept
     real(r8) :: baset_latvary_slope

     !oil palm related new state variables  (Y.Fan)
     !integer  , pointer :: yrop_patch(:)                       ! year of planting
     !integer  , pointer :: idpp_patch(:)                       ! days past planting
     !integer  , pointer :: idpp2_patch(:)                      ! saved from phase 2
     !integer  , pointer :: idop_patch(:)                       ! date of planting (defined in cnveg_state_type)
     integer  , pointer :: np_patch(:)                         ! total number of phytomers having appeared so far
     integer  , pointer :: rankp_patch(:,:)                    ! rank of phytomers from 1=youngest to np=oldest and 0=dead
     integer  , pointer :: plaipeak_patch(:,:)                 ! Flag, 1: max allowed lai per phytomer; 0: not at max
     real(r8) , pointer :: livep_patch(:,:)                    ! Flag, true (=1) if this phytomer is alive
     integer  , pointer :: lfoutday_patch(:,:)                 ! date of leaf/phytomer emergence
     integer  , pointer :: lfoutyr_patch(:,:)                  ! year of leaf/phytomer emergence
     integer  , pointer :: lfdays_patch(:,:)                   ! days past leaf emergence for each phytomer
     real(r8) , pointer :: huileafnp_patch(:,:)                ! hui needed for initiation of successive phytomers
     real(r8) , pointer :: huilfexpnp_patch(:,:)               ! hui needed for leaf expansion of successive phytomers
     real(r8) , pointer :: huigrnnp_patch(:,:)                 ! hui needed for starting grainfill of successive phytomers
     real(r8) , pointer :: huilfmatnp_patch(:,:)               ! hui needed for leaf maturity of successive phytomers
     real(r8) , pointer :: grnmatnp_patch(:,:)                 ! hui needed for grain maturity of successive phytomers
     real(r8) , pointer :: huilfsennp_patch(:,:)               ! hui needed for leaf senescence of successive phytomers
     real(r8) , pointer :: huilfendnp_patch(:,:)               ! hui needed for end of life of successive phytomers
     real(r8) , pointer :: psla_patch(:,:)                      ! specific leaf area (m^2/gC) of each phytomer
     real(r8) , pointer :: plai_patch(:,:)                     ! one-sided leaf area index of each phytomer
     !real(r8) , pointer :: phyllochron2_patch(:)             ! extended phyllochron through maturity
     !real(r8) , pointer :: avgvegsink_patch(:)               ! average vegetative sink size (running mean)
     !real(r8) , pointer :: avgrepsink_patch(:)               ! average reproductive sink size (running mean)
     !real(r8) , pointer :: crop_patch(:)                       !binary flag: crop or not
     !integer ,  pointer :: perennial_patch(:)                  !binary flag for perennial crop phenology (1=perennial, 0=not perennial)
     !integer , pointer :: phytomer_patch(:)                    !total number of phytomers in life time, if >0 use phytomer structure
     !real(r8), pointer :: harvest_flag_patch(:)                !harvest flag for perennial crops 
     real(r8), pointer :: harvest_counter_patch(:)             !harvest counter to tag the phytomer 
     logical , pointer :: prune_patch(:)                       !Flag, if true do pruning
     real(r8), pointer :: bglfr_p_patch(:,:)                   !background litterfall rate for each phytomer 
     real(r8), pointer :: bgtr_p_patch(:,:)                    !background transfer growth rate for each phytomer     
     !real(r8), pointer :: totrepc_patch(:)                     ! (gC/m2) total reproductive carbon, excluding cpool and harvested foodc (added by Y.Fan)
     !real(r8), pointer :: vegsize0_patch(:)                    ! (gC/m2) vegetation size (vegetative C) at the beginning of yield
     

   contains
     ! Public routines
     procedure, public  :: Init               ! Initialize the crop type
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: Restart
     procedure, public  :: ReadNML            ! Read in the crop namelist

     ! NOTE(wjs, 2014-09-29) need to rename this from UpdateAccVars to CropUpdateAccVars
     ! to prevent cryptic error messages with pgi (v. 13.9 on yellowstone)
     ! This is probably related to this bug
     ! <http://www.pgroup.com/userforum/viewtopic.php?t=4285>, which was fixed in pgi 14.7.
     procedure, public  :: CropUpdateAccVars

     procedure, public  :: CropIncrementYear

     ! Private routines
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private, nopass :: checkDates

  end type crop_type

  character(len=*), parameter, private :: baset_map_constant = 'constant'
  character(len=*), parameter, private :: baset_map_latvary  = 'varytropicsbylat'
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------
    
    call this%InitAllocate(bounds)

    if (use_crop) then
       call this%InitHistory(bounds)
       call this%InitCold(bounds)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNML(this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CropType
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'Crop::ReadNML'
    character(len=*), parameter :: nmlname = 'crop'
    !-----------------------------------------------------------------------
    character(len=20) :: baset_mapping
    real(r8) :: baset_latvary_intercept
    real(r8) :: baset_latvary_slope
    namelist /crop/ baset_mapping, baset_latvary_intercept, baset_latvary_slope

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    baset_mapping           = 'constant'
    baset_latvary_intercept = 12._r8
    baset_latvary_slope     = 0.4_r8
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=crop, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (baset_mapping           , mpicom)
    call shr_mpi_bcast (baset_latvary_intercept , mpicom)
    call shr_mpi_bcast (baset_latvary_slope     , mpicom)

    this%baset_mapping           = baset_mapping
    this%baset_latvary_intercept = baset_latvary_intercept
    this%baset_latvary_slope     = baset_latvary_slope
    if (      trim(this%baset_mapping) == baset_map_constant ) then
       if ( masterproc ) write(iulog,*) 'baset mapping for ALL crops are constant'
    else if ( trim(this%baset_mapping) == baset_map_latvary ) then
       if ( masterproc ) write(iulog,*) 'baset mapping for crops vary with latitude'
    else
       call endrun(msg="Bad value for baset_mapping in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=crop)
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------
    
  end subroutine ReadNML

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%nyrs_crop_active_patch(begp:endp)) ; this%nyrs_crop_active_patch(:) = 0
    allocate(this%croplive_patch (begp:endp)) ; this%croplive_patch (:) = .false.
    allocate(this%cropplant_patch(begp:endp)) ; this%cropplant_patch(:) = .false.
    allocate(this%harvdate_patch (begp:endp)) ; this%harvdate_patch (:) = huge(1) 
    allocate(this%fertnitro_patch (begp:endp)) ; this%fertnitro_patch (:) = spval
    allocate(this%gddplant_patch (begp:endp)) ; this%gddplant_patch (:) = spval
    allocate(this%gddtsoi_patch  (begp:endp)) ; this%gddtsoi_patch  (:) = spval
    allocate(this%vf_patch       (begp:endp)) ; this%vf_patch       (:) = 0.0_r8
    allocate(this%cphase_patch   (begp:endp)) ; this%cphase_patch   (:) = 0.0_r8
    allocate(this%latbaset_patch (begp:endp)) ; this%latbaset_patch (:) = spval

    !for phytomer structure (Y.Fan)
    !allocate(this%yrop_patch          (begp:endp))                   ; this%yrop_patch          (:)   = huge(1)
    !allocate(this%idpp_patch          (begp:endp))                   ; this%idpp_patch          (:)   = 0
    !allocate(this%idpp2_patch         (begp:endp))                   ; this%idpp2_patch         (:)   = 0
    !allocate(this%idop_patch          (begp:endp))                   ; this%idop_patch          (:)   = huge(1)
    allocate(this%lfoutday_patch      (begp:endp,1:mxnp))            ; this%lfoutday_patch      (:,:)   = huge(1)
    allocate(this%lfoutyr_patch       (begp:endp,1:mxnp))            ; this%lfoutyr_patch       (:,:)   = huge(1)
    allocate(this%lfdays_patch        (begp:endp,1:mxnp))            ; this%lfdays_patch        (:,:)   = 0
    allocate(this%np_patch            (begp:endp))                   ; this%np_patch            (:)   = 0
    allocate(this%rankp_patch         (begp:endp,1:mxnp))            ; this%rankp_patch         (:,:)   = 0
    allocate(this%plaipeak_patch      (begp:endp,1:mxnp))            ; this%plaipeak_patch      (:,:)   = 0
    allocate(this%livep_patch         (begp:endp,1:mxnp))            ; this%livep_patch         (:,:)   = 0._r8
    allocate(this%plai_patch          (begp:endp,1:mxnp))            ; this%plai_patch          (:,:)   = 0._r8
    allocate(this%psla_patch          (begp:endp,1:mxnp))            ; this%psla_patch           (:,:)   = 0._r8 
    allocate(this%huileafnp_patch     (begp:endp,1:mxnp))            ; this%huileafnp_patch     (:,:)   = spval 
    allocate(this%huilfexpnp_patch    (begp:endp,1:mxnp))            ; this%huilfexpnp_patch    (:,:)   = spval 
    allocate(this%huilfmatnp_patch    (begp:endp,1:mxnp))            ; this%huilfmatnp_patch    (:,:)   = spval 
    allocate(this%huilfsennp_patch    (begp:endp,1:mxnp))            ; this%huilfsennp_patch    (:,:)   = spval 
    allocate(this%huilfendnp_patch    (begp:endp,1:mxnp))            ; this%huilfendnp_patch    (:,:)   = spval 
    allocate(this%huigrnnp_patch      (begp:endp,1:mxnp))            ; this%huigrnnp_patch      (:,:)   = spval
    allocate(this%grnmatnp_patch      (begp:endp,1:mxnp))            ; this%grnmatnp_patch      (:,:)   = spval
    !allocate(this%phyllochron2_patch    (begp:endp))                 ; this%phyllochron2_patch  (:)   = spval 
    !allocate(this%avgvegsink_patch      (begp:endp))                 ; this%avgvegsink_patch    (:)   = 0._r8 
    !allocate(this%avgrepsink_patch      (begp:endp))                 ; this%avgrepsink_patch    (:)   = 0._r8
    !allocate(this%perennial_patch       (begp:endp))                 ; this%perennial_patch     (:)   = 0
    !allocate(this%phytomer_patch        (begp:endp))                 ; this%phytomer_patch      (:)   = 0
    !allocate(this%harvest_flag_patch    (begp:endp))                 ; this%harvest_flag_patch  (:)   = 0._r8
    allocate(this%harvest_counter_patch (begp:endp))                 ; this%harvest_counter_patch (:) = 0._r8
    allocate(this%prune_patch           (begp:endp))                 ; this%prune_patch         (:)   = .false.
    allocate(this%bglfr_p_patch         (begp:endp,1:mxnp))          ; this%bglfr_p_patch       (:,:) = 0._r8
    allocate(this%bgtr_p_patch          (begp:endp,1:mxnp))          ; this%bgtr_p_patch        (:,:) = 0._r8    
    !allocate(this%totrepc_patch        (begp:endp))                  ; this%totrepc_patch         (:)   = spval
    !allocate(this%vegsize0_patch       (begp:endp))                  ; this%vegsize0_patch        (:)   = spval   

   
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    this%fertnitro_patch(begp:endp) = spval
    call hist_addfld1d (fname='FERTNITRO', units='gN/m2/yr', &
         avgflag='A', long_name='Nitrogen fertilizer for each crop', &
         ptr_patch=this%fertnitro_patch, default='inactive')

    this%gddplant_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDPLANT', units='ddays', &
         avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
         ptr_patch=this%gddplant_patch, default='inactive')

    this%gddtsoi_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDTSOI', units='ddays', &
         avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
         ptr_patch=this%gddtsoi_patch, default='inactive')

    this%cphase_patch(begp:endp) = spval
    call hist_addfld1d (fname='CPHASE', units='0-not planted, 1-planted, 2-leaf emerge, 3-grain fill, 4-harvest', &
         avgflag='A', long_name='crop phenology phase', &
         ptr_patch=this%cphase_patch, default='active')

    if ( (trim(this%baset_mapping) == baset_map_latvary) )then
       this%latbaset_patch(begp:endp) = spval
       call hist_addfld1d (fname='LATBASET', units='degree C', &
            avgflag='A', long_name='latitude vary base temperature for gddplant', &
            ptr_patch=this%latbaset_patch, default='inactive')
    end if

    this%livep_patch(begp:endp,1:mxnp) = spval
    call hist_addfld2d (fname='LIVEP', units='none', type2d='phytomer', &
               avgflag='A', long_name='phytomer alive or not', &
               ptr_patch=this%livep_patch, default='inactive')

    this%plai_patch(begp:endp,1:mxnp) = spval
    call hist_addfld2d (fname='PLAI', units='none', type2d='phytomer', &
               avgflag='A', long_name='phytomer LAI', &
               ptr_patch=this%plai_patch, default='inactive')
    
    !this%harvest_flag_patch(begp:endp) = spval
    !call hist_addfld1d (fname='HARVEST_FLAG', units='none', &
            !avgflag='A', long_name='harvest flag', &
            !ptr_patch=this%harvest_flag_patch, default='inactive')

    this%harvest_counter_patch(begp:endp) = spval
    call hist_addfld1d (fname='HARVEST_COUNTER', units='none', &
            avgflag='A', long_name='harvest counter', &
            ptr_patch=this%harvest_counter_patch, default='inactive')

  end subroutine InitHistory

  subroutine InitCold(this, bounds)
    ! !USES:
    use LandunitType, only : lun                
    use landunit_varcon, only : istcrop
    use PatchType, only : patch
    use clm_instur, only : fert_cft
    use pftconMod        , only : pftcon 
    use GridcellType     , only : grc
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g, p, m, ivt ! indices

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

!DLL - added wheat & sugarcane restrictions to base T vary by lat
    do p= bounds%begp,bounds%endp
       g   = patch%gridcell(p)
       ivt = patch%itype(p)

       this%nyrs_crop_active_patch(p) = 0

       if ( grc%latdeg(g) >= 0.0_r8 .and. grc%latdeg(g) <= 30.0_r8) then
          this%latbaset_patch(p)=pftcon%baset(ivt)+12._r8-0.4_r8*grc%latdeg(g)
       else if (grc%latdeg(g) < 0.0_r8 .and. grc%latdeg(g) >= -30.0_r8) then
          this%latbaset_patch(p)=pftcon%baset(ivt)+12._r8+0.4_r8*grc%latdeg(g)
       else
          this%latbaset_patch(p)=pftcon%baset(ivt)
       end if
       if ( trim(this%baset_mapping) == baset_map_constant ) then
          this%latbaset_patch(p) = nan
       end if
    end do
!DLL -- end of mods

    if (use_crop) then
       do p= bounds%begp,bounds%endp
          g = patch%gridcell(p)
          l = patch%landunit(p)
          c = patch%column(p)

          if (lun%itype(l) == istcrop) then
             m = patch%itype(p)
             this%fertnitro_patch(p) = fert_cft(g,m)
          end if
       end do
    end if

  end subroutine InitCold

  !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation. 
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field]. 
    ! To accumulate a field, it must first be defined in subroutine [initAccVars] 
    ! and then accumulated by calls to [updateAccVars].
    !
    ! Should only be called if use_crop is true
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(in) :: this
    type(bounds_type), intent(in) :: bounds  

    !
    ! !LOCAL VARIABLES:
    integer, parameter :: not_used = huge(1)

    !---------------------------------------------------------------------

    call init_accum_field (name='GDDPLANT', units='K', &
         desc='growing degree-days from planting', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='GDDTSOI', units='K', &
         desc='growing degree-days from planting (top two soil layers)', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES:
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    
    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg=" allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif

    nstep = get_nstep()

    call extract_accum_field ('GDDPLANT', rbufslp, nstep) 
    this%gddplant_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('GDDTSOI', rbufslp, nstep) 
    this%gddtsoi_patch(begp:endp)  = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    use PatchType, only : patch
    use pftconMod, only : npcropmin, npcropmax
    !
    ! !ARGUMENTS:
    class(crop_type), intent(inout)  :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer, pointer :: temp1d(:) ! temporary
    integer :: restyear
    integer :: p
    logical :: readvar   ! determine if variable is on initial file

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    if (use_crop) then
       call restartvar(ncid=ncid, flag=flag, varname='nyrs_crop_active', xtype=ncd_int, &
            dim1name='pft', &
            long_name='Number of years this crop patch has been active (0 for non-crop patches)', &
            units='years', &
            interpinic_flag='interp', readvar=readvar, data=this%nyrs_crop_active_patch)
       if (flag == 'read' .and. .not. readvar) then
          ! BACKWARDS_COMPATIBILITY(wjs, 2017-02-17) Old restart files did not have this
          ! patch-level variable. Instead, they had a single scalar tracking the number
          ! of years the crop model ran. Copy this scalar onto all *active* crop patches.

          ! Some arguments in the following restartvar call are irrelevant, because we
          ! only call this for 'read'. I'm simply maintaining the old restartvar call.
          call restartvar(ncid=ncid, flag=flag,  varname='restyear', xtype=ncd_int,  &
               long_name='Number of years prognostic crop ran', units="years", &
               interpinic_flag='copy', readvar=readvar, data=restyear)
          if (readvar) then
             do p = bounds%begp, bounds%endp
                if (patch%itype(p) >= npcropmin .and. patch%itype(p) <= npcropmax .and. &
                     patch%active(p)) then
                   this%nyrs_crop_active_patch(p) = restyear
                end if
             end do
          end if
       end if

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then 
          do p= bounds%begp,bounds%endp
             if (this%croplive_patch(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='croplive', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that crop is alive, but not harvested', &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then 
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%croplive_patch(p) = .true.
             else
                this%croplive_patch(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then 
          do p= bounds%begp,bounds%endp
             if (this%cropplant_patch(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='cropplant', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that crop is planted, but not harvested' , &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then 
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%cropplant_patch(p) = .true.
             else
                this%cropplant_patch(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       call restartvar(ncid=ncid, flag=flag,  varname='harvdate', xtype=ncd_int,  &
            dim1name='pft', long_name='harvest date', units='jday', nvalid_range=(/1,366/), & 
            interpinic_flag='interp', readvar=readvar, data=this%harvdate_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='vf', xtype=ncd_double,  &
            dim1name='pft', long_name='vernalization factor', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%vf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cphase',xtype=ncd_double, &
            dim1name='pft', long_name='crop phenology phase', &
            units='0-not planted, 1-planted, 2-leaf emerge, 3-grain fill, 4-harvest', &
            interpinic_flag='interp', readvar=readvar, data=this%cphase_patch)

       !only create restart variables for oil palm when it is active
       if ( mxnp .gt. 0 ) then

         call restartvar(ncid=ncid, varname='np', xtype=ncd_int,  &
              dim1name='pft',long_name='Total number of phytomers having appeared so far', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%np_patch)
  
         call restartvar(ncid=ncid, varname='rankp', xtype=ncd_int,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='Rank of phytomers from 1=youngest to np=oldest and 0=dead', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%rankp_patch)

         call restartvar(ncid=ncid, varname='psla', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', &
              long_name='phytomer SLA', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%psla_patch)

         call restartvar(ncid=ncid, varname='plai', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', &
              long_name='phytomer LAI', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%plai_patch)
 
         call restartvar(ncid=ncid, varname='plaipeak', xtype=ncd_int,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='Flag, 1: max allowed lai per phytomer; 0: not at max', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%plaipeak_patch)
  
         call restartvar(ncid=ncid, varname='livep', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', long_name='Flag, true if this phytomer is alive', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%livep_patch)
  
         call restartvar(ncid=ncid, varname='lfoutday', xtype=ncd_int,  &
              dim1name='pft',dim2name='phytomer', long_name='Date of leaf/phytomer emergence', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%lfoutday_patch)
  
         call restartvar(ncid=ncid, varname='lfoutyr', xtype=ncd_int,  &
              dim1name='pft',dim2name='phytomer', long_name='Year of leaf/phytomer emergence', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%lfoutyr_patch)
        
         call restartvar(ncid=ncid, varname='lfdays', xtype=ncd_int,  &
              dim1name='pft',dim2name='phytomer', long_name='Days past leaf emergence for each phytomer', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%lfdays_patch)
  
         call restartvar(ncid=ncid, varname='huileafnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer',long_name='hui needed for initiation of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huileafnp_patch)
  
         call restartvar(ncid=ncid, varname='huilfexpnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='hui needed for leaf expansion of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huilfexpnp_patch)
  
         call restartvar(ncid=ncid, varname='huigrnnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='hui needed for starting grainfill of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huigrnnp_patch)
  
         call restartvar(ncid=ncid, varname='huilfmatnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer',long_name='hui needed for leaf maturity of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huilfmatnp_patch)
  
         call restartvar(ncid=ncid, varname='grnmatnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='hui needed for grain maturity of successive phytomers',&
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%grnmatnp_patch)
  
         call restartvar(ncid=ncid, varname='huilfsennp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer', & 
              long_name='hui needed for leaf senescence of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huilfsennp_patch)
  
         call restartvar(ncid=ncid, varname='huilfendnp', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer',long_name='hui needed for end of life of successive phytomers', &
              units='', flag=flag, interpinic_flag='interp', readvar=readvar, data=this%huilfendnp_patch)
  
         !call restartvar(ncid=ncid, varname='harvest_flag', xtype=ncd_double,  &
              !dim1name='pft',long_name='harvest flag',units='unitless', & 
              !flag=flag, interpinic_flag='interp', readvar=readvar, data=this%harvest_flag_patch)
  
         call restartvar(ncid=ncid, varname='harvest_counter', xtype=ncd_double,  &
              dim1name='pft',long_name='harvest counter',units='unitless', &
              flag=flag, interpinic_flag='interp', readvar=readvar, data=this%harvest_counter_patch)
  
         call restartvar(ncid=ncid, varname='bglfr_p', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer',long_name='',units='', &
              flag=flag, interpinic_flag='interp', readvar=readvar, data=this%bglfr_p_patch)
  
         call restartvar(ncid=ncid, varname='bgtr_p', xtype=ncd_double,  &
              dim1name='pft',dim2name='phytomer',long_name='',units='', &
              flag=flag, interpinic_flag='interp', readvar=readvar, data=this%bgtr_p_patch)

        ! call restartvar(ncid=ncid, varname='prune', xtype=ncd_log,  &
        !      dim1name='pft',long_name='Flag if true do pruning', &
        !      flag=flag, interpinic_flag='interp', readvar=readvar,
        !      data=this%prune_patch)
  
         allocate(temp1d(bounds%begp:bounds%endp))
         if (flag == 'write') then
            do p= bounds%begp,bounds%endp
               if (this%prune_patch(p)) then
                  temp1d(p) = 1
               else
                  temp1d(p) = 0
               end if
            end do
         end if
         call restartvar(ncid=ncid, flag=flag,  varname='prune', xtype=ncd_log, &
              dim1name='pft', &
              long_name='Flag, if true do pruning', &
              interpinic_flag='interp', readvar=readvar, data=temp1d)
         if (flag == 'read') then
            do p= bounds%begp,bounds%endp
               if (temp1d(p) == 1) then
                  this%prune_patch(p) = .true.
               else
                  this%prune_patch(p) = .false.
               end if
            end do
         end if
         deallocate(temp1d)

 
       end if 

       if (flag=='read' )then
          call this%checkDates( )  ! Check that restart date is same calendar date (even if year is different)
                                   ! This is so that it properly goes through
                                   ! the crop phases
       end if
    end if

  end subroutine Restart


  !-----------------------------------------------------------------------
  subroutine CropUpdateAccVars(this, bounds, t_ref2m_patch, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables. Should be called every time step.
    ! Should only be called if use_crop is true.
    !
    ! !USES:
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep
    use clm_varpar       , only : nlevsno, nlevgrnd
    use pftconMod        , only : nswheat, nirrig_swheat, pftcon
    use pftconMod        , only : nwwheat, nirrig_wwheat
    use pftconMod        , only : nsugarcane, nirrig_sugarcane
    use ColumnType       , only : col
    use PatchType        , only : patch
    !
    ! !ARGUMENTS:
    implicit none
    class(crop_type)       , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    real(r8)               , intent(in)    :: t_ref2m_patch( bounds%begp:)
    real(r8)               , intent(inout) :: t_soisno_col(bounds%begc:, -nlevsno+1:)
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,g ! indices
    integer :: ivt   ! vegetation type
    integer :: dtime ! timestep size [seconds]
    integer :: nstep ! timestep number
    integer :: ier   ! error status
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    character(len=*), parameter :: subname = 'CropUpdateAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(t_ref2m_patch)  == (/endp/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col)   == (/endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

    dtime = get_step_size()
    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract GDDPLANT
    
    call extract_accum_field ('GDDPLANT', rbufslp, nstep)
    do p = begp,endp
      rbufslp(p) = max(0.0_r8,this%gddplant_patch(p)-rbufslp(p))
    end do
    call update_accum_field  ('GDDPLANT', rbufslp, nstep)
    do p = begp,endp
       if (this%croplive_patch(p)) then ! relative to planting date
          ivt = patch%itype(p)
          if ( (trim(this%baset_mapping) == baset_map_latvary) .and. &
             ((ivt == nswheat) .or. (ivt == nirrig_swheat) .or. &
              (ivt == nsugarcane) .or. (ivt == nirrig_sugarcane)) ) then
             rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
             t_ref2m_patch(p)-(SHR_CONST_TKFRZ + this%latbaset_patch(p)))) &
             * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
             t_ref2m_patch(p)-(SHR_CONST_TKFRZ + pftcon%baset(ivt)))) &
             * dtime/SHR_CONST_CDAY
          end if
          if (ivt == nwwheat .or. ivt == nirrig_wwheat) then
             rbufslp(p) = rbufslp(p) * this%vf_patch(p)
          end if
       else
          rbufslp(p) = accumResetVal
       end if
    end do
    call update_accum_field  ('GDDPLANT', rbufslp, nstep)
    call extract_accum_field ('GDDPLANT', this%gddplant_patch, nstep)

    ! Accumulate and extract GDDTSOI
    ! In agroibis this variable is calculated
    ! to 0.05 m, so here we use the top two soil layers

    do p = begp,endp
       if (this%croplive_patch(p)) then ! relative to planting date
          ivt = patch%itype(p)
          c   = patch%column(p)
          rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
               ((t_soisno_col(c,1)*col%dz(c,1) + &
               t_soisno_col(c,2)*col%dz(c,2))/(col%dz(c,1)+col%dz(c,2))) - &
               (SHR_CONST_TKFRZ + pftcon%baset(ivt)))) * dtime/SHR_CONST_CDAY
          if (ivt == nwwheat .or. ivt == nwwheat) then
             rbufslp(p) = rbufslp(p) * this%vf_patch(p)
          end if
       else
          rbufslp(p) = accumResetVal
       end if
    end do
    call update_accum_field  ('GDDTSOI', rbufslp, nstep)
    call extract_accum_field ('GDDTSOI', this%gddtsoi_patch, nstep)

    deallocate(rbufslp)

  end subroutine CropUpdateAccVars

  !-----------------------------------------------------------------------
  subroutine CropIncrementYear (this, num_pcropp, filter_pcropp)
    !
    ! !DESCRIPTION: 
    ! Increment the crop year, if appropriate
    !
    ! This routine should be called every time step
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date, is_first_step
    !
    ! !ARGUMENTS:
    class(crop_type) :: this
    integer , intent(in) :: num_pcropp       ! number of prog. crop patches in filter
    integer , intent(in) :: filter_pcropp(:) ! filter for prognostic crop patches
    !
    ! !LOCAL VARIABLES:
    integer kyr   ! current year
    integer kmo   ! month of year  (1, ..., 12)
    integer kda   ! day of month   (1, ..., 31)
    integer mcsec ! seconds of day (0, ..., seconds/day)
    integer :: fp, p
    !-----------------------------------------------------------------------

    call get_curr_date (   kyr, kmo, kda, mcsec)
    ! Update nyrs when it's the end of the year (unless it's the very start of the
    ! run). This assumes that, if this patch is active at the end of the year, then it was
    ! active for the whole year.
    if ((kmo == 1 .and. kda == 1 .and. mcsec == 0) .and. .not. is_first_step()) then
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)

          this%nyrs_crop_active_patch(p) = this%nyrs_crop_active_patch(p) + 1
       end do
    end if

  end subroutine CropIncrementYear

  !-----------------------------------------------------------------------
  subroutine checkDates( )
    !
    ! !DESCRIPTION: 
    ! Make sure the dates are compatible. The date given to startup the model
    ! and the date on the restart file must be the same although years can be
    ! different. The dates need to be checked when the restart file is being
    ! read in for a startup or branch case (they are NOT allowed to be different
    ! for a restart case).
    !
    ! For the prognostic crop model the date of planting is tracked and growing
    ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
    ! messes up these bits of saved information.
    !
    ! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
    !
    ! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !-----------------------------------------------------------------------
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( masterproc ) &
            write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
               (stmon_day - stmon_day/100)
          call endrun(msg=' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
               ' and the date on the restart file needs to match (years can be different)'//&
               errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine checkDates

end module CropType

