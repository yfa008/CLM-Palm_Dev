module CNPhenologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CNPhenologyMod
  !
  ! !DESCRIPTION:
  ! Module holding routines used in phenology model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use shr_sys_mod                     , only : shr_sys_flush
  use decompMod                       , only : bounds_type
  use clm_varpar                      , only : maxveg, nlevdecomp_full
  use clm_varctl                      , only : iulog, use_cndv
  use clm_varcon                      , only : tfrz
  use abortutils                      , only : endrun
  use CanopyStateType                 , only : canopystate_type
  use CNDVType                        , only : dgvs_type
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  use CropType                        , only : crop_type
  use pftconMod                       , only : pftcon
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type
  use WaterDiagnosticBulkType                  , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType                  , only : wateratm2lndbulk_type
  use ColumnType                      , only : col                
  use GridcellType                    , only : grc                
  use PatchType                       , only : patch   
  use atm2lndType                     , only : atm2lnd_type             
  use atm2lndType                     , only : atm2lnd_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams           ! Read parameters
  public :: CNPhenologyreadNML   ! Read namelist
  public :: CNPhenologyInit      ! Initialization
  public :: CNPhenology          ! Update
  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: params_type
     real(r8) :: crit_dayl       ! critical day length for senescence
     real(r8) :: ndays_on     	 ! number of days to complete leaf onset
     real(r8) :: ndays_off	 ! number of days to complete leaf offset
     real(r8) :: fstor2tran      ! fraction of storage to move to transfer for each onset
     real(r8) :: crit_onset_fdd  ! critical number of freezing days to set gdd counter
     real(r8) :: crit_onset_swi  ! critical number of days > soilpsi_on for onset
     real(r8) :: soilpsi_on      ! critical soil water potential for leaf onset
     real(r8) :: crit_offset_fdd ! critical number of freezing days to initiate offset
     real(r8) :: crit_offset_swi ! critical number of water stress days to initiate offset
     real(r8) :: soilpsi_off     ! critical soil water potential for leaf offset
     real(r8) :: lwtop   	 ! live wood turnover proportion (annual fraction)
  end type params_type

  type(params_type) :: params_inst

  real(r8) :: dt                            ! radiation time step delta t (seconds)
  real(r8) :: fracday                       ! dtime as a fraction of day
  real(r8) :: crit_dayl                     ! critical daylength for offset (seconds)
  real(r8) :: ndays_on                      ! number of days to complete onset
  real(r8) :: ndays_off                     ! number of days to complete offset
  real(r8) :: fstor2tran                    ! fraction of storage to move to transfer on each onset
  real(r8) :: crit_onset_fdd                ! critical number of freezing days
  real(r8) :: crit_onset_swi                ! water stress days for offset trigger
  real(r8) :: soilpsi_on                    ! water potential for onset trigger (MPa)
  real(r8) :: crit_offset_fdd               ! critical number of freezing degree days to trigger offset
  real(r8) :: crit_offset_swi               ! water stress days for offset trigger
  real(r8) :: soilpsi_off                   ! water potential for offset trigger (MPa)
  real(r8) :: lwtop                         ! live wood turnover proportion (annual fraction)

  ! CropPhenology variables and constants
  real(r8) :: p1d, p1v                      ! photoperiod factor constants for crop vernalization
  real(r8) :: hti                           ! cold hardening index threshold for vernalization
  real(r8) :: tbase                         ! base temperature for vernalization

  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
  integer, parameter :: inNH       = 1      ! Northern Hemisphere
  integer, parameter :: inSH       = 2      ! Southern Hemisphere
  integer, pointer   :: inhemi(:)           ! Hemisphere that patch is in 

  integer, allocatable :: minplantjday(:,:) ! minimum planting julian day
  integer, allocatable :: maxplantjday(:,:) ! maximum planting julian day
  integer              :: jdayyrstart(inSH) ! julian day of start of year

  real(r8), private :: initial_seed_at_planting = 3._r8 ! Initial seed at planting

  !PalmPhenolgy variables
  integer            :: mat1                 !the first mature phytomer (oil palm model)
  integer, parameter :: NOT_Emerged = 9999   ! If leaf not emerged yet
  real(r8), allocatable :: phyllochron2(:)    !phyllochron increases through maturity
  real(r8), allocatable :: huilfexp(:)        !hui needed from leaf initiation to leaf expansion for all phytomers
  real(r8), allocatable :: huilfmat(:)        !hui needed from leaf expansion to leaf maturity for all phytomers
  real(r8), allocatable :: huilfsen(:)        !hui needed from leaf expansion to leaf senescence for all phytomers
  real(r8), allocatable :: huilfend(:)        !hui needed from leaf expansion to leaf end of life for all phytomers

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNPhenologyReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNPhenology
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNPhenologyReadNML'
    character(len=*), parameter :: nmlname = 'cnphenology'
    !-----------------------------------------------------------------------
    namelist /cnphenology/ initial_seed_at_planting

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnphenology, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_seed_at_planting, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnphenology)
       write(iulog,*) ' '
    end if


    !-----------------------------------------------------------------------
    
  end subroutine CNPhenologyReadNML

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNPhenolParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    !
    ! read in parameters
    !   
    tString='crit_dayl'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_dayl=tempr

    tString='ndays_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_on=tempr

    tString='ndays_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_off=tempr

    tString='fstor2tran'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%fstor2tran=tempr

    tString='crit_onset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_onset_fdd=tempr

    tString='crit_onset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_onset_swi=tempr

    tString='soilpsi_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%soilpsi_on=tempr

    tString='crit_offset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_offset_fdd=tempr

    tString='crit_offset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_offset_swi=tempr

    tString='soilpsi_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%soilpsi_off=tempr

    tString='lwtop_ann'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%lwtop=tempr   

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, &
       filter_soilp, num_pcropp, filter_pcropp, &
       doalb, waterdiagnosticbulk_inst, wateratm2lndbulk_inst, temperature_inst, atm2lnd_inst, crop_inst, &
       canopystate_inst, soilstate_inst, dgvs_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, soilbiogeochem_state_inst, &
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       leaf_prof_patch, froot_prof_patch, phase)
    ! !USES:
    use CNSharedParamsMod, only: use_fun
    use pftconMod        , only: mxnp  ! max number of phytomers for oil palm
    !
    ! !DESCRIPTION:
    ! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
    ! 1. grass phenology
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                        , intent(in)    :: num_pcropp      ! number of prog. crop patches in filter
    integer                        , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    logical                        , intent(in)    :: doalb           ! true if time for sfc albedo calc
    type(waterdiagnosticbulk_type)          , intent(in)    :: waterdiagnosticbulk_inst
    type(wateratm2lndbulk_type)          , intent(in)    :: wateratm2lndbulk_inst
    type(temperature_type)         , intent(inout) :: temperature_inst
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(dgvs_type)                , intent(inout) :: dgvs_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_state_type), intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    integer                        , intent(in)    :: phase
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))

    ! each of the following phenology type routines includes a filter
    ! to operate only on the relevant patches


    if ( phase == 1 ) then
       call CNPhenologyClimate(num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
            temperature_inst, cnveg_state_inst, crop_inst)
   
       call CNEvergreenPhenology(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNSeasonDecidPhenology(num_soilp, filter_soilp, &
            temperature_inst, cnveg_state_inst, dgvs_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNStressDecidPhenology(num_soilp, filter_soilp,   &
            soilstate_inst, temperature_inst, atm2lnd_inst, wateratm2lndbulk_inst, cnveg_state_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       if (doalb .and. num_pcropp > 0 ) then

         !palm and other perennial crops use the same perennial flag
         !only set evergreen flag for oil palm, otherwise the it will call the CNEvergreenPhenology
         !PalmPhenology first in the condition because palms are crop and perennial, only call one of them
         if ( mxnp > 0) then
          call PalmPhenology(num_pcropp, filter_pcropp, &
               waterdiagnosticbulk_inst, temperature_inst, crop_inst, canopystate_inst, cnveg_state_inst, &
               cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
               c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
         else
	 call CropPhenology(num_pcropp, filter_pcropp, &
               waterdiagnosticbulk_inst, temperature_inst, soilstate_inst, wateratm2lndbulk_inst, & 
	       crop_inst, canopystate_inst, cnveg_state_inst, &
               cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
               c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
         end if
       end if
    else if ( phase == 2 ) then
       ! the same onset and offset routines are called regardless of
       ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

       call CNOnsetGrowth(num_soilp, filter_soilp, &
            cnveg_state_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNOffsetLitterfall(num_soilp, filter_soilp, &
            num_soilc, filter_soilc, soilbiogeochem_state_inst,&
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, crop_inst)

       call CNBackgroundLitterfall(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, crop_inst)
   
       call CNLivewoodTurnover(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNGrainToProductPools(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, &
            cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       ! gather all patch-level litterfall fluxes to the column for litter C and N inputs

       call CNLitterToColumn(bounds, num_soilc, filter_soilc, &
            cnveg_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, & ! soilbiogeochem_state_inst, &
            crop_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst, &
            leaf_prof_patch(bounds%begp:bounds%endp,1:nlevdecomp_full), & 
            froot_prof_patch(bounds%begp:bounds%endp,1:nlevdecomp_full))
    else
       call endrun( 'bad phase' )
    end if

  end subroutine CNPhenology

  !-----------------------------------------------------------------------
  subroutine CNPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CNPhenology. Must be called after time-manager is
    ! initialized, and after pftcon file is read in.
    !
    ! !USES:
    use clm_time_manager, only: get_step_size
    use clm_varctl      , only: use_crop
    use clm_varcon      , only: secspday
    use pftconMod       , only: mxnp  ! max number of phytomers for oil palm
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    dt      = real( get_step_size(), r8 )
    fracday = dt/secspday

    ! set constants for CNSeasonDecidPhenology 
    ! (critical daylength from Biome-BGC, v4.1.2)
    crit_dayl=params_inst%crit_dayl

    ! Set constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ndays_on=params_inst%ndays_on
    ndays_off=params_inst%ndays_off
    
    ndays_off=1.0_r8

    ! set transfer parameters
    fstor2tran=params_inst%fstor2tran

    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd=params_inst%crit_onset_fdd
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi=params_inst%crit_onset_swi
    soilpsi_on=params_inst%soilpsi_on

    ! offset parameters
    crit_offset_fdd=params_inst%crit_offset_fdd
    crit_offset_swi=params_inst%crit_offset_swi
    soilpsi_off=params_inst%soilpsi_off    

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop=params_inst%lwtop/31536000.0_r8 !annual fraction converted to per second

    ! -----------------------------------------
    ! Call any subroutine specific initialization routines
    ! -----------------------------------------

    if ( use_crop ) call CropPhenologyInit(bounds)

    !Initialize private variables for PalmPhenology (Y.Fan)
    if ( mxnp > 0 )then
       allocate(huilfexp(bounds%begp:bounds%endp))
       allocate(huilfmat(bounds%begp:bounds%endp))
       allocate(huilfsen(bounds%begp:bounds%endp))
       allocate(huilfend(bounds%begp:bounds%endp))
       allocate(phyllochron2(bounds%begp:bounds%endp))
       phyllochron2(bounds%begp:bounds%endp) = nan
       huilfexp(bounds%begp:bounds%endp) = nan
       huilfmat(bounds%begp:bounds%endp) = nan
       huilfsen(bounds%begp:bounds%endp) = nan
       huilfend(bounds%begp:bounds%endp) = nan
    end if

  end subroutine CNPhenologyInit

  !-----------------------------------------------------------------------
  subroutine CNPhenologyClimate (num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
       temperature_inst, cnveg_state_inst, crop_inst)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
    use clm_time_manager , only : get_curr_date, is_first_step
    use clm_varctl      , only: use_crop
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                , intent(in)    :: num_pcropp      ! number of prognostic crops in filter
    integer                , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    type(temperature_type) , intent(inout) :: temperature_inst
    type(cnveg_state_type) , intent(inout) :: cnveg_state_inst
    type(crop_type)        , intent(inout) :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p       ! indices
    integer  :: fp      ! lake filter patch index
    real(r8) :: dayspyr ! days per year (days)
    integer  :: kyr     ! current year
    integer  :: kmo     ! month of year  (1, ..., 12)
    integer  :: kda     ! day of month   (1, ..., 31)
    integer  :: mcsec   ! seconds of day (0, ..., seconds/day)
    real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
    real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
    !-----------------------------------------------------------------------

    associate(                                                & 
         nyrs_crop_active => crop_inst%nyrs_crop_active_patch,   & ! InOut:  [integer (:)  ]  number of years this crop patch has been active
         
         t_ref2m        => temperature_inst%t_ref2m_patch ,   & ! Input:  [real(r8) (:) ]  2m air temperature (K)
         gdd0           => temperature_inst%gdd0_patch    ,   & ! Output: [real(r8) (:) ]  growing deg. days base 0 deg C (ddays)            
         gdd8           => temperature_inst%gdd8_patch    ,   & ! Output: [real(r8) (:) ]     "     "    "    "   8  "  "    "               
         gdd10          => temperature_inst%gdd10_patch   ,   & ! Output: [real(r8) (:) ]     "     "    "    "  10  "  "    "               
         gdd020         => temperature_inst%gdd020_patch  ,   & ! Output: [real(r8) (:) ]  20-yr mean of gdd0 (ddays)                        
         gdd820         => temperature_inst%gdd820_patch  ,   & ! Output: [real(r8) (:) ]  20-yr mean of gdd8 (ddays)                        
         gdd1020        => temperature_inst%gdd1020_patch ,   & ! Output: [real(r8) (:) ]  20-yr mean of gdd10 (ddays)                       
         gdd15          => temperature_inst%gdd15_patch   ,   & ! Output: [real(r8) (:) ]  growing deg. days base 15 deg C (ddays) (Y.Fan)
         gdd1520        => temperature_inst%gdd1520_patch ,   & ! Output: [real(r8) (:) ]  20-yr mean of gdd15 (ddays)  (added by Y.Fan)
         
         tempavg_t2m    => cnveg_state_inst%tempavg_t2m_patch & ! Output: [real(r8) (:) ]  temp. avg 2m air temperature (K)                  
         )

      ! set time steps
      
      dayspyr = get_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         tempavg_t2m(p) = tempavg_t2m(p) + t_ref2m(p) * (fracday/dayspyr)
      end do

      !
      ! The following crop related steps are done here rather than CropPhenology
      ! so that they will be completed each time-step rather than with doalb.
      !
      ! The following lines come from ibis's climate.f + stats.f
      ! gdd SUMMATIONS ARE RELATIVE TO THE PLANTING DATE (see subr. updateAccFlds)
      !!not true: GDD0/GDD8/GDD10 are tracked from April to Sep NH every year. only gddplant and gddtsoi are relatvie to planting date (Y.Fan 2015)

      !if ( use_crop )then ! use this switch so that crop functions can turn on for all pfts, useful for crop experiments (Y.Fan)
      !the num_pcropp filter only counts active crop pfts that have grid weight > 0
      if (num_pcropp > 0) then
         ! get time-related info
         call get_curr_date(kyr, kmo, kda, mcsec)
      end if

!num_soilp filter includes both natural pfts and crops that are active in the column, for testing shrub for crop (Y.Fan 2016)
!      do fp = 1,num_soilp   
!         p = filter_soilp(fp)
      do fp = 1,num_pcropp
         p = filter_pcropp(fp)
         if (kmo == 1 .and. kda == 1 .and. nyrs_crop_active(p) == 0) then ! YR 1:
            gdd020(p)  = 0._r8                             ! set gdd..20 variables to 0
            gdd820(p)  = 0._r8                             ! and crops will not be planted
            gdd1020(p) = 0._r8
            gdd1520(p) = 0._r8
         end if
         if (kmo == 1 .and. kda == 1 .and. mcsec == 0) then        ! <-- END of EVERY YR:
            if (nyrs_crop_active(p) == 1) then                     ! <-- END of YR 1
               gdd020(p)  = gdd0(p)                                ! <-- END of YR 1
               gdd820(p)  = gdd8(p)                                ! <-- END of YR 1
               gdd1020(p) = gdd10(p)                               ! <-- END of YR 1
               gdd1520(p) = gdd15(p)                               ! <-- END of YR 1
            end if                                                 ! <-- END of YR 1
            gdd020(p)  = (yravgm1* gdd020(p)  + gdd0(p))  / yravg  ! gdd..20 must be long term avgs
            gdd820(p)  = (yravgm1* gdd820(p)  + gdd8(p))  / yravg  ! so ignore results for yrs 1 & 2
            gdd1020(p) = (yravgm1* gdd1020(p) + gdd10(p)) / yravg
            gdd1520(p) = (yravgm1* gdd1520(p) + gdd15(p)) / yravg  !GDD base 15 degree for tropical crops (Y.Fan)
         end if
      end do

    end associate

  end subroutine CNPhenologyClimate

  !-----------------------------------------------------------------------
  subroutine CNEvergreenPhenology (num_soilp, filter_soilp , &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)   
       ! cnveg_state_inst) 
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use clm_varcon       , only : secspday
    use clm_time_manager , only : get_days_per_year
    use clm_varctl       , only : CN_evergreen_phenology_opt   
    !
    ! !ARGUMENTS:
    integer           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type), intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst    
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst  
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst     
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst   
    !
    ! !LOCAL VARIABLES:
    real(r8):: dayspyr                    ! Days per year
    integer :: p                          ! indices
    integer :: fp                         ! lake filter patch index
    
    real(r8):: tranr 				      
    real(r8):: t1                         ! temporary variable 
    !-----------------------------------------------------------------------

    associate(                                        & 
         ivt        => patch%itype                    , & ! Input:  [integer  (:) ]  patch vegetation type                                

         evergreen  => pftcon%evergreen             , & ! Input:  binary flag for evergreen leaf habit (0 or 1)     
         leaf_long  => pftcon%leaf_long             , & ! Input:  leaf longevity (yrs)  
         
         woody                               =>    pftcon%woody                                                         , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)     
         
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                             
   		 frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                         
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                         
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                         
   		 livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                  
   		 deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                  
   		 gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage   
   		 leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                              , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                            
   		 frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                             , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer                       
   		 livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                          , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                       
   		 deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                          , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer                       
   		 livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                         , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C transfer                
   		 deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                         , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer   
   		                
   		 leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                         , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                              
   		 frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                         
   		 livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                         
   		 deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                         
   		 livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch                    , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                  
   		 deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch                    , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage               
   		 leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                            , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N transfer                            
   		 frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                           , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N transfer                       
   		 livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                        , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N transfer                       
   		 deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                        , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N transfer                       
   		 livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                       , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N transfer                
   		 deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                       , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer     
   		                 
   		 leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch                    , & ! InOut:  [real(r8) (:)]                                                     
   		 frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch                   , & ! InOut:  [real(r8) (:)]                                                     
   		 livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch                , & ! InOut:  [real(r8) (:)]                                                     
   		 deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch                , & ! InOut:  [real(r8) (:)]                                                     
   		 livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch               , & ! InOut:  [real(r8) (:)]                                                     
   		 deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch               , & ! InOut:  [real(r8) (:)]                                                     
   		 gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch                    , & ! InOut:  [real(r8) (:)]  
   		 leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch                      , & ! InOut:  [real(r8) (:)]                                                    
   		 frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch                    , & ! InOut:  [real(r8) (:)]                                                    
   		 livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch              , & ! InOut:  [real(r8) (:)]                                                    
   		 deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch              , & ! InOut:  [real(r8) (:)]                                                    
   		 livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch            , & ! InOut:  [real(r8) (:)]   
   		                                                    
   		 leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch                  , & ! InOut:  [real(r8) (:)]                                                     
   		 frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch                 , & ! InOut:  [real(r8) (:)]                                                     
   		 livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch              , & ! InOut:  [real(r8) (:)]                                                     
   		 deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch              , & ! InOut:  [real(r8) (:)]                                                     
   		 livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch             , & ! InOut:  [real(r8) (:)]   
   		 deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch             , & ! InOut:  [real(r8) (:)]                                                     
   		 leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch                    , & ! InOut:  [real(r8) (:)]                                                    
   		 frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch                  , & ! InOut:  [real(r8) (:)]                                                    
   		 livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch          , & ! InOut:  [real(r8) (:)]                                                    
   		 deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch          , & ! InOut:  [real(r8) (:)]     
   		           
         bglfr      => cnveg_state_inst%bglfr_patch , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)                  
         bgtr       => cnveg_state_inst%bgtr_patch  , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)             
         lgsf       => cnveg_state_inst%lgsf_patch    & ! Output: [real(r8) (:) ]  long growing season factor [0-1]                  
         )

      dayspyr   = get_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if (evergreen(ivt(p)) == 1._r8) then
            bglfr(p) = 1._r8/(leaf_long(ivt(p)) * dayspyr * secspday)
            bgtr(p)  = 0._r8
            lgsf(p)  = 0._r8
         end if
      end do
               
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (CN_evergreen_phenology_opt == 1) then    
   do fp = 1,num_soilp    
      p = filter_soilp(fp)    
      if (evergreen(ivt(p)) == 1._r8) then    
   
         tranr=0.0002_r8   
         ! set carbon fluxes for shifting storage pools to transfer pools    
         leafc_storage_to_xfer(p)  = tranr * leafc_storage(p)/dt    
         frootc_storage_to_xfer(p) = tranr * frootc_storage(p)/dt    
         if (woody(ivt(p)) == 1.0_r8) then    
            livestemc_storage_to_xfer(p)  = tranr * livestemc_storage(p)/dt    
            deadstemc_storage_to_xfer(p)  = tranr * deadstemc_storage(p)/dt    
            livecrootc_storage_to_xfer(p) = tranr * livecrootc_storage(p)/dt   
            deadcrootc_storage_to_xfer(p) = tranr * deadcrootc_storage(p)/dt   
            gresp_storage_to_xfer(p)      = tranr * gresp_storage(p)/dt        
         end if    

        ! set nitrogen fluxes for shifting storage pools to transfer pools    
        leafn_storage_to_xfer(p)  = tranr * leafn_storage(p)/dt    
        frootn_storage_to_xfer(p) = tranr * frootn_storage(p)/dt   
        if (woody(ivt(p)) == 1.0_r8) then    
            livestemn_storage_to_xfer(p)  = tranr * livestemn_storage(p)/dt    
            deadstemn_storage_to_xfer(p)  = tranr * deadstemn_storage(p)/dt    
            livecrootn_storage_to_xfer(p) = tranr * livecrootn_storage(p)/dt   
            deadcrootn_storage_to_xfer(p) = tranr * deadcrootn_storage(p)/dt   
        end if    
                        
        t1 = 1.0_r8 / dt   
            
        leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)    
        frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)   
            
        leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)    
        frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)   
        if (woody(ivt(p)) == 1.0_r8) then   
            livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)   
            deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)   
            livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)  
            deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)  
                
            livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)   
            deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)   
            livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)  
            deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)  
        end if
                
      end if ! end of if (evergreen(ivt(p)) == 1._r8) then    
     
   end do ! end of pft loop 
   
   end if ! end of if (CN_evergreen_phenology_opt == 1) then    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

  end subroutine CNEvergreenPhenology

  !-----------------------------------------------------------------------
  subroutine CNSeasonDecidPhenology (num_soilp, filter_soilp       , &
       temperature_inst, cnveg_state_inst, dgvs_inst , &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    ! This routine handles the seasonal deciduous phenology code (temperate
    ! deciduous vegetation that has only one growing season per year).
    !
    ! !USES:
    use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
    use clm_varcon      , only: secspday
    use clm_varctl      , only: use_cndv
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(dgvs_type)                , intent(inout) :: dgvs_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p          !indices
    integer :: fp             !lake filter patch index
    real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
    real(r8):: crit_onset_gdd !critical onset growing degree-day sum
    real(r8):: soilt
    !-----------------------------------------------------------------------

    associate(                                                                                                   & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:)   ]  patch vegetation type                                
         dayl                                =>    grc%dayl                                                    , & ! Input:  [real(r8)  (:)   ]  daylength (s)
         prev_dayl                           =>    grc%prev_dayl                                               , & ! Input:  [real(r8)  (:)   ]  daylength from previous time step (s)
         
         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         season_decid                        =>    pftcon%season_decid                                         , & ! Input:  binary flag for seasonal-deciduous leaf habit (0 or 1)
         
         t_soisno                            =>    temperature_inst%t_soisno_col                               , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         
         pftmayexist                         =>    dgvs_inst%pftmayexist_patch                                 , & ! Output: [logical   (:)   ]  exclude seasonal decid patches from tropics           

         annavg_t2m                          =>    cnveg_state_inst%annavg_t2m_patch                           , & ! Input:  [real(r8)  (:)   ]  annual average 2m air temperature (K)             
         dormant_flag                        =>    cnveg_state_inst%dormant_flag_patch                         , & ! Output: [real(r8)  (:)   ]  dormancy flag                                     
         days_active                         =>    cnveg_state_inst%days_active_patch                          , & ! Output: [real(r8)  (:)   ]  number of days since last dormancy                
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Output: [real(r8)  (:)   ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Output: [real(r8)  (:)   ]  onset counter (seconds)                           
         onset_gddflag                       =>    cnveg_state_inst%onset_gddflag_patch                        , & ! Output: [real(r8)  (:)   ]  onset freeze flag                                 
         onset_gdd                           =>    cnveg_state_inst%onset_gdd_patch                            , & ! Output: [real(r8)  (:)   ]  onset growing degree days                         
         offset_flag                         =>    cnveg_state_inst%offset_flag_patch                          , & ! Output: [real(r8)  (:)   ]  offset flag                                       
         offset_counter                      =>    cnveg_state_inst%offset_counter_patch                       , & ! Output: [real(r8)  (:)   ]  offset counter (seconds)                          
         bglfr                               =>    cnveg_state_inst%bglfr_patch                                , & ! Output: [real(r8)  (:)   ]  background litterfall rate (1/s)                  
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Output: [real(r8)  (:)   ]  background transfer growth rate (1/s)             
         lgsf                                =>    cnveg_state_inst%lgsf_patch                                 , & ! Output: [real(r8)  (:)   ]  long growing season factor [0-1]                  
         
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer               
         
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer               

         prev_leafc_to_litter                =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch            , & ! Output: [real(r8)  (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch           , & ! Output: [real(r8)  (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch          , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch        , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch      & ! Output:  [real(r8) (:)   ]                                                    
         )

      ! start patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         if (season_decid(ivt(p)) == 1._r8) then

            ! set background litterfall rate, background transfer rate, and
            ! long growing season factor to 0 for seasonal deciduous types
            bglfr(p) = 0._r8
            bgtr(p) = 0._r8
            lgsf(p) = 0._r8

            ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))

            ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
            if (dayl(g) >= prev_dayl(g)) then
               ws_flag = 1._r8
            else
               ws_flag = 0._r8
            end if

            ! update offset_counter and test for the end of the offset period
            if (offset_flag(p) == 1.0_r8) then
               ! decrement counter for offset period
               offset_counter(p) = offset_counter(p) - dt

               ! if this is the end of the offset_period, reset phenology
               ! flags and indices
               if (offset_counter(p) == 0.0_r8) then
                  ! this code block was originally handled by call cn_offset_cleanup(p)
                  ! inlined during vectorization

                  offset_flag(p) = 0._r8
                  offset_counter(p) = 0._r8
                  dormant_flag(p) = 1._r8
                  days_active(p) = 0._r8
                  if (use_cndv) then
                     pftmayexist(p) = .true.
                  end if

                  ! reset the previous timestep litterfall flux memory
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! update onset_counter and test for the end of the onset period
            if (onset_flag(p) == 1.0_r8) then
               ! decrement counter for onset period
               onset_counter(p) = onset_counter(p) - dt

               ! if this is the end of the onset period, reset phenology
               ! flags and indices
               if (onset_counter(p) == 0.0_r8) then
                  ! this code block was originally handled by call cn_onset_cleanup(p)
                  ! inlined during vectorization

                  onset_flag(p) = 0.0_r8
                  onset_counter(p) = 0.0_r8
                  ! set all transfer growth rates to 0.0
                  leafc_xfer_to_leafc(p)   = 0.0_r8
                  frootc_xfer_to_frootc(p) = 0.0_r8
                  leafn_xfer_to_leafn(p)   = 0.0_r8
                  frootn_xfer_to_frootn(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0.0_r8
                     deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
                     livecrootc_xfer_to_livecrootc(p) = 0.0_r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
                     livestemn_xfer_to_livestemn(p)   = 0.0_r8
                     deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
                     livecrootn_xfer_to_livecrootn(p) = 0.0_r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0.0_r8
                  leafn_xfer(p) = 0.0_r8
                  frootc_xfer(p) = 0.0_r8
                  frootn_xfer(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0.0_r8
                     livestemn_xfer(p) = 0.0_r8
                     deadstemc_xfer(p) = 0.0_r8
                     deadstemn_xfer(p) = 0.0_r8
                     livecrootc_xfer(p) = 0.0_r8
                     livecrootn_xfer(p) = 0.0_r8
                     deadcrootc_xfer(p) = 0.0_r8
                     deadcrootn_xfer(p) = 0.0_r8
                  end if
               end if
            end if

            ! test for switching from dormant period to growth period
            if (dormant_flag(p) == 1.0_r8) then

               ! Test to turn on growing degree-day sum, if off.
               ! switch on the growing degree day sum on the winter solstice

               if (onset_gddflag(p) == 0._r8 .and. ws_flag == 1._r8) then
                  onset_gddflag(p) = 1._r8
                  onset_gdd(p) = 0._r8
               end if

               ! Test to turn off growing degree-day sum, if on.
               ! This test resets the growing degree day sum if it gets past
               ! the summer solstice without reaching the threshold value.
               ! In that case, it will take until the next winter solstice
               ! before the growing degree-day summation starts again.

               if (onset_gddflag(p) == 1._r8 .and. ws_flag == 0._r8) then
                  onset_gddflag(p) = 0._r8
                  onset_gdd(p) = 0._r8
               end if

               ! if the gdd flag is set, and if the soil is above freezing
               ! then accumulate growing degree days for onset trigger

               soilt = t_soisno(c,3)
               if (onset_gddflag(p) == 1.0_r8 .and. soilt > SHR_CONST_TKFRZ) then
                  onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
               end if

               ! set onset_flag if critical growing degree-day sum is exceeded
               if (onset_gdd(p) > crit_onset_gdd) then
                  onset_flag(p) = 1.0_r8
                  dormant_flag(p) = 0.0_r8
                  onset_gddflag(p) = 0.0_r8
                  onset_gdd(p) = 0.0_r8
                  onset_counter(p) = ndays_on * secspday

                  ! move all the storage pools into transfer pools,
                  ! where they will be transfered to displayed growth over the onset period.
                  ! this code was originally handled with call cn_storage_to_xfer(p)
                  ! inlined during vectorization

                  ! set carbon fluxes for shifting storage pools to transfer pools
                  leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
                  frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                     deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                     livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                     deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                     gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
                  end if

                  ! set nitrogen fluxes for shifting storage pools to transfer pools
                  leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
                  frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                     deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                     livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                     deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
                  end if
               end if

               ! test for switching from growth period to offset period
            else if (offset_flag(p) == 0.0_r8) then
               if (use_cndv) then
                  ! If days_active > 355, then remove patch in
                  ! CNDVEstablishment at the end of the year.
                  ! days_active > 355 is a symptom of seasonal decid. patches occurring in
                  ! gridcells where dayl never drops below crit_dayl.
                  ! This results in TLAI>1e4 in a few gridcells.
                  days_active(p) = days_active(p) + fracday
                  if (days_active(p) > 355._r8) pftmayexist(p) = .false.
               end if

               ! only begin to test for offset daylength once past the summer sol
               if (ws_flag == 0._r8 .and. dayl(g) < crit_dayl) then
                  offset_flag(p) = 1._r8
                  offset_counter(p) = ndays_off * secspday
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

         end if ! end if seasonal deciduous

      end do ! end of patch loop

    end associate
 
  end subroutine CNSeasonDecidPhenology

  !-----------------------------------------------------------------------
  subroutine CNStressDecidPhenology (num_soilp, filter_soilp , &
       soilstate_inst, temperature_inst, atm2lnd_inst, wateratm2lndbulk_inst, cnveg_state_inst, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! This routine handles phenology for vegetation types, such as grasses and
    ! tropical drought deciduous trees, that respond to cold and drought stress
    ! signals and that can have multiple growing seasons in a given year.
    ! This routine allows for the possibility that leaves might persist year-round
    ! in the absence of a suitable stress trigger, by switching to an essentially
    ! evergreen habit, but maintaining a deciduous leaf longevity, while waiting
    ! for the next stress trigger.  This is in contrast to the seasonal deciduous
    ! algorithm (for temperate deciduous trees) that forces a single growing season
    ! per year.
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
    use CNSharedParamsMod, only : use_fun
    use clm_varcon       , only : secspday
    use shr_const_mod    , only : SHR_CONST_TKFRZ, SHR_CONST_PI
    use CNSharedParamsMod, only : CNParamsShareInst
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(wateratm2lndbulk_type)             , intent(in)    :: wateratm2lndbulk_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
    integer :: g,c,p           ! indices
    integer :: fp              ! lake filter patch index
    real(r8):: dayspyr         ! days per year
    real(r8):: crit_onset_gdd  ! degree days for onset trigger
    real(r8):: soilt           ! temperature of top soil layer
    real(r8):: psi             ! water stress of top soil layer
    real(r8):: rain_threshold  ! rain threshold for leaf on [mm]
    logical :: additional_onset_condition ! additional condition for leaf onset
    !-----------------------------------------------------------------------

    associate(                                                                                                   & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:)   ]  patch vegetation type                                
         dayl                                =>    grc%dayl                                                    , & ! Input:  [real(r8)  (:)   ]  daylength (s)
         
         prec10                              => wateratm2lndbulk_inst%prec10_patch                                     , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation
         leaf_long                           =>    pftcon%leaf_long                                            , & ! Input:  leaf longevity (yrs)                              
         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         stress_decid                        =>    pftcon%stress_decid                                         , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)
         
         soilpsi                             =>    soilstate_inst%soilpsi_col                                  , & ! Input:  [real(r8)  (:,:) ]  soil water potential in each soil layer (MPa)   
         
         t_soisno                            =>    temperature_inst%t_soisno_col                               , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         dormant_flag                        =>    cnveg_state_inst%dormant_flag_patch                         , & ! Output:  [real(r8) (:)   ]  dormancy flag                                     
         days_active                         =>    cnveg_state_inst%days_active_patch                          , & ! Output:  [real(r8) (:)   ]  number of days since last dormancy                
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Output:  [real(r8) (:)   ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Output:  [real(r8) (:)   ]  onset counter (seconds)                           
         onset_gddflag                       =>    cnveg_state_inst%onset_gddflag_patch                        , & ! Output:  [real(r8) (:)   ]  onset freeze flag                                 
         onset_fdd                           =>    cnveg_state_inst%onset_fdd_patch                            , & ! Output:  [real(r8) (:)   ]  onset freezing degree days counter                
         onset_gdd                           =>    cnveg_state_inst%onset_gdd_patch                            , & ! Output:  [real(r8) (:)   ]  onset growing degree days                         
         onset_swi                           =>    cnveg_state_inst%onset_swi_patch                            , & ! Output:  [real(r8) (:)   ]  onset soil water index                            
         offset_flag                         =>    cnveg_state_inst%offset_flag_patch                          , & ! Output:  [real(r8) (:)   ]  offset flag                                       
         offset_counter                      =>    cnveg_state_inst%offset_counter_patch                       , & ! Output:  [real(r8) (:)   ]  offset counter (seconds)                          
         offset_fdd                          =>    cnveg_state_inst%offset_fdd_patch                           , & ! Output:  [real(r8) (:)   ]  offset freezing degree days counter               
         offset_swi                          =>    cnveg_state_inst%offset_swi_patch                           , & ! Output:  [real(r8) (:)   ]  offset soil water index                           
         lgsf                                =>    cnveg_state_inst%lgsf_patch                                 , & ! Output:  [real(r8) (:)   ]  long growing season factor [0-1]                  
         bglfr                               =>    cnveg_state_inst%bglfr_patch                                , & ! Output:  [real(r8) (:)   ]  background litterfall rate (1/s)                  
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Output:  [real(r8) (:)   ]  background transfer growth rate (1/s)             
         annavg_t2m                          =>    cnveg_state_inst%annavg_t2m_patch                           , & ! Output:  [real(r8) (:)   ]  annual average 2m air temperature (K)             
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                          , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C 
     
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                         , & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer               
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer               
         
         prev_leafc_to_litter                =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch            , & ! Output:  [real(r8) (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch           , & ! Output:  [real(r8) (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch          , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch        , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch      & ! Output:  [real(r8) (:)   ]                                                    
         )

      ! set time steps
      dayspyr = get_days_per_year()

      ! specify rain threshold for leaf onset
      rain_threshold = 20._r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         if (stress_decid(ivt(p)) == 1._r8) then
            soilt = t_soisno(c,3)
            psi = soilpsi(c,3)

            ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))


            ! update offset_counter and test for the end of the offset period
            if (offset_flag(p) == 1._r8) then
               ! decrement counter for offset period
               offset_counter(p) = offset_counter(p) - dt

               ! if this is the end of the offset_period, reset phenology
               ! flags and indices
               if (offset_counter(p) == 0._r8) then
                  ! this code block was originally handled by call cn_offset_cleanup(p)
                  ! inlined during vectorization
                  offset_flag(p) = 0._r8
                  offset_counter(p) = 0._r8
                  dormant_flag(p) = 1._r8
                  days_active(p) = 0._r8

                  ! reset the previous timestep litterfall flux memory
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! update onset_counter and test for the end of the onset period
            if (onset_flag(p) == 1.0_r8) then
               ! decrement counter for onset period
               onset_counter(p) = onset_counter(p) - dt

               ! if this is the end of the onset period, reset phenology
               ! flags and indices
               if (onset_counter(p) == 0.0_r8) then
                  ! this code block was originally handled by call cn_onset_cleanup(p)
                  ! inlined during vectorization
                  onset_flag(p) = 0._r8
                  onset_counter(p) = 0._r8
                  ! set all transfer growth rates to 0.0
                  leafc_xfer_to_leafc(p)   = 0._r8
                  frootc_xfer_to_frootc(p) = 0._r8
                  leafn_xfer_to_leafn(p)   = 0._r8
                  frootn_xfer_to_frootn(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0._r8
                     deadstemc_xfer_to_deadstemc(p)   = 0._r8
                     livecrootc_xfer_to_livecrootc(p) = 0._r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0._r8
                     livestemn_xfer_to_livestemn(p)   = 0._r8
                     deadstemn_xfer_to_deadstemn(p)   = 0._r8
                     livecrootn_xfer_to_livecrootn(p) = 0._r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0._r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0._r8
                  leafn_xfer(p) = 0._r8
                  frootc_xfer(p) = 0._r8
                  frootn_xfer(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0._r8
                     livestemn_xfer(p) = 0._r8
                     deadstemc_xfer(p) = 0._r8
                     deadstemn_xfer(p) = 0._r8
                     livecrootc_xfer(p) = 0._r8
                     livecrootn_xfer(p) = 0._r8
                     deadcrootc_xfer(p) = 0._r8
                     deadcrootn_xfer(p) = 0._r8
                  end if
               end if
            end if

            ! test for switching from dormant period to growth period
            if (dormant_flag(p) == 1._r8) then

               ! keep track of the number of freezing degree days in this
               ! dormancy period (only if the freeze flag has not previously been set
               ! for this dormancy period

               if (onset_gddflag(p) == 0._r8 .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday

               ! if the number of freezing degree days exceeds a critical value,
               ! then onset will require both wet soils and a critical soil
               ! temperature sum.  If this case is triggered, reset any previously
               ! accumulated value in onset_swi, so that onset now depends on
               ! the accumulated soil water index following the freeze trigger

               if (onset_fdd(p) > crit_onset_fdd) then
                  onset_gddflag(p) = 1._r8
                  onset_fdd(p) = 0._r8
                  onset_swi(p) = 0._r8
               end if

               ! if the freeze flag is set, and if the soil is above freezing
               ! then accumulate growing degree days for onset trigger

               if (onset_gddflag(p) == 1._r8 .and. soilt > SHR_CONST_TKFRZ) then
                  onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
               end if

              ! if soils are wet, accumulate soil water index for onset trigger
               additional_onset_condition = .true.
               if(CNParamsShareInst%constrain_stress_deciduous_onset) then
                  ! if additional constraint condition not met,  set to false
                  if ((prec10(p) * (3600.0_r8*10.0_r8*24.0_r8)) < rain_threshold) then
                     additional_onset_condition = .false.
                  endif
               endif

               if (psi >= soilpsi_on) then
                  onset_swi(p) = onset_swi(p) + fracday
               endif

               ! if critical soil water index is exceeded, set onset_flag, and
               ! then test for soil temperature criteria

               ! Adding in Kyla's rainfall trigger when fun on. RF. prec10 (mm/s) needs to be higher than 8mm over 10 days. 

               if (onset_swi(p) > crit_onset_swi.and. additional_onset_condition)  then
                  onset_flag(p) = 1._r8
              
                  ! only check soil temperature criteria if freeze flag set since
                  ! beginning of last dormancy.  If freeze flag set and growing
                  ! degree day sum (since freeze trigger) is lower than critical
                  ! value, then override the onset_flag set from soil water.

                  if (onset_gddflag(p) == 1._r8 .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0._r8
               end if

               ! only allow onset if dayl > 6hrs
               if (onset_flag(p) == 1._r8 .and. dayl(g) <= secspqtrday) then
                  onset_flag(p) = 0._r8
               end if

               ! if this is the beginning of the onset period
               ! then reset the phenology flags and indices

               if (onset_flag(p) == 1._r8) then
                  dormant_flag(p) = 0._r8
                  days_active(p) = 0._r8
                  onset_gddflag(p) = 0._r8
                  onset_fdd(p) = 0._r8
                  onset_gdd(p) = 0._r8
                  onset_swi(p) = 0._r8
                  onset_counter(p) = ndays_on * secspday

                  ! call subroutine to move all the storage pools into transfer pools,
                  ! where they will be transfered to displayed growth over the onset period.
                  ! this code was originally handled with call cn_storage_to_xfer(p)
                  ! inlined during vectorization

                  ! set carbon fluxes for shifting storage pools to transfer pools
                  leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
                  frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                     deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                     livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                     deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                     gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
                  end if

                  ! set nitrogen fluxes for shifting storage pools to transfer pools
                  leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
                  frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                     deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                     livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                     deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
                  end if
               end if

               ! test for switching from growth period to offset period
            else if (offset_flag(p) == 0._r8) then

               ! if soil water potential lower than critical value, accumulate
               ! as stress in offset soil water index

               if (psi <= soilpsi_off) then
                  offset_swi(p) = offset_swi(p) + fracday

                  ! if the offset soil water index exceeds critical value, and
                  ! if this is not the middle of a previously initiated onset period,
                  ! then set flag to start the offset period and reset index variables

                  if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8

                  ! if soil water potential higher than critical value, reduce the
                  ! offset water stress index.  By this mechanism, there must be a
                  ! sustained period of water stress to initiate offset.

               else if (psi >= soilpsi_on) then
                  offset_swi(p) = offset_swi(p) - fracday
                  offset_swi(p) = max(offset_swi(p),0._r8)
               end if

               ! decrease freezing day accumulator for warm soil
               if (offset_fdd(p) > 0._r8 .and. soilt > SHR_CONST_TKFRZ) then
                  offset_fdd(p) = offset_fdd(p) - fracday
                  offset_fdd(p) = max(0._r8, offset_fdd(p))
               end if

               ! increase freezing day accumulator for cold soil
               if (soilt <= SHR_CONST_TKFRZ) then
                  offset_fdd(p) = offset_fdd(p) + fracday

                  ! if freezing degree day sum is greater than critical value, initiate offset
                  if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
               end if

               ! force offset if daylength is < 6 hrs
               if (dayl(g) <= secspqtrday) then
                  offset_flag(p) = 1._r8
               end if

               ! if this is the beginning of the offset period
               ! then reset flags and indices
               if (offset_flag(p) == 1._r8) then
                  offset_fdd(p) = 0._r8
                  offset_swi(p) = 0._r8
                  offset_counter(p) = ndays_off * secspday
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! keep track of number of days since last dormancy for control on
            ! fraction of new growth to send to storage for next growing season

            if (dormant_flag(p) == 0.0_r8) then
               days_active(p) = days_active(p) + fracday
            end if

            ! calculate long growing season factor (lgsf)
            ! only begin to calculate a lgsf greater than 0.0 once the number
            ! of days active exceeds days/year.
            lgsf(p) = max(min(3.0_r8*(days_active(p)-leaf_long(ivt(p))*dayspyr )/dayspyr, 1._r8),0._r8)
            ! RosieF. 5 Nov 2015.  Changed this such that the increase in leaf turnover is faster after
            ! trees enter the 'fake evergreen' state. Otherwise, they have a whole year of 
            ! cheating, with less litterfall than they should have, resulting in very high LAI. 
            ! Further, the 'fake evergreen' state (where lgsf>0) is entered at the end of a single leaf lifespan
            ! and not a whole year. The '3' is arbitrary, given that this entire system is quite abstract. 


            ! set background litterfall rate, when not in the phenological offset period
            if (offset_flag(p) == 1._r8) then
               bglfr(p) = 0._r8
            else
               ! calculate the background litterfall rate (bglfr)
               ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

               bglfr(p) = (1._r8/(leaf_long(ivt(p))*dayspyr*secspday))*lgsf(p)
            end if

            ! set background transfer rate when active but not in the phenological onset period
            if (onset_flag(p) == 1._r8) then
               bgtr(p) = 0._r8
            else
               ! the background transfer rate is calculated as the rate that would result
               ! in complete turnover of the storage pools in one year at steady state,
               ! once lgsf has reached 1.0 (after 730 days active).

               bgtr(p) = (1._r8/(dayspyr*secspday))*lgsf(p)

               ! set carbon fluxes for shifting storage pools to transfer pools

               ! reduced the amount of stored carbon flowing to display pool by only counting the delta
               ! between leafc and leafc_store in the flux. RosieF, Nov5 2015. 
               leafc_storage_to_xfer(p)  = max(0.0_r8,(leafc_storage(p)-leafc(p))) * bgtr(p)
               frootc_storage_to_xfer(p) = max(0.0_r8,(frootc_storage(p)-frootc(p))) * bgtr(p)
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
                  deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
                  livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
                  deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
                  gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
               frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
                  deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
                  livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
                  deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
               end if
            end if

         end if ! end if stress deciduous

      end do ! end of patch loop

    end associate

  end subroutine CNStressDecidPhenology

  !-----------------------------------------------------------------------
  subroutine PalmPhenology(num_pcropp, filter_pcropp                     , &
       waterdiagnosticbulk_inst, temperature_inst, crop_inst, canopystate_inst, cnveg_state_inst , &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)

    ! !DESCRIPTION:
    ! !New code specifically for the phenology of palms which have a sub-canopy phytomer structure
    ! !Developed by Y. Fan (CLM-Palm, Y.Fan et al. (2015), doi:10.5194/gmd-8-3785-2015)

    ! !USES:
    use clm_time_manager , only : get_curr_date, get_curr_calday, get_days_per_year, get_rad_step_size
    use pftconMod        , only : noilpalm, nirrig_oilpalm, mxnp    ! !,              &
    ! !                                  gddmin, lfemerg, mxgdd, minplanttemp, planttemp, &
    ! !                                  lfmat, lfsen, lfexp, grnfill, grnmat, mxmat, &
    ! !                                  phyllochron, mxlivenp, laimx, transplant
    use clm_varcon       , only : spval, secspday
    use clm_varctl       , only : use_fertilizer
    use clm_varctl       , only : use_c13, use_c14
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                        , intent(in)    :: filter_pcropp(:) ! filter for prognostic crop patches
    type(waterdiagnosticbulk_type)          , intent(in)    :: waterdiagnosticbulk_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    !
    ! LOCAL VARAIBLES:
    integer kyr       ! current year
    integer kmo       ! month of year  (1, ..., 12)
    integer kda       ! day of month   (1, ..., 31)
    integer mcsec     ! seconds of day (0, ..., seconds/day)
    integer jday      ! julian day of the year
    integer fp,p      ! patch indices
    integer c         ! column indices
    integer g         ! gridcell indices
    integer h         ! hemisphere indices
    real(r8) :: dtrad ! radiation time step delta t (seconds)
    real(r8) crmcorn  ! comparitive relative maturity for corn
    real(r8) ndays_on ! number of days to fertilize
    integer n,m,i     ! phytomer indices
    !integer nyrs      ! number of years prognostic crop has run (now is nyrs_crop_active calculated in Crop_Type)
    real(r8) dayspyr  ! days per year
    integer phyllday  ! average days per phyllochron
    integer buddays   ! days for leaf storage (bud) growth period
    real(r8) gddperday ! average GDD per day (base 15 degree)
    integer np0, np1, np2 !temporary phytomer number
    !integer update_rank(1:mxnp) !phytomer rank update

    !------------------------------------------------------------------------

    associate(                                                                   &
         ivt               =>    patch%itype                                     , & ! Input:  [integer  (:) ]  patch vegetation type
         phytomer          =>    pftcon%phytomer                             , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
         perennial         =>    pftcon%perennial                            , & ! Input:  [integer (:)]  binary flag for perennial crop phenology (1=perennial, 0=not perennial) (added by Y.Fan)
         mxlivenp          =>    pftcon%mxlivenp                             , & ! Input:  [integer (:)]   maximum number of live phytomer
         transplant        =>    pftcon%transplant                           , & ! Input:  [integer (:)] 
         pprodharv10       =>    pftcon%pprodharv10                          , & ! Input:  harvest mortality proportion of deadstem to 10-yr pool
         deadwdcn          =>    pftcon%deadwdcn                             , & ! Input:  [real(r8) (:)] live wood C:N (gC/gN)
         slatop            =>    pftcon%slatop                               , & ! Input:  [real(r8) (:)] specific leaf area at top of canopy, projected area basis [m^2/gC]
         leaf_long         =>    pftcon%leaf_long                              , & ! Input:  leaf longevity (yrs)
         leafcn            =>    pftcon%leafcn                                 , & ! Input:  leaf C:N (gC/gN)
         manunitro         =>    pftcon%manunitro                              , & ! Input:  max manure to be applied in total (kgN/m2)
         mxmat             =>    pftcon%mxmat                                  , & ! Input:
         minplanttemp      =>    pftcon%minplanttemp                           , & ! Input:
         planttemp         =>    pftcon%planttemp                              , & ! Input:
         gddmin            =>    pftcon%gddmin                                 , & ! Input:
         mxgdd             =>    pftcon%mxgdd                                  , & ! Input:
         hybgdd            =>    pftcon%hybgdd                                 , & ! Input:
         phyllochron       =>    pftcon%phyllochron                            , & ! Input:
         lfemerg           =>    pftcon%lfemerg                                , & ! Input:
         lfexp             =>    pftcon%lfexp                                 , & ! Input:
         lfsen             =>    pftcon%lfsen                                 , & ! Input:
         lfmat             =>    pftcon%lfmat                                 , & ! Input:
         grnmat            =>    pftcon%grnmat                                , & ! Input:
         grnfill           =>    pftcon%grnfill                               , & ! Input:
         laimx             =>    pftcon%laimx                                  , & ! Input:
         t_ref2m_min       =>    temperature_inst%t_ref2m_min_patch            , & ! Input:  [real(r8) (:) ]  daily minimum of average 2 m height surface air temperature (K)
         t10               =>    temperature_inst%t_a10_patch                  , & ! Input:  [real(r8) (:) ]  10-day running mean of the 2 m temperature (K)
         a5tmin            =>    temperature_inst%t_a5min_patch                , & ! Input:  [real(r8) (:) ]  5-day running mean of min 2-m temperature
         a10tmin           =>    temperature_inst%t_a10min_patch               , & ! Input:  [real(r8) (:) ]  10-day running mean of min 2-m temperature
         gdd020            =>    temperature_inst%gdd020_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd0
         gdd820            =>    temperature_inst%gdd820_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd8
         gdd1020           =>    temperature_inst%gdd1020_patch                , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd10
         gdd1520           =>    temperature_inst%gdd1520_patch                , & ! Input:  [real(r8) (:)]  20 yr mean of gdd15
         gdd15             =>    temperature_inst%gdd15_patch                  , & ! Input:  [real(r8) (:)]  growing deg. days base 15 deg C (ddays) (Y.Fan)

         idpp              =>    cnveg_state_inst%idpp_patch                             , & ! Output: [integer (:)]  days past planting
         yrop              =>    cnveg_state_inst%yrop_patch                             , & ! Output: [integer (:)]  year of planting (added Y.Fan)
         lfoutday          =>    crop_inst%lfoutday_patch                                , & ! Output: [integer (:,:)]  date of leaf/phytomer emergence (Y.Fan)
         lfoutyr           =>    crop_inst%lfoutyr_patch                                 , & ! Output: [integer (:,:)]  year of leaf/phytomer emergence (Y.Fan)
         lfdays            =>    crop_inst%lfdays_patch                                  , & ! Output: [integer (:,:)]  days past leaf emergence for each phytomer (Y.Fan)
      !  harvyr            =>    crop_inst%harvyr_patch                                  , & ! Output: [integer (:)]  harvest year (added Y.Fan)
         harvdate          =>    crop_inst%harvdate_patch                                , & ! Output: [integer (:)]  harvest date
         huileafnp         =>    crop_inst%huileafnp_patch                               , & ! Output: [real(r8) (:,:)]  hui needed for initiation of successive phytomers
         huilfexpnp        =>    crop_inst%huilfexpnp_patch                              , & ! Output: [real(r8) (:,:)]  hui needed for leaf expansion of successive phytomers
         huilfmatnp        =>    crop_inst%huilfmatnp_patch                              , & ! Output: [real(r8) (:,:)]  hui needed for leaf maturity of successive phytomers^M
         huilfsennp        =>    crop_inst%huilfsennp_patch                              , & ! Output: [real(r8) (:,:)]  hui needed for leaf senescence of successive phytomers^M
         huilfendnp        =>    crop_inst%huilfendnp_patch                              , & ! Output: [real(r8) (:,:)]  hui needed for end of life of successive phytomers
         huigrnnp          =>    crop_inst%huigrnnp_patch                                , & ! Output: [real(r8) (:,:)]  hui needed for starting grainfill of successive phytomers
         grnmatnp          =>    crop_inst%grnmatnp_patch                                , & ! Output: [real(r8) (:,:)]  hui needed for grain maturity of successive phytomers
         np                =>    crop_inst%np_patch                                      , & ! Output: [integer (:)]   total number of phytomers having appeared so far
         rankp             =>    crop_inst%rankp_patch                                   , & ! Output: [integer (:,:)]  rank of phytomers from 1=youngest to np=oldest and 0=dead
         livep             =>    crop_inst%livep_patch                                   , & ! Output: [real(r8) (:,:)]  Flag, 1 if this phytomer is alive
 
         fertnitro         =>    crop_inst%fertnitro_patch                     , & ! Input:  [real(r8) (:) ]  fertilizer nitrogen
         hui               =>    crop_inst%gddplant_patch                      , & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)
         leafout           =>    crop_inst%gddtsoi_patch                       , & ! Input:  [real(r8) (:) ]  gdd from top soil layer temperature
         croplive          =>    crop_inst%croplive_patch                      , & ! Output: [logical  (:) ]  Flag, true if planted, not harvested
         cropplant         =>    crop_inst%cropplant_patch                     , & ! Output: [logical  (:) ]  Flag, true if crop may be planted
         vf                =>    crop_inst%vf_patch                            , & ! Output: [real(r8) (:) ]  vernalization factor
         peaklai           =>    cnveg_state_inst%peaklai_patch                , & ! Output: [integer  (:) ] 1: max allowed lai; 0: not at max
         tlai              =>    canopystate_inst%tlai_patch                   , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow

         idop              =>    cnveg_state_inst%idop_patch                   , & ! Output: [integer  (:) ]  date of planting
         gddmaturity       =>    cnveg_state_inst%gddmaturity_patch            , & ! Output: [real(r8) (:) ]  gdd needed to harvest
         huileaf           =>    cnveg_state_inst%huileaf_patch                , & ! Output: [real(r8) (:) ]  heat unit index needed from planting to leaf emergence
         huigrain          =>    cnveg_state_inst%huigrain_patch               , & ! Output: [real(r8) (:) ]  same to reach vegetative maturity
         cumvd             =>    cnveg_state_inst%cumvd_patch                  , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?
         hdidx             =>    cnveg_state_inst%hdidx_patch                  , & ! Output: [real(r8) (:) ]  cold hardening index?
         bglfr             =>    cnveg_state_inst%bglfr_patch                  , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)
         bgtr              =>    cnveg_state_inst%bgtr_patch                   , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)
         lgsf              =>    cnveg_state_inst%lgsf_patch                   , & ! Output: [real(r8) (:) ]  long growing season factor [0-1]
         onset_flag        =>    cnveg_state_inst%onset_flag_patch             , & ! Output: [real(r8) (:) ]  onset flag
         offset_flag       =>    cnveg_state_inst%offset_flag_patch            , & ! Output: [real(r8) (:) ]  offset flag
         onset_counter     =>    cnveg_state_inst%onset_counter_patch          , & ! Output: [real(r8) (:) ]  onset counter
         offset_counter    =>    cnveg_state_inst%offset_counter_patch         , & ! Output: [real(r8) (:) ]  offset counter

         leafc_xfer        =>    cnveg_carbonstate_inst%leafc_xfer_patch       , & ! Output: [real(r8) (:) ]  (gC/m2)   leaf C transfer

         crop_seedc_to_leaf =>   cnveg_carbonflux_inst%crop_seedc_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gC/m2/s) seed source to leaf

         fert_counter      =>    cnveg_nitrogenflux_inst%fert_counter_patch    , & ! Output: [real(r8) (:) ]  >0 fertilize; <=0 not (seconds)
         leafn_xfer        =>    cnveg_nitrogenstate_inst%leafn_xfer_patch     , & ! Output: [real(r8) (:) ]  (gN/m2)   leaf N transfer
         crop_seedn_to_leaf =>   cnveg_nitrogenflux_inst%crop_seedn_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gN/m2/s) seed source to leaf
         cphase            =>    crop_inst%cphase_patch                        , & ! Output: [real(r8) (:)]   phenology phase
         fert              =>    cnveg_nitrogenflux_inst%fert_patch            , & ! Output: [real(r8) (:) ]  (gN/m2/s) fertilizer applied each timestep
         bglfr_p                 =>    crop_inst%bglfr_p_patch                                , & ! Output: [real(r8) (:,:)]  background litterfall rate for each phytomer (Y.Fan)
         bgtr_p                  =>    crop_inst%bgtr_p_patch                                 , & ! Output: [real(r8) (:,:)]  background transfer rate for each phytomer (Y.Fan)
         harvest_flag            =>    crop_inst%harvest_flag_patch                           , & ! Output: [real(r8) (:)]  harvest flag (added by Y.Fan)
         harvest_counter         =>    crop_inst%harvest_counter_patch                        , & ! Output: [real(r8) (:)]  harvest counter to tag phytomer
         prune                   =>    crop_inst%prune_patch                                  , & ! Output: [real(r8) (:)]  leaf pruning flag
        ! pleafc                  =>    cnveg_carbonstate_inst%pleafc_patch                           , & ! InOut:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C
         pgrainc                 =>    cnveg_carbonstate_inst%pgrainc_patch                          , & ! Input:  [real(r8) (:,:)] (gC/m2) grain C on each phytomer
         deadstemc_xfer          =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                   , & ! Output:  [real(r8) (:)] (gC/m2) live stem C transfer
         pleafc_storage          =>    cnveg_carbonstate_inst%pleafc_storage_patch                   , & ! Input:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C transfer
         pleafc_storage_to_xfer  =>    cnveg_carbonflux_inst%pleafc_storage_to_xfer_patch            , & ! InOut:  [real(r8) (:,:)]^M
         leafc_storage_to_xfer   =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch             , & ! InOut:  [real(r8) (:)]
         pleafc_xfer             =>    cnveg_carbonstate_inst%pleafc_xfer_patch                      , & ! Output: [real(r8) (:,:)]  (gC/m2) phytomer leaf C transfer
         pleafn_xfer             =>    cnveg_nitrogenstate_inst%pleafn_xfer_patch                    , & ! Output: [real(r8) (:,:)]  (gN/m2) phytomer leaf N transfer
         deadstemn_xfer          =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                 , & ! Output: [real(r8) (:)]  (gN/m2) live stem N transfer
         pleafn_storage          =>    cnveg_nitrogenstate_inst%pleafn_storage_patch                 , & ! Input:  [real(r8) (:,:)]  (gN/m2) phytomer leaf N transfer
         pleafn_storage_to_xfer  =>    cnveg_nitrogenflux_inst%pleafn_storage_to_xfer_patch          , & ! InOut:  [real(r8) (:,:)]^M
         leafn_storage_to_xfer   =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch           , & ! InOut:  [real(r8) (:)]
         crop_seedc_to_deadstem   =>   cnveg_carbonflux_inst%crop_seedc_to_deadstem_patch      , & ! Output: [real(r8) (:)]  (gC/m2/s) CN dynamic landcover fluxes: seed source to PFT-level; initialized as 0
         crop_seedn_to_deadstem   =>   cnveg_nitrogenflux_inst%crop_seedn_to_deadstem_patch    & ! Output: [real(r8) (:)]  (gN/m2/s) seed source to PFT-level
   )


   ! get time info
   dayspyr = get_days_per_year()
   jday    = get_curr_calday()
   call get_curr_date(kyr, kmo, kda, mcsec)
   dtrad   = real( get_rad_step_size(), r8 )
  
   if (use_fertilizer) then
    ndays_on = 20._r8 ! number of days to fertilize
   else
    ndays_on = 0._r8 ! number of days to fertilize
   end if
  
   do fp = 1, num_pcropp
      p = filter_pcropp(fp)
      c = patch%column(p)
      g = patch%gridcell(p)
      h = inhemi(p)
  
      if ( phytomer(ivt(p)) > 0) then
  
         ! background litterfall and transfer rates
         bglfr(p) = 0._r8 ! this value changes later in a crop's life cycle
         bgtr(p)  = 0._r8
  
         ! in order to allow a crop to be planted only once each year
         ! initialize cropplant = .false., but hold it = .true. through the end of the year
         if ( jday == jdayyrstart(h) .and. mcsec == 0 )then
             if (.not. croplive(p))  then
                cropplant(p) = .false.
                idop(p)      = NOT_Planted
             end if
         end if
  
         ! Phase 1: Planting to leaf emergence
  
         if ( (.not. croplive(p)) .and. (.not. cropplant(p)) ) then
  
             if ((t10(p) /= spval .and. a10tmin(p) /= spval .and. &
                 t10(p)     > planttemp(ivt(p))             .and. &
                 a10tmin(p) > minplanttemp(ivt(p))          .and. &
                 jday       >= minplantjday(ivt(p),h)       .and. &
                 jday       <= maxplantjday(ivt(p),h)       .and. &
                 gdd15(p) /= spval .and. gdd15(p) >= gddmin(ivt(p))) &
                 .or.  &
                (jday == maxplantjday(ivt(p),h) .and. &
                gdd15(p) /= spval .and. gdd15(p) > 0._r8 )    ) then
                !impose limit on growing season length and gddmin needed or hit the max planting julian day -- go ahead and plant
                !gdd1520(p) is used for long-term crop spinup, use current year gdd15 for this study
  
                croplive(p)  = .true.
                cropplant(p) = .true.
                idop(p)      = jday
                yrop(p)      = kyr
                harvdate(p)  = NOT_Harvested
  
                ! At planting, each crop is assigned 1g leaf C/m2 to be transferred to the leaves upon leaf emergence.
                ! Similar as evergreen PFTs which are initialized with leafc = 3 gC/m2
                ! It could be set as a higher value for crops transplanted from nursery (Y.Fan)
                ! oil palm transplanted when seedlings between 10 and 20 months (1m tall) after the germinated seed stage.
                if (transplant(ivt(p)) > 0._r8) then  ! transplant means initial LAI of seedling
                    leafc_xfer(p)  = transplant(ivt(p))/slatop(ivt(p))
                else
                    leafc_xfer(p)  = initial_seed_at_planting  ! now the initial seed is a variable 
                end if
                leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! An equivalent amount of seed leaf N is assigned with onset
                crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                !for woody crops like palm, assign 10% seeding C to deadstem at transplanting
                if (pftcon%woody(ivt(p)) == 1._r8) then 
                   deadstemc_xfer(p) = 0.1_r8 * leafc_xfer(p)
                   deadstemn_xfer(p) = deadstemc_xfer(p) / deadwdcn(ivt(p))
                   crop_seedc_to_deadstem(p) = deadstemc_xfer(p)/dt
                   crop_seedn_to_deadstem(p) = deadstemn_xfer(p)/dt
                end if
  
             end if
  
             onset_counter(p) = 0._r8 !onset_counter initialized to zero when .not. croplive
             offset_counter(p) = 0._r8 !offset_counter initialized to zero when .not. croplive
             harvest_counter(p) = 0._r8
             lfoutday(p,:) = NOT_Emerged ! all phytomers initialized to NOT_Emerged when .not. croplive
             lfoutyr(p,:) = NOT_Emerged
             lfdays(p,:) = 0
             livep(p,:) = 0._r8       ! all phytomers initialized to zero (not alive) when .not. croplive
             rankp(p,:) = 0
             np(p) = 0
         end if ! palm not live nor planted
  
         ! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline or stabilize)
         ! Harvest: export fruits and prune old leaves
  
         onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
         offset_flag(p) = 0._r8 ! carbon and nitrogen transfers
         harvest_flag(p) = 0._r8
  
         ! calculate threshold for attaining leaf emergence(initiation)
         ! initiation of the first phytomer, assume immediately after transplanting if transplant(ivt(p)) > 0
         ! lfemerg has a missing value set to 0; zeros will be read as NAs; 
         ! to avoid conflict, here lfemerg is set to 0.001 for transplanting
         ! (Y.Fan 2022.08).
         huileaf(p) = lfemerg(ivt(p))
    
         ! pre-expansion phase: from initiation to leaf expansion
         huilfexp(p) = huileaf(p) + lfexp(ivt(p))
    
         ! post-expansion phases:
         ! from leaf expansion to leaf maturity
         huilfmat(p) = huilfexp(p) + lfmat(ivt(p))
         ! from leaf expansion to start of leaf senescence
         huilfsen(p) = huilfexp(p) + lfsen(ivt(p))
         ! from leaf expansion to beginning of fruit-fill
         huigrain(p) = huilfexp(p) + grnfill(ivt(p))
         ! from leaf expansion to fruit maturity and harvest
         gddmaturity(p)= huilfexp(p) + grnmat(ivt(p))
         ! from leaf expansion to end of leaf senescence (die)
         huilfend(p) = huilfexp(p) + mxgdd(ivt(p))
    
         if (croplive(p)) then
             cphase(p) = 1._r8

             ! days past planting may determine harvest
             ! new counter for perennial crops (Y.Fan)
             if (perennial(ivt(p)) == 1) then
                idpp(p) = int(dayspyr)*(kyr-yrop(p)) + jday - idop(p)
             else  !only for annual or biannual crops
                if (jday >= idop(p)) then
                   idpp(p) = jday - idop(p)
                else
                   idpp(p) = int(dayspyr) + jday - idop(p)
                end if
             end if

             !gddperday = max(10._r8, max(gdd15(p)/jday, gdd1520(p)/dayspyr)) !do not use gdd1520 for long term crop spinup^M
             gddperday = max(10._r8, gdd15(p)/jday)
             phyllday = int(phyllochron(ivt(p))/gddperday)

             ! =======================
             ! enter onset for one time step:
             ! transfer seed carbon to leaf emergence; the amount of transfer (1g leaf C) already allocated at planting

             ! onset_counter initialized to zero when .not. croplive
             ! offset_counter relevant only at time step of final rotation
             onset_counter(p) = onset_counter(p) - dt
             if (leafout(p) >= huileaf(p) .and. hui(p) < huilfmat(p)) then
               cphase(p) = 2._r8
               if (abs(onset_counter(p)) > 1.e-6_r8) then
                  onset_flag(p)    = 1._r8
                  onset_counter(p) = dt

                  !for transplanting from nursery: the number of initial phytomers depends on initial C amount
                  if (transplant(ivt(p)) > 0._r8) then

                     np(p) = floor(transplant(ivt(p))/(0.1_r8*laimx(ivt(p))/2)) !assume 10% of laimx for a mature leaf of seedling^M
                     np(p) = min(mxlivenp(ivt(p)), np(p))
                     rankp(p,1:np(p)) = (/np(p):1:-1/) !already expanded phytomers

                                 !leaf expansion date for already existing leaves at transplanting^M
                                 !assign the current date to the youngest expanded leaf^M
                     lfoutyr(p,:) = kyr
                     lfoutday(p,1:np(p)) = (/ ((jday-(i-1)*phyllday), i=np(p), 1, -1) /)
                     where (lfoutday(p,:) < 0)
                         lfoutyr(p,:) = kyr - (int(abs(lfoutday(p,:))/dayspyr) + 1)
                                       lfoutday(p,:) = lfoutday(p,:) + dayspyr*(int(abs(lfoutday(p,:))/dayspyr) + 1)
                     endwhere

                     !assign initial C to already expanded leaves at transplanting
                     np1 = int((huilfmat(p)-huilfexp(p))/phyllochron(ivt(p))) !number of expanded immature phytomers
                     if (np(p) > np1) then
                         !pleafc_xfer(p,1:(np(p)-1)) = leafc_xfer(p)/np(p)
                         !pleafc_xfer(p,np(p)) = leafc_xfer(p) - sum(pleafc_xfer(p,1:(np(p)-1)))
                         pleafc_xfer(p,(np(p)-np1+1):np(p)) = 0.1_r8*laimx(ivt(p))/slatop(ivt(p))*(/np1:1:-1/)/(np1+1)  !linear decrease in leaf size till the youngest phytomer
                         pleafc_xfer(p,1:(np(p)-np1)) = 0.1_r8*laimx(ivt(p))/slatop(ivt(p))*(/1:(np(p)-np1):1/)/(np(p)-np1+1) 
                         !the oldest phytomer of a seedling is smaller than the middle ones due to limited resource at the beginning
                         pleafc_xfer(p,1:np(p)) = pleafc_xfer(p,1:np(p)) + (leafc_xfer(p) - sum(pleafc_xfer(p,:)))/np(p) !distribute the remaining leaf C
                     else if (np(p) > 1) then
                         pleafc_xfer(p,1:np(p)) = 0.1_r8*laimx(ivt(p))/slatop(ivt(p))*(/np(p):1:-1/)/(np(p)+1)  
                         !np plus 1 results in exactly a half of np(p) when summing the ratio series
                         pleafc_xfer(p,1:np(p)) = pleafc_xfer(p,1:np(p)) + (leafc_xfer(p) - sum(pleafc_xfer(p,:)))/np(p)
                     else
                         pleafc_xfer(p,1) = leafc_xfer(p)
                     end if
                     pleafn_xfer(p,:) = pleafc_xfer(p,:) / leafcn(ivt(p))

                     !add non-expanded phytomers but already initiated before transplanting
                     !do not assign initial C to them
                     np0 = int((huilfexp(p)-huileaf(p))/phyllochron(ivt(p)))
                     if (np0 >= 1) then
                         rankp(p,(np(p)+1):(np(p)+np0)) = (/-1:(-np0):-1/)       !unexpanded phytomers have negative ranks
                         np(p) = np(p) + np0
                     end if
                     livep(p,1:np(p)) = 1._r8 
                     !estimate the expansion/maturity/harvest time for each existing phytomer
                     !assume the youngest phytomer is initiated at planting (not expanded yet)
                     huileafnp(p,1:np(p)) = (/ ((huileaf(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     huilfexpnp(p,1:np(p)) = (/ ((huilfexp(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     huilfmatnp(p,1:np(p)) = (/ ((huilfmat(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     huigrnnp(p,1:np(p)) = (/ ((huigrain(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     grnmatnp(p,1:np(p)) = (/ ((gddmaturity(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     huilfsennp(p,1:np(p)) = (/ ((huilfsen(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /)
                     huilfendnp(p,1:np(p)) = (/ ((huilfend(p)-(i-1)*phyllochron(ivt(p))), i=np(p), 1, -1) /) 
                  else !from seed
                     np(p) = 1  !first phytomer initiated after seed planting
                     livep(p,1) = 1._r8
                     rankp(p,1) = -1
                     pleafc_xfer(p,1) = leafc_xfer(p)
                     pleafn_xfer(p,1) = pleafc_xfer(p,1) / leafcn(ivt(p))
                     huileafnp(p,1) = huileaf(p)
                     huilfexpnp(p,1) = huilfexp(p)
                     huilfmatnp(p,1) = huilfmat(p)
                     huigrnnp(p,1) = huigrain(p)
                     grnmatnp(p,1) = gddmaturity(p)
                     huilfsennp(p,1) = huilfsen(p)
                     huilfendnp(p,1) = huilfend(p)
                  end if
               else
                  ! this ensures no re-entry to onset of phase2
                  ! b/c onset_counter(p) = onset_counter(p) - dt
                  ! at every time step
                  onset_counter(p) = dt
               end if
             end if !end if leaf emergence to leaf maturity period

             ! determine initiation and maturity thresholds of successive phytomers
             ! for oil palm, phyllochron increases (max 1.5 times) through age until 10 years old
             if (np(p) > 0) then
                phyllochron2(p) = phyllochron(ivt(p)) * &
                               min(1.5_r8, 1._r8 + (1.5_r8 - 1._r8) * &
                               min(1._r8, real(idpp(p))/(10._r8*dayspyr)))
                huileafnp(p, np(p)+1) = huileafnp(p, np(p)) + phyllochron2(p)

             else    ! plant never emerged from the ground
                croplive(p) = .false.
                crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                leafc_xfer(p) = 0._r8  ! revert planting transfers
                leafn_xfer(p) = 0._r8
                crop_seedc_to_deadstem(p) = crop_seedc_to_deadstem(p) - deadstemc_xfer(p)/dt
                crop_seedn_to_deadstem(p) = crop_seedn_to_deadstem(p) - deadstemn_xfer(p)/dt
                deadstemc_xfer(p) = 0._r8  ! revert planting transfers
                deadstemn_xfer(p) = 0._r8
             end if

             !turn on successive phytomers after each leaf initiation, later turn off by ageing or pruning
             if (np(p) > 0 .and. hui(p) >= huileafnp(p, np(p)+1)) then
                np(p) = np(p) + 1
                livep(p,np(p)) = 1._r8
                rankp(p,np(p)) = minval(rankp(p,:)) - 1
                huilfexpnp(p,np(p)) = huileafnp(p,np(p)) + lfexp(ivt(p))  !pre-expansion
                huilfmatnp(p,np(p)) = huilfexpnp(p,np(p)) + lfmat(ivt(p)) !post-expansion
                huigrnnp(p,np(p)) = huilfexpnp(p,np(p)) + grnfill(ivt(p))
                grnmatnp(p,np(p)) = huilfexpnp(p,np(p)) + grnmat(ivt(p))
                huilfsennp(p,np(p)) = huilfexpnp(p,np(p)) + lfsen(ivt(p))
                huilfendnp(p,np(p)) = huilfexpnp(p,np(p)) + mxgdd(ivt(p))
             end if

             ! =======================
             ! enter phase 2 from leaf initiation to leaf expansion:
             ! switch from storage growth (bud & "spear" leaf stage) to photosynthetic active LAI growth
             np2 = sum(maxloc(rankp(p,:), mask=rankp(p,:) < 0)) !record the index of current largest unexpanded phytomer (rankp=-1)
             if (livep(p,np2) == 1._r8 .and. rankp(p,np2) < 0 .and. &
                hui(p) >= huilfexpnp(p,np2)) then
                !the largest "spear" leaf expands and its rank changes from -1 to +1
                rankp(p,np2) = 1
                !record leaf out date (only once)
                lfoutday(p,np2) = jday
                lfoutyr(p,np2) = kyr

                !update the rank of other live phytomers
                where (rankp(p,1:(np2-1)) > 0) rankp(p,1:(np2-1)) = rankp(p,1:(np2-1)) + 1 !expanded
                rankp(p,(np2+1):np(p)) = rankp(p,(np2+1):np(p)) + 1  !unexpanded
             end if

             !count days past each leaf/phytomer expansion
             where (rankp(p,:) > 0 .and. lfoutday(p,:) /= NOT_Emerged) &
                lfdays(p,:) = int(dayspyr)*(kyr-lfoutyr(p,:)) + jday - lfoutday(p,:)


             ! =======================
             ! enter phase 2 for each phytomer: from leaf initiation to leaf expansion
             ! switch from storage growth to photosynthetic active LAI growth
             ! the phytomer storage pool represents displayed but non-photosynthetic leafc
             ! this is very different from that of deciduous PFTs
             if (idpp(p) < mxmat(ivt(p))) then !from planting to crop rotation
                 where (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:))
                     ! the background transfer rate gives a complete turnover of the storage pool
                     ! to displayed pool from leaf expansion until maturity.
                     ! bgtr flux merges with current time growth allocation flux
                     !bgtr_p(p,:) = 2._r8/((huilfmat(p) - huilfexp(p)) / gddperday * secspday)
                     !above equation gives nearly ~90% turnover in nonlinear decreasing rate

                     bgtr_p(p,:) = 1.0_r8 /((huilfmatnp(p,:) - hui(p)) / gddperday * secspday) !100% turnover
                     !if multiplier > 2.0_r8, it gives steeper nonlinear decreasing rate; 1.0_r8 gives a linear rate
                     !not exceed the max rate 1/dt
                     where (bgtr_p(p,:) > 1.0_r8/dt) bgtr_p(p,:) = 1.0_r8/dt

                     ! set carbon fluxes for shifting storage pools to transfer pools

                     pleafc_storage_to_xfer(p,:)  = pleafc_storage(p,:) * bgtr_p(p,:)
                     !do N flux at the same time with the same rate to avoid using CN ratio
                     pleafn_storage_to_xfer(p,:)  = pleafn_storage(p,:) * bgtr_p(p,:)

                 elsewhere (hui(p) >= huilfmatnp(p,:) .and. hui(p) < huilfsennp(p,:)) !transfer all remaining storage, if any
                     pleafc_storage_to_xfer(p,:)  = pleafc_storage(p,:)/dt
                     pleafn_storage_to_xfer(p,:)  = pleafn_storage(p,:)/dt
                     bgtr_p(p,:) = 0._r8

                     !litter fall: for both big-leaf canopy and multilayer canopy
                     !if CN ratio fixed, use background litterfall for updating photosynthetic LAI
                     !if dynamic CN scheme, only retranslocate leaf N gradually during senescence (LAI no change)
                     !but use pruning to move leaf C to litter
                 elsewhere (hui(p) >= huilfsennp(p,:) ) !start of senescence
                     livep(p,:) = 0._r8

                     !a constant bglfr_p value gives ~90% turnover within the senescence period (non-linear decreasing rate)
                     bglfr_p(p,:) = 2._r8/((mxgdd(ivt(p)) - lfsen(ivt(p))) / gddperday * secspday)

                     !a non-constant value as a function of hui gives complete turnover 100%
                     !senescence rate is the inverse of the remaining time until the end of a phytomer life cycle
                     !bglfr_p(p,:) = 1.0_r8 /((huilfendnp(p,:) - hui(p)) / gddperday * secspday)
                     !1.0_r8 gives a linear rate
                     !if multiply by 2 or 3, 4.. in the above equation, it gives a non-linear decreasing rate with 100% turnover
                     !not exceed the max rate 1/dt
                     where (bglfr_p(p,:) > 1.0_r8/dt) bglfr_p(p,:) = 1.0_r8/dt

                 end where

                 !senescence fluxes finish: when leaf is pruned
                 !or when leaf C:N ratio reaches leaf litter C:N ratio (dynamic N scheme)
                 where (rankp(p,:) <= 0 ) bglfr_p(p,:) = 0._r8

                 !summarize for all leaves to update the leafc pool
                 leafc_storage_to_xfer(p) = sum(pleafc_storage_to_xfer(p,:))
                 leafn_storage_to_xfer(p) = sum(pleafn_storage_to_xfer(p,:))
                 !turn on bgtr (only for leaves)
                 if (maxval(bgtr_p(p,:)) > 0._r8) bgtr(p)  = 1.e-6_r8

             end if !end if phytomer level phase 2

             ! =======================
             ! enter harvests:
             !enforce a threshold gddmin for oil palm to start first fruiting (for the whole plant)
             huigrain(p) = max(huigrain(p), gddmin(ivt(p)))

             ! - transfer phytomer grainc to crop yield
             mat1 = sum(minloc(grnmatnp(p,:), mask= huigrnnp(p,:) >= huigrain(p))) !the first mature fruit
             if (hui(p) >= grnmatnp(p,mat1) .and. idpp(p) < mxmat(ivt(p))) then
                 cphase(p) = 3._r8
                 n= mat1 + int(harvest_counter(p)) !this ensure each phytomer is harvested once
                 if (tlai(p) <= 0._r8) then ! plant died during growing season
                   if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday  !reset xsmrpool when died
                   croplive(p) = .false.
                   cphase(p) = 4._r8
                   offset_flag(p) = 1._r8
                   offset_counter(p) = dt
                 else if (hui(p) >= grnmatnp(p,n)) then
                   if (pgrainc(p,n) > 0._r8) then
                      harvest_flag(p) = 1._r8 !(if fruit exists, then do harvest and update harvest_counter afterwards)
                   else
                      harvest_counter(p) = harvest_counter(p) + 1 !(otherwise move on for the next phytomer without doing harvest)
                   endif
                 end if

                 !apply fertilizer at every 12 harvest ~half year (Y.Fan)
                 !ensure fert_counter is renewed only for one time step
                 if (idpp(p) >= int(6*dayspyr) .and. mod(idpp(p),180) ==0 .and. fert_counter(p) <= 0._r8) then
                    fert_counter(p)  = ndays_on * secspday
                    fert(p) = fertnitro(ivt(p)) * 1000._r8 / fert_counter(p)  !gN/m^2/s
                 end if

                 ! continue fertilizer application before final rotation;
                 ! fert_counter resumes every 6 month
                 if (fert_counter(p) <= 0._r8) then
                    fert(p) = 0._r8
                 else ! continue same fert application every timestep
                    fert_counter(p) = fert_counter(p) - dt
                 end if

             ! =======================
             ! final harvest/rotation
             ! - send xsmrpool to the atmosphere
             else if (idpp(p) >= mxmat(ivt(p))) then
                 if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                 croplive(p) = .false.
                 cphase(p) = 4._r8
                 if (tlai(p) > 0._r8) then ! plant had emerged before rotation
                     offset_flag(p) = 1._r8
                     offset_counter(p) = dt
                 else                      ! plant never emerged from the ground
                     crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                     crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                     leafc_xfer(p) = 0._r8  ! revert planting transfers
                     leafn_xfer(p) = 0._r8
                     crop_seedc_to_deadstem(p) = crop_seedc_to_deadstem(p) - deadstemc_xfer(p)/dt
                     crop_seedn_to_deadstem(p) = crop_seedn_to_deadstem(p) - deadstemn_xfer(p)/dt
                     deadstemc_xfer(p) = 0._r8  ! revert planting transfers
                     deadstemn_xfer(p) = 0._r8
                 end if
             end if

             !Continue froot background litterfall and livestem turnover year around;
             ! set background litterfall rate, when not in the phenological offset/rotation
             if (offset_flag(p) == 1._r8) then !must be zero during offset to avoid recalculating of frootc_to_litter
                 bglfr(p) = 0._r8
             else
                 !bglfr is related to leaf turnover rate, used for both froot and livestem turnover
                 !correction for longer growing phase of leaf storage (bud) growth
                 buddays = lfexp(ivt(p))/ gddperday
                 !mxgdd: from leaf expansion to end of leaf senescence (die)
                 bglfr(p) = 1._r8/((mxgdd(ivt(p)) / gddperday + buddays)*secspday)
                 !bglfr(p) = 1._r8/((leaf_long(ivt(p))*dayspyr + buddays)*secspday)
             end if

             ! =======================
             !do pruning:
             !senescent leaves at the bottom layer are pruned to maintain a max number of alive phytomers
             !adjust mxlivenp, mxgdd and phyllochron to determine pruning frequency
             !add 10 for senescent leaves that still on the tree (rankp>=1)
             !if (count(mask=rankp(p,:) >0) > (mxlivenp(ivt(p))+10)) then
             if (harvest_flag(p) /= 1._r8) then !only do pruning at harvest
                 prune(p) = .false.
             else if (count(mask=rankp(p,:) >0) > mxlivenp(ivt(p))) then
                 prune(p) = .true.
             else
                 prune(p) = .false.
             end if

         else   ! crop not live
            ! next 2 lines conserve mass if leaf*_xfer > 0 due to interpinic
            ! We subtract from any existing value in crop_seedc_to_leaf /
            ! crop_seedn_to_leaf in the unlikely event that we enter this block of
            ! code in the same time step where the planting transfer originally
            ! occurred.
            crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
            crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
            onset_counter(p) = 0._r8
            leafc_xfer(p) = 0._r8
            leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
            if (use_c13) then
               c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif
            if (use_c14) then
               c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif

            ! revert transplanting transfers for juvenile plants
            crop_seedc_to_deadstem(p) = crop_seedc_to_deadstem(p) - deadstemc_xfer(p)/dt
            crop_seedn_to_deadstem(p) = crop_seedn_to_deadstem(p) - deadstemn_xfer(p)/dt
            deadstemc_xfer(p) = 0._r8  ! revert planting transfers
            deadstemn_xfer(p) = 0._r8
         end if ! croplive
  
      end if  !end if palm phenology

    end do ! prognostic crops loop

    end associate
  end subroutine PalmPhenology


  !-----------------------------------------------------------------------       
  subroutine CropPhenology(num_pcropp, filter_pcropp                     , &
       waterdiagnosticbulk_inst, temperature_inst, soilstate_inst, wateratm2lndbulk_inst, &
	   crop_inst, canopystate_inst, cnveg_state_inst , &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)

    ! !DESCRIPTION:
    ! Code from AgroIBIS to determine crop phenology and code from CN to
    ! handle CN fluxes during the phenological onset                       & offset periods.
    
    ! !USES:
    use clm_time_manager , only : get_curr_date, get_curr_calday, get_days_per_year, get_rad_step_size
    use pftconMod        , only : ntmp_corn, nswheat, nwwheat, ntmp_soybean
    use pftconMod        , only : nirrig_tmp_corn, nirrig_swheat, nirrig_wwheat, nirrig_tmp_soybean
    use pftconMod        , only : ntrp_corn, nsugarcane, ntrp_soybean, ncotton, nrice
    use pftconMod        , only : nirrig_trp_corn, nirrig_sugarcane, nirrig_trp_soybean
    use pftconMod        , only : nirrig_cotton, nirrig_rice
    use pftconMod        , only : noilpalm, nirrig_oilpalm, ncocoa, nirrig_cocoa
    use clm_varcon       , only : spval, secspday
    use clm_varctl       , only : use_fertilizer 
    use clm_varctl       , only : use_c13, use_c14
    use clm_varcon       , only : c13ratio, c14ratio
    use CNSharedParamsMod, only : use_fun    
    use shr_const_mod    , only : SHR_CONST_TKFRZ, SHR_CONST_PI
    use CNSharedParamsMod, only : CNParamsShareInst
    
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                        , intent(in)    :: filter_pcropp(:) ! filter for prognostic crop patches
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(soilstate_type)           , intent(in)    :: soilstate_inst 
    type(wateratm2lndbulk_type)    , intent(in)    :: wateratm2lndbulk_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    !
    ! LOCAL VARAIBLES:
    integer kyr       ! current year
    integer kmo       ! month of year  (1, ..., 12)
    integer kda       ! day of month   (1, ..., 31)
    integer mcsec     ! seconds of day (0, ..., seconds/day)
    integer jday      ! julian day of the year
    integer fp,p      ! patch indices
    integer c         ! column indices
    integer g         ! gridcell indices
    integer h         ! hemisphere indices
    !integer idpp      ! number of days past planting
    real(r8) :: dtrad ! radiation time step delta t (seconds)
    real(r8) dayspyr  ! days per year
    real(r8) crmcorn  ! comparitive relative maturity for corn
    real(r8) ndays_on ! number of days to fertilize
    real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
    real(r8):: crit_onset_gdd  ! degree days for onset trigger
    real(r8):: soilt           ! temperature of top soil layer
    real(r8):: psi             ! water stress of top soil layer
    real(r8):: rain_threshold  ! rain threshold for leaf on [mm]
    logical :: additional_onset_condition ! additional condition for leaf onset
    !------------------------------------------------------------------------

    associate(                                                                   & 
         ivt               =>    patch%itype                                   , & ! Input:  [integer  (:) ]  patch vegetation type      
         nyrs_crop_active  =>    crop_inst%nyrs_crop_active_patch              ,& ! InOut: [integer (:)  ]  number of years this crop patch has been active         
         leaf_long         =>    pftcon%leaf_long                              , & ! Input:  leaf longevity (yrs)                              
         leafcn            =>    pftcon%leafcn                                 , & ! Input:  leaf C:N (gC/gN)                                  
         manunitro         =>    pftcon%manunitro                              , & ! Input:  max manure to be applied in total (kgN/m2)
         mxmat             =>    pftcon%mxmat                                  , & ! Input:  
         minplanttemp      =>    pftcon%minplanttemp                           , & ! Input:  
         planttemp         =>    pftcon%planttemp                              , & ! Input:  
         gddmin            =>    pftcon%gddmin                                 , & ! Input:  
         hybgdd            =>    pftcon%hybgdd                                 , & ! Input:  
         lfemerg           =>    pftcon%lfemerg                                , & ! Input:  
         grnfill           =>    pftcon%grnfill                                , & ! Input: 
         deadwdcn          =>    pftcon%deadwdcn                               , & ! Input:  dead wood C:N (gC/gN)   
         t_ref2m_min       =>    temperature_inst%t_ref2m_min_patch            , & ! Input:  [real(r8) (:) ]  daily minimum of average 2 m height surface air temperature (K)
         t10               =>    temperature_inst%t_a10_patch                  , & ! Input:  [real(r8) (:) ]  10-day running mean of the 2 m temperature (K)    
         a5tmin            =>    temperature_inst%t_a5min_patch                , & ! Input:  [real(r8) (:) ]  5-day running mean of min 2-m temperature         
         a10tmin           =>    temperature_inst%t_a10min_patch               , & ! Input:  [real(r8) (:) ]  10-day running mean of min 2-m temperature        
         gdd020            =>    temperature_inst%gdd020_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd0                                
         gdd820            =>    temperature_inst%gdd820_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd8                                
         gdd1020           =>    temperature_inst%gdd1020_patch                , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd10                               
         gdd1520           =>    temperature_inst%gdd1520_patch                , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd15
         fertnitro         =>    crop_inst%fertnitro_patch                     , & ! Input:  [real(r8) (:) ]  fertilizer nitrogen
         hui               =>    crop_inst%gddplant_patch                      , & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)                    
         leafout           =>    crop_inst%gddtsoi_patch                       , & ! Input:  [real(r8) (:) ]  gdd from top soil layer temperature              
         harvdate          =>    crop_inst%harvdate_patch                      , & ! Output: [integer  (:) ]  harvest date                                       
         croplive          =>    crop_inst%croplive_patch                      , & ! Output: [logical  (:) ]  Flag, true if planted, not harvested               
         cropplant         =>    crop_inst%cropplant_patch                     , & ! Output: [logical  (:) ]  Flag, true if crop may be planted                  
         vf                =>    crop_inst%vf_patch                            , & ! Output: [real(r8) (:) ]  vernalization factor                              
         peaklai           =>    cnveg_state_inst%peaklai_patch                , & ! Output: [integer  (:) ] 1: max allowed lai; 0: not at max                  
         tlai              =>    canopystate_inst%tlai_patch                   , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     
         
         idop              =>    cnveg_state_inst%idop_patch                   , & ! Output: [integer  (:) ]  date of planting                                   
         gddmaturity       =>    cnveg_state_inst%gddmaturity_patch            , & ! Output: [real(r8) (:) ]  gdd needed to harvest                             
         huileaf           =>    cnveg_state_inst%huileaf_patch                , & ! Output: [real(r8) (:) ]  heat unit index needed from planting to leaf emergence
         huigrain          =>    cnveg_state_inst%huigrain_patch               , & ! Output: [real(r8) (:) ]  same to reach vegetative maturity                 
         cumvd             =>    cnveg_state_inst%cumvd_patch                  , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?             
         hdidx             =>    cnveg_state_inst%hdidx_patch                  , & ! Output: [real(r8) (:) ]  cold hardening index?                             
         bglfr             =>    cnveg_state_inst%bglfr_patch                  , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)                  
         bgtr              =>    cnveg_state_inst%bgtr_patch                   , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)             
         lgsf              =>    cnveg_state_inst%lgsf_patch                   , & ! Output: [real(r8) (:) ]  long growing season factor [0-1]                  
         onset_flag        =>    cnveg_state_inst%onset_flag_patch             , & ! Output: [real(r8) (:) ]  onset flag                                        
         offset_flag       =>    cnveg_state_inst%offset_flag_patch            , & ! Output: [real(r8) (:) ]  offset flag                                       
         onset_counter     =>    cnveg_state_inst%onset_counter_patch          , & ! Output: [real(r8) (:) ]  onset counter                                     
         offset_counter    =>    cnveg_state_inst%offset_counter_patch         , & ! Output: [real(r8) (:) ]  offset counter                                    

         perennial         =>    pftcon%perennial                              , & ! Input:  [integer (:)]  binary flag for perennial crop phenology (1=perennial, 0=not perennial) (added by Y.Fan)
         idpp              =>    cnveg_state_inst%idpp_patch                   , & ! Output: [integer (:)]  days past planting (Y.Fan)
         gddmaturity2      =>    cnveg_state_inst%gddmaturity2_patch           , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest since previous harvest (Y.Fan)
         huigrain2         =>    cnveg_state_inst%huigrain2_patch              , & ! Input:  [real(r8) (:)]  gdd needed from last harvest to start of next grainfill (Y.Fan)
         harvest_flag      =>    crop_inst%harvest_flag_patch                 , & ! Output: [real(r8) (:)]  harvest flag (added by Y.Fan)
         grainc            =>    cnveg_carbonstate_inst%grainc_patch          , & ! Input:  [real(r8) (:) ]  (gC/m2) grain C
         leafc_xfer        =>    cnveg_carbonstate_inst%leafc_xfer_patch       , & ! Output: [real(r8) (:) ]  (gC/m2)   leaf C transfer                           

         crop_seedc_to_leaf =>   cnveg_carbonflux_inst%crop_seedc_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gC/m2/s) seed source to leaf

         fert_counter      =>    cnveg_nitrogenflux_inst%fert_counter_patch    , & ! Output: [real(r8) (:) ]  >0 fertilize; <=0 not (seconds)                   
         leafn_xfer        =>    cnveg_nitrogenstate_inst%leafn_xfer_patch     , & ! Output: [real(r8) (:) ]  (gN/m2)   leaf N transfer                           
         crop_seedn_to_leaf =>   cnveg_nitrogenflux_inst%crop_seedn_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gN/m2/s) seed source to leaf
         cphase            =>    crop_inst%cphase_patch                        , & ! Output: [real(r8) (:)]   phenology phase
	 dayl              =>    grc%dayl                                      , & ! Input:  [real(r8)  (:)   ]  daylength (s)
	 max_dayl          =>    grc%max_dayl                                  , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)
	 min_dayl          =>    grc%min_dayl                                  , & ! Input:  [real(r8) (:)   ]  minimum daylength for this grid cell (s)
	 prev_dayl         =>    grc%prev_dayl                                 , & ! Input:  [real(r8)  (:)   ]  daylength from previous time step (s)
         
         prec10            =>    wateratm2lndbulk_inst%prec10_patch               , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation
	 prec60            =>    wateratm2lndbulk_inst%prec60_patch               , & ! Input:  [real(r8) (:)     ]  60-day running mean of tot. precipitation
	 
	 semi_decid        =>    pftcon%semi_decid                                , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)
	 clearcut_yr       =>    pftcon%clearcut_yr                                , & ! Input: [integer (:)]  clearcut yr for semi-deciduous 
         
         soilpsi           =>    soilstate_inst%soilpsi_col                       , & ! Input:  [real(r8)  (:,:) ]  soil water potential in each soil layer (MPa)   
         
         t_soisno          =>    temperature_inst%t_soisno_col                    , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         dormant_flag      =>    cnveg_state_inst%dormant_flag_patch              , & ! Output:  [real(r8) (:)   ]  dormancy flag                                     
         days_active       =>    cnveg_state_inst%days_active_patch               , & ! Output:  [real(r8) (:)   ]  number of days since last dormancy                
         onset_gddflag     =>    cnveg_state_inst%onset_gddflag_patch             , & ! Output:  [real(r8) (:)   ]  onset freeze flag                                 
         onset_fdd         =>    cnveg_state_inst%onset_fdd_patch                 , & ! Output:  [real(r8) (:)   ]  onset freezing degree days counter                
         onset_gdd         =>    cnveg_state_inst%onset_gdd_patch                 , & ! Output:  [real(r8) (:)   ]  onset growing degree days                         
         onset_swi         =>    cnveg_state_inst%onset_swi_patch                 , & ! Output:  [real(r8) (:)   ]  onset soil water index                                                                 
         offset_fdd        =>    cnveg_state_inst%offset_fdd_patch                , & ! Output:  [real(r8) (:)   ]  offset freezing degree days counter               
         offset_swi        =>    cnveg_state_inst%offset_swi_patch                , & ! Output:  [real(r8) (:)   ]  offset soil water index                           
         annavg_t2m        =>    cnveg_state_inst%annavg_t2m_patch                           , & ! Output:  [real(r8) (:)   ]  annual average 2m air temperature (K)             
         leafc             =>    cnveg_carbonstate_inst%leafc_patch                          , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C      
         frootc            =>    cnveg_carbonstate_inst%frootc_patch                         , & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C
         leafc_storage     =>    cnveg_carbonstate_inst%leafc_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage                            
         frootc_storage    =>    cnveg_carbonstate_inst%frootc_storage_patch                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage                       
         gresp_storage     =>    cnveg_carbonstate_inst%gresp_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage              
         frootc_xfer       =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer
         leafn_storage     =>    cnveg_nitrogenstate_inst%leafn_storage_patch                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage                            
         frootn_storage    =>    cnveg_nitrogenstate_inst%frootn_storage_patch               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage                                                  
         frootn_xfer       =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer                      
         prev_leafc_to_litter     =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch     , & ! Output:  [real(r8) (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter    =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch    , & ! Output:  [real(r8) (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc      =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_xfer_to_frootc    =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         leafc_storage_to_xfer    =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_storage_to_xfer   =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         gresp_storage_to_xfer    =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch        , & ! Output:  [real(r8) (:)   ]
	 
	 woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
	 
	 livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer
	 livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer               
         
	 livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage              
	 livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage 
	 
	 livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:)   ]
	 livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:)   ]                                                  
                                                
	 livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch    ,  & ! Output:  [real(r8) (:)   ]        
         fert              =>    cnveg_nitrogenflux_inst%fert_patch              & ! Output: [real(r8) (:) ]  (gN/m2/s) fertilizer applied each timestep 
         )

      ! get time info
      dayspyr = get_days_per_year()
      jday    = get_curr_calday()
      call get_curr_date(kyr, kmo, kda, mcsec)
      dtrad   = real( get_rad_step_size(), r8 )
	  
      ! specify rain threshold for leaf onset
      rain_threshold = 20._r8

      if (use_fertilizer) then
       ndays_on = 20._r8 ! number of days to fertilize
      else
       ndays_on = 0._r8 ! number of days to fertilize
      end if

      do fp = 1, num_pcropp
         p = filter_pcropp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)
         h = inhemi(p)
	 
	 ! threshold value of soilpsi off
	 !if (prec10(p) * (3600.0_r8*10.0_r8*24.0_r8) == 0._r8 .and. &
	     !prec60(p) * (3600.0_r8*60.0_r8*24.0_r8) < 30._r8) then
            !soilpsi_off = -2.0_r8  	   
	 !else 
	    !soilpsi_off = -2.5_r8  !-0.8_r8	    
         !end if	

         ! background litterfall and transfer rates; long growing season factor

         bglfr(p) = 0._r8 ! this value changes later in a crop's life cycle
         bgtr(p)  = 0._r8
         lgsf(p)  = 0._r8

         ! ---------------------------------
         ! from AgroIBIS subroutine planting
         ! ---------------------------------

         ! in order to allow a crop to be planted only once each year
         ! initialize cropplant = .false., but hold it = .true. through the end of the year

         ! initialize other variables that are calculated for crops
         ! on an annual basis in cropresidue subroutine

         if ( jday == jdayyrstart(h) .and. mcsec == 0 )then

            ! make sure variables aren't changed at beginning of the year
            ! for a crop that is currently planted, such as
            ! WINTER TEMPERATE CEREAL = winter (wheat + barley + rye)
            ! represented here by the winter wheat pft

            if (.not. croplive(p))  then
               cropplant(p) = .false.
               idop(p)      = NOT_Planted

               ! keep next for continuous, annual winter temperate cereal crop;
               ! if we removed elseif,
               ! winter cereal grown continuously would amount to a cereal/fallow
               ! rotation because cereal would only be planted every other year

            else if (croplive(p) .and. (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat)) then
               cropplant(p) = .false.
               !           else ! not possible to have croplive and ivt==cornORsoy? (slevis)
            end if

         end if

         if ( (.not. croplive(p)) .and. (.not. cropplant(p)) ) then

            ! gdd needed for * chosen crop and a likely hybrid (for that region) *
            ! to reach full physiological maturity

            ! based on accumulated seasonal average growing degree days from
            ! April 1 - Sept 30 (inclusive)
            ! for corn and soybeans in the United States -
            ! decided upon by what the typical average growing season length is
            ! and the gdd needed to reach maturity in those regions

            ! first choice is used for spring temperate cereal and/or soybeans and maize

            ! slevis: ibis reads xinpdate in io.f from control.crops.nc variable name 'plantdate'
            !         According to Chris Kucharik, the dataset of
            !         xinpdate was generated from a previous model run at 0.5 deg resolution

            ! winter temperate cereal : use gdd0 as a limit to plant winter cereal

            if (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat) then

               ! add check to only plant winter cereal after other crops (soybean, maize)
               ! have been harvested

               ! *** remember order of planting is crucial - in terms of which crops you want
               ! to be grown in what order ***

               ! in this case, corn or soybeans are assumed to be planted before
               ! cereal would be in any particular year that both patches are allowed
               ! to grow in the same grid cell (e.g., double-cropping)

               ! slevis: harvdate below needs cropplant(p) above to be cropplant(p,ivt(p))
               !         where ivt(p) has rotated to winter cereal because
               !         cropplant through the end of the year for a harvested crop.
               !         Also harvdate(p) should be harvdate(p,ivt(p)) and should be
               !         updated on Jan 1st instead of at harvest (slevis)
               if (a5tmin(p)             /= spval                  .and. &
                    a5tmin(p)             <= minplanttemp(ivt(p))   .and. &
                    jday                  >= minplantjday(ivt(p),h) .and. &
                    (gdd020(p)            /= spval                  .and. &
                    gdd020(p)             >= gddmin(ivt(p)))) then

                  cumvd(p)       = 0._r8
                  hdidx(p)       = 0._r8
                  vf(p)          = 0._r8
                  croplive(p)    = .true.
                  cropplant(p)   = .true.
                  idop(p)        = jday
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = initial_seed_at_planting
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
                  if (use_c13) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
                     endif
                  endif
                  if (use_c14) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
                     endif
                  endif

                  ! latest possible date to plant winter cereal and after all other 
                  ! crops were harvested for that year

               else if (jday       >=  maxplantjday(ivt(p),h) .and. &
                    gdd020(p)  /= spval                   .and. &
                    gdd020(p)  >= gddmin(ivt(p))) then

                  cumvd(p)       = 0._r8
                  hdidx(p)       = 0._r8
                  vf(p)          = 0._r8
                  croplive(p)    = .true.
                  cropplant(p)   = .true.
                  idop(p)        = jday
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = initial_seed_at_planting
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
                  if (use_c13) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
                     endif
                  endif
                  if (use_c14) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
                     endif
                  endif
               else
                  gddmaturity(p) = 0._r8
               end if

            else ! not winter cereal... slevis: added distinction between NH and SH
               ! slevis: The idea is that jday will equal idop sooner or later in the year
               !         while the gdd part is either true or false for the year.
               if (t10(p) /= spval.and. a10tmin(p) /= spval   .and. &
                    t10(p)     > planttemp(ivt(p))             .and. &
                    a10tmin(p) > minplanttemp(ivt(p))          .and. &
                    jday       >= minplantjday(ivt(p),h)       .and. &
                    jday       <= maxplantjday(ivt(p),h)       .and. &
                    t10(p) /= spval .and. a10tmin(p) /= spval  .and. &
                    gdd820(p) /= spval                         .and. &
                    gdd820(p) >= gddmin(ivt(p))) then

                  ! impose limit on growing season length needed
                  ! for crop maturity - for cold weather constraints
                  croplive(p)  = .true.
                  cropplant(p) = .true.
                  idop(p)      = jday
                  harvdate(p)  = NOT_Harvested

                  ! go a specified amount of time before/after
                  ! climatological date
                  if (ivt(p) == ntmp_soybean .or. ivt(p) == nirrig_tmp_soybean .or. &
                       ivt(p) == ntrp_soybean .or. ivt(p) == nirrig_trp_soybean) then
                     gddmaturity(p) = min(gdd1020(p), hybgdd(ivt(p)))
                  end if
                  if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
                      ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
                      ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane) then
                     gddmaturity(p) = max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                     gddmaturity(p) = max(950._r8, min(gddmaturity(p)+150._r8, 1850._r8))
                  end if
                  if (ivt(p) == nswheat .or. ivt(p) == nirrig_swheat .or. &
                      ivt(p) == ncotton .or. ivt(p) == nirrig_cotton .or. &
                      ivt(p) == nrice   .or. ivt(p) == nirrig_rice) then
                     gddmaturity(p) = min(gdd020(p), hybgdd(ivt(p)))
                  end if
		  
		          if (ivt(p)==noilpalm .or. ivt(p) == nirrig_oilpalm) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*2.5_r8, hybgdd(ivt(p))))
		  
		          if (ivt(p)==ncocoa .or. ivt(p) == nirrig_cocoa) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*6.0_r8, hybgdd(ivt(p))))

                  if (ivt(p)==noilpalm .or. ivt(p) == nirrig_oilpalm) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*2.5_r8, hybgdd(ivt(p))))
                  !2.5 years to vegetative maturity /first harvest (Y.Fan)
		  
		  if (ivt(p)==ncocoa .or. ivt(p) == nirrig_cocoa) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*6.0_r8, hybgdd(ivt(p))))

                  leafc_xfer(p)  = initial_seed_at_planting
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
                  if (use_c13) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
                     endif
                  endif
                  if (use_c14) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
                     endif
                  endif


                  ! If hit the max planting julian day -- go ahead and plant
               else if (jday == maxplantjday(ivt(p),h) .and. gdd820(p) > 0._r8 .and. &
                    gdd820(p) /= spval ) then
                  croplive(p)  = .true.
                  cropplant(p) = .true.
                  idop(p)      = jday
                  harvdate(p)  = NOT_Harvested

                  if (ivt(p) == ntmp_soybean .or. ivt(p) == nirrig_tmp_soybean .or. &
                      ivt(p) == ntrp_soybean .or. ivt(p) == nirrig_trp_soybean) then
                     gddmaturity(p) = min(gdd1020(p), hybgdd(ivt(p)))
                  end if
                  if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
                      ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
                      ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane) then
                     gddmaturity(p) = max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                  end if
                  if (ivt(p) == nswheat .or. ivt(p) == nirrig_swheat .or. &
                      ivt(p) == ncotton .or. ivt(p) == nirrig_cotton .or. &
                      ivt(p) == nrice   .or. ivt(p) == nirrig_rice) then
                     gddmaturity(p) = min(gdd020(p), hybgdd(ivt(p)))
                  end if
		  
		  if (ivt(p)==noilpalm .or. ivt(p) == nirrig_oilpalm) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*2.5_r8, hybgdd(ivt(p))))
				  !2.5 years to vegetative maturity /first harvest (Y.Fan)
				  
		  if (ivt(p)==ncocoa .or. ivt(p) == nirrig_cocoa) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*6.0_r8, hybgdd(ivt(p))))
				  !6.0 years to vegetative maturity /first harvest (A.Ali)

                  if (ivt(p)==noilpalm .or. ivt(p) == nirrig_oilpalm) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*2.5_r8, hybgdd(ivt(p))))
                  !2.5 years to vegetative maturity /first harvest (Y.Fan)
		  
		  if (ivt(p)==ncocoa .or. ivt(p) == nirrig_cocoa) gddmaturity(p)=max(10000._r8, min(gdd1520(p)*6.0_r8, hybgdd(ivt(p))))

                  leafc_xfer(p)  = initial_seed_at_planting
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
                  if (use_c13) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
                     endif
                  endif
                  if (use_c14) then
                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
                     else
                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
                     endif
                  endif

               else
                  gddmaturity(p) = 0._r8
               end if
            end if ! crop patch distinction

            ! crop phenology (gdd thresholds) controlled by gdd needed for
            ! maturity (physiological) which is based on the average gdd
            ! accumulation and hybrids in United States from April 1 - Sept 30

            ! calculate threshold from phase 1 to phase 2:
            ! threshold for attaining leaf emergence (based on fraction of
            ! gdd(i) -- climatological average)
            ! Hayhoe and Dwyer, 1990, Can. J. Soil Sci 70:493-497
            ! Carlson and Gage, 1989, Agric. For. Met., 45: 313-324
            ! J.T. Ritchie, 1991: Modeling Plant and Soil systems

            huileaf(p) = lfemerg(ivt(p)) * gddmaturity(p) ! 3-7% in cereal

            ! calculate threshhold from phase 2 to phase 3:
            ! from leaf emergence to beginning of grain-fill period
            ! this hypothetically occurs at the end of tassling, not the beginning
            ! tassel initiation typically begins at 0.5-0.55 * gddmaturity

            ! calculate linear relationship between huigrain fraction and relative
            ! maturity rating for maize

            if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
                ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
                ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane) then
               ! the following estimation of crmcorn from gddmaturity is based on a linear
               ! regression using data from Pioneer-brand corn hybrids (Kucharik, 2003,
               ! Earth Interactions 7:1-33: fig. 2)
               crmcorn = max(73._r8, min(135._r8, (gddmaturity(p)+ 53.683_r8)/13.882_r8))

               ! the following adjustment of grnfill based on crmcorn is based on a tuning
               ! of Agro-IBIS to give reasonable results for max LAI and the seasonal
               ! progression of LAI growth (pers. comm. C. Kucharik June 10, 2010)
               huigrain(p) = -0.002_r8  * (crmcorn - 73._r8) + grnfill(ivt(p))

               huigrain(p) = min(max(huigrain(p), grnfill(ivt(p))-0.1_r8), grnfill(ivt(p)))
               huigrain(p) = huigrain(p) * gddmaturity(p)     ! Cabelguenne et al. 1999

               !oil palm begins to produce 2.5 years after planted, later grainfill takes place every year
               !the PFT value grnfill(ivt(p)) is only used here for the first grainfill after planting
            else
               huigrain(p) = grnfill(ivt(p)) * gddmaturity(p)
               !initialize huigrain2 for perennial crops, update annually or monthly in the harvest cycle(Y.Fan)
               huigrain2(p) = huigrain(p)
            end if
            onset_counter(p) = 0._r8 !onset_counter initialized to zero when .not. croplive (Y.Fan)
            offset_counter(p) = 0._r8 !offset_counter initialized to zero when .not. croplive (Y.Fan)
            if (semi_decid(ivt(p)) == 1._r8) onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
            if (semi_decid(ivt(p)) == 1._r8) offset_flag(p) = 0._r8 ! carbon and nitrogen transfers 

         end if ! crop not live nor planted

         ! ----------------------------------
         ! from AgroIBIS subroutine phenocrop
         ! ----------------------------------

         ! all of the phenology changes are based on the total number of gdd needed
         ! to change to the next phase - based on fractions of the total gdd typical
         ! for  that region based on the April 1 - Sept 30 window of development

         ! crop phenology (gdd thresholds) controlled by gdd needed for
         ! maturity (physiological) which is based on the average gdd
         ! accumulation and hybrids in United States from April 1 - Sept 30

         ! Phase 1: Planting to leaf emergence (now in CNAllocation)
         ! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline)
         ! Harvest: if gdd past grain fill initiation exceeds limit
         ! or number of days past planting reaches a maximum, the crop has
         ! reached physiological maturity and plant is harvested;
         ! crop could be live or dead at this stage - these limits
         ! could lead to reaching physiological maturity or determining
         ! a harvest date for a crop killed by an early frost (see next comments)
         ! --- --- ---
         ! keeping comments without the code (slevis):
         ! if minimum temperature, t_ref2m_min <= freeze kill threshold, tkill
         ! for 3 consecutive days and lai is above a minimum,
         ! plant will be damaged/killed. This function is more for spring freeze events
         ! or for early fall freeze events

         ! spring temperate cereal is affected by this, winter cereal kill function
         ! is determined in crops.f - is a more elaborate function of
         ! cold hardening of the plant

         ! currently simulates too many grid cells killed by freezing temperatures

         ! removed on March 12 2002 - C. Kucharik
         ! until it can be a bit more refined, or used at a smaller scale.
         ! we really have no way of validating this routine
         ! too difficult to implement on 0.5 degree scale grid cells
         ! --- --- ---

         if (semi_decid(ivt(p)) /= 1._r8) onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
         if (semi_decid(ivt(p)) /= 1._r8) offset_flag(p) = 0._r8 ! carbon and nitrogen transfers
         !onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
         !offset_flag(p) = 0._r8 ! carbon and nitrogen transfers
         harvest_flag(p) = 0._r8 ! annual harvest flag for perennial crops (Y.Fan)

         if (croplive(p)) then
            cphase(p) = 1._r8

            ! call vernalization if winter temperate cereal planted, living, and the
            ! vernalization factor is not 1;
            ! vf affects the calculation of gddtsoi & gddplant

            if (t_ref2m_min(p) < 1.e30_r8 .and. vf(p) /= 1._r8 .and. &
               (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat)) then
               call vernalization(p, &
                    canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, &
                    crop_inst)
            end if

            ! days past planting may determine harvest
            ! add new counter for perennial crops (Y.Fan)
            if (perennial(ivt(p)) == 1) then
               if (int(dayspyr)*(nyrs_crop_active(p)-1) <= mxmat(ivt(p))) then !nyrs=1 for the start of first year
                  idpp(p) = int(dayspyr)*(nyrs_crop_active(p)-1) + jday - idop(p)
               else !mod is to get the reminder for each rotation; minus 1 to account for the offset
                  !because crops are allowed to be planted only once a year, controlled by flag cropplant(p)
                  idpp(p) = int(dayspyr)*(mod((nyrs_crop_active(p)-1), nint(mxmat(ivt(p))/dayspyr)) - 1) + jday - idop(p)
               end if
            else  !only for annual or biannual crops
               if (jday >= idop(p)) then
                  idpp(p) = jday - idop(p)
               else
                  idpp(p) = int(dayspyr) + jday - idop(p)
               end if
            end if
    
            ! Ashehad Ali added codes from stress_deciduous phenology here
   
            if (semi_decid(ivt(p)) == 1._r8) then
               soilt = t_soisno(c,3)
               psi = soilpsi(c,3)
               
               ! onset gdd sum from Biome-BGC, v4.1.2
               crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))
               
                  
               ! update offset_counter and test for the end of the offset period
               if (offset_flag(p) == 1._r8) then
                  ! decrement counter for offset period
                  offset_counter(p) = offset_counter(p) - dt
               
                  ! if this is the end of the offset_period, reset phenology
                  ! flags and indices				   
                  if (offset_counter(p) == dt) then
               	  ! this code block was originally handled by call cn_offset_cleanup(p)
               	  ! inlined during vectorization
               	  offset_flag(p) = 0._r8
               	  offset_counter(p) = 0._r8
               	  dormant_flag(p) = 1._r8
               	  days_active(p) = 0._r8
               	  
               	  ! reset the previous timestep litterfall flux memory
               	  prev_leafc_to_litter(p) = 0._r8
               	  prev_frootc_to_litter(p) = 0._r8
                  end if
               end if
               
               ! update onset_counter and test for the end of the onset period
               if (onset_flag(p) == 1.0_r8) then
                  ! decrement counter for onset period
                  onset_counter(p) = onset_counter(p) - dt
               
                  ! if this is the end of the onset period, reset phenology
                  ! flags and indices				   
                  if (onset_counter(p) == dt) then
               	  ! this code block was originally handled by call cn_onset_cleanup(p)
               	  ! inlined during vectorization
               	  onset_flag(p) = 0._r8
               	  onset_counter(p) = 0._r8
               	  ! set all transfer growth rates to 0.0
               	  leafc_xfer_to_leafc(p)   = 0._r8
               	  frootc_xfer_to_frootc(p) = 0._r8
               	  leafn_xfer_to_leafn(p)   = 0._r8
               	  frootn_xfer_to_frootn(p) = 0._r8
                 if (woody(ivt(p)) == 1.0_r8) then
               		 livestemc_xfer_to_livestemc(p)   = 0._r8
               		 deadstemc_xfer_to_deadstemc(p)   = 0._r8
               		 livecrootc_xfer_to_livecrootc(p) = 0._r8
               		 deadcrootc_xfer_to_deadcrootc(p) = 0._r8
               		 livestemn_xfer_to_livestemn(p)   = 0._r8
               		 deadstemn_xfer_to_deadstemn(p)   = 0._r8
               		 livecrootn_xfer_to_livecrootn(p) = 0._r8
               		 deadcrootn_xfer_to_deadcrootn(p) = 0._r8
               	  end if                  
               	  ! set transfer pools to 0.0
               	  leafc_xfer(p) = 0._r8
               	  leafn_xfer(p) = 0._r8
               	  frootc_xfer(p) = 0._r8
               	  frootn_xfer(p) = 0._r8
                 if (woody(ivt(p)) == 1.0_r8) then
               		 livestemc_xfer(p) = 0._r8
               		 livestemn_xfer(p) = 0._r8
               		 deadstemc_xfer(p) = 0._r8
               		 deadstemn_xfer(p) = 0._r8
               		 livecrootc_xfer(p) = 0._r8
               		 livecrootn_xfer(p) = 0._r8
               		 deadcrootc_xfer(p) = 0._r8
               		 deadcrootn_xfer(p) = 0._r8
               	  end if
               	  
                  end if
               end if		
               
               ! test for switching from dormant period to growth period
               if (dormant_flag(p) == 1._r8) then
               
                  ! keep track of the number of freezing degree days in this
                  ! dormancy period (only if the freeze flag has not previously been set
                  ! for this dormancy period
               
                  if (onset_gddflag(p) == 0._r8 .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday
               
                  ! if the number of freezing degree days exceeds a critical value,
                  ! then onset will require both wet soils and a critical soil
                  ! temperature sum.  If this case is triggered, reset any previously
                  ! accumulated value in onset_swi, so that onset now depends on
                  ! the accumulated soil water index following the freeze trigger
               
                  if (onset_fdd(p) > crit_onset_fdd) then
               	  onset_gddflag(p) = 1._r8
               	  onset_fdd(p) = 0._r8
               	  onset_swi(p) = 0._r8
                  end if
               
                  ! if the freeze flag is set, and if the soil is above freezing
                  ! then accumulate growing degree days for onset trigger
               
                  if (onset_gddflag(p) == 1._r8 .and. soilt > SHR_CONST_TKFRZ) then
               	  onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
                  end if
               
                 ! if soils are wet, accumulate soil water index for onset trigger
                  additional_onset_condition = .true.
                  if(CNParamsShareInst%constrain_stress_deciduous_onset) then
               	   !if additional constraint condition not met,  set to false
               	  if ((prec10(p) * (3600.0_r8*10.0_r8*24.0_r8)) < rain_threshold) then
               		 additional_onset_condition = .false.
               	  endif
                  endif
               
                  if (psi >= soilpsi_on) then
               	  onset_swi(p) = onset_swi(p) + fracday
                  endif
               
                  ! if critical soil water index is exceeded, set onset_flag, and
                  ! then test for soil temperature criteria
               
                  ! Adding in Kyla's rainfall trigger when fun on. RF. prec10 (mm/s) needs to be higher than 8mm over 10 days. 
               
                  if (onset_swi(p) > crit_onset_swi .and. additional_onset_condition)  then
               	  onset_flag(p) = 1._r8
                 
               	  ! only check soil temperature criteria if freeze flag set since
               	  ! beginning of last dormancy.  If freeze flag set and growing
               	  ! degree day sum (since freeze trigger) is lower than critical
               	  ! value, then override the onset_flag set from soil water.
               
               	  if (onset_gddflag(p) == 1._r8 .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0._r8
                  end if
               
                  ! only allow onset if dayl > 6hrs
                  if (onset_flag(p) == 1._r8 .and. dayl(g) <= secspqtrday) then
               	  onset_flag(p) = 0._r8
                  end if
               
                  ! if this is the beginning of the onset period
                  ! then reset the phenology flags and indices
               
                   if (onset_flag(p) == 1._r8) then
               	  dormant_flag(p) = 0._r8
               	  days_active(p) = 0._r8
               	  onset_gddflag(p) = 0._r8
               	  onset_fdd(p) = 0._r8
               	  onset_gdd(p) = 0._r8
               	  onset_swi(p) = 0._r8
               	  onset_counter(p) = ndays_on * secspday
               	  cphase(p) = 2.0_r8
               	  
               	  fert_counter(p)  = ndays_on * secspday
                             if (ndays_on .gt. 0) then
                                fert(p) = (manunitro(ivt(p)) * 1000._r8 + fertnitro(p))/ fert_counter(p)
                             else
                                fert(p) = 0._r8
                             end if				  
               	  
               	  ! call subroutine to move all the storage pools into transfer pools,
               	  ! where they will be transfered to displayed growth over the onset period.
               	  ! this code was originally handled with call cn_storage_to_xfer(p)
               	  ! inlined during vectorization
               
               	  ! set carbon fluxes for shifting storage pools to transfer pools					  
               	  					  
               	  leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt 
               	  frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
                 if (woody(ivt(p)) == 1.0_r8) then
               		 livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt 
               		 deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt 
               		 livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt 
               		 deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt 
               		 gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt  
               	  end if                  
               
               	  ! set nitrogen fluxes for shifting storage pools to transfer pools
               	  leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt 
               	  frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
                 if (woody(ivt(p)) == 1.0_r8) then
               		 livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
               		 deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
               		 livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
               		 deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               	  end if                  
                   end if
               
                  ! test for switching from growth period to offset period
               else if (offset_flag(p) == 0._r8) then
               
                  ! if soil water potential lower than critical value, accumulate
                  ! as stress in offset soil water index
                  
                  if (psi <= -4._r8) then
               	  offset_swi(p) = offset_swi(p) + fracday
               
               	  ! if the offset soil water index exceeds critical value, and
               	  ! if this is not the middle of a previously initiated onset period,
               	  ! then set flag to start the offset period and reset index variables
               
               	  if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
               
               	  ! if soil water potential higher than critical value, reduce the
               	  ! offset water stress index.  By this mechanism, there must be a
               	  ! sustained period of water stress to initiate offset.
               
                  else if (psi >= soilpsi_on) then
               	  offset_swi(p) = offset_swi(p) - fracday
               	  offset_swi(p) = max(offset_swi(p),0._r8)
                  end if
               
                  ! decrease freezing day accumulator for warm soil
                  if (offset_fdd(p) > 0._r8 .and. soilt > SHR_CONST_TKFRZ) then
               	  offset_fdd(p) = offset_fdd(p) - fracday
               	  offset_fdd(p) = max(0._r8, offset_fdd(p))
                  end if
               
                  ! increase freezing day accumulator for cold soil
                  if (soilt <= SHR_CONST_TKFRZ) then
               	  offset_fdd(p) = offset_fdd(p) + fracday
               
               	  ! if freezing degree day sum is greater than critical value, initiate offset
               	  if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
                  end if
               
                  ! force offset if daylength is < 6 hrs
                  if (dayl(g) <= secspqtrday) then
               	  offset_flag(p) = 1._r8
                  end if
                  
                  call get_curr_date(kyr, kmo, kda, mcsec)
                  
                 ! force offset if dayl condition is met  
                  
                  if (kyr >= clearcut_yr(ivt(p)) + 4) then				   
               	   if (max_dayl(g) <= 12.5_r8 * 3600._r8 .and. &
               	       dayl(g) >= min_dayl(g) + 0.028_r8 * 3600._r8 .and. & 
                              dayl(g) >= prev_dayl(g) .and. &
                              dayl(g) <= min_dayl(g) + 0.028_r8 * 3600._r8 * 1.1_r8) then
               		   offset_flag(p) = 1._r8		   
               	   end if					   
               	   if (max_dayl(g) > 12.5_r8 * 3600._r8 .and. &					    
                              dayl(g) >= prev_dayl(g) .and. &
                              dayl(g) <= min_dayl(g) + 0.02_r8 * 3600._r8) then
               	       offset_flag(p) = 1._r8		   
               	   end if
                  end if			   
                   
                  write(12,*) kyr, max_dayl(g), dayl(g), min_dayl(g), prev_dayl(g)		
                  
                  if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p) .and. idpp(p) < mxmat(ivt(p))) then
                      cphase(p) = 2._r8	   
                    else if (hui(p) >= gddmaturity(p) .and. idpp(p) < mxmat(ivt(p)) .and. &
                      perennial(ivt(p)) == 1) then
                      gddmaturity2(p) = huigrain2(p) + (1._r8 - grnfill(ivt(p))) * gddmaturity(p)/12._r8 !add a new parameter grnfill2 (Y.Fan 2022.09.13)		
                     ! part of leaves are pruned off at each harvest (at one-time step)
                      if (hui(p) >= gddmaturity2(p)) then
                        !the following ensures monthly harvests for perennial crops
                        !Assume a percentage 0.5 of the initial grnfill=around 1 year period between subsequent grnfills
                        huigrain2(p) = huigrain2(p) + 0.5_r8 * grnfill(ivt(p)) * gddmaturity(p)/12._r8 !to add a new parameter grnfill2 (Y.Fan 2022.09.13)		   
                        !apply fertilizer after each grainfill stage before next grainfill (Y.Fan)
                        !fert continue until the end of ndays_on because this loop is entered one time step only
                        !and fert_counter(p) = fert_counter(p) - dt (see later)
                        fert_counter(p)  = ndays_on * secspday
                        fert(p) = manunitro(ivt(p)) * 1000._r8 / fert_counter(p)
                        if (tlai(p) <= 0._r8) then ! plant never emerged or died
                           croplive(p) = .false.
                           if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                              crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                              crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                              leafc_xfer(p) = 0._r8  ! revert planting transfers
                              leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))     

                              !Ashehad added this (not needed when
                              !crop_seedc_to_deadstemc is not set at
                              !the time of planting, YFan 2023)
                              !  if (woody(ivt(p)) == 1._r8) then
                              !    deadstemc_xfer(p) = 0._r8
                              !    deadstemn_xfer(p) = deadstemc_xfer(p) / deadwdcn(ivt(p))
                              !  end if
                        
                         else if (grainc(p) > 0._r8) then !only harvest when there is positive grainc accumulated during grainfill
                            harvest_flag(p) = 1._r8                                             
                         end if
                 
                      else
                           harvest_flag(p) = 0._r8
                           ! continue background litterfall year around because perennials keeps alive after harvest (Y.Fan)
                           bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
                      end if
                      
                   ! ashehad added this for high latitudes
                   else if (hui(p) >= gddmaturity(p) .and. idpp(p) < mxmat(ivt(p)) .and. &
                      perennial(ivt(p)) == 1 .and. max_dayl(g) > 12.5_r8 * 3600._r8 .and. & 
                      kmo > 4 .and. kmo < 12) then
                      gddmaturity2(p) = huigrain2(p) + (1._r8 - grnfill(ivt(p))) * gddmaturity(p)/12._r8 !add a new parameter grnfill2 (Y.Fan 2022.09.13)		
                     ! part of leaves are pruned off at each harvest (at one-time step)
                      if (hui(p) >= gddmaturity2(p)) then
                        !the following ensures monthly harvests for perennial crops
                        !Assume a percentage 0.5 of the initial grnfill=around 1 year period between subsequent grnfills
                        huigrain2(p) = huigrain2(p) + 0.5_r8 * grnfill(ivt(p)) * gddmaturity(p)/12._r8 !add a new parameter grnfill2 (Y.Fan 2022.09.13)		   
                        !apply fertilizer after each grainfill stage before next grainfill (Y.Fan)
                        !fert continue until the end of ndays_on because this loop is entered one time step only
                        !and fert_counter(p) = fert_counter(p) - dt (see later)
                        fert_counter(p)  = ndays_on * secspday
                        fert(p) = manunitro(ivt(p)) * 1000._r8 / fert_counter(p)
                        if (tlai(p) <= 0._r8) then ! plant never emerged or died
                           croplive(p) = .false.
                           if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                              crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                              crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                              leafc_xfer(p) = 0._r8  ! revert planting transfers
                              leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))    
                  
                              !Ashehad added this (not needed when
                              !crop_seedc_to_deadstemc is not set at
                              !the time of planting, YFan 2023)
                              !  if (woody(ivt(p)) == 1._r8) then
                              !    deadstemc_xfer(p) = 0._r8
                              !    deadstemn_xfer(p) = deadstemc_xfer(p) / deadwdcn(ivt(p))
                              !  end if
                        
                        else if (grainc(p) > 0._r8) then !only harvest when there is positive grainc accumulated during grainfill
                           harvest_flag(p) = 1._r8                                             
                        end if
                 
                      else
                              harvest_flag(p) = 0._r8
                              ! continue background litterfall year around because perennials keeps alive after harvest (Y.Fan)
                              bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
                      end if 
                   ! end high latitudes
               
                  else if (hui(p) >= gddmaturity(p) .or. idpp(p) >= mxmat(ivt(p))) then
                      if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                         croplive(p) = .false.     ! no re-entry in greater if-block
                         cphase(p) = 4._r8
                         if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                            offset_flag(p) = 1._r8
                            offset_counter(p) = dt
                         else                      ! plant never emerged from the ground
                        ! Revert planting transfers; this will replenish the crop seed deficit.
                        ! We subtract from any existing value in crop_seedc_to_leaf /
                        ! crop_seedn_to_leaf in the unlikely event that we enter this block of
                        ! code in the same time step where the planting transfer originally
                        ! occurred.
                        crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                        crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                        leafc_xfer(p) = 0._r8
                        leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
                        !Ashehad added this (not needed when
                        !crop_seedc_to_deadstemc is not set at
                        !the time of planting, YFan 2023) 
                        !  if (woody(ivt(p)) == 1._r8) then
                        !     deadstemc_xfer(p) = 0._r8
                        !     deadstemn_xfer(p) = deadstemc_xfer(p) / deadwdcn(ivt(p))
                        !  end if
                       end if
                     
                  else if (hui(p) >= huigrain(p)) then
                          cphase(p) = 3._r8
                          bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
                  end if
               
                 ! continue fertilizer application while in phase 2;
                 ! assumes that onset of phase 2 took one time step only
               
                 if (fert_counter(p) <= 0._r8) then
                    fert(p) = 0._r8
                 else ! continue same fert application every timestep
                    fert_counter(p) = fert_counter(p) - dtrad
                 end if
               
               ! if this is the beginning of the offset period
                  ! then reset flags and indices
                  if (offset_flag(p) == 1._r8) then
               	  offset_fdd(p) = 0._r8
               	  offset_swi(p) = 0._r8
               	  offset_counter(p) = ndays_off * secspday * 1.0_r8  
               	  prev_leafc_to_litter(p) = 0._r8
               	  prev_frootc_to_litter(p) = 0._r8
                  end if
               end if
               
               ! keep track of number of days since last dormancy for control on
               ! fraction of new growth to send to storage for next growing season
               
               if (dormant_flag(p) == 0.0_r8) then
                  days_active(p) = days_active(p) + fracday
               end if
               
               ! calculate long growing season factor (lgsf)
               ! only begin to calculate a lgsf greater than 0.0 once the number
               ! of days active exceeds days/year.
               lgsf(p) = max(min(3.0_r8*(days_active(p)-leaf_long(ivt(p))*dayspyr )/dayspyr, 1._r8),0._r8)
               ! RosieF. 5 Nov 2015.  Changed this such that the increase in leaf turnover is faster after
               ! trees enter the 'fake evergreen' state. Otherwise, they have a whole year of 
               ! cheating, with less litterfall than they should have, resulting in very high LAI. 
               ! Further, the 'fake evergreen' state (where lgsf>0) is entered at the end of a single leaf lifespan
               ! and not a whole year. The '3' is arbitrary, given that this entire system is quite abstract. 
               
               
               ! set background litterfall rate, when not in the phenological offset period
               ! 1.5_r8 and 0.5_r8 as in Land paper			
               
               call get_curr_date(kyr, kmo, kda, mcsec)
               
               if (offset_flag(p) == 1._r8) then
                  bglfr(p) = 0._r8 
               else
                  bglfr(p) = 1._r8/(leaf_long(ivt(p)) * dayspyr * secspday)   
               end if
               
               ! set background transfer rate when active but not in the phenological onset period
               if (onset_flag(p) == 1._r8) then
                  bgtr(p) = 0._r8
               else
                  ! the background transfer rate is calculated as the rate that would result
                  ! in complete turnover of the storage pools in one year at steady state,
                  ! once lgsf has reached 1.0 (after 730 days active).
               
                  bgtr(p) = (1._r8/(dayspyr*secspday))*lgsf(p)
               
                  ! set carbon fluxes for shifting storage pools to transfer pools
               
                  ! reduced the amount of stored carbon flowing to display pool by only counting the delta
                  ! between leafc and leafc_store in the flux. RosieF, Nov5 2015. 
                  leafc_storage_to_xfer(p)  = max(0.0_r8,(leafc_storage(p)-leafc(p))) * bgtr(p)
                  frootc_storage_to_xfer(p) = max(0.0_r8,(frootc_storage(p)-frootc(p))) * bgtr(p)
                  if (woody(ivt(p)) == 1.0_r8) then
               	  livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
               	  deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
               	  livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
               	  deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
               	  gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
                  end if
               
                  ! set nitrogen fluxes for shifting storage pools to transfer pools
                  leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
                  frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
                  if (woody(ivt(p)) == 1.0_r8) then
               	  livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
               	  deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
               	  livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
               	  deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
                  end if               
               end if

           end if ! end if semi deciduous		
   
           ! continue with other non semi-decid plants
           if (semi_decid(ivt(p)) /= 1._r8)  then
    
    
            ! onset_counter initialized to zero when .not. croplive
            ! offset_counter relevant only at time step of harvest

            onset_counter(p) = onset_counter(p) - dt

            ! enter phase 2 onset for one time step:
            ! transfer seed carbon to leaf emergence

            if (peaklai(p) >= 1) then
               hui(p) = max(hui(p),huigrain(p))
            endif

            if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p) .and. idpp(p) < mxmat(ivt(p))) then
               cphase(p) = 2._r8
               if (abs(onset_counter(p)) > 1.e-6_r8) then
                  onset_flag(p)    = 1._r8
                  onset_counter(p) = dt
                    fert_counter(p)  = ndays_on * secspday
                    if (ndays_on .gt. 0) then
                       fert(p) = (manunitro(ivt(p)) * 1000._r8 + fertnitro(p))/ fert_counter(p)
                    else
                       fert(p) = 0._r8
                    end if
               else
                  ! this ensures no re-entry to onset of phase2
                  ! b/c onset_counter(p) = onset_counter(p) - dt
                  ! at every time step

                  onset_counter(p) = dt
               end if

               ! enter harvest for one time step:
               ! - transfer live biomass to litter and to crop yield
               ! - send xsmrpool to the atmosphere
               ! if onset and harvest needed to last longer than one timestep
               ! the onset_counter would change from dt and you'd need to make
               ! changes to the offset subroutine below

            !New harvest for perennial plants to keep them alive after several harvests.
            !A new harvest flag is added. (Y.Fan 06.06.2014)
            else if (hui(p) >= gddmaturity(p) .and. idpp(p) < mxmat(ivt(p)) .and. &
                     perennial(ivt(p)) == 1) then
               !gddmaturity2(p) = huigrain2(p) + (1._r8 - grnfill(ivt(p))) * gddmaturity(p) !gdd till next harvests every year
               gddmaturity2(p) = huigrain2(p) + (1._r8 - grnfill(ivt(p))) * gddmaturity(p)/12._r8 !gdd till next harvests every month
               !(to do) add a new parameter grnfill2 and set grnfill2(ivt(p)) * gddmaturity to be around one month of GDD for monthly harvests (Y.Fan 2022.09)
               ! gddmaturity2(p) = huigrain2(p) + grnfill2(ivt(p)) * gddmaturity(p)
               !part of leaves are pruned off at each harvest (at one-time step)
               if (hui(p) >= gddmaturity2(p)) then
                  !the following ensures a one-time harvest each year for perennial crops
                  !Assume a percentage 0.5 of the initial grnfill=around 1 year period between subsequent grnfills
                  !in the future, add another PFT variable grnfill2 for perennial crops
                  !huigrain2(p) = huigrain2(p) + 0.5_r8 * grnfill(ivt(p)) * gddmaturity(p) !for yearly harvests
                  huigrain2(p) = huigrain2(p) + 0.5_r8 * grnfill(ivt(p)) * gddmaturity(p)/12._r8 !for monthly harvests
                 ! huigrain2(p) = huigrain2(p) + grnfill2(ivt(p)) * gddmaturity(p) !to add a new parameter grnfill2 (Y.Fan 2022.09.13)
                  !apply fertilizer after each grainfill stage before next grainfill (Y.Fan)
                  !fert continue until the end of ndays_on because this loop is entered one time step only
                  !and fert_counter(p) = fert_counter(p) - dt (see later)
                  fert_counter(p)  = ndays_on * secspday
                  fert(p) = manunitro(ivt(p)) * 1000._r8 / fert_counter(p)
                  if (tlai(p) <= 0._r8) then ! plant never emerged or died
                     croplive(p) = .false.
                     if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                     crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                     crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                     leafc_xfer(p) = 0._r8  ! revert planting transfers
                     leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
                  else if (grainc(p) > 0._r8) then !only harvest when there is positive grainc accumulated during grainfill
                     harvest_flag(p) = 1._r8
                  end if
               else
                  harvest_flag(p) = 0._r8
                 ! continue background litterfall year around because perennials keeps alive after harvest (Y.Fan)
                  bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
               end if
   
            ! trigger full litterfall (offset) when maximum age is reached and start rotation (added by Y.Fan)
            else if (hui(p) >= gddmaturity(p) .and. idpp(p) >= mxmat(ivt(p)) .and. &
                     perennial(ivt(p)) == 1) then
               if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
               croplive(p) = .false.     !
               if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                  offset_flag(p) = 1._r8
                  offset_counter(p) = dt
               else                      ! plant never emerged from the ground
                  crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                  leafc_xfer(p) = 0._r8  ! revert planting transfers
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
               end if

            else if (hui(p) >= gddmaturity(p) .or. idpp(p) >= mxmat(ivt(p))) then
               if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
               croplive(p) = .false.     ! no re-entry in greater if-block
               cphase(p) = 4._r8
               if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                  offset_flag(p) = 1._r8
                  offset_counter(p) = dt
               else                      ! plant never emerged from the ground
                  ! Revert planting transfers; this will replenish the crop seed deficit.
                  ! We subtract from any existing value in crop_seedc_to_leaf /
                  ! crop_seedn_to_leaf in the unlikely event that we enter this block of
                  ! code in the same time step where the planting transfer originally
                  ! occurred.
                  crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                  leafc_xfer(p) = 0._r8
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
                  if (use_c13) then
                     c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
                  endif
                  if (use_c14) then
                     c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
                  endif

               end if

               ! enter phase 3 while previous criteria fail and next is true;
               ! in terms of order, phase 3 occurs before harvest, but when
               ! harvest *can* occur, we want it to have first priority.
               ! AgroIBIS uses a complex formula for lai decline.
               ! Use CN's simple formula at least as a place holder (slevis)

            else if (hui(p) >= huigrain(p)) then
               cphase(p) = 3._r8
               bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
            end if

            ! continue fertilizer application while in phase 2;
            ! assumes that onset of phase 2 took one time step only

              if (fert_counter(p) <= 0._r8) then
                 fert(p) = 0._r8
              else ! continue same fert application every timestep
                 fert_counter(p) = fert_counter(p) - dtrad
              end if

	   end if ! end other non-semi decid plants

         else   ! crop not live
            ! next 2 lines conserve mass if leaf*_xfer > 0 due to interpinic.
            ! We subtract from any existing value in crop_seedc_to_leaf /
            ! crop_seedn_to_leaf in the unlikely event that we enter this block of
            ! code in the same time step where the planting transfer originally
            ! occurred.
    
            crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
            crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
            onset_counter(p) = 0._r8
            leafc_xfer(p) = 0._r8
            leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
   
            !Ashehad added this (double check if it is needed, YFan) 
            if (woody(ivt(p)) == 1._r8) then
               deadstemc_xfer(p) = 0._r8
               deadstemn_xfer(p) = deadstemc_xfer(p) / deadwdcn(ivt(p))
            end if
    
            if (use_c13) then
               c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif
            if (use_c14) then
               c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif
         end if ! croplive

      end do ! prognostic crops loop

    end associate

  end subroutine CropPhenology
  
  !-----------------------------------------------------------------------
  subroutine CropPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CropPhenology. Must be called after time-manager is
    ! initialized, and after pftcon file is read in.
    !
    ! !USES:
    use pftconMod       , only: npcropmin, npcropmax
    use clm_time_manager, only: get_calday
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! LOCAL VARAIBLES:
    integer           :: p,g,n,i                     ! indices
    !------------------------------------------------------------------------

    allocate( inhemi(bounds%begp:bounds%endp) )

    allocate( minplantjday(0:maxveg,inSH)) ! minimum planting julian day
    allocate( maxplantjday(0:maxveg,inSH)) ! minimum planting julian day

    ! Julian day for the start of the year (mid-winter)
    jdayyrstart(inNH) =   1
    jdayyrstart(inSH) = 182

    ! Convert planting dates into julian day
    minplantjday(:,:) = huge(1)
    maxplantjday(:,:) = huge(1)
    do n = npcropmin, npcropmax
       if (pftcon%is_pft_known_to_model(n)) then
          minplantjday(n, inNH) = int( get_calday( pftcon%mnNHplantdate(n), 0 ) )
          maxplantjday(n, inNH) = int( get_calday( pftcon%mxNHplantdate(n), 0 ) )

          minplantjday(n, inSH) = int( get_calday( pftcon%mnSHplantdate(n), 0 ) )
          maxplantjday(n, inSH) = int( get_calday( pftcon%mxSHplantdate(n), 0 ) )
       end if
    end do

    ! Figure out what hemisphere each PATCH is in
    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)
       ! Northern hemisphere
       if ( grc%latdeg(g) > 0.0_r8 )then
          inhemi(p) = inNH
       else
          inhemi(p) = inSH
       end if
    end do

    !
    ! Constants for Crop vernalization
    !
    ! photoperiod factor calculation
    ! genetic constant - can be modified

    p1d = 0.004_r8  ! average for genotypes from Ritchey, 1991.
    ! Modeling plant & soil systems: Wheat phasic developmt
    p1v = 0.003_r8  ! average for genotypes from Ritchey, 1991.

    hti   = 1._r8
    tbase = 0._r8

  end subroutine CropPhenologyInit

  !-----------------------------------------------------------------------
  subroutine vernalization(p, &
       canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, crop_inst)
    !
    ! !DESCRIPTION:
    !
    ! * * * only call for winter temperate cereal * * *
    !
    ! subroutine calculates vernalization and photoperiod effects on
    ! gdd accumulation in winter temperate cereal varieties. Thermal time accumulation
    ! is reduced in 1st period until plant is fully vernalized. During this
    ! time of emergence to spikelet formation, photoperiod can also have a
    ! drastic effect on plant development.
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: p    ! PATCH index running over
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(cnveg_state_type) , intent(inout) :: cnveg_state_inst
    type(crop_type)        , intent(inout) :: crop_inst
    !
    ! LOCAL VARAIBLES:
    real(r8) tcrown                     ! ?
    real(r8) vd, vd1, vd2               ! vernalization dependence
    real(r8) tkil                       ! Freeze kill threshold
    integer  c,g                        ! indices
    !------------------------------------------------------------------------

    associate(                                               & 
         tlai        => canopystate_inst%tlai_patch        , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     

         t_ref2m     => temperature_inst%t_ref2m_patch     , & ! Input:  [real(r8) (:) ]  2 m height surface air temperature (K)            
         t_ref2m_min => temperature_inst%t_ref2m_min_patch , & ! Input:  [real(r8) (:) ] daily minimum of average 2 m height surface air temperature (K)
         t_ref2m_max => temperature_inst%t_ref2m_max_patch , & ! Input:  [real(r8) (:) ] daily maximum of average 2 m height surface air temperature (K)

         snow_depth  => waterdiagnosticbulk_inst%snow_depth_col     , & ! Input:  [real(r8) (:) ]  snow height (m)                                   

         hdidx       => cnveg_state_inst%hdidx_patch       , & ! Output: [real(r8) (:) ]  cold hardening index?                             
         cumvd       => cnveg_state_inst%cumvd_patch       , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?             
         gddmaturity => cnveg_state_inst%gddmaturity_patch , & ! Output: [real(r8) (:) ]  gdd needed to harvest                             
         huigrain    => cnveg_state_inst%huigrain_patch    , & ! Output: [real(r8) (:) ]  heat unit index needed to reach vegetative maturity

         vf          => crop_inst%vf_patch                   & ! Output: [real(r8) (:) ]  vernalization factor for cereal
         )

      c = patch%column(p)

      ! for all equations - temperatures must be in degrees (C)
      ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
      ! snow depth in centimeters
      
      if (t_ref2m(p) < tfrz) then !slevis: t_ref2m inst of td=daily avg (K)
         tcrown = 2._r8 + (t_ref2m(p) - tfrz) * (0.4_r8 + 0.0018_r8 * &
              (min(snow_depth(c)*100._r8, 15._r8) - 15._r8)**2)
      else !slevis: snow_depth inst of adsnod=daily average (m)
         tcrown = t_ref2m(p) - tfrz
      end if

      ! vernalization factor calculation
      ! if vf(p) = 1.  then plant is fully vernalized - and thermal time
      ! accumulation in phase 1 will be unaffected
      ! refers to gddtsoi & gddplant, defined in the accumulation routines (slevis)
      ! reset vf, cumvd, and hdidx to 0 at planting of crop (slevis)

      if (t_ref2m_max(p) > tfrz) then
         if (t_ref2m_min(p) <= tfrz+15._r8) then
            vd1      = 1.4_r8 - 0.0778_r8 * tcrown
            vd2      = 0.5_r8 + 13.44_r8 / ((t_ref2m_max(p)-t_ref2m_min(p)+3._r8)**2) * tcrown
            vd       = max(0._r8, min(1._r8, vd1, vd2))
            cumvd(p) = cumvd(p) + vd
         end if

         if (cumvd(p) < 10._r8 .and. t_ref2m_max(p) > tfrz+30._r8) then
            cumvd(p) = cumvd(p) - 0.5_r8 * (t_ref2m_max(p) - tfrz - 30._r8)
         end if
         cumvd(p) = max(0._r8, cumvd(p))       ! must be > 0

         vf(p) = 1._r8 - p1v * (50._r8 - cumvd(p))
         vf(p) = max(0._r8, min(vf(p), 1._r8)) ! must be between 0 - 1
      end if

      ! calculate cold hardening of plant
      ! determines for winter cereal varieties whether the plant has completed
      ! a period of cold hardening to protect it from freezing temperatures. If
      ! not, then exposure could result in death or killing of plants.

      ! there are two distinct phases of hardening

      if (t_ref2m_min(p) <= tfrz-3._r8 .or. hdidx(p) /= 0._r8) then
         if (hdidx(p) >= hti) then   ! done with phase 1
            hdidx(p) = hdidx(p) + 0.083_r8
            hdidx(p) = min(hdidx(p), hti*2._r8)
         end if

         if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
            hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            hdidx(p) = max(0._r8, hdidx(p))
         end if

      else if (tcrown >= tbase-1._r8) then
         if (tcrown <= tbase+8._r8) then
            hdidx(p) = hdidx(p) + 0.1_r8 - (tcrown-tbase+3.5_r8)**2 / 506._r8
            if (hdidx(p) >= hti .and. tcrown <= tbase + 0._r8) then
               hdidx(p) = hdidx(p) + 0.083_r8
               hdidx(p) = min(hdidx(p), hti*2._r8)
            end if
         end if

         if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
            hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            hdidx(p) = max(0._r8, hdidx(p))
         end if
      end if

      ! calculate what the cereal killing temperature
      ! there is a linear inverse relationship between
      ! hardening of the plant and the killing temperature or
      ! threshold that the plant can withstand
      ! when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C

      ! will have to develop some type of relationship that reduces LAI and
      ! biomass pools in response to cold damaged crop

      if (t_ref2m_min(p) <= tfrz - 6._r8) then
         tkil = (tbase - 6._r8) - 6._r8 * hdidx(p)
         if (tkil >= tcrown) then
            if ((0.95_r8 - 0.02_r8 * (tcrown - tkil)**2) >= 0.02_r8) then
               write (iulog,*)  'crop damaged by cold temperatures at p,c =', p,c
            else if (tlai(p) > 0._r8) then ! slevis: kill if past phase1
               gddmaturity(p) = 0._r8      !         by forcing through
               huigrain(p)    = 0._r8      !         harvest
               write (iulog,*)  '95% of crop killed by cold temperatures at p,c =', p,c
            end if
         end if
      end if

    end associate 

  end subroutine vernalization

  !-----------------------------------------------------------------------
  subroutine CNOnsetGrowth (num_soilp, filter_soilp, &
       cnveg_state_inst, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of stored C and N from transfer pools to display
    ! pools during the phenological onset period.
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)             , intent(in)    :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter patch index
    real(r8):: t1           ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                             & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:) ]  patch vegetation type                                

         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Input:  [real(r8)  (:) ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Input:  [real(r8)  (:) ]  onset days counter                                
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Input:  [real(r8)  (:) ]  background transfer growth rate (1/s)             
         phytomer                            =>    pftcon%phytomer                                             , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
	 perennial                           =>    pftcon%perennial                                             , & ! Input:  [integer (:)]   
         pleafc_xfer                         =>    cnveg_carbonstate_inst%pleafc_xfer_patch                          , & ! Output: [real(r8) (:,:)]  (gC/m2) phytomer leaf C transfer
         pleafc_xfer_to_pleafc               =>    cnveg_carbonflux_inst%pleafc_xfer_to_pleafc_patch                 , & ! Output: [real(r8) (:,:)]  added by Y.Fan
         pleafn_xfer                         =>    cnveg_nitrogenstate_inst%pleafn_xfer_patch                        , & ! Output: [real(r8) (:,:)]  (gN/m2) phytomer leaf N transfer
         pleafn_xfer_to_pleafn               =>    cnveg_nitrogenflux_inst%pleafn_xfer_to_pleafn_patch               , & ! Output: [real(r8) (:,:)]  added by Y.Fan
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Input:  [real(r8)  (:) ]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Input:  [real(r8)  (:) ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Input:  [real(r8)  (:) ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead coarse root C transfer               
         
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Input:  [real(r8)  (:) ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Input:  [real(r8)  (:) ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Input:  [real(r8)  (:) ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Input:  [real(r8)  (:) ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead coarse root N transfer               
         
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:) ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:) ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:) ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:) ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:) ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:) ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:) ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:) ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:) ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:) ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:) ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch   & ! Output:  [real(r8) (:) ]                                                    
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes during onset period
         if (onset_flag(p) == 1._r8) then

            ! The transfer rate is a linearly decreasing function of time,
            ! going to zero on the last timestep of the onset period

            if (onset_counter(p) == dt) then
               t1 = 1.0_r8 / dt
            else
               t1 = 2.0_r8 / (onset_counter(p))
            end if
            leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)
            frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)
            leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)
            frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
               deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
               livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
               deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
               livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
               deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
               livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
               deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)
            end if

            !phytomer-based structure (Y.Fan):
            !For oil palm, use initial C transfer to represent transplanting
            !from nursery to field (incl. leafc_xfer and livestemc_xfer)
            !depending on the amount of initial leafc_xfer, first phytomers get
            !allocation not exceeding the maximum per leaf
            if (phytomer(ivt(p)) > 0) then
               pleafc_xfer_to_pleafc(p,:) = t1 * pleafc_xfer(p,:)
               pleafn_xfer_to_pleafn(p,:) = t1 * pleafn_xfer(p,:)
               !below two lines added for oil palm, so that it can avoid using the
               !woody flag (Y.Fan 2022)
               livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
               livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
            end if

         end if ! end if onset period

         ! calculate the background rate of transfer growth (used for stress
         ! deciduous algorithm). In this case, all of the mass in the transfer
         ! pools should be moved to displayed growth in each timestep.

         !for palm phytomer post-expansion period: 
         !(bgtr only used for leaves) move transfer pools to displayed growth
         !the leaf xfer pool is refilled and depleted at every timestep (Y.Fan)
         if (bgtr(p) > 0._r8 .and. phytomer(ivt(p)) > 0) then
            pleafc_xfer_to_pleafc(p,:) = pleafc_xfer(p,:) / dt
            pleafn_xfer_to_pleafn(p,:) = pleafn_xfer(p,:) / dt
            leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
            leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
         else if (bgtr(p) > 0._r8) then
         !if (bgtr(p) > 0._r8) then
            leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
            frootc_xfer_to_frootc(p) = frootc_xfer(p) / dt
            leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
            frootn_xfer_to_frootn(p) = frootn_xfer(p) / dt
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
               deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
               livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
               deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
               livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
               deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
               livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
               deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
            end if
         end if ! end if bgtr

      end do ! end patch loop

    end associate

  end subroutine CNOnsetGrowth

  !-----------------------------------------------------------------------
  subroutine CNOffsetLitterfall (num_soilp, filter_soilp, &
       num_soilc, filter_soilc, soilbiogeochem_state_inst,&
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, crop_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools during the phenological offset period.
    !
    ! !USES:
    use pftconMod        , only : npcropmin !, mxlivenp, pprodharv10
    use CNSharedParamsMod, only : use_fun
    use clm_varctl       , only : CNratio_floating   
    use dynHarvestMod , only: CNHarvestPftToColumn !for woody crop types 
    !
    ! !ARGUMENTS:
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
    type(crop_type)               , intent(inout) :: crop_inst
    type(cnveg_carbonstate_type)  , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type), intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, n, m         ! indices
    integer :: fp           ! lake filter patch index
    real(r8):: t1           ! temporary variable
    real(r8):: denom        ! temporary variable for divisor
    real(r8) :: ntovr_leaf  
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    !-----------------------------------------------------------------------

    associate(                                                                           & 
         ivt                   =>    patch%itype                                       , & ! Input:  [integer  (:) ]  patch vegetation type                                
         leafcn                =>    pftcon%leafcn                                     , & ! Input:  leaf C:N (gC/gN)                                  
         lflitcn               =>    pftcon%lflitcn                                    , & ! Input:  leaf litter C:N (gC/gN)                           
         frootcn               =>    pftcon%frootcn                                    , & ! Input:  fine root C:N (gC/gN)                             
         graincn               =>    pftcon%graincn                                    , & ! Input:  grain C:N (gC/gN)
	 woody                 =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)                                 
         phytomer              =>    pftcon%phytomer                                   , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
         perennial             =>    pftcon%perennial                              , & ! Input:  [integer (:)]  binary flag for perennial crop phenology (1=perennial, 0=not perennial) (added by Y.Fan)
         mxlivenp              =>    pftcon%mxlivenp                                   , & ! Input:  [integer (:)]   maximum number of live phytomers 
         pprodharv10           =>    pftcon%pprodharv10                             , & ! Input:  harvest mortality proportion of deadstem to 10-yr pool
         leaf_long             =>    pftcon%leaf_long                            , & ! Input:  [real(r8) (:)]  leaf longevity (yrs)
         prune                 =>    crop_inst%prune_patch                                  , & ! Input:  [real(r8) (:)]  flag for pruning
         np                    =>    crop_inst%np_patch                                      , & ! Input:  [integer (:)]   total number of phytomers having appeared so far
         rankp                 =>    crop_inst%rankp_patch                                   , & ! Input:  [integer (:,:)]  rank of phytomers from 1=youngest to np=oldest and 0=dead
         livep                 =>    crop_inst%livep_patch                                   , & ! Input:  [real(r8) (:,:)]  Flag, 1 if this phytomer is alive
         lfdays                =>    crop_inst%lfdays_patch                                  , & ! Input:  [integer (:,:)]  days past leaf emergence for each phytomer
         pgrainc               =>    cnveg_carbonstate_inst%pgrainc_patch                                 , & ! Input:  [real(r8) (:,:)]  (gC/m2) phytomer grain C
         pgrainn               =>    cnveg_nitrogenstate_inst%pgrainn_patch                                 , & ! Input:  [real(r8) (:,:)]  (gN/m2) phytomer grain N
         pleafc                =>    cnveg_carbonstate_inst%pleafc_patch                                  , & ! Input:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C
         pleafn                =>    cnveg_nitrogenstate_inst%pleafn_patch                                  , & ! Input:  [real(r8) (:,:)]  (gN/m2) phytomer leaf N
         pleafc_storage        =>    cnveg_carbonstate_inst%pleafc_storage_patch                          , & ! InOut:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C stroage (unexpanded)
         pleafn_storage        =>    cnveg_nitrogenstate_inst%pleafn_storage_patch                          , & ! InOut:  [real(r8) (:,:)]  (gN/m2) phytomer leaf N storage
         pleafc_xfer           =>    cnveg_carbonstate_inst%pleafc_xfer_patch                             , & ! InOut:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C transfer
         pleafn_xfer           =>    cnveg_nitrogenstate_inst%pleafn_xfer_patch                             , & ! InOut:  [real(r8) (:,:)]
         pleafc_to_litter      =>    cnveg_carbonflux_inst%pleafc_to_litter_patch                        , & ! InOut:  [real(r8) (:,:)]  phytomer leaf C litterfall (gC/m2/s)
         pleafc_storage_to_litter      =>    cnveg_carbonflux_inst%pleafc_storage_to_litter_patch                , & ! InOut:  [real(r8) (:,:)]  phytomer leaf C storage to litter at final rotation (gC/m2/s)
         pleafc_xfer_to_litter         =>    cnveg_carbonflux_inst%pleafc_xfer_to_litter_patch                   , & ! InOut:  [real(r8) (:,:)]
         pgrainc_to_food               =>    cnveg_carbonflux_inst%pgrainc_to_food_patch                         , & ! InOut:  [real(r8) (:,:)]  phytomer grain C to food (gC/m2/s)
         npool_to_pgrainn              =>    cnveg_nitrogenflux_inst%npool_to_pgrainn_patch                        , & ! Input:  [real(r8) (:,:)]  allocation to phytomer grain N (gN/m2/s)^M
         npool_to_pleafn               =>    cnveg_nitrogenflux_inst%npool_to_pleafn_patch                         , & ! Input:  [real(r8) (:,:)]  allocation to phytomer leaf N (gN/m2/s)
         npool_to_pleafn_storage       =>    cnveg_nitrogenflux_inst%npool_to_pleafn_storage_patch                 , & ! Input:  [real(r8) (:,:)]
         pleafn_to_litter              =>    cnveg_nitrogenflux_inst%pleafn_to_litter_patch                        , & ! InOut:  [real(r8) (:,:)]  phytomer leaf N litterfall (gN/m2/s)
         pleafn_storage_to_litter      =>    cnveg_nitrogenflux_inst%pleafn_storage_to_litter_patch                , & ! InOut:  [real(r8) (:,:)]  phytomer leaf N storage to litter at final rotation (gN/m2/s)
         pleafn_xfer_to_litter         =>    cnveg_nitrogenflux_inst%pleafn_xfer_to_litter_patch                   , & ! InOut:  [real(r8) (:,:)]
         pgrainn_to_food               =>    cnveg_nitrogenflux_inst%pgrainn_to_food_patch                         , & ! InOut:  [real(r8) (:,:)]  phytomer grain N to food (gN/m2/s)
         pleafc_xfer_to_pleafc =>    cnveg_carbonflux_inst%pleafc_xfer_to_pleafc_patch       , & ! Output: [real(r8) (:,:)]  added by Y.Fan
         pleafn_xfer_to_pleafn =>    cnveg_nitrogenflux_inst%pleafn_xfer_to_pleafn_patch     , & ! Output: [real(r8) (:,:)]  added by Y.Fan
         cpool_to_pleafc       =>    cnveg_carbonflux_inst%cpool_to_pleafc_patch                         , & ! Input:  [real(r8) (:,:)]  allocation to phytomer leaf C (gC/m2/s)
         cpool_to_pleafc_storage =>    cnveg_carbonflux_inst%cpool_to_pleafc_storage_patch                 , & ! Input:  [real(r8) (:,:)]^M
         cpool_to_pgrainc      =>    cnveg_carbonflux_inst%cpool_to_pgrainc_patch                        , & ! Input:  [real(r8) (:,:)]  allocation to phytomer grain C (gC/m2/s)
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch     , & ! Input: [real(r8) (:)]  (gC/m2) leaf C storage
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch        , & ! Input: [real(r8) (:)]  (gC/m2) leaf C transfer
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch     , & ! Input: [real(r8) (:)]  (gC/m2) growth respiration storage
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch        , & ! Input: [real(r8) (:)]  (gC/m2) growth respiration transfer
         xsmrpool                            =>    cnveg_carbonstate_inst%xsmrpool_patch        , & ! Input: [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand
         hrv_leafc_to_litter                 =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_frootc_to_litter                =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_livestemc_to_litter             =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch , & ! Output: [real(r8) (:)]
         wood_harvestc                       =>    cnveg_carbonflux_inst%wood_harvestc_patch            , & ! Output: [real(r8) (:)]
         hrv_livecrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_patch , & ! Output: [real(r8) (:)]
         hrv_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch , & ! Output: [real(r8) (:)]^M
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input: [real(r8) (:)]  (gN/m2) leaf N storage
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch         , & ! Input: [real(r8) (:)]  (gN/m2) leaf N transfer
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch        , & ! Input: [real(r8) (:)] (gN/m2) plant pool of retranslocated N
         hrv_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_frootn_to_litter                =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch , & ! Output: [real(r8) (:)]
         wood_harvestn                       =>    cnveg_nitrogenflux_inst%wood_harvestn_patch                    , & ! Output: [real(r8) (:)]
         hrv_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_retransn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch , & ! Output: [real(r8) (:)]
         hrv_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch   , & ! Output: [real(r8) (:)]
         harvest_flag          =>    crop_inst%harvest_flag_patch                     , & ! Input:  [real(r8) (:)]  harvest flag (added by Y.Fan)
         harvest_counter       =>    crop_inst%harvest_counter_patch                  , & ! Input:  [real(r8) (:)]  harvest counter to tag phytomer
         offset_flag           =>    cnveg_state_inst%offset_flag_patch                , & ! Input:  [real(r8) (:) ]  offset flag                                       
         offset_counter        =>    cnveg_state_inst%offset_counter_patch             , & ! Input:  [real(r8) (:) ]  offset days counter                               
         leafc                 =>    cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C                                    
         frootc                =>    cnveg_carbonstate_inst%frootc_patch               , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C                               
         grainc                =>    cnveg_carbonstate_inst%grainc_patch               , & ! Input:  [real(r8) (:) ]  (gC/m2) grain C                                   
         cropseedc_deficit     =>    cnveg_carbonstate_inst%cropseedc_deficit_patch    , & ! Input:  [real(r8) (:) ]  (gC/m2) crop seed C deficit
         livestemc             =>    cnveg_carbonstate_inst%livestemc_patch            , & ! Input:  [real(r8) (:) ]  (gC/m2) livestem C                           
         deadstemc             =>    cnveg_carbonstate_inst%deadstemc_patch            , & ! Input: [real(r8) (:)]  (gC/m2) dead stem C 
         cropseedn_deficit     =>    cnveg_nitrogenstate_inst%cropseedn_deficit_patch  , & ! Input:  [real(r8) (:) ]  (gC/m2) crop seed N deficit
         livestemn             =>    cnveg_nitrogenstate_inst%livestemn_patch          , & ! Input:  [real(r8) (:) ]  (gN/m2) livestem N
         deadstemn             =>    cnveg_nitrogenstate_inst%deadstemn_patch          , & ! Input: [real(r8) (:)]  (gN/m2) dead stem N
         cpool_to_grainc       =>    cnveg_carbonflux_inst%cpool_to_grainc_patch       , & ! Input:  [real(r8) (:) ]  allocation to grain C (gC/m2/s)                   
         npool_to_grainn       =>    cnveg_nitrogenflux_inst%npool_to_grainn_patch     , & ! Input: [real(r8) (:)  ]  allocation to grain N (gN/m2/s)
         grainn                =>    cnveg_nitrogenstate_inst%grainn_patch             , & ! Input: [real(r8) (:)  ]  (kgN/m2) grain N
         cpool_to_livestemc    =>    cnveg_carbonflux_inst%cpool_to_livestemc_patch    , & ! Input:  [real(r8) (:) ]  allocation to live stem C (gC/m2/s)               
         cpool_to_leafc        =>    cnveg_carbonflux_inst%cpool_to_leafc_patch        , & ! Input:  [real(r8) (:) ]  allocation to leaf C (gC/m2/s)                    
         cpool_to_frootc       =>    cnveg_carbonflux_inst%cpool_to_frootc_patch       , & ! Input:  [real(r8) (:) ]  allocation to fine root C (gC/m2/s)               
         prev_leafc_to_litter  =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch  , & ! Output: [real(r8) (:) ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch , & ! Output: [real(r8) (:) ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_to_litter       =>    cnveg_carbonflux_inst%leafc_to_litter_patch       , & ! Output: [real(r8) (:) ]  leaf C litterfall (gC/m2/s)                       
         frootc_to_litter      =>    cnveg_carbonflux_inst%frootc_to_litter_patch      , & ! Output: [real(r8) (:) ]  fine root C litterfall (gC/m2/s)                  
         livestemc_to_litter   =>    cnveg_carbonflux_inst%livestemc_to_litter_patch   , & ! Output: [real(r8) (:) ]  live stem C litterfall (gC/m2/s)                  
         grainc_to_food        =>    cnveg_carbonflux_inst%grainc_to_food_patch        , & ! Output: [real(r8) (:) ]  grain C to food (gC/m2/s)                         
         grainc_to_seed        =>    cnveg_carbonflux_inst%grainc_to_seed_patch        , & ! Output: [real(r8) (:) ]  grain C to seed (gC/m2/s)
         leafn                 =>    cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N      
         frootn                =>    cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:) ]  (gN/m2) fine root N                        
         livestemn_to_litter   =>    cnveg_nitrogenflux_inst%livestemn_to_litter_patch , & ! Output: [real(r8) (:) ]  livestem N to litter (gN/m2/s)                    
         grainn_to_food        =>    cnveg_nitrogenflux_inst%grainn_to_food_patch      , & ! Output: [real(r8) (:) ]  grain N to food (gN/m2/s)                         
         grainn_to_seed        =>    cnveg_nitrogenflux_inst%grainn_to_seed_patch      , & ! Output: [real(r8) (:) ]  grain N to seed (gN/m2/s)
         leafn_to_litter       =>    cnveg_nitrogenflux_inst%leafn_to_litter_patch     , & ! Output: [real(r8) (:) ]  leaf N litterfall (gN/m2/s)                       
         leafn_to_retransn     =>    cnveg_nitrogenflux_inst%leafn_to_retransn_patch   , & ! Input: [real(r8) (:) ]  leaf N to retranslocated N pool (gN/m2/s)         
         free_retransn_to_npool=>    cnveg_nitrogenflux_inst%free_retransn_to_npool_patch  , & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         paid_retransn_to_npool=>    cnveg_nitrogenflux_inst%retransn_to_npool_patch, & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         frootn_to_litter      =>    cnveg_nitrogenflux_inst%frootn_to_litter_patch    , & ! Output: [real(r8) (:) ]  fine root N litterfall (gN/m2/s)                  
         leafc_to_litter_fun   =>    cnveg_carbonflux_inst%leafc_to_litter_fun_patch   , & ! Output:  [real(r8) (:) ]  leaf C litterfall used by FUN (gC/m2/s)
         leafcn_offset         =>    cnveg_state_inst%leafcn_offset_patch               & ! Output:  [real(r8) (:) ]  Leaf C:N used by FUN
         )

      ! The litterfall transfer rate starts at 0.0 and increases linearly
      ! over time, with displayed growth going to 0.0 on the last day of litterfall
      
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !incorporate a one-time harvest: annually for perennial crops; seasonally for phytomer phenology (Y.Fan)
         !leaf flux is treated according to leaf property (evergreen or deciduous) and management (prune or not)
         !grain flux goes to food
         if (harvest_flag(p) == 1._r8) then
           if (phytomer(ivt(p)) > 0) then
              t1 = 1.0_r8 / dt
              n= mat1 + int(harvest_counter(p))
              pgrainc_to_food(p,n) = t1 * pgrainc(p,n)  + cpool_to_pgrainc(p,n)
              pgrainn_to_food(p,n) = t1 * pgrainn(p,n)  + npool_to_pgrainn(p,n)
              grainc_to_food(p) = pgrainc_to_food(p,n)
              grainn_to_food(p) = pgrainn_to_food(p,n)
   
              harvest_counter(p) = harvest_counter(p) + 1
             !do pruning at the time of harvest
             !adjust mxlivenp and leaf_long/mxgdd so as to prune nearly every 6 month
             !  if (nlevcan > 1 .and. prune(p) .and. livep(p,n) == 0._r8) then
             if (prune(p)) then
               where (livep(p,:) == 0._r8)
                  pleafc_to_litter(p,:) = pleafc(p,:) /dt  !all remaining leafc moves to litter
                  pleafn_to_litter(p,:) = pleafn(p,:) /dt  !no N retranslocation when pruning
                  rankp(p,:) = 0  !update rankp only after finishing pruning and C/N update
               end where
               leafc_to_litter(p) = sum(pleafc_to_litter(p,:))
               leafn_to_litter(p) = sum(pleafn_to_litter(p,:))
             end if
   
           else if (perennial(ivt(p)) == 1) then
              t1 = 1.0_r8 / dt
              grainc_to_food(p) = t1 * grainc(p)  + cpool_to_grainc(p)
              grainn_to_food(p) = t1 * grainn(p)  + npool_to_grainn(p)
           end if
         end if
   
         ! calculate fluxes during offset period (a one time step for crops; multiple days for trees)
         ! for oil palm offset means the final harvest/rotation at every 25-30 years (Y.Fan)
         if (offset_flag(p) == 1._r8) then
            ! adding harvest/rotation stem C/N into litter pool will dramatically increase LITTERC_HR which will drain out SMINN very soon!
            ! adding deadstem (accounting for 90% of total pft C after 25 years) to product pool instead of litter pool to avoid using up soil mineral N
            ! otherwise have to use fertilization continuously to maintain sufficient SMINN!  (Y.Fan 2016)
            if (phytomer(ivt(p)) > 0) then
               t1 = 1.0_r8 / dt
   
               ! clear-cut carbon fluxes, remove all displayed/storage/transfer pools
               ! the following hrv fluxes will not be added to litter pool unless using function NHarvestPftToColumn, otherwise will have C balance error
               ! summarize into litter in CNLitterToColumn mod by calling NHarvestPftToColumn
               ! displayed pools (similar to wood harvest fluxes in dynHarvestMod)
               hrv_leafc_to_litter(p)               = leafc(p)          * t1
               hrv_leafn_to_litter(p)               = leafn(p)          * t1
               hrv_frootc_to_litter(p)              = frootc(p)         * t1
               hrv_frootn_to_litter(p)              = frootn(p)         * t1
               hrv_livestemc_to_litter(p)           = livestemc(p)      * t1
               hrv_livestemn_to_litter(p)           = livestemn(p)      * t1
   
               ! storage pools
               hrv_leafc_storage_to_litter(p)       = leafc_storage(p)  * t1
               hrv_leafn_storage_to_litter(p)       = leafn_storage(p)  * t1
               hrv_gresp_storage_to_litter(p)       = gresp_storage(p)  * t1
   
               ! transfer pools
               hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)     * t1
               hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)     * t1
               hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)     * t1
   
               ! N retransn pool also move to litter
               hrv_retransn_to_litter(p)            = retransn(p)       * t1
   
               ! updated in CNSummary and included in C balance check
               hrv_xsmrpool_to_atm(p)               = xsmrpool(p)       * t1
   
               !wood product pool (prod10c/prod100c updates even without fpftdyn), included in totcolc for C balance
               !this pool is large and cannot be added into litter, otherwise will drain out SMINN soon! (Y.Fan 2016)
               wood_harvestc(p)                     = deadstemc(p)     * t1 
               wood_harvestn(p)                     = deadstemn(p)     * t1 

               !summarize hrv fluxes to column level in similar way to natural
               !tree PFTs with the following function from dynHarvestMod.F90
               call CNHarvestPftToColumn(num_soilc, filter_soilc, &
                    soilbiogeochem_state_inst,cnveg_carbonflux_inst,cnveg_nitrogenflux_inst)
               !wood harvest function is designed for natural tree PFTs but can be
               !applied to perennial woody crops like oil palm (Y.Fan 2022)

               ! clear-up phytomer pools too
               pgrainc_to_food(p,:)           = t1 * pgrainc(p,:)  + cpool_to_pgrainc(p,:)
               pgrainn_to_food(p,:)           = t1 * pgrainn(p,:)  + npool_to_pgrainn(p,:)
               pleafc_to_litter(p,:)          = t1 * pleafc(p,:) + cpool_to_pleafc(p,:)
               pleafn_to_litter(p,:)          = t1 * pleafn(p,:) + npool_to_pleafn(p,:)
               pleafc_storage_to_litter(p,:)  = t1 * pleafc_storage(p,:) + cpool_to_pleafc_storage(p,:)
               pleafn_storage_to_litter(p,:)  = t1 * pleafn_storage(p,:) + npool_to_pleafn_storage(p,:)
               pleafc_xfer_to_litter(p,:)     = t1 * pleafc_xfer(p,:)
               pleafn_xfer_to_litter(p,:)     = t1 * pleafn_xfer(p,:)


               !at final rotation remove all ripe or inmature fruits (YFan,2022)
               !pgrainc_to_food now includes inmature fruits
               grainc_to_seed(p) = t1 * min(-cropseedc_deficit(p), grainc(p))
               grainn_to_seed(p) = t1 * min(-cropseedn_deficit(p), grainn(p))

               grainc_to_food(p) = sum(pgrainc_to_food(p,:)) - grainc_to_seed(p)
               grainn_to_food(p) = sum(pgrainn_to_food(p,:)) - grainn_to_seed(p)
               !grain C/N all saved in food export pool (cropprodc)
               !grainc_to_food(p) = t1 *grainc(p)  + cpool_to_grainc(p) - grainc_to_seed(p)
               !grainn_to_food(p) = t1 *grainn(p)  + npool_to_grainn(p) - grainn_to_seed(p)

            !for annual crops and trees
            else if (offset_counter(p) == dt) then
            !if (offset_counter(p) == dt) then
               t1 = 1.0_r8 / dt
               leafc_to_litter(p)  = t1 * leafc(p)  + cpool_to_leafc(p)
               frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
               ! this assumes that offset_counter == dt for crops
               ! if this were ever changed, we'd need to add code to the "else"
               if (ivt(p) >= npcropmin) then
                  ! Replenish the seed deficits from grain, if there is enough
                  ! available grain. (If there is not enough available grain, the seed
                  ! deficits will accumulate until there is eventually enough grain to
                  ! replenish them.)
                  grainc_to_seed(p) = t1 * min(-cropseedc_deficit(p), grainc(p))
                  grainn_to_seed(p) = t1 * min(-cropseedn_deficit(p), grainn(p))
                  ! Send the remaining grain to the food product pool
                  grainc_to_food(p) = t1 * grainc(p)  + cpool_to_grainc(p) - grainc_to_seed(p)
                  grainn_to_food(p) = t1 * grainn(p)  + npool_to_grainn(p) - grainn_to_seed(p)

                  livestemc_to_litter(p) = t1 * livestemc(p)  + cpool_to_livestemc(p)		  
		  
		  ! Ashehad added this
		  if (woody(ivt(p)) == 1.0_r8) then
		      wood_harvestc(p)                     = deadstemc(p)     * t1 
                      wood_harvestn(p)                     = deadstemn(p)     * t1 
		  end if
		  
               end if	      
            else
               t1 = dt * 2.0_r8 / (offset_counter(p) * offset_counter(p))
               leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
               frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))

            end if
            
            if ( use_fun ) then
               if(leafc_to_litter(p)*dt.gt.leafc(p))then
                   leafc_to_litter(p) = leafc(p)/dt + cpool_to_leafc(p)
               endif
               if(frootc_to_litter(p)*dt.gt.frootc(p))then
                   frootc_to_litter(p) = frootc(p)/dt + cpool_to_frootc(p)
               endif
            end if
            
            
            if ( use_fun ) then
               leafc_to_litter_fun(p)      =  leafc_to_litter(p)
               leafn_to_retransn(p)        =  paid_retransn_to_npool(p) + free_retransn_to_npool(p)
               if (leafn(p).gt.0._r8) then
                  if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
                      leafcn_offset(p)     =  leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
                  else
                      leafcn_offset(p)     =  leafc(p)/leafn(p)
                  end if
               else
                  leafcn_offset(p)         =  leafcn(ivt(p))
               end if
               leafn_to_litter(p)          =  leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
               leafn_to_litter(p)          =  max(leafn_to_litter(p),0._r8)
               
               denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
               if ( denom /= 0.0_r8 ) then
                  fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
               else if ( leafn_to_litter(p) == 0.0_r8 ) then
                  fr_leafn_to_litter =  0.0_r8
               else
                  fr_leafn_to_litter =  1.0_r8
               end if

            else
               if (CNratio_floating .eqv. .true.) then    
                  fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
               end if
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
               leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

            end if    

            ! calculate fine root N litterfall (no retranslocation of fine root N)
            frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))
            
            if (CNratio_floating .eqv. .true.) then    
               if (leafc(p) == 0.0_r8) then    
                  ntovr_leaf = 0.0_r8    
               else    
                  ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))   
               end if   
           
               leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
               leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
               if (frootc(p) == 0.0_r8) then    
                   frootn_to_litter(p) = 0.0_r8    
                else    
                   frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))   
                end if   
            end if  
            
            if ( use_fun ) then
               if(frootn_to_litter(p)*dt.gt.frootn(p))then
                   frootn_to_litter(p) = frootn(p)/dt
               endif    
            end if
            
            if (ivt(p) >= npcropmin) then
               ! NOTE(slevis, 2014-12) results in -ve livestemn and -ve totpftn
               !X! livestemn_to_litter(p) = livestemc_to_litter(p) / livewdcn(ivt(p))
               ! NOTE(slevis, 2014-12) Beth Drewniak suggested this instead
               livestemn_to_litter(p) = livestemn(p) / dt
            end if

            ! save the current litterfall fluxes
            prev_leafc_to_litter(p)  = leafc_to_litter(p)
            prev_frootc_to_litter(p) = frootc_to_litter(p)

         end if ! end if offset period

      end do ! end patch loop

    end associate 

  end subroutine CNOffsetLitterfall

  !-----------------------------------------------------------------------
  subroutine CNBackgroundLitterfall (num_soilp, filter_soilp, &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, crop_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools as the result of background litter fall.
    !
    ! !USES:
    use CNSharedParamsMod   , only : use_fun
    use clm_varctl          , only : CNratio_floating    
    ! !ARGUMENTS:
    implicit none
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
    type(crop_type)                   , intent(inout) :: crop_inst
    type(cnveg_carbonstate_type)  , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type), intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter patch index
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    real(r8) :: ntovr_leaf  
    real(r8) :: denom       
    !-----------------------------------------------------------------------
    associate(                                                                     & 
         ivt               =>    patch%itype                                     , & ! Input:  [integer  (:) ]  patch vegetation type                                
         leafcn            =>    pftcon%leafcn                                   , & ! Input:  leaf C:N (gC/gN)                                  
         lflitcn           =>    pftcon%lflitcn                                  , & ! Input:  leaf litter C:N (gC/gN)                           
         frootcn           =>    pftcon%frootcn                                  , & ! Input:  fine root C:N (gC/gN)                             
         phytomer          =>    pftcon%phytomer                                 , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
         leaf_long         =>    pftcon%leaf_long                            , & ! Input:  [real(r8) (:)]  leaf longevity (yrs)
         bglfr             =>    cnveg_state_inst%bglfr_patch                    , & ! Input:  [real(r8) (:) ]  background litterfall rate (1/s)                  
         offset_flag       =>    cnveg_state_inst%offset_flag_patch                      , & ! Input:  [real(r8) (:)]  offset flag for perennial crop ratation (Y.Fan)
         bglfr_p           =>    crop_inst%bglfr_p_patch                        , & ! InOut:  [real(r8) (:,:)]  background litterfall rate for each phytomer (Y.Fan)
         prune                 =>    crop_inst%prune_patch                            , & ! Input:  [real(r8) (:)]  flag for pruning
         np                    =>    crop_inst%np_patch                                      , & ! Input:  [integer (:)]   total number of phytomers having appeared so far
         rankp                 =>    crop_inst%rankp_patch                                   , & ! Input:  [integer (:,:)]  rank of phytomers from 1=youngest to np=oldest and 0=dead
         livep                 =>    crop_inst%livep_patch                                   , & ! Input:  [real(r8) (:,:)]  Flag, 1 if this phytomer is alive
         lfdays                =>    crop_inst%lfdays_patch                                  , & ! Input:  [integer (:,:)]  days past leaf emergence for each phytomer
         pleafc                =>    cnveg_carbonstate_inst%pleafc_patch                     , & ! Input:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C
         pleafn                =>    cnveg_nitrogenstate_inst%pleafn_patch                   , & ! Input:  [real(r8) (:,:)]  (gN/m2) phytomer leaf N
         pleafc_to_litter      =>    cnveg_carbonflux_inst%pleafc_to_litter_patch            , & ! InOut:  [real(r8) (:,:)]  phytomer leaf C litterfall (gC/m2/s)
         pleafn_to_litter      =>    cnveg_nitrogenflux_inst%pleafn_to_litter_patch          , & ! InOut:  [real(r8) (:,:)]  phytomer leaf N litterfall (gN/m2/s)
         pleafn_to_retransn    =>    cnveg_nitrogenflux_inst%pleafn_to_retransn_patch        , & ! Input:  [real(r8) (:,:)] phytomer leaf N restranslocation (gN/m2/s)
         !leafc_senescent       =>    cnveg_carbonstate_inst%leafc_senescent_patch            , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C saved for pruning (added by Y.Fan)
         !leafn_senescent       =>    cnveg_nitrogenstate_inst%leafn_senescent_patch          , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N saved for pruning (added by Y.Fan)

         leafc             =>    cnveg_carbonstate_inst%leafc_patch              , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C                                    
         frootc            =>    cnveg_carbonstate_inst%frootc_patch             , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C                               
         
         leafc_to_litter   =>    cnveg_carbonflux_inst%leafc_to_litter_patch     , & ! Output: [real(r8) (:) ]                                                    
         frootc_to_litter  =>    cnveg_carbonflux_inst%frootc_to_litter_patch    , & ! Output: [real(r8) (:) ]                                                    
         leafn             =>    cnveg_nitrogenstate_inst%leafn_patch            , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N  
         frootn            =>    cnveg_nitrogenstate_inst%frootn_patch           , & ! Input:  [real(r8) (:) ]  (gN/m2) fine root N 
         leafn_to_litter   =>    cnveg_nitrogenflux_inst%leafn_to_litter_patch   , & ! Output: [real(r8) (:) ]                                                    
         leafn_to_retransn =>    cnveg_nitrogenflux_inst%leafn_to_retransn_patch , & ! Output: [real(r8) (:) ]                                                    
         frootn_to_litter  =>    cnveg_nitrogenflux_inst%frootn_to_litter_patch  , & ! Output: [real(r8) (:) ]                                                    
         leafc_to_litter_fun   => cnveg_carbonflux_inst%leafc_to_litter_fun_patch, & ! Output:  [real(r8) (:) ] leaf C litterfall used by FUN (gC/m2/s)
         leafcn_offset         => cnveg_state_inst%leafcn_offset_patch           , & ! Output:  [real(r8) (:) ] Leaf C:N used by FUN
         free_retransn_to_npool=>    cnveg_nitrogenflux_inst%free_retransn_to_npool_patch  , & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         paid_retransn_to_npool=>    cnveg_nitrogenflux_inst%retransn_to_npool_patch   & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

        ! for palm phenology, use phytomer-specific senescence function
        ! if prune = True, the remaining CN flux is calculated in CNOffsetLitterfall
         if (.not. prune(p) .and. offset_flag(p) /= 1._r8 .and. phytomer(ivt(p)) > 0) then
            !update leaf C to reduce photosynthetic LAI of senescent phytomer
            !leaf N is also reduced via retranslocation
            where (bglfr_p(p,:) > 0._r8 .and. pleafc(p,:) > 0._r8 )
                pleafc_to_litter(p,:) = bglfr_p(p,:) * pleafc(p,:)
                pleafn_to_litter(p,:) = pleafc_to_litter(p,:) / lflitcn(ivt(p))
                pleafn_to_retransn(p,:) = bglfr_p(p,:) * pleafn(p,:) - pleafn_to_litter(p,:) !avoild using leaf CN ratio
            endwhere
!            !for dynamic N scheme, possible further N downscaling during leaf senescence
!            !additional N retranslocation to decrease leaf N concentration until C:N ratio reaches that of leaf litter
!            !LAI has more control on photosynthetic capacity than N downscaling because vcmax can only reduce a half by N removal
!                where (bglfr_p(p,:) > 0._r8 .and. pleafn(p,:) > pleafc(p,:)/lflitcn(ivt(p)) )
!                    !pleafn_to_retransn(p,:) = bglfr_p(p,:) * (pleafn(p,:) - pleafc(p,:)/lflitcn(ivt(p)))
!                    !this equation is equivalent to the above pleafn_to_retransn
!                    pleafn_to_retransn(p,:) = 1.2_r8 * bglfr_p(p,:) * (pleafn(p,:) - pleafc(p,:)/lflitcn(ivt(p)))
!                endwhere
!            end if
!            !Since SUM is a nonelemental function, it is evaluated fully for all items of a variable included in sum().
!            !do not use sum in where/elsewhere clause, use mask condition instead
!            !leafc_senescent(p) = leafc_senescent(p) + sum(pleafc_to_litter(p,:), mask=(livep(p,:) == 0._r8 .and. rankp(p,:) > 0))*dt
!            !leafn_senescent(p) = leafn_senescent(p) + sum(pleafn_to_litter(p,:), mask=(livep(p,:) == 0._r8 .and. rankp(p,:) > 0))*dt
!            leafc_senescent(p) = leafc_senescent(p) + sum(pleafc_to_litter(p,:))*dt
!            leafn_senescent(p) = leafn_senescent(p) + sum(pleafn_to_litter(p,:))*dt
!            !at pruning, move senescent pools to litter in CNLitterToColum

            leafc_to_litter(p) = sum(pleafc_to_litter(p,:)) 
            leafn_to_litter(p) = sum(pleafn_to_litter(p,:)) 
            leafn_to_retransn(p) = sum(pleafn_to_retransn(p,:)) !use explicit N profile for each layer
         end if

         ! only calculate these fluxes if the background litterfall rate is non-zero
         if (bglfr(p) > 0._r8) then

          if (phytomer(ivt(p)) > 0) then
            !bglfr rate for fine root, assume the same as leaf turnover rate, like other PFTs
            frootc_to_litter(p) = bglfr(p) * frootc(p)
            frootn_to_litter(p) = bglfr(p) * frootn(p)
          else !for trees and other crops

            ! units for bglfr are already 1/s
            leafc_to_litter(p)  = bglfr(p) * leafc(p)
            frootc_to_litter(p) = bglfr(p) * frootc(p)
            if ( use_fun ) then
               leafc_to_litter_fun(p)     = leafc_to_litter(p)
               leafn_to_retransn(p)       = paid_retransn_to_npool(p) + free_retransn_to_npool(p)
               if (leafn(p).gt.0._r8) then
                  if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
                     leafcn_offset(p)     = leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
                  else
                     leafcn_offset(p)     = leafc(p)/leafn(p)
                  end if
               else
                  leafcn_offset(p)        = leafcn(ivt(p))
               end if
               leafn_to_litter(p)         = leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
               leafn_to_litter(p)         = max(leafn_to_litter(p),0._r8)

               denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
               if ( denom /= 0.0_r8 ) then
                  fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
               else if ( leafn_to_litter(p) == 0.0_r8 ) then
                  fr_leafn_to_litter =  0.0_r8
               else
                  fr_leafn_to_litter =  1.0_r8
               end if


            else
               if (CNratio_floating .eqv. .true.) then    
                  fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
               end if
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
               leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

            end if    

            ! calculate fine root N litterfall (no retranslocation of fine root N)
            frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))
            
            if (CNratio_floating .eqv. .true.) then    
               if (leafc(p) == 0.0_r8) then    
                  ntovr_leaf = 0.0_r8    
               else    
                  ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))   
               end if   
           
               leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
               leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
               if (frootc(p) == 0.0_r8) then    
                   frootn_to_litter(p) = 0.0_r8    
                else    
                   frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))   
                end if   
            end if    

            if ( use_fun ) then
               if(frootn_to_litter(p)*dt.gt.frootn(p))then
                    frootn_to_litter(p) = frootn(p)/dt
               endif
            end if

          end if !condition for palm specific bglfr

         end if

      end do

    end associate 

  end subroutine CNBackgroundLitterfall

  !-----------------------------------------------------------------------
  subroutine CNLivewoodTurnover (num_soilp, filter_soilp, &
       cnveg_state_inst,cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from live wood to
    ! dead wood pools, for stem and coarse root.
    !
    use CNSharedParamsMod, only: use_fun
    use clm_varctl          , only : CNratio_floating 
    use pftconMod        , only : npcropmin   
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter patch index
    real(r8):: ctovr        ! temporary variable for carbon turnover
    real(r8):: ntovr        ! temporary variable for nitrogen turnover
    !-----------------------------------------------------------------------

    associate(                                                                                   & 
         ivt                      =>    patch%itype                                              , & ! Input:  [integer  (:) ]  patch vegetation type                                
         bglfr                    =>    cnveg_state_inst%bglfr_patch                          , & ! Input:  [real(r8) (:) ]  background litterfall rate (1/s)      
         phytomer                 =>    pftcon%phytomer                                       , & ! Input:  [integer (:)]   total number of phytomers in life time (>0 use the new PhytomerPhenology) (added by Y.Fan) 
         woody                    =>    pftcon%woody                                           , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         livewdcn                 =>    pftcon%livewdcn                                        , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN) 
         deadwdcn                 =>    pftcon%deadwdcn                                        , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)       
         livestemc                =>    cnveg_carbonstate_inst%livestemc_patch                 , & ! Input:  [real(r8) (:) ]  (gC/m2) live stem C                               
         livecrootc               =>    cnveg_carbonstate_inst%livecrootc_patch                , & ! Input:  [real(r8) (:) ]  (gC/m2) live coarse root C                        
         livestemn                =>    cnveg_nitrogenstate_inst%livestemn_patch               , & ! Input:  [real(r8) (:) ]  (gN/m2) live stem N                               
         livecrootn               =>    cnveg_nitrogenstate_inst%livecrootn_patch              , & ! Input:  [real(r8) (:) ]  (gN/m2) live coarse root N                        
         
         livestemc_to_deadstemc   =>    cnveg_carbonflux_inst%livestemc_to_deadstemc_patch     , & ! Output: [real(r8) (:) ]                                                    
         livecrootc_to_deadcrootc =>    cnveg_carbonflux_inst%livecrootc_to_deadcrootc_patch   , & ! Output: [real(r8) (:) ]                                                    
         livestemn_to_deadstemn   =>    cnveg_nitrogenflux_inst%livestemn_to_deadstemn_patch   , & ! Output: [real(r8) (:) ]                                                    
         livestemn_to_retransn    =>    cnveg_nitrogenflux_inst%livestemn_to_retransn_patch    , & ! Output: [real(r8) (:) ]                                                    
         livecrootn_to_deadcrootn =>    cnveg_nitrogenflux_inst%livecrootn_to_deadcrootn_patch , & ! Output: [real(r8) (:) ]                                                    
         livecrootn_to_retransn   =>    cnveg_nitrogenflux_inst%livecrootn_to_retransn_patch     & ! Output: [real(r8) (:) ]                                                    
         )



      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes for woody types	 
         if (woody(ivt(p)) > 0._r8) then

           ! live stem to dead stem turnover

            ctovr = livestemc(p) * lwtop
            !ntovr = ctovr / livewdcn(ivt(p))  !ashehad turned this off
	    ntovr = livestemn(p) * lwtop
            livestemc_to_deadstemc(p) = ctovr
            livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
            
            if (CNratio_floating .eqv. .true.) then    
               if (livestemc(p) == 0.0_r8) then    
                   ntovr = 0.0_r8    
                else    
                   ntovr = ctovr * (livestemn(p) / livestemc(p))   
                end if   

                livestemn_to_deadstemn(p) = 0.5_r8 * ntovr   ! assuming 50% goes to deadstemn 
            end if    
            
            livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)    

            ! live coarse root to dead coarse root turnover

            ctovr = livecrootc(p) * lwtop
            !ntovr = ctovr / livewdcn(ivt(p))  !ashehad turned this off
	    ntovr = livecrootn(p) * lwtop
            livecrootc_to_deadcrootc(p) = ctovr
            livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
            
            if (CNratio_floating .eqv. .true.) then    
              if (livecrootc(p) == 0.0_r8) then    
                  ntovr = 0.0_r8    
               else    
                  ntovr = ctovr * (livecrootn(p) / livecrootc(p))   
               end if   

               livecrootn_to_deadcrootn(p) = 0.5_r8 * ntovr   ! assuming 50% goes to deadstemn 
            end if    
            
            livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)
            if(use_fun)then
               !TURNED OFF FLUXES TO CORRECT N ACCUMULATION ISSUE. RF. Oct 2015. 
               livecrootn_to_retransn(p) = 0.0_r8
               livestemn_to_retransn(p)  = 0.0_r8
            endif

          end if

          !for oil palm, use the same turnover rate for leaf and froot except
          !that 6% of stem is assumed alive according to van Kraalingen et al.1989
          !oil palm's trunk below the lowest leaf contains only 6% actively
          !respiring tissue compared to a 100% active upper trunk section
          !it is equivalent to consider that most stem C (94%) accumulated
          !during the life of a leaf attached to the stem section becomes dead.
          !only the stem C being accumulated together with alive leaves are
          !considered 100% active (Y.Fan)
          if (phytomer(ivt(p)) > 0) then
            livestemc_to_deadstemc(p) = bglfr(p) * livestemc(p)*(1._r8 - 0.06_r8)
            livestemn_to_deadstemn(p) = livestemc_to_deadstemc(p) / deadwdcn(ivt(p))
            livestemn_to_retransn(p) = bglfr(p) * livestemn(p)*(1._r8 - 0.06_r8) - livestemn_to_deadstemn(p)
            livecrootc_to_deadcrootc(p) = 0.0_r8
            livecrootn_to_deadcrootn(p) = 0.0_r8
          end if

      end do

    end associate

  end subroutine CNLivewoodTurnover

  !-----------------------------------------------------------------------
  subroutine CNGrainToProductPools(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! If using prognostic crop along with use_grainproduct, then move the patch-level
    ! grain-to-food fluxes into the column-level grain-to-cropprod fluxes
    !
    ! !USES:
    use clm_varctl    , only : use_crop
    use clm_varctl    , only : use_grainproduct
    use subgridAveMod , only : p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)             , intent(in)    :: bounds
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                       , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                       , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p

    character(len=*), parameter :: subname = 'CNGrainToProductPools'
    !-----------------------------------------------------------------------

    ! Explicitly checking use_crop is probably unnecessary here (because presumably
    ! use_grainproduct is only true if use_crop is true), but we do it for safety because
    ! the grain*_to_food_patch fluxes are not set if use_crop is false.
    if (use_crop .and. use_grainproduct) then
       do fp = 1, num_soilp
          p = filter_soilp(fp)
          cnveg_carbonflux_inst%grainc_to_cropprodc_patch(p) = &
               cnveg_carbonflux_inst%grainc_to_food_patch(p)
          cnveg_nitrogenflux_inst%grainn_to_cropprodn_patch(p) = &
               cnveg_nitrogenflux_inst%grainn_to_food_patch(p)
       end do

       call p2c (bounds, num_soilc, filter_soilc, &
            cnveg_carbonflux_inst%grainc_to_cropprodc_patch(bounds%begp:bounds%endp), &
            cnveg_carbonflux_inst%grainc_to_cropprodc_col(bounds%begc:bounds%endc))

       call p2c (bounds, num_soilc, filter_soilc, &
            cnveg_nitrogenflux_inst%grainn_to_cropprodn_patch(bounds%begp:bounds%endp), &
            cnveg_nitrogenflux_inst%grainn_to_cropprodn_col(bounds%begc:bounds%endc))
    end if

    ! No else clause: if use_grainproduct is false, then the grain*_to_cropprod fluxes
    ! will remain at their initial value (0).

  end subroutine CNGrainToProductPools

  !-----------------------------------------------------------------------
  subroutine CNLitterToColumn (bounds, num_soilc, filter_soilc,         &
       cnveg_state_inst,cnveg_carbonflux_inst,cnveg_nitrogenflux_inst, &
       crop_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst, & 
       leaf_prof_patch, froot_prof_patch)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all patch-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : max_patch_per_col, nlevdecomp
    use pftconMod  , only : npcropmin
    use clm_varctl , only : use_grainproduct
    !use dynHarvestMod , only: CNHarvestPftToColumn !need to make this subroutine public
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_state_type)          , intent(in)    :: cnveg_state_inst
    type(crop_type)                 , intent(inout) :: crop_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    !type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    !class(cnharvest_type) , intent(inout) :: cnharvest

    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j       ! indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))

    associate(                                                                                & 
         leaf_prof                 => leaf_prof_patch                                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                => froot_prof_patch                                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         !phytomer                  =>    pftcon%phytomer                                    , & ! Input:  [integer (:)]   total number of phytomers in life time (>0 use the new PhytomerPhenology) (added by Y.Fan)
         !leafc_senescent           =>    cnveg_carbonstate_inst%leafc_senescent_patch       , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C saved for pruning (added by Y.Fan)
         !leafn_senescent           =>    cnveg_nitrogenstate_inst%leafn_senescent_patch     , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N saved for pruning (added by Y.Fan)
         !offset_flag               =>    cnveg_state_inst%offset_flag_patch                 , & ! Input:  [real(r8) (:)]  offset flag for perennial crop ratation (Y.Fan)
         !prune                     =>    crop_inst%prune_patch                            , & ! Input:  [real(r8) (:)]  flag for pruning (Y.Fan)
         woody                     =>    pftcon%woody                                      , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
	 perennial                 =>    pftcon%perennial                                      , & ! Input:  [real(r8) (:)]  
         deadstemc_to_litter       =>    cnveg_carbonflux_inst%deadstemc_to_litter_patch     , & ! InOut:  [real(r8) (:)]  dead stem C litterfall (gC/m2/s) (Y.Fan)
         deadstemn_to_litter       =>    cnveg_nitrogenflux_inst%deadstemn_to_litter_patch   , & ! InOut:  [real(r8) (:)]  deadstem N to litter (gN/m2/s) (Y.Fan)
         livecrootc_to_litter      =>    cnveg_carbonflux_inst%livecrootc_to_litter_patch    , & ! InOut:  [real(r8) (:)]  live coarse root C litterfall (gC/m2/s) (Y.Fan)
         deadcrootc_to_litter      =>    cnveg_carbonflux_inst%deadcrootc_to_litter_patch    , & ! InOut:  [real(r8) (:)]  dead coarse root C litterfall (gC/m2/s) (Y.Fan)
         livecrootn_to_litter      =>    cnveg_nitrogenflux_inst%livecrootn_to_litter_patch  , & ! InOut:  [real(r8) (:)]  live coarse root N litter (gN/m2/s) (Y.Fan)
         deadcrootn_to_litter      =>    cnveg_nitrogenflux_inst%deadcrootn_to_litter_patch  , & ! InOut:  [real(r8) (:)]  dead coarse root N litter (gN/m2/s) (Y.Fan)

         ivt                       => patch%itype                                             , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                     => patch%wtcol                                             , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this patch (0-1)
	 
	 wood_harvestc             =>    cnveg_carbonflux_inst%wood_harvestc_patch            , & ! Output: [real(r8) (:)]
	 wood_harvestn             =>    cnveg_nitrogenflux_inst%wood_harvestn_patch          , & ! Output: [real(r8) (:)]    

         lf_flab                   => pftcon%lf_flab                                        , & ! Input:  leaf litter labile fraction                       
         lf_fcel                   => pftcon%lf_fcel                                        , & ! Input:  leaf litter cellulose fraction                    
         lf_flig                   => pftcon%lf_flig                                        , & ! Input:  leaf litter lignin fraction                       
         fr_flab                   => pftcon%fr_flab                                        , & ! Input:  fine root litter labile fraction                  
         fr_fcel                   => pftcon%fr_fcel                                        , & ! Input:  fine root litter cellulose fraction               
         fr_flig                   => pftcon%fr_flig                                        , & ! Input:  fine root litter lignin fraction                  

         leafc_to_litter           => cnveg_carbonflux_inst%leafc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]  leaf C litterfall (gC/m2/s)                       
         frootc_to_litter          => cnveg_carbonflux_inst%frootc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)                  
         livestemc_to_litter       => cnveg_carbonflux_inst%livestemc_to_litter_patch       , & ! Input:  [real(r8) (:)   ]  live stem C litterfall (gC/m2/s)                  
         grainc_to_food            => cnveg_carbonflux_inst%grainc_to_food_patch            , & ! Input:  [real(r8) (:)   ]  grain C to food (gC/m2/s)                         
         phenology_c_to_litr_met_c => cnveg_carbonflux_inst%phenology_c_to_litr_met_c_col   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
         phenology_c_to_litr_cel_c => cnveg_carbonflux_inst%phenology_c_to_litr_cel_c_col   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
         phenology_c_to_litr_lig_c => cnveg_carbonflux_inst%phenology_c_to_litr_lig_c_col   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)

         livestemn_to_litter       => cnveg_nitrogenflux_inst%livestemn_to_litter_patch     , & ! Input:  [real(r8) (:)   ]  livestem N to litter (gN/m2/s)                    
         grainn_to_food            => cnveg_nitrogenflux_inst%grainn_to_food_patch          , & ! Input:  [real(r8) (:)   ]  grain N to food (gN/m2/s)                         
         leafn_to_litter           => cnveg_nitrogenflux_inst%leafn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]  leaf N litterfall (gN/m2/s)                       
         frootn_to_litter          => cnveg_nitrogenflux_inst%frootn_to_litter_patch        , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)                  
         phenology_n_to_litr_met_n => cnveg_nitrogenflux_inst%phenology_n_to_litr_met_n_col , & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
         phenology_n_to_litr_cel_n => cnveg_nitrogenflux_inst%phenology_n_to_litr_cel_n_col , & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
         phenology_n_to_litr_lig_n => cnveg_nitrogenflux_inst%phenology_n_to_litr_lig_n_col   & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)
         )
    
      do j = 1, nlevdecomp
         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if ( pi <=  col%npatches(c) ) then
                  p = col%patchi(c) + pi - 1
                  if (patch%active(p)) then

                   !****now simplify the leaf litterfall routine for oil palm,
                   !all leaf litter fluxes are treated in the background
                   !litterfall and offset litterfall subroutines (2022.10)
  
                   ! if (phytomer(ivt(p)) > 0) then
                   !    ! leaf litter carbon fluxes
                   !    !for oil palm, only move senescnet pool to litter at the time of pruning
                   !    if (prune(p) .or. offset_flag(p) == 1._r8) then
                   !       !do not write two conditions together like (phytomer(ivt(p)) > 0 .and. prune(p)),
                   !       !otherwise the else block below will execute and leafc_to_litter will move to litter pool at every time step and cause cbalance error
                   !       !for palm structure, leafc_to_litter(p) only meant for real-time updating leafc in CNCStateUpdate1Mod for respiration cost
                   !       !here leafc_to_litter(p) only comes together with leafc_senescent(p) for one time step at each pruning
                   !       phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                   !           + (leafc_to_litter(p) + leafc_senescent(p)/dt) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   !       phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                   !           + (leafc_to_litter(p) + leafc_senescent(p)/dt) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   !       phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                   !           + (leafc_to_litter(p) + leafc_senescent(p)/dt) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                   !       ! leaf litter nitrogen fluxes
                   !       phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                   !           + (leafn_to_litter(p) + leafn_senescent(p)/dt) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   !       phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                   !           + (leafn_to_litter(p) + leafn_senescent(p)/dt) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   !       phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                   !           + (leafn_to_litter(p) + leafn_senescent(p)/dt) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                   !       ! clear up the senescent pools
                   !       leafc_senescent(p) = 0._r8
                   !       leafn_senescent(p) = 0._r8

                   !       !at final rotation (don't worry about harvest fluxes
                   !       !of some inmmature palm fruits at the final rotation, 14.03.2022
                   !       if (offset_flag(p) == 1._r8) then
                   !          call CNHarvestPftToColumn(num_soilc, filter_soilc, &
                   !               soilbiogeochem_state_inst,cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
                   !          !must need to summarize hrv fluxes to column level, use the function from pftdynMod
                   !       end if
                   !       !the CNHarvestPftToColumn function will summarize both leaf/froot/livestem pools to litter materials met_c/cel_c/lig_c &
                   !       !and also summarize wood product pools (hrv_deadstemc_to_prod10c/hrv_deadstemc_to_prod100c) to column level

                   !       ! deadstem is added to 10-year wood product pool instead (Y.Fan 2016)
                   !       ! adding deadstem into litter pool will dramatically increase LITTERC_HR which will drain out SMINN very soon!
                   !   ! end if

                   ! else

                     ! leaf litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! leaf litter nitrogen fluxes
                     phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                          + leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                          + leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                          + leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   ! end if

                     ! fine root litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! fine root litter nitrogen fluxes
                     phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                          + frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                          + frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                          + frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! agroibis puts crop stem litter together with leaf litter
                     ! so I've used the leaf lf_f* parameters instead of making
                     ! new ones for now (slevis)
                     ! also for simplicity I've put "food" into the litter pools

                     if (ivt(p) >= npcropmin) then ! add livestemc to litter
                        ! stem litter carbon fluxes
                        phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                             + livestemc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                             + livestemc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                             + livestemc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                        ! stem litter nitrogen fluxes
                        phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                             + livestemn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                             + livestemn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                             + livestemn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                        if(woody(ivt(p)) == 1.0_r8) then   
                        !below will be zero for non-woody crop types (Y.Fan)
                            !dead stem
                            phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                              + deadstemc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                              + deadstemc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                              + deadstemc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           !live coarse root
                            phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                              + livecrootc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                              + livecrootc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                              + livecrootc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           !dead coarse root
                            phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                              + deadcrootc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                              + deadcrootc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                            phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                              + deadcrootc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                           !Nitrogen fluxes: 
                           !dead stem
                           phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                             + deadstemn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                             + deadstemn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                             + deadstemn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                          !live coarse root
                           phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                             + livecrootn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                             + livecrootn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                             + livecrootn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                          !dead coarse root
                           phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                             + deadcrootn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                             + deadcrootn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                           phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                             + deadcrootn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        end if

                        if (.not. use_grainproduct) then
                         ! grain litter carbon fluxes
                         phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                              + grainc_to_food(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                         phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                              + grainc_to_food(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                         phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                              + grainc_to_food(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
 
                         ! grain litter nitrogen fluxes
                         phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                              + grainn_to_food(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                         phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                              + grainn_to_food(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                         phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                              + grainn_to_food(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        end if


                     end if
                  end if
               end if

            end do

         end do
      end do
      
     end associate 

  end subroutine CNLitterToColumn

end module CNPhenologyMod
