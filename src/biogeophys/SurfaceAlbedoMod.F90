module SurfaceAlbedoMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs surface albedo calculations
  !
  ! !PUBLIC TYPES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use landunit_varcon   , only : istsoil, istcrop, istdlak
  use clm_varcon        , only : grlnd, namep
  use clm_varpar        , only : numrad, nlevcan, nlevsno, nlevcan
  use clm_varctl        , only : fsurdat, iulog, use_snicar_frc, use_SSRE
  use pftconMod         , only : pftcon
  use SnowSnicarMod     , only : sno_nbr_aer, SNICAR_RT, DO_SNO_AER, DO_SNO_OC
  use AerosolMod        , only : aerosol_type
  use CanopyStateType   , only : canopystate_type
  use LakeStateType     , only : lakestate_type
  use SurfaceAlbedoType , only : surfalb_type
  use TemperatureType   , only : temperature_type
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use GridcellType      , only : grc                
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : patch                
  !
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedo_readnl
  public :: SurfaceAlbedoInitTimeConst
  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilAlbedo    ! Determine ground surface albedo
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer
  private :: Multilayer    ! Multilayer radiative transfer scheme for sub-PFT palm structure (Y.Fan 2014)
  private :: Gfunc_all     ! The new G function for multilayer canopy (Y.Fan 2014)
  private :: Beta_Dist     ! Leaf angle Beta distribution used for G function (Y.Fan 2014)

  !
  ! !PUBLIC DATA MEMBERS:
  ! The CLM default albice values are too high.
  ! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
  ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: albice(numrad) = (/ 0.80_r8, 0.55_r8 /)

  ! namelist default setting for inputting alblakwi
  real(r8), public  :: lake_melt_icealb(numrad) = (/ 0.10_r8, 0.10_r8/)


  ! !PRIVATE DATA MEMBERS:
 ! real(r8), allocatable :: gfunc_solar(:)   ! G function for direct beam
 ! real(r8), allocatable :: gfunc_sky(:,:)   ! G function for diffuse sky radiation
 ! real(r8), allocatable :: lad(:)   		 ! leaf angle distribution probability density function

  ! albedo frozen lakes by waveband (1=vis, 2=nir) 
  ! unclear what the reference is for this
  real(r8), private :: alblak(numrad) = (/0.60_r8, 0.40_r8/)

  ! albedo of melting lakes due to puddling, open water, or white ice
  ! From D. Mironov (2010) Boreal Env. Research
  ! To revert albedo of melting lakes to the cold snow-free value, set
  ! lake_melt_icealb namelist to 0.60, 0.40 like alblak above.
  real(r8), private :: alblakwi(numrad)            

  ! Coefficient for calculating ice "fraction" for lake surface albedo
  ! From D. Mironov (2010) Boreal Env. Research
  real(r8), parameter :: calb = 95.6_r8   

  !
  ! !PRIVATE DATA MEMBERS:
  logical, private :: snowveg_affects_radiation = .true. ! Whether snow on the vegetation canopy affects the radiation/albedo calculations

  !
  ! !PRIVATE DATA FUNCTIONS:
  real(r8), allocatable, private :: albsat(:,:) ! wet soil albedo by color class and waveband (1=vis,2=nir)
  real(r8), allocatable, private :: albdry(:,:) ! dry soil albedo by color class and waveband (1=vis,2=nir)
  integer , allocatable, private :: isoicol(:)  ! column soil color class

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedo_readnl( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for SurfaceAlbedo
    !
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = "surfacealbedo_inparm"

    character(len=*), parameter :: subname = 'SurfaceAlbedo_readnl'
    !-----------------------------------------------------------------------

    namelist /surfacealbedo_inparm/ snowveg_affects_radiation

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=surfacealbedo_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(snowveg_affects_radiation, mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) nmlname, ' settings'
       write(iulog,nml=surfacealbedo_inparm)
       write(iulog,*)
    end if

  end subroutine SurfaceAlbedo_readnl


  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedoInitTimeConst(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module time constant variables
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use fileutils  , only : getfil
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    use spmdMod    , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g          ! indices
    integer            :: mxsoil_color ! maximum number of soil color classes
    type(file_desc_t)  :: ncid         ! netcdf id
    character(len=256) :: locfn        ! local filename
    integer            :: ier          ! error status
    logical            :: readvar 
    integer  ,pointer  :: soic2d (:)   ! read in - soil color 
    !---------------------------------------------------------------------

    ! Allocate module variable for soil color

    allocate(isoicol(bounds%begc:bounds%endc)) 

    ! Determine soil color and number of soil color classes 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, readvar=readvar)
    if ( .not. readvar ) then
       call endrun(msg=' ERROR: mxsoil_color NOT on surfdata file '//errMsg(sourcefile, __LINE__))
    end if

    allocate(soic2d(bounds%begg:bounds%endg)) 
    call ncd_io(ncid=ncid, varname='SOIL_COLOR', flag='read', data=soic2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: SOIL_COLOR NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       isoicol(c) = soic2d(g)
    end do
    deallocate(soic2d)

    call ncd_pio_closefile(ncid)

    ! Determine saturated and dry soil albedos for n color classes and 
    ! numrad wavebands (1=vis, 2=nir)

    allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'allocation error for albsat, albdry'
       call endrun(msg=errMsg(sourcefile, __LINE__)) 
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil colo data .....'
    end if
    
    if (mxsoil_color == 8) then
       albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
       albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
    else if (mxsoil_color == 20) then
       albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
            0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
       albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
            0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
            0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
            0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
    else
       write(iulog,*)'maximum color class = ',mxsoil_color,' is not supported'
       call endrun(msg=errMsg(sourcefile, __LINE__)) 
    end if

    ! Set alblakwi
    alblakwi(:) = lake_melt_icealb(:)

  end subroutine SurfaceAlbedoInitTimeConst

  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedo(bounds,nc,  &
        num_nourbanc, filter_nourbanc, &
        num_nourbanp, filter_nourbanp, &
        num_urbanc   , filter_urbanc,  &
        num_urbanp   , filter_urbanp,  &
        nextsw_cday  , declinp1,       &
        clm_fates,                     &
        atm2lnd_inst,       & !add atm2lnd_inst for multilayer radiative transfer model (Y.fan) 
        aerosol_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
        lakestate_inst, temperature_inst, surfalb_inst)
    !
    ! !DESCRIPTION:
    ! Surface albedo and two-stream fluxes
    ! Surface albedos. Also fluxes (per unit incoming direct and diffuse
    ! radiation) reflected, transmitted, and absorbed by vegetation.
    ! Calculate sunlit and shaded fluxes as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy to calculate APAR profile
     !
    ! The calling sequence is:
    ! -> SurfaceAlbedo:  albedos for next time step
    !    -> SoilAlbedo:  soil/lake/glacier/wetland albedos
    !    -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
    !    -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
    !    -> TwoStream:   absorbed, reflected, transmitted solar fluxes (vis dir,vis dif, nir dir, nir dif)
    !
    ! Note that this is called with the "inactive_and_active" version of the filters, because
    ! the variables computed here are needed over inactive points that might later become
    ! active (due to landuse change). Thus, this routine cannot depend on variables that are
    ! only computed over active points.
    !
    ! !USES:
    use shr_orb_mod
    use clm_time_manager   , only : get_nstep
    use abortutils         , only : endrun
    use clm_varctl         , only : subgridflag, use_snicar_frc, use_fates
    use CLMFatesInterfaceMod, only : hlm_fates_interface_type
    use pftconMod          , only : pftcon, mxnp, noilpalm, nirrig_oilpalm  !  added by Y.Fan 2016
    use clm_varpar         , only : nlevcan, radiative_transfer       ! added by Y.Fan 2016
    use clm_varctl         , only : use_cn   !Y.Fan
    use atm2lndType        , only : atm2lnd_type !Y.Fan

    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)            :: bounds             ! bounds
    integer                , intent(in)            :: nc                 ! clump index
    integer                , intent(in)            :: num_nourbanc       ! number of columns in non-urban filter
    integer                , intent(in)            :: filter_nourbanc(:) ! column filter for non-urban points
    integer                , intent(in)            :: num_nourbanp       ! number of patches in non-urban filter
    integer                , intent(in)            :: filter_nourbanp(:) ! patch filter for non-urban points
    integer                , intent(in)            :: num_urbanc         ! number of columns in urban filter
    integer                , intent(in)            :: filter_urbanc(:)   ! column filter for urban points
    integer                , intent(in)            :: num_urbanp         ! number of patches in urban filter
    integer                , intent(in)            :: filter_urbanp(:)   ! patch filter for rban points
    real(r8)               , intent(in)            :: nextsw_cday        ! calendar day at Greenwich (1.00, ..., days/year)
    real(r8)               , intent(in)            :: declinp1           ! declination angle (radians) for next time step
    type(hlm_fates_interface_type), intent(inout)  :: clm_fates
    type(aerosol_type)     , intent(in)            :: aerosol_inst
    type(canopystate_type) , intent(in)            :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)            :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)            :: waterdiagnosticbulk_inst
    type(lakestate_type)   , intent(in)            :: lakestate_inst
    type(temperature_type) , intent(in)            :: temperature_inst
    type(surfalb_type)     , intent(inout)         :: surfalb_inst
    type(atm2lnd_type)     , intent(in)            :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer, allocatable :: ix(:)                                                         ! indices of expanded phytomers with positive LAI
    integer  :: ixn                                                                       ! the variable size of array ix (Y.Fan)
    integer  :: nump                                                                      ! the number of phytomers per layer
    integer  :: i                                                                         ! index for layers [idx]
    integer  :: aer                                                                       ! index for sno_nbr_aer
    real(r8) :: extkn                                                                     ! nitrogen allocation coefficient
    integer  :: fp,fc,g,c,p,iv                                                            ! indices
    integer  :: ib                                                                        ! band index
    integer  :: ic                                                                        ! 0=unit incoming direct; 1=unit incoming diffuse
    real(r8) :: dinc                                                                      ! lai+sai increment for canopy layer
    real(r8) :: dincmax                                                                   ! maximum lai+sai increment for canopy layer
    real(r8) :: dincmax_sum                                                               ! cumulative sum of maximum lai+sai increment for canopy layer
    real(r8) :: laisum                                                                    ! sum of canopy layer lai for error check
    real(r8) :: saisum                                                                    ! sum of canopy layer sai for error check
    integer  :: flg_slr                                                                   ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer  :: flg_snw_ice                                                               ! flag for SNICAR (=1 when called from CLM, =2 when called from sea-ice)
    integer  :: num_vegsol                                                                ! number of vegetated patches where coszen>0
    integer  :: num_novegsol                                                              ! number of vegetated patches where coszen>0
    integer  :: filter_vegsol   (bounds%endp-bounds%begp+1)                               ! patch filter where vegetated and coszen>0
    integer  :: filter_novegsol (bounds%endp-bounds%begp+1)                               ! patch filter where vegetated and coszen>0
    real(r8) :: wl              (bounds%begp:bounds%endp)                                 ! fraction of LAI+SAI that is LAI
    real(r8) :: ws              (bounds%begp:bounds%endp)                                 ! fraction of LAI+SAI that is SAI
    real(r8) :: blai(bounds%begp:bounds%endp)                                             ! lai buried by snow: tlai - elai
    real(r8) :: bsai(bounds%begp:bounds%endp)                                             ! sai buried by snow: tsai - esai
    real(r8) :: coszen_gcell    (bounds%begg:bounds%endg)                                 ! cosine solar zenith angle for next time step (grc)
    real(r8) :: coszen_patch    (bounds%begp:bounds%endp)                                 ! cosine solar zenith angle for next time step (patch)
    real(r8) :: rho(bounds%begp:bounds%endp,numrad)                                       ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8) :: tau(bounds%begp:bounds%endp,numrad)                                       ! leaf/stem tran weighted by fraction LAI and SAI
    real(r8) :: albsfc          (bounds%begc:bounds%endc,numrad)                          ! albedo of surface underneath snow (col,bnd) 
    real(r8) :: albsnd(bounds%begc:bounds%endc,numrad)                                    ! snow albedo (direct)
    real(r8) :: albsni(bounds%begc:bounds%endc,numrad)                                    ! snow albedo (diffuse)
    real(r8) :: albsnd_pur      (bounds%begc:bounds%endc,numrad)                          ! direct pure snow albedo (radiative forcing)
    real(r8) :: albsni_pur      (bounds%begc:bounds%endc,numrad)                          ! diffuse pure snow albedo (radiative forcing)
    real(r8) :: albsnd_bc       (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without BC (radiative forcing)
    real(r8) :: albsni_bc       (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without BC (radiative forcing)
    real(r8) :: albsnd_oc       (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without OC (radiative forcing)
    real(r8) :: albsni_oc       (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without OC (radiative forcing)
    real(r8) :: albsnd_dst      (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without dust (radiative forcing)
    real(r8) :: albsni_dst      (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without dust (radiative forcing)
    real(r8) :: flx_absd_snw    (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! flux absorption factor for just snow (direct) [frc]
    real(r8) :: flx_absi_snw    (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! flux absorption factor for just snow (diffuse) [frc]
    real(r8) :: foo_snw         (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! dummy array for forcing calls
    real(r8) :: h2osno_liq      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! liquid snow content (col,lyr) [kg m-2]
    real(r8) :: h2osno_ice      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! ice content in snow (col,lyr) [kg m-2]
    integer  :: snw_rds_in      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! snow grain size sent to SNICAR (col,lyr) [microns]
    real(r8) :: mss_cnc_aer_in_frc_pur (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_bc  (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for BC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_oc  (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for OC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_dst (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for dust forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_fdb     (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of all aerosol species for feedback calculation (col,lyr,aer) [kg kg-1]
    real(r8), parameter :: mpe = 1.e-06_r8                                                ! prevents overflow for division by zero
    integer , parameter :: nband =numrad                                                  ! number of solar radiation waveband classes
    !-----------------------------------------------------------------------

   associate(&
          rhol          =>    pftcon%rhol                         , & ! Input:  leaf reflectance: 1=vis, 2=nir        
          rhos          =>    pftcon%rhos                         , & ! Input:  stem reflectance: 1=vis, 2=nir        
          taul          =>    pftcon%taul                         , & ! Input:  leaf transmittance: 1=vis, 2=nir      
          taus          =>    pftcon%taus                         , & ! Input:  stem transmittance: 1=vis, 2=nir      

          tlai          =>    canopystate_inst%tlai_patch         , & ! Input:  [real(r8)  (:)   ]  one-sided leaf area index, no burying by snow 
          tsai          =>    canopystate_inst%tsai_patch         , & ! Input:  [real(r8)  (:)   ]  one-sided stem area index, no burying by snow
          elai          =>    canopystate_inst%elai_patch         , & ! Input:  [real(r8)  (:)   ]  one-sided leaf area index with burying by snow
          esai          =>    canopystate_inst%esai_patch         , & ! Input:  [real(r8)  (:)   ]  one-sided stem area index with burying by snow

          frac_sno      =>    waterdiagnosticbulk_inst%frac_sno_col        , & ! Input:  [real(r8)  (:)   ]  fraction of ground covered by snow (0 to 1)
          h2osno        =>    waterstatebulk_inst%h2osno_col          , & ! Input:  [real(r8)  (:)   ]  snow water (mm H2O)                     
          h2osoi_liq    =>    waterstatebulk_inst%h2osoi_liq_col      , & ! Input:  [real(r8)  (:,:) ]  liquid water content (col,lyr) [kg/m2]
          h2osoi_ice    =>    waterstatebulk_inst%h2osoi_ice_col      , & ! Input:  [real(r8)  (:,:) ]  ice lens content (col,lyr) [kg/m2]    
          snw_rds       =>    waterdiagnosticbulk_inst%snw_rds_col         , & ! Input:  [real(r8)  (:,:) ]  snow grain radius (col,lyr) [microns] 
          
          ivt           =>    patch%itype                         , & ! Input:  [integer   (:)   ]  pft vegetation type
          !np           =>    pps%np                              , & ! Input:  [integer   (:)   ]  total number of phytomers having appeared so far
          !livep        =>    pps%livep                           , & ! Input:  [real(r8)  (:,:) ]  Flag, true if this phytomer is alive
          plai          =>    canopystate_inst%plai_patch         , & ! Input:  [real(r8)  (:,:) ]  one-sided leaf area index of each phytomer
          rankp         =>    canopystate_inst%rankp_patch        , & ! Input:  [integer   (:,:) ]  rank of phytomers from 1=youngest to np=oldest and 0=dead
          phytomer      =>    pftcon%phytomer                     , & ! Input:  [integer   (:)   ]  total number of phytomers in life time (added by Y.Fan)
          mxlivenp      =>    pftcon%mxlivenp                     , & ! Input: 

          mss_cnc_bcphi =>    aerosol_inst%mss_cnc_bcphi_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophilic BC (col,lyr) [kg/kg]
          mss_cnc_bcpho =>    aerosol_inst%mss_cnc_bcpho_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophobic BC (col,lyr) [kg/kg]
          mss_cnc_ocphi =>    aerosol_inst%mss_cnc_ocphi_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophilic OC (col,lyr) [kg/kg]
          mss_cnc_ocpho =>    aerosol_inst%mss_cnc_ocpho_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophobic OC (col,lyr) [kg/kg]
          mss_cnc_dst1  =>    aerosol_inst%mss_cnc_dst1_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
          mss_cnc_dst2  =>    aerosol_inst%mss_cnc_dst2_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
          mss_cnc_dst3  =>    aerosol_inst%mss_cnc_dst3_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
          mss_cnc_dst4  =>    aerosol_inst%mss_cnc_dst4_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]

          fsun_z        =>    surfalb_inst%fsun_z_patch           , & ! Output:  [real(r8) (:,:) ]  sunlit fraction of canopy layer       
          tlai_z        =>    surfalb_inst%tlai_z_patch           , & ! Output:  [real(r8) (:,:) ]  tlai increment for canopy layer       
          tsai_z        =>    surfalb_inst%tsai_z_patch           , & ! Output:  [real(r8) (:,:) ]  tsai increment for canopy layer       
          vcmaxcintsun  =>    surfalb_inst%vcmaxcintsun_patch     , & ! Output:  [real(r8) (:)   ]  leaf to canopy scaling coefficient, sunlit leaf vcmax
          vcmaxcintsha  =>    surfalb_inst%vcmaxcintsha_patch     , & ! Output:  [real(r8) (:)   ]  leaf to canopy scaling coefficient, shaded leaf vcmax
          ncan          =>    surfalb_inst%ncan_patch             , & ! Output:  [integer  (:)   ]  number of canopy layers                  
          nrad          =>    surfalb_inst%nrad_patch             , & ! Output:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          coszen_col    =>    surfalb_inst%coszen_col             , & ! Output:  [real(r8) (:)   ]  cosine of solar zenith angle            
          albgrd        =>    surfalb_inst%albgrd_col             , & ! Output:  [real(r8) (:,:) ]  ground albedo (direct)                
          albgri        =>    surfalb_inst%albgri_col             , & ! Output:  [real(r8) (:,:) ]  ground albedo (diffuse)               
          albsod        =>    surfalb_inst%albsod_col             , & ! Output:  [real(r8) (:,:) ]  direct-beam soil albedo (col,bnd) [frc]
          albsoi        =>    surfalb_inst%albsoi_col             , & ! Output:  [real(r8) (:,:) ]  diffuse soil albedo (col,bnd) [frc]   
          albgrd_pur    =>    surfalb_inst%albgrd_pur_col         , & ! Output:  [real(r8) (:,:) ]  pure snow ground albedo (direct)      
          albgri_pur    =>    surfalb_inst%albgri_pur_col         , & ! Output:  [real(r8) (:,:) ]  pure snow ground albedo (diffuse)     
          albgrd_bc     =>    surfalb_inst%albgrd_bc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without BC (direct)     
          albgri_bc     =>    surfalb_inst%albgri_bc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without BC (diffuse)    
          albgrd_oc     =>    surfalb_inst%albgrd_oc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without OC (direct)     
          albgri_oc     =>    surfalb_inst%albgri_oc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without OC (diffuse)    
          albgrd_dst    =>    surfalb_inst%albgrd_dst_col         , & ! Output:  [real(r8) (:,:) ]  ground albedo without dust (direct)   
          albgri_dst    =>    surfalb_inst%albgri_dst_col         , & ! Output:  [real(r8) (:,:) ]  ground albedo without dust (diffuse)  
          albsnd_hst    =>    surfalb_inst%albsnd_hst_col         , & ! Output:  [real(r8) (:,:) ]  snow albedo, direct, for history files (col,bnd) [frc]
          albsni_hst    =>    surfalb_inst%albsni_hst_col         , & ! Output:  [real(r8) (:,:) ]  snow ground albedo, diffuse, for history files (col,bnd) [frc]
          albd          =>    surfalb_inst%albd_patch             , & ! Output:  [real(r8) (:,:) ]  surface albedo (direct)
          albi          =>    surfalb_inst%albi_patch             , & ! Output:  [real(r8) (:,:) ]  surface albedo (diffuse)
          albdSF        =>    surfalb_inst%albdSF_patch           , & ! Output:  [real(r8) (:,:) ]  diagnostic snow-free surface albedo (direct)
          albiSF        =>    surfalb_inst%albiSF_patch           , & ! Output:  [real(r8) (:,:) ]  diagnostic snow-free surface albedo (diffuse)
          fabd          =>    surfalb_inst%fabd_patch             , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
          fabd_sun      =>    surfalb_inst%fabd_sun_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
          fabd_sha      =>    surfalb_inst%fabd_sha_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
          fabi          =>    surfalb_inst%fabi_patch             , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
          fabi_sun      =>    surfalb_inst%fabi_sun_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha      =>    surfalb_inst%fabi_sha_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
          ftdd          =>    surfalb_inst%ftdd_patch             , & ! Output:  [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
          ftid          =>    surfalb_inst%ftid_patch             , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
          ftii          =>    surfalb_inst%ftii_patch             , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flux
          flx_absdv     =>    surfalb_inst%flx_absdv_col          , & ! Output:  [real(r8) (:,:) ]  direct flux absorption factor (col,lyr): VIS [frc]
          flx_absdn     =>    surfalb_inst%flx_absdn_col          , & ! Output:  [real(r8) (:,:) ]  direct flux absorption factor (col,lyr): NIR [frc]
          flx_absiv     =>    surfalb_inst%flx_absiv_col          , & ! Output:  [real(r8) (:,:) ]  diffuse flux absorption factor (col,lyr): VIS [frc]
          flx_absin     =>    surfalb_inst%flx_absin_col          , & ! Output:  [real(r8) (:,:) ]  diffuse flux absorption factor (col,lyr): NIR [frc]
          fabd_sun_z    =>    surfalb_inst%fabd_sun_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z    =>    surfalb_inst%fabd_sha_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z    =>    surfalb_inst%fabi_sun_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z    =>    surfalb_inst%fabi_sha_z_patch         & ! Output:  [real(r8) (:,:) ]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          )

    ! Cosine solar zenith angle for next time step

    do g = bounds%begg,bounds%endg
       coszen_gcell(g) = shr_orb_cosz (nextsw_cday, grc%lat(g), grc%lon(g), declinp1)
    end do
    do c = bounds%begc,bounds%endc
       g = col%gridcell(c)
       coszen_col(c) = coszen_gcell(g)
    end do
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       g = patch%gridcell(p)
          coszen_patch(p) = coszen_gcell(g)
    end do

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          albsod(c,ib)     = 0._r8
          albsoi(c,ib)     = 0._r8
          albgrd(c,ib)     = 0._r8
          albgri(c,ib)     = 0._r8
          albgrd_pur(c,ib) = 0._r8
          albgri_pur(c,ib) = 0._r8
          albgrd_bc(c,ib)  = 0._r8
          albgri_bc(c,ib)  = 0._r8
          albgrd_oc(c,ib)  = 0._r8
          albgri_oc(c,ib)  = 0._r8
          albgrd_dst(c,ib) = 0._r8
          albgri_dst(c,ib) = 0._r8
          do i=-nlevsno+1,1,1
             flx_absdv(c,i) = 0._r8
             flx_absdn(c,i) = 0._r8
             flx_absiv(c,i) = 0._r8
             flx_absin(c,i) = 0._r8
          enddo
       end do

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          albd(p,ib) = 1._r8
          albi(p,ib) = 1._r8
          if (use_SSRE) then
             albdSF(p,ib) = 1._r8
             albiSF(p,ib) = 1._r8
          end if
          fabd(p,ib) = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib) = 0._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 0._r8
       end do

    end do  ! end of numrad loop

    ! SoilAlbedo called before SNICAR_RT
    ! so that reflectance of soil beneath snow column is known
    ! ahead of time for snow RT calculation.

    ! Snow albedos
    ! Note that snow albedo routine will only compute nonzero snow albedos
    ! where h2osno> 0 and coszen > 0

    ! Ground surface albedos
    ! Note that ground albedo routine will only compute nonzero snow albedos
    ! where coszen > 0

    call SoilAlbedo(bounds, &
         num_nourbanc, filter_nourbanc, &
         coszen_col(bounds%begc:bounds%endc), &
         albsnd(bounds%begc:bounds%endc, :), &
         albsni(bounds%begc:bounds%endc, :), &  
         lakestate_inst, temperature_inst, waterstatebulk_inst, surfalb_inst)

    ! set variables to pass to SNICAR.

    flg_snw_ice = 1   ! calling from CLM, not CSIM
    do c=bounds%begc,bounds%endc
       albsfc(c,:)     = albsoi(c,:)
       h2osno_liq(c,:) = h2osoi_liq(c,-nlevsno+1:0)
       h2osno_ice(c,:) = h2osoi_ice(c,-nlevsno+1:0)
       snw_rds_in(c,:) = nint(snw_rds(c,:))
    end do

    ! zero aerosol input arrays
    do aer = 1, sno_nbr_aer
       do i = -nlevsno+1, 0
          do c = bounds%begc, bounds%endc
             mss_cnc_aer_in_frc_pur(c,i,aer) = 0._r8
             mss_cnc_aer_in_frc_bc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_oc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_dst(c,i,aer) = 0._r8
             mss_cnc_aer_in_fdb(c,i,aer)     = 0._r8
          end do
       end do
    end do

    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if (DO_SNO_AER) then
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)

       ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
       !  1) Knowledge of their optical properties is primitive
       !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow,
       !     it has a negligible darkening effect.
       if (DO_SNO_OC) then
          mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
    endif

    ! If radiative forcing is being calculated, first estimate clean-snow albedo

    if (use_snicar_frc) then
       ! 1. BC input array:
       !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       ! BC FORCING CALCULATIONS
          flg_slr = 1; ! direct-beam
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsnd_bc(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)

          flg_slr = 2; ! diffuse
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsni_bc(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)

       ! 2. OC input array:
       !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)

          ! OC FORCING CALCULATIONS
             flg_slr = 1; ! direct-beam
          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                         coszen_col(bounds%begc:bounds%endc), &
                         flg_slr, &
                         h2osno_liq(bounds%begc:bounds%endc, :), &
                         h2osno_ice(bounds%begc:bounds%endc, :), &
                         snw_rds_in(bounds%begc:bounds%endc, :), &
                         mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                         albsfc(bounds%begc:bounds%endc, :), &
                         albsnd_oc(bounds%begc:bounds%endc, :), &
                  foo_snw(bounds%begc:bounds%endc, :, :), &
                  waterstatebulk_inst, waterdiagnosticbulk_inst)

             flg_slr = 2; ! diffuse
          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                         coszen_col(bounds%begc:bounds%endc), &
                         flg_slr, &
                         h2osno_liq(bounds%begc:bounds%endc, :), &
                         h2osno_ice(bounds%begc:bounds%endc, :), &
                         snw_rds_in(bounds%begc:bounds%endc, :), &
                         mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                         albsfc(bounds%begc:bounds%endc, :), &
                         albsni_oc(bounds%begc:bounds%endc, :), &
                  foo_snw(bounds%begc:bounds%endc, :, :), &
                  waterstatebulk_inst, waterdiagnosticbulk_inst)
       endif

       ! 3. DUST input array:
       ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
       mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       ! DUST FORCING CALCULATIONS
          flg_slr = 1; ! direct-beam
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsnd_dst(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)

          flg_slr = 2; ! diffuse
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsni_dst(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)

       ! 4. ALL AEROSOL FORCING CALCULATION
       ! (pure snow albedo)
          flg_slr = 1; ! direct-beam
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsnd_pur(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)

          flg_slr = 2; ! diffuse
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsni_pur(bounds%begc:bounds%endc, :), &
               foo_snw(bounds%begc:bounds%endc, :, :), &
               waterstatebulk_inst, waterdiagnosticbulk_inst)
    end if

    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
       flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                   coszen_col(bounds%begc:bounds%endc), &
                   flg_slr, &
                   h2osno_liq(bounds%begc:bounds%endc, :), &
                   h2osno_ice(bounds%begc:bounds%endc, :), &
                   snw_rds_in(bounds%begc:bounds%endc, :), &
                   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                   albsfc(bounds%begc:bounds%endc, :), &
                   albsnd(bounds%begc:bounds%endc, :), &
            flx_absd_snw(bounds%begc:bounds%endc, :, :), &
            waterstatebulk_inst, waterdiagnosticbulk_inst)

       flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                   coszen_col(bounds%begc:bounds%endc), &
                   flg_slr, &
                   h2osno_liq(bounds%begc:bounds%endc, :), &
                   h2osno_ice(bounds%begc:bounds%endc, :), &
                   snw_rds_in(bounds%begc:bounds%endc, :), &
                   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                   albsfc(bounds%begc:bounds%endc, :), &
                   albsni(bounds%begc:bounds%endc, :), &
            flx_absi_snw(bounds%begc:bounds%endc, :, :), &
            waterstatebulk_inst, waterdiagnosticbulk_inst)

    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
             if (coszen_col(c) > 0._r8) then
             ! ground albedo was originally computed in SoilAlbedo, but is now computed here
             ! because the order of SoilAlbedo and SNICAR_RT was switched for SNICAR.
             albgrd(c,ib) = albsod(c,ib)*(1._r8-frac_sno(c)) + albsnd(c,ib)*frac_sno(c)
             albgri(c,ib) = albsoi(c,ib)*(1._r8-frac_sno(c)) + albsni(c,ib)*frac_sno(c)

             ! albedos for radiative forcing calculations:
             if (use_snicar_frc) then
                ! BC forcing albedo
                albgrd_bc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_bc(c,ib)*frac_sno(c)
                albgri_bc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_bc(c,ib)*frac_sno(c)

                if (DO_SNO_OC) then
                   ! OC forcing albedo
                   albgrd_oc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_oc(c,ib)*frac_sno(c)
                   albgri_oc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_oc(c,ib)*frac_sno(c)
                endif

                ! dust forcing albedo
                albgrd_dst(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_dst(c,ib)*frac_sno(c)
                albgri_dst(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_dst(c,ib)*frac_sno(c)

                ! pure snow albedo for all-aerosol radiative forcing
                albgrd_pur(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_pur(c,ib)*frac_sno(c)
                albgri_pur(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_pur(c,ib)*frac_sno(c)
             end if

             ! also in this loop (but optionally in a different loop for vectorized code)
             !  weight snow layer radiative absorption factors based on snow fraction and soil albedo
             !  (NEEDED FOR ENERGY CONSERVATION)
             do i = -nlevsno+1,1,1
              if (subgridflag == 0 .or. lun%itype(col%landunit(c)) == istdlak) then 
                 if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                endif
             else
                if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)
                endif
             endif
             enddo
          endif
       enddo
    enddo

       ! For diagnostics, set snow albedo to spval over non-snow non-urban points 
       ! so that it is not averaged in history buffer (OPTIONAL)
       ! TODO - this is set to 0 not spval - seems wrong since it will be averaged in

    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
             if ((coszen_col(c) > 0._r8) .and. (h2osno(c) > 0._r8)) then
             albsnd_hst(c,ib) = albsnd(c,ib)
             albsni_hst(c,ib) = albsni(c,ib)
          else
             albsnd_hst(c,ib) = 0._r8
             albsni_hst(c,ib) = 0._r8
          endif
       enddo
    enddo

    ! Create solar-vegetated filter for the following calculations

    num_vegsol = 0
    num_novegsol = 0
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
          if (coszen_patch(p) > 0._r8) then
             if ((lun%itype(patch%landunit(p)) == istsoil .or.  &
                  lun%itype(patch%landunit(p)) == istcrop     ) &
                 .and. (elai(p) + esai(p)) > 0._r8) then
                    num_vegsol = num_vegsol + 1
                    filter_vegsol(num_vegsol) = p
             else
                num_novegsol = num_novegsol + 1
                filter_novegsol(num_novegsol) = p
             end if
          end if
    end do

    ! Weight reflectance/transmittance by lai and sai
       ! Only perform on vegetated patches where coszen > 0

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       wl(p) = elai(p) / max( elai(p)+esai(p), mpe )
       ws(p) = esai(p) / max( elai(p)+esai(p), mpe )
    end do

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          rho(p,ib) = max( rhol(patch%itype(p),ib)*wl(p) + rhos(patch%itype(p),ib)*ws(p), mpe )
          tau(p,ib) = max( taul(patch%itype(p),ib)*wl(p) + taus(patch%itype(p),ib)*ws(p), mpe )
       end do
    end do

    ! Diagnose number of canopy layers for radiative transfer, in increments of dincmax.
    ! Add to number of layers so long as cumulative leaf+stem area does not exceed total
    ! leaf+stem area. Then add any remaining leaf+stem area to next layer and exit the loop.
    ! Do this first for elai and esai (not buried by snow) and then for the part of the
    ! canopy that is buried by snow.
    ! ------------------
    ! tlai_z = leaf area increment for a layer
    ! tsai_z = stem area increment for a layer
    ! nrad   = number of canopy layers above snow
    ! ncan   = total number of canopy layers
    !
    ! tlai_z summed from 1 to nrad = elai
    ! tlai_z summed from 1 to ncan = tlai

    ! tsai_z summed from 1 to nrad = esai
    ! tsai_z summed from 1 to ncan = tsai
    ! ------------------
    !
    ! Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol"
    ! because layering is needed for all time steps regardless of radiation
    !
    ! Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered by
    ! nlevcan = 1

    dincmax = 0.25_r8
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)


       if (nlevcan == 1) then
          nrad(p) = 1
          ncan(p) = 1
          tlai_z(p,1) = elai(p)
          tsai_z(p,1) = esai(p)

        !--for the oil palm PFT use dynamic multilayer structure (Y.Fan)--!
        !reference: Fan, Y. et al. A sub-canopy structure for
        !simulating oil palm in the Community Land Model (CLM-Palm):
        !phenology, allocation and yield. Geoscientific Model Development
        !8, 3785–3800 (2015). doi:10.5194/gmd-8-3785-2015

        else if (phytomer(ivt(p)) > 0) then
            if (use_cn) then
                !use pack function to gather indices of phytomers (range 1:mxnp) with positive LAI
                !all dead or too small phytomers are filtered from radiation and photosynthesis calculations
                !the last element of ix is the index of top youngest expanded phytomer
                ixn = count(mask=(rankp(p,:) > 0 .and. plai(p,:) > 1.e-03_r8))
                allocate (ix(1:ixn))
                ix = pack((/1:mxnp:1/), mask=(rankp(p,:) > 0 .and. plai(p,:) > 1.e-03_r8))
                !0.001 is the minimum phytomer LAI but 0.01 might be a safer threshold for both TwoStream and Multilayer
                if (elai(p)+esai(p) == 0._r8 .or. ixn == 0) then
                    nrad(p) = 0
                    !below condition not needed: this will make a jump in GPP when change from one layer to multilayer suddenly
                    !            else if (elai(p)+esai(p) <= 0.25_r8 .or. ixn == 1) then !seedlings one layer
                    !                nrad(p) = 1
                    !                tlai_z(p,1) = elai(p)
                    !                tsai_z(p,1) = esai(p)
                !for oil palm, every 8 phytomers girth the stem at the same height
                else !either each phytomer is a layer; or several phytomers form a layer
                    nump = max(1, NINT(real(mxlivenp(ivt(p)))/real(nlevcan))) !number of phytomers per layer
                    dincmax_sum = 0._r8
                    do iv = 1, nlevcan
                        nrad(p) = iv
                        if (ixn <= nump) then
                            tlai_z(p,1) = tlai(p)
                            tsai_z(p,1) = tsai(p)
                            exit
                        else if ((ixn - nump*iv) > 0 .and. iv < nlevcan) then
                            tlai_z(p,iv) = sum(plai(p,ix((ixn+1-nump*iv):(ixn-nump*(iv-1)))))
                            tsai_z(p,iv) = tlai_z(p,iv) * tsai(p) / max(tlai(p), mpe)
                        else !all the rest leaves go to the last layer
                            tlai_z(p,iv) = tlai(p) - sum(tlai_z(p,1:(iv-1)))
                            tsai_z(p,iv) = tsai(p) - sum(tsai_z(p,1:(iv-1)))
                            exit
                        end if
                        dincmax_sum = dincmax_sum + tlai_z(p,iv)
                        if ((elai(p)-dincmax_sum) < 1.e-03_r8) then
                            if (iv ==1) then
                                tlai_z(p,iv) = elai(p)
                                tsai_z(p,iv) = esai(p)
                            else
                                tlai_z(p,iv) = elai(p) - sum(tlai_z(p,1:(iv-1)))
                                tsai_z(p,iv) = esai(p) - sum(tsai_z(p,1:(iv-1)))
                            end if
                            exit
                        end if
                    end do
                end if
                deallocate(ix)
            !SP mode readin tlai_z/tsai_z from SatellitePhenologyMod
            else
                if (elai(p)+esai(p) == 0._r8) then
                    nrad(p) = 0
                else
                    dincmax_sum = 0._r8
                    do iv = 1, nlevcan
                        nrad(p) = iv
                        dincmax_sum = dincmax_sum + tlai_z(p,iv)
                        if ((elai(p)-dincmax_sum) <= 1.e-03_r8) then
                            if (iv ==1) then
                                tlai_z(p,iv) = elai(p)
                                tsai_z(p,iv) = esai(p)
                            else
                                tlai_z(p,iv) = elai(p) - sum(tlai_z(p,1:(iv-1)))
                                tsai_z(p,iv) = esai(p) - sum(tsai_z(p,1:(iv-1)))
                            end if
                            exit
                        end if
                        if (tlai_z(p,iv) > 0._r8 .and. tlai_z(p,iv) < 1.e-03_r8) then
                            tlai_z(p,iv) = 1.e-03_r8 !enforce minimum lai on live/active layers
                            tsai_z(p,iv) = 0.05_r8*tlai_z(p,iv)
                        end if
                    end do
                end if
                !although palm structure only applicable in tropical climate,
                !the above loops keep formative to consider snowcover
            end if
        !---other PFTs use normal multilayer with uniform layer thickness--!
       else if (nlevcan > 1) then
          if (elai(p)+esai(p) == 0._r8) then
             nrad(p) = 0
          else
             dincmax_sum = 0._r8
             do iv = 1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((elai(p)+esai(p))-dincmax_sum) > 1.e-06_r8) then
                   nrad(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                else
                   nrad(p) = iv
                   dinc = dincmax - (dincmax_sum - (elai(p)+esai(p)))
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                   exit
                end if
             end do



             ! Mimumum of 4 canopy layers

             if (nrad(p) < 4) then
                nrad(p) = 4
                do iv = 1, nrad(p)
                   tlai_z(p,iv) = elai(p) / nrad(p)
                   tsai_z(p,iv) = esai(p) / nrad(p)
                end do
             end if
          end if
       end if

       ! Error check: make sure cumulative of increments does not exceed total

       laisum = 0._r8
       saisum = 0._r8
       do iv = 1, nrad(p)
          laisum = laisum + tlai_z(p,iv)
          saisum = saisum + tsai_z(p,iv)
       end do
       !only do check when nrad(p) > 0 (Y.Fan)
       ! if (nrad(p) > 0 .and. (abs(laisum-elai(p)) > 1.e-06_r8 .or. abs(saisum-esai(p)) > 1.e-06_r8)) then

       if (abs(laisum-elai(p)) > 1.e-06_r8 .or. abs(saisum-esai(p)) > 1.e-06_r8) then
          write (iulog,*) 'multi-layer canopy error 01 in SurfaceAlbedo: ',&
               nrad(p),elai(p),laisum,esai(p),saisum
          call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
       end if

       ! Repeat to find canopy layers buried by snow

       if (nlevcan > 1) then
          blai(p) = tlai(p) - elai(p)
          bsai(p) = tsai(p) - esai(p)
          if (blai(p)+bsai(p) == 0._r8) then
             ncan(p) = nrad(p)
          else
             dincmax_sum = 0._r8
             do iv = nrad(p)+1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((blai(p)+bsai(p))-dincmax_sum) > 1.e-06_r8) then
                   ncan(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                else
                   ncan(p) = iv
                   dinc = dincmax - (dincmax_sum - (blai(p)+bsai(p)))
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                   exit
                end if
             end do
          end if

          ! Error check: make sure cumulative of increments does not exceed total

          laisum = 0._r8
          saisum = 0._r8
          do iv = 1, ncan(p)
             laisum = laisum + tlai_z(p,iv)
             saisum = saisum + tsai_z(p,iv)
          end do
          !only do check when nrad(p) > 0 (Y.Fan)
          !if (nrad(p) > 0 .and. (abs(laisum-tlai(p)) > 1.e-06_r8 .or. abs(saisum-tsai(p)) > 1.e-06_r8)) then
          if (abs(laisum-tlai(p)) > 1.e-06_r8 .or. abs(saisum-tsai(p)) > 1.e-06_r8) then
             write (iulog,*) 'multi-layer canopy error 02 in SurfaceAlbedo: ',nrad(p),ncan(p)
             write (iulog,*) tlai(p),elai(p),blai(p),laisum,tsai(p),esai(p),bsai(p),saisum
             call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

    end do

    ! Zero fluxes for active canopy layers

    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       do iv = 1, nrad(p)
          fabd_sun_z(p,iv) = 0._r8
          fabd_sha_z(p,iv) = 0._r8
          fabi_sun_z(p,iv) = 0._r8
          fabi_sha_z(p,iv) = 0._r8
          fsun_z(p,iv) = 0._r8
       end do
    end do

    ! Default leaf to canopy scaling coefficients, used when coszen <= 0.
    ! This is the leaf nitrogen profile integrated over the full canopy.
    ! Integrate exp(-kn*x) over x=0 to x=elai and assign to shaded canopy,
    ! because sunlit fraction is 0. Canopy scaling coefficients are set in
    ! TwoStream for coszen > 0. So kn must be set here and in TwoStream.

    extkn = 0.30_r8
    if (phytomer(ivt(p)) > 0) then !for palm PFT no need Kn (Y.Fan 2015)

       extkn = 0.01_r8 ! use 0.01_r8 to turn off N scaling for two-stream


    end if  !Kn for multilayer canopy model is set by kn_leaf in CanopyFluxesMod
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       if (nlevcan == 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn
          if (elai(p) > 0._r8) then
             vcmaxcintsha(p) = vcmaxcintsha(p) / elai(p)
          else
             vcmaxcintsha(p) = 0._r8
          end if
       !Kn for multilayer canopy model is set by kn_leaf in CanopyFluxesMod
       else if (nlevcan > 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = 0._r8
       end if
    end do

    ! Calculate surface albedos and fluxes
    ! Only perform on vegetated pfts where coszen > 0

    if (use_fates) then
          
       call clm_fates%wrap_canopy_radiation(bounds, nc, &
            num_vegsol, filter_vegsol, &
            coszen_patch(bounds%begp:bounds%endp), surfalb_inst)

    else if (radiative_transfer == 1) then
       !Multilayer radiative transfer scheme (Norman) (Y.Fan)
       !use a G_Function for multilayer canopy
       call Gfunc_all (bounds, filter_vegsol, num_vegsol, &
            coszen_patch(bounds%begp:bounds%endp), &
            surfalb_inst)
       call Multilayer (bounds, filter_vegsol, num_vegsol, &
            coszen_patch(bounds%begp:bounds%endp), &
            rho(bounds%begp:bounds%endp, :), &
            tau(bounds%begp:bounds%endp, :), &
            atm2lnd_inst, temperature_inst, waterdiagnosticbulk_inst, surfalb_inst)
    else

       call TwoStream (bounds, filter_vegsol, num_vegsol, &
            coszen_patch(bounds%begp:bounds%endp), &
            rho(bounds%begp:bounds%endp, :), &
            tau(bounds%begp:bounds%endp, :), &
            canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, surfalb_inst)
       ! Run TwoStream again just to calculate the Snow Free (SF) albedo's
       if (use_SSRE) then
          if ( nlevcan > 1 )then
             call endrun( 'ERROR: use_ssre option was NOT developed with allowance for multi-layer canopy: '// &
                          'nlevcan can ONLY be 1 in when use_ssre is on')
          end if
          call TwoStream (bounds, filter_vegsol, num_vegsol, &
               coszen_patch(bounds%begp:bounds%endp), &
               rho(bounds%begp:bounds%endp, :), &
               tau(bounds%begp:bounds%endp, :), &
               canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, surfalb_inst, &
               SFonly=.true.)
       end if
       
    endif

    ! Determine values for non-vegetated patches where coszen > 0

    do ib = 1,numrad
       do fp = 1,num_novegsol
          p = filter_novegsol(fp)
          c = patch%column(p)
          fabd(p,ib)     = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib)     = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib)     = 1._r8
          ftid(p,ib)     = 0._r8
          ftii(p,ib)     = 1._r8
          albd(p,ib)     = albgrd(c,ib)
          albi(p,ib)     = albgri(c,ib)
          if (use_SSRE) then
             albdSF(p,ib)    = albsod(c,ib)
             albiSF(p,ib)    = albsoi(c,ib)
          end if
       end do
    end do

     end associate

   end subroutine SurfaceAlbedo

   !-----------------------------------------------------------------------
   subroutine SoilAlbedo (bounds, &
        num_nourbanc, filter_nourbanc, &
        coszen, albsnd, albsni, &
        lakestate_inst, temperature_inst, waterstatebulk_inst, surfalb_inst)
     !
     ! !DESCRIPTION:
     ! Determine ground surface albedo, accounting for snow
     !
     ! !USES:
    use clm_varpar      , only : numrad
    use clm_varcon      , only : tfrz
    use landunit_varcon , only : istice_mec, istdlak
    use LakeCon         , only : lakepuddling
    !
    ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds             
     integer , intent(in) :: num_nourbanc               ! number of columns in non-urban points in column filter
     integer , intent(in) :: filter_nourbanc(:)          ! column filter for non-urban points
     real(r8), intent(in) :: coszen( bounds%begc: )      ! cos solar zenith angle next time step [col]
     real(r8), intent(in) :: albsnd( bounds%begc: , 1: ) ! snow albedo (direct) [col, numrad]
     real(r8), intent(in) :: albsni( bounds%begc: , 1: ) ! snow albedo (diffuse) [col, numrad]
     type(temperature_type) , intent(in)    :: temperature_inst
     type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
     type(lakestate_type)   , intent(in)    :: lakestate_inst
     type(surfalb_type)     , intent(inout) :: surfalb_inst
     !
     ! !LOCAL VARIABLES:
     !
    integer, parameter :: nband =numrad ! number of solar radiation waveband classes
    integer  :: fc            ! non-urban filter column index
    integer  :: c,l           ! indices
    integer  :: ib            ! waveband number (1=vis, 2=nir)
    real(r8) :: inc           ! soil water correction factor for soil albedo
    integer  :: soilcol       ! soilcolor
    real(r8) :: sicefr        ! Lake surface ice fraction (based on D. Mironov 2010)
    !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endc/)),         errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(albsnd) == (/bounds%endc, numrad/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(albsni) == (/bounds%endc, numrad/)), errMsg(sourcefile, __LINE__))

   associate(&
          snl          => col%snl                         , & ! Input:  [integer  (:)   ]  number of snow layers                    

          t_grnd       => temperature_inst%t_grnd_col     , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             

          h2osoi_vol   => waterstatebulk_inst%h2osoi_vol_col  , & ! Input:  [real(r8) (:,:) ]  volumetric soil water [m3/m3]         

          lake_icefrac => lakestate_inst%lake_icefrac_col , & ! Input:  [real(r8) (:,:) ]  mass fraction of lake layer that is frozen
          
          albgrd       => surfalb_inst%albgrd_col         , & ! Output: [real(r8) (:,:) ]  ground albedo (direct)                
          albgri       => surfalb_inst%albgri_col         , & ! Output: [real(r8) (:,:) ]  ground albedo (diffuse)               
          albsod       => surfalb_inst%albsod_col         , & ! Output: [real(r8) (:,:) ]  soil albedo (direct)                  
          albsoi       => surfalb_inst%albsoi_col           & ! Output: [real(r8) (:,:) ]  soil albedo (diffuse)                 
   )

    ! Compute soil albedos

    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          if (coszen(c) > 0._r8) then
             l = col%landunit(c)

             if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)  then ! soil
                inc    = max(0.11_r8-0.40_r8*h2osoi_vol(c,1), 0._r8)
                soilcol = isoicol(c)
                ! changed from local variable to clm_type:
                !albsod = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                !albsoi = albsod
                albsod(c,ib) = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                albsoi(c,ib) = albsod(c,ib)
             else if (lun%itype(l) == istice_mec)  then  ! land ice
                ! changed from local variable to clm_type:
                !albsod = albice(ib)
                !albsoi = albsod
                albsod(c,ib) = albice(ib)
                albsoi(c,ib) = albsod(c,ib)
             ! unfrozen lake, wetland
             else if (t_grnd(c) > tfrz .or. (lakepuddling .and. lun%itype(l) == istdlak .and. t_grnd(c) == tfrz .and. &
                      lake_icefrac(c,1) < 1._r8 .and. lake_icefrac(c,2) > 0._r8) ) then

                albsod(c,ib) = 0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8)
                ! This expression is apparently from BATS according to Yongjiu Dai.

                ! The diffuse albedo should be an average over the whole sky of an angular-dependent direct expression.
                ! The expression above may have been derived to encompass both (e.g. Henderson-Sellers 1986),
                ! but I'll assume it applies more appropriately to the direct form for now.

                ! ZMS: Attn EK, currently restoring this for wetlands even though it is wrong in order to try to get
                ! bfb baseline comparison when no lakes are present. I'm assuming wetlands will be phased out anyway.
                if (lun%itype(l) == istdlak) then
                   albsoi(c,ib) = 0.10_r8
                else
                   albsoi(c,ib) = albsod(c,ib)
                end if

             else                                     ! frozen lake, wetland
                ! Introduce crude surface frozen fraction according to D. Mironov (2010)
                ! Attn EK: This formulation is probably just as good for "wetlands" if they are not phased out.
                ! Tenatively I'm restricting this to lakes because I haven't tested it for wetlands. But if anything
                ! the albedo should be lower when melting over frozen ground than a solid frozen lake.
                !
                if (lun%itype(l) == istdlak .and. .not. lakepuddling .and. snl(c) == 0) then
                    ! Need to reference snow layers here because t_grnd could be over snow or ice
                                      ! but we really want the ice surface temperature with no snow
                   sicefr = 1._r8 - exp(-calb * (tfrz - t_grnd(c))/tfrz)
                   albsod(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), &
                                  0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8))
                   albsoi(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), 0.10_r8)
                   ! Make sure this is no less than the open water albedo above.
                   ! Setting lake_melt_icealb(:) = alblak(:) in namelist reverts the melting albedo to the cold
                   ! snow-free value.
                else
                   albsod(c,ib) = alblak(ib)
                   albsoi(c,ib) = albsod(c,ib)
                end if
             end if

             ! Weighting is done in SurfaceAlbedo, after the call to SNICAR_RT
             ! This had to be done, because SoilAlbedo is called before SNICAR_RT, so at
             ! this point, snow albedo is not yet known.
          end if
       end do
    end do

    end associate
   end subroutine SoilAlbedo

   !-----------------------------------------------------------------------
   subroutine TwoStream (bounds, &
        filter_vegsol, num_vegsol, &
        coszen, rho, tau, &
        canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, surfalb_inst, &
        SFonly)
     !
     ! !DESCRIPTION:
     ! Two-stream fluxes for canopy radiative transfer
     ! Use two-stream approximation of Dickinson (1983) Adv Geophysics
     ! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
     ! to calculate fluxes absorbed by vegetation, reflected by vegetation,
     ! and transmitted through vegetation for unit incoming direct or diffuse
     ! flux given an underlying surface with known albedo.
     ! Calculate sunlit and shaded fluxes as described by
     ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
     ! a multi-layer canopy to calculate APAR profile
     !
     ! !USES:
     use clm_varpar, only : numrad, nlevcan
     use clm_varcon, only : omegas, tfrz, betads, betais
     use clm_varctl, only : iulog
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds           
     integer                , intent(in)    :: filter_vegsol (:)        ! filter for vegetated patches with coszen>0
     integer                , intent(in)    :: num_vegsol               ! number of vegetated patches where coszen>0
     real(r8), intent(in)  :: coszen( bounds%begp: )   ! cosine solar zenith angle for next time step [pft]
     real(r8), intent(in)  :: rho( bounds%begp: , 1: ) ! leaf/stem refl weighted by fraction LAI and SAI [pft, numrad]
     real(r8), intent(in)  :: tau( bounds%begp: , 1: ) ! leaf/stem tran weighted by fraction LAI and SAI [pft, numrad]
     type(canopystate_type) , intent(in)    :: canopystate_inst
     type(temperature_type) , intent(in)    :: temperature_inst
     type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
     type(surfalb_type)     , intent(inout) :: surfalb_inst
     logical, optional      , intent(in)    :: SFonly                              ! If should just calculate the Snow Free albedos
     !
     ! !LOCAL VARIABLES:
     integer  :: fp,p,c,iv        ! array indices
     integer  :: ib               ! waveband number
     real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
     real(r8) :: asu              ! single scattering albedo
     real(r8) :: chil(bounds%begp:bounds%endp)    ! -0.4 <= xl <= 0.6
     real(r8) :: gdir(bounds%begp:bounds%endp)    ! leaf projection in solar direction (0 to 1)
     real(r8) :: twostext(bounds%begp:bounds%endp)! optical depth of direct beam per unit leaf area
     real(r8) :: avmu(bounds%begp:bounds%endp)    ! average diffuse optical depth
     real(r8) :: omega(bounds%begp:bounds%endp,numrad)   ! fraction of intercepted radiation that is scattered (0 to 1)
     real(r8) :: omegal           ! omega for leaves
     real(r8) :: betai            ! upscatter parameter for diffuse radiation
     real(r8) :: betail           ! betai for leaves
     real(r8) :: betad            ! upscatter parameter for direct beam radiation
     real(r8) :: betadl           ! betad for leaves
     real(r8) :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
     real(r8) :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
     real(r8) :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
     real(r8) :: phi1,phi2,sigma                                   ! temporary
     real(r8) :: temp1                                             ! temporary
     real(r8) :: temp0    (bounds%begp:bounds%endp)                ! temporary
     real(r8) :: temp2(bounds%begp:bounds%endp)                    ! temporary
     real(r8) :: t1                                                ! temporary
     real(r8) :: a1,a2                                             ! parameter for sunlit/shaded leaf radiation absorption
     real(r8) :: v,dv,u,du                                         ! temporary for flux derivatives
     real(r8) :: dh2,dh3,dh5,dh6,dh7,dh8,dh9,dh10                  ! temporary for flux derivatives
     real(r8) :: da1,da2                                           ! temporary for flux derivatives
     real(r8) :: d_ftid,d_ftii                                     ! ftid, ftii derivative with respect to lai+sai
     real(r8) :: d_fabd,d_fabi                                     ! fabd, fabi derivative with respect to lai+sai
     real(r8) :: d_fabd_sun,d_fabd_sha                             ! fabd_sun, fabd_sha derivative with respect to lai+sai
     real(r8) :: d_fabi_sun,d_fabi_sha                             ! fabi_sun, fabi_sha derivative with respect to lai+sai
     real(r8) :: laisum                                            ! cumulative lai+sai for canopy layer (at middle of layer)
     real(r8) :: extkb                                             ! direct beam extinction coefficient
     real(r8) :: extkn                                             ! nitrogen allocation coefficient
     logical  :: lSFonly                                           ! Local version of SFonly (Snow Free) flag
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endp/)),         errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(rho)    == (/bounds%endp, numrad/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(tau)    == (/bounds%endp, numrad/)), errMsg(sourcefile, __LINE__))

     if ( present(SFonly) )then
        lSFonly = SFonly
     else
        lSFonly = .false.
     end if

   associate(&
          xl           =>    pftcon%xl                           , & ! Input:  ecophys const - leaf/stem orientation index
          ivt          =>    patch%itype                         , & ! Input:  [integer (:)]  pft vegetation type
          phytomer     =>    pftcon%phytomer                     , & ! Input:  [integer (:)]  number of phytomers for palm PFT (Y.Fan)

          t_veg        =>    temperature_inst%t_veg_patch        , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)         

          fwet         =>    waterdiagnosticbulk_inst%fwet_patch          , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1) 
          fcansno      =>    waterdiagnosticbulk_inst%fcansno_patch       , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is snow-covered (0 to 1) 

          elai         =>    canopystate_inst%elai_patch         , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
          esai         =>    canopystate_inst%esai_patch         , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow

          tlai_z       =>    surfalb_inst%tlai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tlai increment for canopy layer       
          tsai_z       =>    surfalb_inst%tsai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tsai increment for canopy layer       
          nrad         =>    surfalb_inst%nrad_patch             , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          albgrd       =>    surfalb_inst%albgrd_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (direct) (column-level) 
          albgri       =>    surfalb_inst%albgri_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (diffuse)(column-level) 

          ! For non-Snow Free
          fsun_z       =>    surfalb_inst%fsun_z_patch           , & ! Output: [real(r8) (:,:) ]  sunlit fraction of canopy layer       
          vcmaxcintsun =>    surfalb_inst%vcmaxcintsun_patch     , & ! Output: [real(r8) (:)   ]  leaf to canopy scaling coefficient, sunlit leaf vcmax
          vcmaxcintsha =>    surfalb_inst%vcmaxcintsha_patch     , & ! Output: [real(r8) (:)   ]  leaf to canopy scaling coefficient, shaded leaf vcmax
          fabd_sun_z   =>    surfalb_inst%fabd_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z   =>    surfalb_inst%fabd_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z   =>    surfalb_inst%fabi_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z   =>    surfalb_inst%fabi_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          albd         =>    surfalb_inst%albd_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (direct)               
          albi         =>    surfalb_inst%albi_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (diffuse)              
          fabd         =>    surfalb_inst%fabd_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
          fabd_sun     =>    surfalb_inst%fabd_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
          fabd_sha     =>    surfalb_inst%fabd_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
          fabi         =>    surfalb_inst%fabi_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
          fabi_sun     =>    surfalb_inst%fabi_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha     =>    surfalb_inst%fabi_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
          ftdd         =>    surfalb_inst%ftdd_patch             , & ! Output: [real(r8) (:,:) ]  down direct flux below canopy per unit direct flx
          ftid         =>    surfalb_inst%ftid_patch             , & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flx
          ftii         =>    surfalb_inst%ftii_patch             , & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flx

          ! Needed for SF Snow free case
          albsod       =>    surfalb_inst%albsod_col             , & ! Input: [real(r8)  (:,:) ]  soil albedo (direct)
          albsoi       =>    surfalb_inst%albsoi_col             , & ! Input: [real(r8)  (:,:) ]  soil albedo (diffuse)
          albdSF       =>    surfalb_inst%albdSF_patch           , & ! Output: [real(r8) (:,:) ]  Snow Free surface albedo (direct)
          albiSF       =>    surfalb_inst%albiSF_patch             & ! Output: [real(r8) (:,:) ]  Snow Free surface albedo (diffuse)
   )

    ! Calculate two-stream parameters that are independent of waveband:
    ! chil, gdir, twostext, avmu, and temp0 and temp2 (used for asu)

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)

       ! note that the following limit only acts on cosz values > 0 and less than
       ! 0.001, not on values cosz = 0, since these zero have already been filtered
       ! out in filter_vegsol
       cosz = max(0.001_r8, coszen(p))

       chil(p) = min( max(xl(patch%itype(p)), -0.4_r8), 0.6_r8 )
       if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8
       phi1 = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
       phi2 = 0.877_r8 * (1._r8-2._r8*phi1)
       gdir(p) = phi1 + phi2*cosz
       twostext(p) = gdir(p)/cosz
       avmu(p) = ( 1._r8 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
       ! Restrict this calculation of temp0. We have seen cases where small temp0
       ! can cause unrealistic single scattering albedo (asu) associated with the
       ! log calculation in temp2 below, thereby eventually causing a negative soil albedo
       ! See bugzilla bug 2431: http://bugs.cgd.ucar.edu/show_bug.cgi?id=2431
       temp0(p) = max(gdir(p) + phi2*cosz,1.e-6_r8)
       temp1 = phi1*cosz
       temp2(p) = ( 1._r8 - temp1/temp0(p) * log((temp1+temp0(p))/temp1) )
    end do

   ! Loop over all wavebands to calculate for the full canopy the scattered fluxes
   ! reflected upward and transmitted downward by the canopy and the flux absorbed by the
   ! canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given
   ! an underlying surface of known albedo.
   !
   ! Output:
   ! ------------------
   ! Direct beam fluxes
   ! ------------------
   ! albd       - Upward scattered flux above canopy (per unit direct beam flux)
   ! ftid       - Downward scattered flux below canopy (per unit direct beam flux)
   ! ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
   ! fabd       - Flux absorbed by canopy (per unit direct beam flux)
   ! fabd_sun   - Sunlit portion of fabd
   ! fabd_sha   - Shaded portion of fabd
   ! fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
   ! ------------------
   ! Diffuse fluxes
   ! ------------------
   ! albi       - Upward scattered flux above canopy (per unit diffuse flux)
   ! ftii       - Downward scattered flux below canopy (per unit diffuse flux)
   ! fabi       - Flux absorbed by canopy (per unit diffuse flux)
   ! fabi_sun   - Sunlit portion of fabi
   ! fabi_sha   - Shaded portion of fabi
   ! fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          c = patch%column(p)

          ! Calculate two-stream parameters omega, betad, and betai.
          ! Omega, betad, betai are adjusted for snow. Values for omega*betad
          ! and omega*betai are calculated and then divided by the new omega
          ! because the product omega*betai, omega*betad is used in solution.
          ! Also, the transmittances and reflectances (tau, rho) are linear
          ! weights of leaf and stem values.

          omegal = rho(p,ib) + tau(p,ib)
          asu = 0.5_r8*omegal*gdir(p)/temp0(p) *temp2(p)
          betadl = (1._r8+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu
          betail = 0.5_r8 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) &
                 * ((1._r8+chil(p))/2._r8)**2) / omegal

          if ( lSFonly .or. ( (.not. snowveg_affects_radiation) .and. (t_veg(p) > tfrz) ) ) then
             ! Keep omega, betad, and betai as they are (for Snow free case or
             ! when there is no snow
             tmp0 = omegal
             tmp1 = betadl
             tmp2 = betail
          else
             ! Adjust omega, betad, and betai for intercepted snow
             if (snowveg_affects_radiation) then
                tmp0 =   (1._r8-fcansno(p))*omegal        + fcansno(p)*omegas(ib)
                tmp1 = ( (1._r8-fcansno(p))*omegal*betadl + fcansno(p)*omegas(ib)*betads ) / tmp0
                tmp2 = ( (1._r8-fcansno(p))*omegal*betail + fcansno(p)*omegas(ib)*betais ) / tmp0
             else
                tmp0 =   (1._r8-fwet(p))*omegal        + fwet(p)*omegas(ib)
                tmp1 = ( (1._r8-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
                tmp2 = ( (1._r8-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
             end if
          end if  ! end Snow free

          omega(p,ib) = tmp0
          betad = tmp1
          betai = tmp2

          ! Common terms

          b = 1._r8 - omega(p,ib) + omega(p,ib)*betai
          c1 = omega(p,ib)*betai
          tmp0 = avmu(p)*twostext(p)
          d = tmp0 * omega(p,ib)*betad
          f = tmp0 * omega(p,ib)*(1._r8-betad)
          tmp1 = b*b - c1*c1
          h = sqrt(tmp1) / avmu(p)
          sigma = tmp0*tmp0 - tmp1
          p1 = b + avmu(p)*h
          p2 = b - avmu(p)*h
          p3 = b + tmp0
          p4 = b - tmp0

          ! Absorbed, reflected, transmitted fluxes per unit incoming radiation
          ! for full canopy

          t1 = min(h*(elai(p)+esai(p)), 40._r8)
          s1 = exp(-t1)
          t1 = min(twostext(p)*(elai(p)+esai(p)), 40._r8)
          s2 = exp(-t1)

          ! Direct beam
          if ( .not. lSFonly )then
             u1 = b - c1/albgrd(c,ib)
             u2 = b - c1*albgrd(c,ib)
             u3 = f + c1*albgrd(c,ib)
          else
             ! Snow Free (SF) only 
             ! albsod instead of albgrd here:
             u1 = b - c1/albsod(c,ib)
             u2 = b - c1*albsod(c,ib)
             u3 = f + c1*albsod(c,ib)
          end if
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
          if ( .not. lSFonly )then
            albd(p,ib) = h1/sigma + h2 + h3
            ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1
            ftdd(p,ib) = s2
            fabd(p,ib) = 1._r8 - albd(p,ib) - (1._r8-albgrd(c,ib))*ftdd(p,ib) - (1._r8-albgri(c,ib))*ftid(p,ib)
          else
            albdSF(p,ib) = h1/sigma + h2 + h3
          end if
          

          a1 = h1 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h2         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h3         * (1._r8 - s2/s1) / (twostext(p) - h)

          a2 = h4 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h5         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h6         * (1._r8 - s2/s1) / (twostext(p) - h)
          if ( .not. lSFonly )then
            fabd_sun(p,ib) = (1._r8 - omega(p,ib)) * ( 1._r8 - s2 + 1._r8 / avmu(p) * (a1 + a2) )
            fabd_sha(p,ib) = fabd(p,ib) - fabd_sun(p,ib)
          end if

          ! Diffuse
          if ( .not. lSFonly )then
            u1 = b - c1/albgri(c,ib)
            u2 = b - c1*albgri(c,ib)
          else
             ! Snow Free (SF) only 
             ! albsoi instead of albgri here:
            u1 = b - c1/albsoi(c,ib)
            u2 = b - c1*albsoi(c,ib)
          end if
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

  
          ! Final Snow Free albedo
          if ( lSFonly )then
            albiSF(p,ib) = h7 + h8
          else
            ! For non snow Free case, adjustments continue
            albi(p,ib) = h7 + h8
            ftii(p,ib) = h9*s1 + h10/s1
            fabi(p,ib) = 1._r8 - albi(p,ib) - (1._r8-albgri(c,ib))*ftii(p,ib)
            a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
            a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)

            fabi_sun(p,ib) = (1._r8 - omega(p,ib)) / avmu(p) * (a1 + a2)
            fabi_sha(p,ib) = fabi(p,ib) - fabi_sun(p,ib)
  
            ! Repeat two-stream calculations for each canopy layer to calculate derivatives.
            ! tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are
            ! calculated at the center of the layer. Derivatives are needed only for the
            ! visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
            ! Derivatives are calculated first per unit lai+sai and then normalized for sunlit
            ! or shaded fraction of canopy layer.
  
            ! Sun/shade big leaf code uses only one layer, with canopy integrated values from above
            ! and also canopy-integrated scaling coefficients
  
            if (ib == 1) then
               if (nlevcan == 1) then
  
                  ! sunlit fraction of canopy
                  fsun_z(p,1) = (1._r8 - s2) / t1
  
                  ! absorbed PAR (per unit sun/shade lai+sai)
                  laisum = elai(p)+esai(p)
                  fabd_sun_z(p,1) = fabd_sun(p,ib) / (fsun_z(p,1)*laisum)
                  fabi_sun_z(p,1) = fabi_sun(p,ib) / (fsun_z(p,1)*laisum)
                  fabd_sha_z(p,1) = fabd_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)
                  fabi_sha_z(p,1) = fabi_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)
  
                  ! leaf to canopy scaling coefficients
                  extkn = 0.30_r8
                  if (phytomer(ivt(p)) > 0) then !for palm PFT no need Kn (Y.Fan 2015)
                     extkn = 0.01_r8 ! use 0.01_r8 to turn off N scaling
                  end if
                  extkb = twostext(p)
                  vcmaxcintsun(p) = (1._r8 - exp(-(extkn+extkb)*elai(p))) / (extkn + extkb)
                  vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn - vcmaxcintsun(p)
                  if (elai(p)  >  0._r8) then
                    vcmaxcintsun(p) = vcmaxcintsun(p) / (fsun_z(p,1)*elai(p))
                    vcmaxcintsha(p) = vcmaxcintsha(p) / ((1._r8 - fsun_z(p,1))*elai(p))
                  else
                    vcmaxcintsun(p) = 0._r8
                    vcmaxcintsha(p) = 0._r8
                  end if
  
               else if (nlevcan > 1)then
                  do iv = 1, nrad(p)
  
                     ! Cumulative lai+sai at center of layer
  
                     if (iv == 1) then
                        laisum = 0.5_r8 * (tlai_z(p,iv)+tsai_z(p,iv))
                     else
                        laisum = laisum + 0.5_r8 * ((tlai_z(p,iv-1)+tsai_z(p,iv-1))+(tlai_z(p,iv)+tsai_z(p,iv)))
                     end if
  
                     ! Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction
     
                     t1 = min(h*laisum, 40._r8)
                     s1 = exp(-t1)
                     t1 = min(twostext(p)*laisum, 40._r8)
                     s2 = exp(-t1)
                     fsun_z(p,iv) = s2
  
                     ! ===============
                     ! Direct beam
                     ! ===============
  
                     ! Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai
  
                     u1 = b - c1/albgrd(c,ib)
                     u2 = b - c1*albgrd(c,ib)
                     u3 = f + c1*albgrd(c,ib)
  
                     ! Derivatives for h2, h3, h5, h6 and a1, a2
  
                     v = d1
                     dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1
  
                     u = tmp6 * tmp2 / s1 - p2 * tmp7
                     du = h * tmp6 * tmp2 / s1 + twostext(p) * p2 * tmp7
                     dh2 = (v * du - u * dv) / (v * v)
  
                     u = -tmp6 * tmp3 * s1 + p1 * tmp7
                     du = h * tmp6 * tmp3 * s1 - twostext(p) * p1 * tmp7
                     dh3 = (v * du - u * dv) / (v * v)
  
                     v = d2
                     dv = h * tmp4 / s1 + h * tmp5 * s1
     
                     u = -h4/sigma * tmp4 / s1 - tmp9
                     du = -h * h4/sigma * tmp4 / s1 + twostext(p) * tmp9
                     dh5 = (v * du - u * dv) / (v * v)
  
                     u = h4/sigma * tmp5 * s1 + tmp9
                     du = -h * h4/sigma * tmp5 * s1 - twostext(p) * tmp9
                     dh6 = (v * du - u * dv) / (v * v)
  
                     da1 = h1/sigma * s2*s2 + h2 * s2*s1 + h3 * s2/s1 &
                         + (1._r8 - s2*s1) / (twostext(p) + h) * dh2 &
                         + (1._r8 - s2/s1) / (twostext(p) - h) * dh3
                     da2 = h4/sigma * s2*s2 + h5 * s2*s1 + h6 * s2/s1 &
                         + (1._r8 - s2*s1) / (twostext(p) + h) * dh5 &
                         + (1._r8 - s2/s1) / (twostext(p) - h) * dh6
  
                     ! Flux derivatives
     
                     d_ftid = -twostext(p)*h4/sigma*s2 - h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1
                     d_fabd = -(dh2+dh3) + (1._r8-albgrd(c,ib))*twostext(p)*s2 - (1._r8-albgri(c,ib))*d_ftid
                     d_fabd_sun = (1._r8 - omega(p,ib)) * (twostext(p)*s2 + 1._r8 / avmu(p) * (da1 + da2))
                     d_fabd_sha = d_fabd - d_fabd_sun
  
                     fabd_sun_z(p,iv) = max(d_fabd_sun, 0._r8)
                     fabd_sha_z(p,iv) = max(d_fabd_sha, 0._r8)
  
                     ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                     ! to normalize derivatives by sunlit or shaded fraction to get
                     ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
  
                     fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / fsun_z(p,iv)
                     fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))
  
                     ! ===============
                     ! Diffuse
                     ! ===============
  
                     ! Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai
  
                     u1 = b - c1/albgri(c,ib)
                     u2 = b - c1*albgri(c,ib)

                     a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
                     a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)
     
                     ! Derivatives for h7, h8, h9, h10 and a1, a2
  
                     v = d1
                     dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1
     
                     u = c1 * tmp2 / s1
                     du = h * c1 * tmp2 / s1
                     dh7 = (v * du - u * dv) / (v * v)
  
                     u = -c1 * tmp3 * s1
                     du = h * c1 * tmp3 * s1
                     dh8 = (v * du - u * dv) / (v * v)
  
                     v = d2
                     dv = h * tmp4 / s1 + h * tmp5 * s1
  
                     u = tmp4 / s1
                     du = h * tmp4 / s1
                     dh9 = (v * du - u * dv) / (v * v)
  
                     u = -tmp5 * s1
                     du = h * tmp5 * s1
                     dh10 = (v * du - u * dv) / (v * v)
  
                     da1 = h7*s2*s1 +  h8*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh7 + (1._r8-s2/s1)/(twostext(p)-h)*dh8
                     da2 = h9*s2*s1 + h10*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh9 + (1._r8-s2/s1)/(twostext(p)-h)*dh10
  
                     ! Flux derivatives
  
                     d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1
                     d_fabi = -(dh7+dh8) - (1._r8-albgri(c,ib))*d_ftii
                     d_fabi_sun = (1._r8 - omega(p,ib)) / avmu(p) * (da1 + da2)
                     d_fabi_sha = d_fabi - d_fabi_sun
  
                     fabi_sun_z(p,iv) = max(d_fabi_sun, 0._r8)
                     fabi_sha_z(p,iv) = max(d_fabi_sha, 0._r8)
  
                     ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                     ! to normalize derivatives by sunlit or shaded fraction to get
                     ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
  
                     fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / fsun_z(p,iv)
                     fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))
  
                  end do ! end of iv loop
               end if ! nlevcan
            end if   ! first band
          end if  ! NOT lSFonly

       end do   ! end of pft loop
    end do   ! end of radiation band loop

     end associate 
  end subroutine TwoStream

   !-----------------------------------------------------------------------
   subroutine Multilayer (bounds, filter_vegsol, num_vegsol, coszen, rho, tau, &
              atm2lnd_inst, temperature_inst, waterdiagnosticbulk_inst, surfalb_inst)
     !
     ! !DESCRIPTION: (added by Yuanchao Fan 01.12.2014)
     ! Multilayer canopy radiative transfer based on the scheme of Norman.
     ! Norman, J.M. 1979. Modeling the complete crop canopy.
     ! in Modification of the Aerial Environment of Crops.
     ! B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
     ! Translated and modified from CANOAK model from C to Fortran (DENNIS BALDOCCHI 2008)

     ! This subroutine is developed by Y.Fan and is unpublished. Please cite:
     ! Fan, Y. Modeling oil palm monoculture and its associated impacts on
     ! land-atmosphere carbon, water and energy fluxes in Indonesia. PhD Thesis.(University of
     ! Göttingen, 2016). 
     ! Please contact Y.Fan before using this code (yfansunny@gmail.com)
   

     !
     ! !USES:
     !use clm_atmlnd, only : clm_a2l
     use atm2lndType, only : atm2lnd_type
     use clm_varpar, only : numrad, nlevcan
     use clm_varcon, only : omegas, tfrz, betads, betais
     use pftconMod , only : pftcon !use pointer for clumping,phytomer !Y.Fa n
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds           ! bounds
     integer , intent(in)  :: filter_vegsol(:)         ! filter for vegetated pfts with coszen>0
     integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0
     real(r8), intent(in)  :: coszen( bounds%begp: )   ! cosine solar zenith angle for next time step [pft]
     real(r8), intent(in)  :: rho( bounds%begp: , 1: ) ! leaf/stem refl weighted by fraction LAI and SAI [pft, numrad]
     real(r8), intent(in)  :: tau( bounds%begp: , 1: ) ! leaf/stem tran weighted by fraction LAI and SAI [pft, numrad]
     type(temperature_type) , intent(in)    :: temperature_inst
     type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
     !type(waterstate_type)  , intent(in)    :: waterstate_inst
     type(surfalb_type)     , intent(inout) :: surfalb_inst
     type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
     !
     ! !LOCAL VARIABLES:
     REAL, PARAMETER :: PI = 3.1415926_r8
     REAL, PARAMETER :: PI180 = 0.017453292_r8   ! PI divided by 180, radians per degree
     real(r8), parameter :: mpe = 1.e-06_r8 ! prevents overflow for division by zero
     integer  :: g,fp,p,c,iv,iz,IREP        ! array indices
     integer  :: ib                 ! waveband number
     real(r8) :: cosz               ! 0.001 <= coszen <= 1.000
	 real(r8) :: rho1, tau1         ! reflectance and transmittance of vegetation elements
     real(r8) :: zen, d_zen, d_lai, laisum
     real(r8) :: fraction_beam, DOWN, UP, fabd_nir_z, fabi_nir_z
    ! real(r8) :: beam(bounds%begp:bounds%endp,1:numrad,1:(nlevcan+1))	 !
    ! real(r8) :: down_dif(bounds%begp:bounds%endp,1:numrad,1:(nlevcan+1))	 !
    ! real(r8) :: up_dif(bounds%begp:bounds%endp,1:numrad,1:(nlevcan+1))	 !
     real(r8) :: UPDOWN(bounds%begp:bounds%endp,1:numrad,1:(nlevcan+1))	 !
     real(r8) :: ext_dif(bounds%begp:bounds%endp,1:nlevcan)	 	 !
     real(r8) :: ext_beam(bounds%begp:bounds%endp,1:nlevcan)	 !
     real(r8) :: PEN_dif(bounds%begp:bounds%endp,1:nlevcan)	 	 !
     real(r8) :: PEN_beam(bounds%begp:bounds%endp,0:nlevcan)	 !
     real(r8) :: refl_dif(bounds%begp:bounds%endp,1:numrad,1:nlevcan)	 !
     real(r8) :: refl_beam(bounds%begp:bounds%endp,1:numrad,1:nlevcan)	 !
     real(r8) :: tran_dif(bounds%begp:bounds%endp,1:numrad,1:nlevcan) 	 !
     real(r8) :: tran_beam(bounds%begp:bounds%endp,1:numrad,1:nlevcan)	 !
     real(r8) :: abso_dif(bounds%begp:bounds%endp,1:numrad,1:nlevcan) 	 !
     real(r8) :: abso_beam(bounds%begp:bounds%endp,1:numrad,1:nlevcan)	 !
     !real(r8) :: gfunc_solar(1:nlevcan)
     !real(r8) :: gfunc_sky(1:nlevcan,1:9)
    !!the following will cause error "forrtl: severe (151): allocatable array is already allocated"
    !!should allocate and initiate in initSurfAlbMod.F90
    !   allocate(gfunc_solar(1:nlevcan))
    !   allocate(gfunc_sky(1:nlevcan,1:9))
    !   allocate(lad(1:9))
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endp/)),         errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(rho)    == (/bounds%endp, numrad/)), errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(tau)    == (/bounds%endp, numrad/)), errMsg(__FILE__, __LINE__))

   associate(&
          forc_solad   =>    atm2lnd_inst%forc_solad_grc         , & ! Input:  [real(r8) (:,:)]  direct beam radiation (W/m**2)
          forc_solai   =>    atm2lnd_inst%forc_solai_grc         , & ! Input:  [real(r8) (:,:)]  diffuse radiation (W/m**2)

          albgrd       =>    surfalb_inst%albgrd_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (direct) (column-level)
          albgri       =>    surfalb_inst%albgri_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (diffuse)(column-level)
          t_veg        =>    temperature_inst%t_veg_patch        , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)
          fwet         =>    waterdiagnosticbulk_inst%fwet_patch          , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)
          albd         =>    surfalb_inst%albd_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (direct)
          albi         =>    surfalb_inst%albi_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (diffuse)
          fabd         =>    surfalb_inst%fabd_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
          fabd_sun     =>    surfalb_inst%fabd_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
          fabd_sha     =>    surfalb_inst%fabd_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
          fabi         =>    surfalb_inst%fabi_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
          fabi_sun     =>    surfalb_inst%fabi_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha     =>    surfalb_inst%fabi_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
          ftdd         =>    surfalb_inst%ftdd_patch             , & ! Output: [real(r8) (:,:) ]  down direct flux below canopy per unit direct flx
          ftid         =>    surfalb_inst%ftid_patch             , & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flx
          ftii         =>    surfalb_inst%ftii_patch             , & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flx
          n_iter       =>    surfalb_inst%n_iter_patch           , & ! Output: [real(r8) (:,:)]  number of iterations for solving diffuse flux
          balerr       =>    surfalb_inst%balerr_patch           , & ! Output: [real(r8) (:,:)]  remaining energy balance error after iteration
          tlai_z       =>    surfalb_inst%tlai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tlai increment for canopy layer
          tsai_z       =>    surfalb_inst%tsai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tsai increment for canopy layer
          nrad         =>    surfalb_inst%nrad_patch             , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          fabd_sun_z   =>    surfalb_inst%fabd_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z   =>    surfalb_inst%fabd_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z   =>    surfalb_inst%fabi_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z   =>    surfalb_inst%fabi_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          beam         =>    surfalb_inst%beam_patch             , & ! Output: [real(r8) (:,:,:)] Direct beam flux above each layer (incident on the horizontal)
          up_dif       =>    surfalb_inst%up_dif_patch           , & ! Output: [real(r8) (:,:,:)] Upward scattered flux above each layer (per unit total radiation)
          down_dif     =>    surfalb_inst%down_dif_patch         , & ! Output: [real(r8) (:,:,:)] Downward scattered flux above each layer (per unit total radiation)
          gfunc_sky    =>    surfalb_inst%gfunc_sky_patch        , & ! Input:  [real(r8) (:,:,:)] the G function for diffuse sky radiation
          gfunc_solar  =>    surfalb_inst%gfunc_solar_patch      , & ! Input:  [real(r8) (:,:)] the G function for direct beam radiation
          fsun_z       =>    surfalb_inst%fsun_z_patch           , & ! Output: [real(r8) (:,:) ]  sunlit fraction of canopy layer
          phytomer     =>    pftcon%phytomer                     , & ! Input:  [integer (:)] total number of phytomers in life time (added by Y.Fan)
          clumping     =>    pftcon%clumping                     , & ! Input: canopy clumping index for multilayer radiative transfer (Y.Fan) 
          ivt          =>    patch%itype                          & ! Input:  [integer (:)]  pft vegetation type
   )

    ! Calculate extinction coefficients that are independent of waveband:
    do fp = 1,num_vegsol
       p = filter_vegsol(fp)

       ! note that the following limit only acts on cosz values > 0 and less than 0.001,
       ! not on values cosz = 0, since these zero have already been filtered out in filter_vegsol
       cosz = max(0.001_r8, coszen(p))
       !cosz = max(mpe, coszen(p))

       laisum = 0._r8
       PEN_beam(p,0) = 1._r8 !initial value above canopy
       do iv = 1, nrad(p) !from canopy top to bottom

          ! Cumulative lai+sai
          d_lai  = tlai_z(p,iv)+tsai_z(p,iv)
          laisum = laisum + d_lai

        !!Diffuse radiation
          !Calculate probability of penetration of diffuse radiation for an isotropic sky
          PEN_dif(p,iv) = 0._r8
          do iz=1, 9 !9 zenith angles: 5 to 85 degree
            zen = ((iz - 1) * 10._r8 + 5._r8) * PI180
            d_zen = 10._r8 * PI180 !delta radian

            !probability of diffuse photon transfer through a single layer within a single sky sector
            ext_dif(p,iv) = exp(-d_lai*clumping(ivt(p))*gfunc_sky(p,iv,iz)/cos(zen))
            !Integrated probability of diffuse sky radiation penetration through the hemisphere for each layer
            PEN_dif(p,iv) = PEN_dif(p,iv) + 2._r8*ext_dif(p,iv)*sin(zen)*cos(zen)*d_zen
            PEN_dif(p,iv) = min(1._r8, max(0._r8, PEN_dif(p,iv)))
            !negative PEN_dif values may occur at bottom layers with small d_lai when sun elevation is low (2015.11)
          end do

        !!Direct beam
          !probability of direct beam penetration through a layer
          ext_beam(p,iv) = exp(-d_lai*clumping(ivt(p))*gfunc_solar(p,iv)/cosz)
          !Integrated probability of beam penetration (unintercepted) from top layer to iv
          !PEN_beam(p,iv) = exp(-laisum*clumping(ivt(p))*gfunc_solar(p,iv)/cosz) !Canoak method: wrong
          !gfunc_solar is only per layer, have to integrate ext_beam from canoopy top to the current layer
          PEN_beam(p,iv) = PEN_beam(p,iv-1)*ext_beam(p,iv) !PEN_beam decreases with increasing depth
       end do
    end do

   ! Loop over all wavebands to calculate for each canopy layer the scattered and absorbed fluxes
   ! using the Norman Cupid multilayer radiative transfer model
   !
   ! Output:
   ! ------------------
   ! Direct beam fluxes
   ! ------------------
   ! beam       - Direct beam flux above each layer (incident on the horizontal)
   ! refl_beam  - Direct beam reflectance of a canopy layer
   ! tran_beam  - Direct beam transmittance of a canopy layer

   ! ------------------
   ! Diffuse fluxes
   ! ------------------
   ! up_dif     - Upward scattered flux above each layer (per unit total radiation)
   ! down_dif   - Downward scattered flux above each layer (per unit total radiation)
   ! refl_dif   - Diffuse radiation reflectance of a canopy layer
   ! tran_dif   - Diffuse radiation transmittance of a canopy layer

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          c = patch%column(p)
          g = patch%gridcell(p)

        cosz = max(0.001_r8, coszen(p))
       !cosz = max(mpe, coszen(p))

      !Night time flux (including sunset when sun elevation < 5 degree)
!    if (cosz < 0.1) then
       !do nothing
!    else !calculate radiative fluxes only at daytime

       !avoid division by zero at sunset/evening time
       !sometimes there are negative radiation values in input forcing file when the sun is at or below horizon
!       if (forc_solad(g,ib) <= 0._r8) then
!          forc_solad(g,ib) = 0.001_r8
!       end if
       fraction_beam = forc_solad(g,ib) / max(mpe,(forc_solad(g,ib)+forc_solai(g,ib)))
       fraction_beam = min(0.99_r8, max(0.01_r8, fraction_beam)) !min threshold 1% beam, otherwise NetCDF I/O error


       !Adjust direct/diffuse reflectance and transmittance for intercepted snow on the canopy
       !snow adjustment not needed for tropical PFTs
       !currently use snow omegas and betads(=betais=0.5) to estimate snow reflectance and transmittance
       !in the future parametrize rho_s and tau_s values for snow

       if (t_veg(p) > tfrz) then                             !no snow
          rho1 = rho(p,ib)
          tau1 = tau(p,ib)
       else !transmittance has to go through both snow layer and leaf layer
          rho1 =  (1._r8-fwet(p))*rho(p,ib) + fwet(p)*omegas(ib)*betads
          tau1 =  (1._r8-fwet(p))*tau(p,ib) + fwet(p)*omegas(ib)*(1._r8-betads)*tau(p,ib)
       end if


    if (nrad(p) >= 1) then !when there is canopy

    ! initiate scattering using the technique of NORMAN (1979).
    ! scattering is computed using an iterative technique.

       !top boundary: downward diffuse radiation from above canopy
       down_dif(p,ib,1) = 1._r8 - fraction_beam  !per unit incoming total radiation
!       down_dif(p,ib,1) = max(mpe,forc_solai(g,ib))

       do iv = 1, nrad(p) !from canopy top (iv=1) to bottom (iv=nrad)

          !Direct beam flux above each layer (incident on the horizontal)
          beam(p,ib,1) = fraction_beam !incident beam, top of canopy (per unit incoming total)
!          beam(p,ib,1) = max(mpe,forc_solad(g,ib))
          beam(p,ib,iv+1) = beam(p,ib,iv) * ext_beam(p,iv)
        !or!
        !  beam(p,iv+1) = PEN_beam(p,iv)*beam(p,1)

          !!Contribute of direct beam to scattered diffuse fluxes
          !Direct beam reflectance of a canopy layer
          !intercepted and reflected still within the same layer iv
          refl_beam(p,ib,iv) = (1._r8 - ext_beam(p,iv)) * rho1
          !beam reflectance of soil/snow ground: refl_beam(nrad(p)+1) = albgrd(c,ib)

          !Direct beam transmittance of a canopy layer
          !intercepted by layer iv but transmitted downward to lower layer iv+1
          tran_beam(p,ib,iv) = (1._r8 - ext_beam(p,iv)) * tau1
          !no transmittance below soil/snow ground: tran_beam(nrad(p)+1) = 0._r8

          !Direct beam absorptance of a canopy layer (the portion of intercepted beam after refl and trans)
          abso_beam(p,ib,iv) = (1._r8 - ext_beam(p,iv)) *(1._r8 - tau1 - rho1)
                               !* & gfunc_solar(p,iv)/cosz
          !PSUN (for absorption) should be the radiation incident on the mean leaf normal (or leaf area projected in sun direction)? (no)
          !beam absorptance by soil/snow ground: abso_beam(nrad(p)+1) = (1._r8-albgrd(c,ib))
        
          !Diffuse radiation reflectance of a canopy layer, related to rho of vegetation elements
          !intercepted and reflected by layer iv
          refl_dif(p,ib,iv) = (1._r8 - PEN_dif(p,iv)) * rho1
          !Diffuse radiation transmittance of a canopy layer, related to tau of vegetation elements
          !intercepted and transmitted by layer iv, plus the portion unintercepted
          tran_dif(p,ib,iv) = (1._r8 - PEN_dif(p,iv)) * tau1 + PEN_dif(p,iv)
          !Diffuse radiation absorptance of a canopy layer
          !intercepted and absorbed by layer iv
          abso_dif(p,ib,iv) = (1._r8 - PEN_dif(p,iv)) * (1._r8 - tau1 - rho1)
          !the above three terms refl_dif, tran_dif, abso_dif sum to 1

        !Following Eq. 1-5 for display and derivation only
        !!Diffuse: upward and downward at each layer
        !  up_dif(p,ib,iv) = down_dif(p,ib,iv)*refl_dif(p,ib,iv) + up_dif(p,ib,iv+1)*tran_dif(p,ib,iv) + beam(p,ib,iv)*refl_beam(p,ib,iv)        !Eq.1

        !  down_dif(p,ib,iv+1) = up_dif(p,ib,iv+1)*refl_dif(p,ib,iv) + down_dif(p,ib,iv)*tran_dif(p,ib,iv) + beam(p,ib,iv)*tran_beam(p,ib,iv)    !Eq.2

        !!solve above Eq.1 & 2 for upward and downward diffuse, assuming beam(p,ib,iv)=0
        !Define UPDOWN(p,ib,iv)=up_dif(p,ib,iv)/down_dif(p,ib,iv), ratio of up/down diffuse radiation for layer iv
        !divide Eq.1 by down_dif(p,ib,iv) and Eq.2 by down_dif(p,ib,iv+1), multiply the two equations, one gets:
        !  UPDOWN(p,ib,iv) = ((tran_dif(p,ib,iv)*tran_dif(p,ib,iv) - refl_dif(p,ib,iv)*refl_dif(p,ib,iv))*UPDOWN(p,ib,iv+1) + &               !Eq.3
        !               refl_dif(p,ib,iv)) / (1._r8 - refl_dif(p,ib,iv)*UPDOWN(p,ib,iv+1))

        !!one also gets the ratio of downward fluxes in successive layers from Eq.2 (without beam):
        !  down_dif(p,ib,iv+1) = down_dif(p,ib,iv)*tran_dif(p,ib,iv) / (1._r8 - UPDOWN(p,ib,iv+1)*refl_dif(p,ib,iv))		                  !Eq.4

       end do

       !Starting from soil layer, calculate UPDOWN successively using Eq.3
       UPDOWN(p,ib,nrad(p)+1) = albgri(c,ib) !diffuse reflectance of soil/snow ground
       do iv = nrad(p), 1, -1 !from canopy bottom to  top
          UPDOWN(p,ib,iv) = UPDOWN(p,ib,iv+1) * tran_dif(p,ib,iv)*tran_dif(p,ib,iv) / &   !same as Eq.3
                            max(mpe,(1._r8 - UPDOWN(p,ib,iv+1)*refl_dif(p,ib,iv))) + refl_dif(p,ib,iv)
       end do

       !!Then solve downward fluxes with Eq.4 (when assuming no beam)
       !start with downward flux from the sky: down_dif(p,ib,1)
       do iv = 1, nrad(p)
          down_dif(p,ib,iv+1) = down_dif(p,ib,iv)*tran_dif(p,ib,iv)/max(mpe,(1._r8 - UPDOWN(p,ib,iv+1)*refl_dif(p,ib,iv)))      !based on Eq.4
                    ! + beam(p,ib,iv)*tran_beam(p,ib,iv)                  !based on Eq.2 and Eq.4
        !Then solve upward fluxes
          up_dif(p,ib,iv) = UPDOWN(p,ib,iv)*down_dif(p,ib,iv)     !based on Eq.3
                    ! + beam(p,ib,iv)*refl_beam(p,ib,iv)     !based on Eq.1 and Eq.3
       end do

       !lower boundary: upward radiation from soil/snow ground
       up_dif(p,ib,nrad(p)+1) = albgri(c,ib)*down_dif(p,ib,nrad(p)+1)
                        !  + albgrd(c,ib)*beam(p,ib,nrad(p)+1)
       !soil albedo is set the same for direct and diffuse (see line 946)

    !!--------------------------------------------------------------------------
    !Now one has a first approximation for all the up and down fluxes using the assumed UPDOWN ratio
    !substitute these fluxes to the right sides of Eq.1 and 2 to estimate left sides iteratively when direct beam is on
    !finally get better estimates of up/down fluxes when the iteration converges
    
    !!Iterative calculation of upward diffuse and downward beam + diffuse
       n_iter(p,ib) = 0._r8
       IREP =1
       do while (IREP==1 .and. n_iter(p,ib)<=100._r8) !maximum iteration 100
          IREP = 0
          n_iter(p,ib) = n_iter(p,ib) + 1._r8

          !loop for downward diffuse
          do iv = 1, nrad(p) !from canopy top to bottom
             DOWN = up_dif(p,ib,iv+1)*refl_dif(p,ib,iv) + down_dif(p,ib,iv)*tran_dif(p,ib,iv) + beam(p,ib,iv)*tran_beam(p,ib,iv)
             !even if only one layer does not converge, iteration continue
             if (abs(DOWN - down_dif(p,ib,iv+1)) > 0.001_r8) IREP = 1
             down_dif(p,ib,iv+1) = DOWN
             !down_dif(p,ib,iv+1) = min(1._r8, max(0._r8, DOWN)) !no negative values occur here
          end do
          !since down_dif all adjusted, up_dif at soil needs to adjust again
          up_dif(p,ib,nrad(p)+1) = albgri(c,ib)*down_dif(p,ib,nrad(p)+1) + albgrd(c,ib)*beam(p,ib,nrad(p)+1)

          !loop for upward diffuse
          do iv = nrad(p), 1, -1 !from canopy bottom to top
             UP = down_dif(p,ib,iv)*refl_dif(p,ib,iv) + up_dif(p,ib,iv+1)*tran_dif(p,ib,iv) + beam(p,ib,iv)*refl_beam(p,ib,iv)
             !even if only one layer does not converge, iteration continue
             if (abs(UP - up_dif(p,ib,iv)) > 0.001_r8) IREP = 1
             up_dif(p,ib,iv) = UP
             !up_dif(p,ib,iv) = min(1._r8, max(0._r8, UP)) !no negative values occur here
          end do
       end do

       !since up_dif all adjusted, down_dif may need to be adjusted again
       !there might be remaining energy balance error in the above iteration
       !Check the error threshold later (line 1942-1952)
       !and Repeat the above two loops when necessary

       !to do: report iter in history file to check iteration times


    !!--------------------------------------------------------------
      !derive absorbed fluxes by each canopy layer and
      !then derive canopy level fabd/fabi, ftdd/ftid/ftii, albd/albi

      if (ib == 1) then !visible band
        !only visible band is used to derive APAR per layer to drive photosynthesis and energy exchanges
        fabd_sun(p,ib) = 0._r8
        fabd_sha(p,ib) = 0._r8
        fabi_sun(p,ib) = 0._r8
        fabi_sha(p,ib) = 0._r8
        do iv = 1, nrad(p)
        !Sunlit fraction of each canopy layer
        !Note: PEN_beam is similar to 's2' in two-stream but calculation of laisum at a layer is different
        !for tree canopies, leaf elements are assumed a turbid medium, laisum is cumulative lai+sai at center of layer
        !but for oil palm, a frond is assumed a thin layer, so laisum is for each whole layer
        !  if (phytomer(ivt(p)) > 0) then !for oil palm, assume a very thin layer
        !     d_lai  = tlai_z(p,iv)+tsai_z(p,iv)
        !     laisum = laisum + d_lai
        !    ! fsun_z(p,1) = 1._r8 !the top layer is fully sunlit (not really! because leaves have self-shading/umbral)
        !     fsun_z(p,iv) = max(0._r8, min(1._r8, clumping(ivt(p))*PEN_beam(p,iv)))
        !  else !for other tree PFTs, assume turbid medium
             ! Cumulative lai+sai at center of layer
        !     if (iv == 1) then
        !        d_lai  = 0.5_r8 * (tlai_z(p,iv)+tsai_z(p,iv))
        !        laisum = d_lai
        !     else
        !	    d_lai  = 0.5_r8 * ((tlai_z(p,iv-1)+tsai_z(p,iv-1))+(tlai_z(p,iv)+tsai_z(p,iv)))
        !        laisum = laisum + d_lai
        !     end if
        !     fsun_z(p,iv) = max(0._r8, min(1._r8, clumping(ivt(p))*PEN_beam(p,iv)))
        !  end if

        !  fsun_z(p,iv) = min(1._r8-mpe, max(mpe, PEN_beam(p,iv)))  !most models use this equation
!          fsun_z(p,iv) = min(1._r8-mpe, max(mpe, clumping(ivt(p))*PEN_beam(p,iv))) !Canoak method, (Leuning et al. Spitters et al.)
         !above fsun_z is the integrated probability of sun beam penetration at a layer, which is not umbral
          !Only if gfunc_solar is integrated value from canopy top, then divide the intercepted beam ratio (1-PEN_beam)
          !by the mean projection (G/cosz) will yield the cumulative sunlit leaf area from top to layer iv:
          !laisum_sun_z(p,iv) = cosz*(1._r8 - PEN_beam(p,iv))/(clumping(ivt(p))*gfunc_solar(p,iv))
          !since gfunc_solar is only per layer, the following equation is used to calculate fsun_z
        !fsun_z calculated by the following equation (Wilson 1967) has a different meaning: the fraction of sunlit leaf area
        !this should be used to evaluate fluxes for each layer if the revised diffuse method is used (Y.Fan 2015.07)
          d_lai  = tlai_z(p,iv)+tsai_z(p,iv)
          fsun_z(p,iv) = cosz*(PEN_beam(p,iv-1)- PEN_beam(p,iv))/(clumping(ivt(p))*gfunc_solar(p,iv)) / d_lai
          fsun_z(p,iv) = min(1._r8-mpe, max(mpe, fsun_z(p,iv)))

         !!*****Direct*****
           !absorbed direct beam flux by sunlit leaves per layer (incident along the mean leaf normal)
!	   if (cosz > 0.1_r8) then !avoid division by small cosz to give very high APAR at sunset
!             fabd_sun_z(p,iv) = (beam(p,ib,1)*gfunc_solar(p,iv)/cosz) *(1._r8 - tau1 - rho1) !Canoak method: element absorptance, per unit sunlit LAI
!          !above already normalized per unit sunlit LAI in the sun direction (expected unit by other modules)
!          !(beam(p,ib,1) already per unit sunlit area (W/m2), later fsun_z (=PEN_beam) will account for intercepted beam (by sunlit leaf area)
!          else !only keep diffuse light when sun is at or below horizon: sun elevation < 5 degree
!             fabd_sun_z(p,iv) = 0._r8
!          end if
          fabd_sun_z(p,iv) = beam(p,ib,iv) * abso_beam(p,ib,iv) !canopy absorption for the current sunlit LAI (need to normalize later)
          fabd_sha_z(p,iv) = 0._r8  !shaded leaves receive NO direct beam

          !need to normalize per unit sunlit(LAI+SAI) and per unit shaded(LAI+SAI)
          !for canopy layers (to be compatible with other modules)
          fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / (fsun_z(p,iv)*d_lai)
          fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / ((1._r8 - fsun_z(p,iv))*d_lai)


         !!*****Diffuse*****
           !absorbed diffuse radiation (received on top and bottom of a layer) by sunlit and shaded leaves
!          fabi_sun_z(p,iv) = (down_dif(p,ib,iv) + up_dif(p,ib,iv+1)) *(1._r8 - tau1 - rho1)
!		  !Canoak method use vegetation element absorptance, already per unit sunlit LAI (but the absorption may be wrong)
!          fabi_sha_z(p,iv) = (down_dif(p,ib,iv) + up_dif(p,ib,iv+1)) *(1._r8 - tau1 - rho1)
          !revised diffuse radiation
          fabi_sun_z(p,iv) = fsun_z(p,iv) * (down_dif(p,ib,iv) + up_dif(p,ib,iv+1))* abso_dif(p,ib,iv)
          fabi_sha_z(p,iv) = (1._r8 - fsun_z(p,iv)) * (down_dif(p,ib,iv) + up_dif(p,ib,iv+1))* abso_dif(p,ib,iv)
          !normalize per unit sunlit(LAI+SAI) and per unit shaded(LAI+SAI)
          fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / (fsun_z(p,iv)*d_lai)
          fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / ((1._r8 - fsun_z(p,iv))*d_lai)


          !Above are absolute APARsun and APARsha per unit incoming total solar radiation
          !need to normalize per unit direct beam flux and per unit diffuse flux
          fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / beam(p,ib,1) !fraction_beam
          fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / beam(p,ib,1) !fraction_beam
          fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / down_dif(p,ib,1) !(1._r8 - fraction_beam)
          fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / down_dif(p,ib,1) !(1._r8 - fraction_beam)


          !absorbed by sunlit/shaded canopy: from unit LAI per layer to canopy
          fabd_sun(p,ib) = fabd_sun(p,ib) + fabd_sun_z(p,iv)*fsun_z(p,iv)*d_lai
          fabd_sha(p,ib) = 0._r8 !only sunlit leaf absorb direct beam
          fabi_sun(p,ib) = fabi_sun(p,ib) + fabi_sun_z(p,iv)*fsun_z(p,iv)*d_lai
          fabi_sha(p,ib) = fabi_sha(p,ib) + fabi_sha_z(p,iv)*(1._r8 - fsun_z(p,iv))*d_lai


        end do

        !Direct beam flux absorbed by whole canopy (per unit direct beam flux)
        fabd(p,ib) = fabd_sun(p,ib) + fabd_sha(p,ib)

        !Diffuse flux absorbed by whole canopy (per unit diffuse flux)
        fabi(p,ib) = fabi_sun(p,ib) + fabi_sha(p,ib)

        !Average sunlit/shaded fraction for the whole canopy
        !or fsun_z(p,1) when nlevcan =1
        !  fsun(p) = cosz*(1._r8 - PEN_beam(p,iv))/(laisum*clumping(ivt(p))*gfunc_solar(p,iv))

      else !NIR bands: drive energy exchanges
        !no need to differentiate sun/shade for NIR
        !fabd_sun/fabd_sha only used for visible bands
        fabd(p,ib) = 0._r8
        fabi(p,ib) = 0._r8
        do iv = 1, nrad(p)
!          if (cosz > 0.1_r8) then
!             fabd_nir_z = (beam(p,ib,1)*gfunc_solar(p,iv)/cosz) *(1._r8 - tau1 - rho1)  !Canoak method
!          else
!             fabd_nir_z = 0._r8
!          end if
!          fabi_nir_z = (down_dif(p,ib,iv) + up_dif(p,ib,iv+1)) *(1._r8 - tau1 - rho1) !Canoak may be wrong: diffuse has different penetration function
          fabd_nir_z = beam(p,ib,iv) * abso_beam(p,ib,iv)
          fabi_nir_z = (down_dif(p,ib,iv) + up_dif(p,ib,iv+1))* abso_dif(p,ib,iv) !no need to normalize by d_lai, only used to calculate fabi total
          !absorbed by whole canopy
          fabd(p,ib) = fabd(p,ib) + fabd_nir_z
          fabi(p,ib) = fabi(p,ib) + fabi_nir_z                     !layer to canopy (revised method)
!          fabd(p,ib) = fabd(p,ib) + fabd_nir_z *fsun_z(p,iv)*d_lai !leaf to canopy (only sunlit leaf absorb direct beam)
!          fabi(p,ib) = fabi(p,ib) + fabi_nir_z *d_lai			    !Canoak method
        end do
        !normalize to get per unit direct beam flux and per unit diffuse flux
        fabd(p,ib) = fabd(p,ib) / beam(p,ib,1) !fraction_beam
        fabi(p,ib) = fabi(p,ib) / down_dif(p,ib,1) !(1._r8 - fraction_beam)

      end if


    !!--------------------------------------------------------------
    !derive canopy level albd/albi, ftdd/ftid/ftii to fit with other modules in the CLM structure
    !cannot follow the two-stream equations because a portion of direct beam is converted to diffuse

      !Direct
       !Unintercepted (so called Transmitted in two-stream) direct beam flux below canopy (per unit direct beam flux)
       ftdd(p,ib) = beam(p,ib,nrad(p)+1)/beam(p,ib,1) !=PEN_beam(p,nrad(p)), already normalized

       !Downward scattered beam flux below canopy (per unit direct beam flux)
       !beam(nrad(p))*tran_beam(nrad(p))/beam(1) is only from the bottom canopy layer
       !scattered beam flux from all layers are merged to diffuse flux up_dif/down_dif in this multilayer code
       !so set ftid to zero. this portion will be included in ftii ultimately
       ftid(p,ib) = 0._r8

       !Surface albedo: upward scattered beam flux above canopy (per unit direct beam flux)
       !albd(p,ib) = 0._r8  !reflected beam already converted to diffuse and included in up_dif(p,iv)
       albd(p,ib) = 1._r8 - fabd(p,ib) - (1._r8-albgrd(c,ib))*ftdd(p,ib) !last item is the fraction absorbed by ground
       !ensure energy balance and avoid error

     !Diffuse
       !Downward scattered diffuse flux below canopy (per unit diffuse flux)
       ftii(p,ib) = down_dif(p,ib,nrad(p)+1)/down_dif(p,ib,1)

       !Surface albedo: upward scattered diffuse flux above canopy (per unit diffuse flux)
       !albi(p,ib) = up_dif(p,ib,1)/down_dif(p,ib,1) !!here up_dif includes reflected beam; this option has minor balerr!
       albi(p,ib) = 1._r8 - fabi(p,ib) - (1._r8-albgri(c,ib))*ftii(p,ib)

       !There might be remaining error in diffuse energy balance after certain iterations of solving diffuse fluxes
       !calculate the error as (for diagnosis)
       balerr(p,ib) = up_dif(p,ib,1) + fabd(p,ib)*beam(p,ib,1) + fabi(p,ib)*down_dif(p,ib,1) + &
              down_dif(p,ib,nrad(p)+1)*(1._r8-albgri(c,ib)) + beam(p,ib,nrad(p)+1)*(1._r8-albgrd(c,ib)) - &
              (forc_solad(g,ib)+forc_solai(g,ib))
    !   balerr(p,ib) =	balerr(p,ib) / (1._r8 - fraction_beam)	!proportional error
    !If balerr > 1%, repeat the iteration process in lines 1755-1781
    !   if (balerr > 0.01_r8) go to 100

       !Add the remaining error (<1%) to albi for energy balance
    !   albi(p,ib) =	albi(p,ib) + balerr(p,ib)

       !absolute error (the BalanceCheckMod calculate ERRSOL as follows)
       !balerr(p,ib) = forc_solad(g,ib)*fabd(p,ib) + forc_solai(g,ib)*fabi(p,ib) + &
       !               forc_solad(g,ib)*ftdd(p,ib)*(1._r8-albgrd(c,ib)) + &
       !               forc_solai(g,ib)*ftii(p,ib)*(1._r8-albgri(c,ib)) + &
       !               albd(p,ib)*forc_solad(g,ib) + albi(p,ib)*forc_solai(g,ib) - &
       !              (forc_solad(g,ib)+forc_solai(g,ib))

    else !bare ground or leaf not out yet
       fabd(p,ib) = 0._r8
       ftdd(p,ib) = 1._r8
       ftid(p,ib) = 0._r8
       albd(p,ib) = albgrd(c,ib)

       fabi(p,ib) = 0._r8
       ftii(p,ib) = 1._r8
       albi(p,ib) = albgri(c,ib)

    end if

 !          end if !end evening/day time
       end do   ! end of pft loop
    end do   ! end of radiation band loop (VIS/NIR)


    end associate
   end subroutine Multilayer

   !-----------------------------------------------------------------------
   subroutine Gfunc_all (bounds, filter_vegsol, num_vegsol, coszen, surfalb_inst)
     !
     ! !DESCRIPTION: (added by Y.Fan 01.12.2014)
     !!!This subroutine computes the G Function according to Wang et al. 2007
     !Comparison of leaf angle distribution functions: effects on extinction coefficient and fraction of sunlit foliage.
     !Agricultural and Forest Meteorology, 143(1), 106-122.

     ! This subroutine is developed by Y.Fan and is unpublished. Please cite:
     ! Fan, Y. Modeling oil palm monoculture and its associated impacts on
     ! land-atmosphere carbon, water and energy fluxes in Indonesia. PhD
     ! Thesis.(University of
     ! Göttingen, 2016).
     ! Please contact Y.Fan before using this code (yfansunny@gmail.com)


     ! !USES:
     use clm_varpar, only : nlevcan, fLAD
     !use pftvarcon , only: laimx, mxlivenp, theta, phytomer
     use pftconMod  , only : pftcon
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds           ! bounds
     integer , intent(in)  :: filter_vegsol(:)         ! filter for vegetated pfts with coszen>0
     integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0
     real(r8), intent(in)  :: coszen( bounds%begp: )   ! cosine solar zenith angle for next time step [pft]
    ! real(r8), intent(out) :: gfunc_solar(1:nlevcan)   !
    ! real(r8), intent(out) :: gfunc_sky(1:nlevcan,1:9) !
     type(surfalb_type)     , intent(inout) :: surfalb_inst

     ! !LOCAL VARIABLES:
     REAL, PARAMETER :: PI = 3.1415926_r8
     REAL, PARAMETER :: PI180 = 0.017453292_r8   ! PI divided by 180, radians per degree
     real(r8), parameter :: mpe = 1.e-06_r8 ! prevents overflow for division by zero
     integer  :: fp,p,iz,iv,il,ii       ! array indices
     integer  :: ixm    !index of the layer where the leaf inclination is close to mean
    ! integer  :: z_leaf
     real(r8) :: cosz, zen, laisum, z_ii
     real(r8) :: hth, ath, phi
     real(r8) :: Gfunc, z_mean, z_std
	 real(r8) :: z_dens(1:9)
     real(r8) :: diff(1:nlevcan)

     real(r8) :: chil(bounds%begp:bounds%endp)    ! -0.4 <= xl <= 0.6
     real(r8) :: phi1,phi2
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))

    associate(&
          ncan          =>    surfalb_inst%ncan_patch             , & ! Input:  [integer  (:)   ]  number of canopy layers
          nrad          =>    surfalb_inst%nrad_patch             , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          phytomer      =>    pftcon%phytomer        , & ! Input:  [integer (:)]  total number of phytomers in life time (added by Y.Fan)
          laimx         =>    pftcon%laimx           , & ! Input:
          mxlivenp      =>    pftcon%mxlivenp        , & ! Input: 
          theta         =>    pftcon%theta           , & ! Input:
          ivt           =>    patch%itype            , & ! Input:  [integer (:)]  pft vegetation type
          xl            =>    pftcon%xl              , & ! Input:  [real(r8) (:)]  ecophys const - leaf/stem orientation index
          tlai_z        =>    surfalb_inst%tlai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tlai increment for canopy layer
          tsai_z        =>    surfalb_inst%tsai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tsai increment for canopy layer
          gfunc_sky     =>    surfalb_inst%gfunc_sky_patch        , & ! Output:  [real(r8) (:,:,:)] the G function for diffuse sky radiation
          gfunc_solar   =>    surfalb_inst%gfunc_solar_patch      , & ! Output:  [real(r8) (:,:)] the G function for direct beam radiation
          lad           =>    surfalb_inst%lad_patch                & ! InOut: [real(r8) (:,:) ]  leaf angle distribution function for each canopy layer
!     z_mean                 =>    pps%z_mean             , & ! Output: [real(r8) (:,:)]  mean leaf inclination angle (LAD) per canopy layer
!     z_std                  =>    pps%z_std              , & ! Output: [real(r8) (:,:)]  standard deviation of LAD per canopy layer
!     z_dens                 =>    pps%z_dens             , & ! Output: [real(r8) (:,:)]  leaf angle distribution probability density per canopy layer
    )


    do fp = 1,num_vegsol
        p = filter_vegsol(fp)
      !two options of leaf angle distribution function (fLAD):
      !0= uniform (use xl index); 1 = dynamic (beta distribution) (Y.Fan 2015)
      if (fLAD == 0) then !!for xl leaf angle distribution
        !G function using xl index, canopy uniform leaf angle distribution like two-stream
        chil(p) = min( max(xl(ivt(p)), -0.4_r8), 0.6_r8 )
        if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8
        phi1 = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
        phi2 = 0.877_r8 * (1._r8-2._r8*phi1)
        !!direct beam: at a given solar direction
        cosz = max(0.001_r8, coszen(p))
        Gfunc = phi1 + phi2*cosz
        gfunc_solar(p,:) = max(mpe, Gfunc) !!all layers use the same G function
        !!diffuse: for different sky sectors
        do iz = 1, 9  !also 9 angle classes
          zen = ((iz - 1) * 10._r8 + 5._r8) * PI180
          cosz = max(0.001_r8, cos(zen))
          Gfunc = phi1 + phi2*cosz
          gfunc_sky(p,:,iz) = max(mpe, Gfunc)
        end do

      else !for leaf angle Beta distribution
        !find index of the layer where the leaf inclination is close to mean
        !z_mean always range from 90 degree at the top to 0 at the bottom
        !and the average z_mean of all layers equal to parameter theta
        !here theta is the mean inclination for the whole canopy
        do ixm = 1, nrad(p)-1
            diff(ixm) = abs( (sum(90._r8-(90._r8-theta(ivt(p)))/ixm*(/1:ixm:1/)) + sum(theta(ivt(p)) - &
                theta(ivt(p))/(nrad(p)-ixm)*(/1:(nrad(p)-ixm):1/)) ) / nrad(p) - theta(ivt(p)))
        end do
        ixm = sum(minloc(diff(1:(nrad(p)-1))))  !the layer with minimal difference from the mean inclination

        laisum = 0._r8
        do iv = 1, nrad(p) !from canopy top to bottom
		  !leaf inclination probability density function per canopy layer follows the Gaussian or Beta distribution (Y.Fan)
            ! z_leaf = (10 * il - 5)*PI180 !leaf angle interval centred at z_leaf
            ! z_dens(iv,il) = 1._r8/sqrt(2._r8*PI*z_var)*  &   !Gaussian distribution
			!        exp(-((z_leaf-z_mean)**2/(2._r8*z_var)))

          ! Cumulative lai+sai
          laisum = laisum +(tlai_z(p,iv)+tsai_z(p,iv))
          !Calculate the leaf angle probability density, LAD
		  !the fraction of leaf area within a leaf inclination angles interval

          if (phytomer(ivt(p)) > 0) then
		     !when canopy is still growing and for the top young fronds (45 degree)
!            if (laisum <= 1._r8) then
!               z_mean = 45*PI180
!            else
!!              z_mean = 45*PI180*((iv - nrad(p)/2._r8)/(nrad(p)/2._r8))**2  !mean leaf inclination per layer
!               z_mean = maxincl(ivt(p))*PI180*abs((iv - 0.8_r8*nrad(p))/(0.8_r8*nrad(p)))
!              !z_std = 1._r8/(2._r8*PI)*(1._r8+(iv /mxlivenp(ivt(p)))**2)  ! variance = z_std**2
!!              z_mean = PI/2._r8*((laisum - laimx(ivt(p))*mxlivenp(ivt(p))/2._r8)/ &
!!		             (laimx(ivt(p))*mxlivenp(ivt(p))/2._r8))**2  !mean leaf inclination per layer
!              !z_std = 1._r8/(2._r8*PI)*(1._r8+(laisum/(mxlivenp(ivt(p))*laimx(ivt(p))))**2)

             if (iv <= ixm) then
                z_mean = (90._r8-(90._r8-theta(ivt(p)))/ixm*iv) * PI180
             else
                z_mean = (theta(ivt(p)) - theta(ivt(p))/(nrad(p)-ixm)*(iv-ixm)) * PI180
             end if
          !other PFTs follow Canoak and GOEL AND STREBEL (1984)
          else if (laisum < 2.6) then
             z_mean = 45*PI180
          else
             z_mean = (78.99 - 14.84 * laisum + 0.244 * laisum * laisum)*PI180
          end if

          z_mean = max(10*PI180, min(80*PI180, z_mean))  !within 10-80 degree
          z_std = 18*PI180 !standard deviation is fixed for Beta distribution

          call Beta_Dist(z_mean,z_std,z_dens)
          lad(p,iv,:) = z_dens
        end do

	!-----------------------------------
		!for direct beam
        iz = 1  ! a single solar direction for direct beam
        cosz = max(0.001_r8, coszen(p))
        !cosz = max(mpe, coszen(p))
        zen = acos(cosz) !inverse of cos() to obtain solar zenith angle
        zen = min(max(zen, mpe), (PI/2._r8-mpe)) !0< zen <90
        !Compute the G function for each canopy layer

        do iv = 1, nrad(p) !from canopy top to bottom
          Gfunc = 0._r8
          do il = 1, 9   !9 leaf normal zenith angle classes: 5 to 85 degree
             hth = 0._r8
             do ii = (10*il-10), (10*il-1), 1  !integral between intervals, increment of 1 degree
               z_ii = (ii+0.5_r8) * PI180 !intervals: 0.5~9.5, 10.5~19.5,..., 80.5~89.5
			   if (abs(1._r8/(tan(zen)*tan(z_ii))) > 1._r8) then
			     ath = cos(zen)*cos(z_ii)
			   else
			     phi = acos(1._r8/(tan(zen)*tan(z_ii)))
			     ath = cos(zen)*cos(z_ii)*(1._r8 + 2._r8/PI*(tan(phi)-phi))
			   end if !ath is the G function at a given leaf inclination angle
               hth = hth + ath * 1._r8*PI180
             end do
             hth = hth /(10._r8*PI180)  !hth is the average G function for a 10 degree interval
             Gfunc = Gfunc + hth * lad(p,iv,il)
          end do
          gfunc_solar(p,iv) = max(mpe, Gfunc)
        end do

	!-----------------------------------
      !for diffuse sky radiation
      !diffuse extinction function will integrate over all sky sector
       do iz = 1, 9  !also 9 angle classes
        zen = ((iz - 1) * 10._r8 + 5._r8) * PI180

        !Compute the G function for each canopy layer
        do iv = 1, nrad(p) !from canopy top to bottom
          Gfunc = 0._r8
          do il = 1, 9 !!leaf inclination angle: 5 to 85 degree
             hth = 0._r8
             do ii = (10*il-10), (10*il-1), 1  !integral between leaf angle intervals, increment of 1 degree
               z_ii = (ii+0.5_r8) * PI180 !intervals: 0.5~9.5, 10.5~19.5,..., 80.5~89.5
			   if (abs(1._r8/(tan(zen)*tan(z_ii))) > 1._r8) then
			     ath = cos(zen)*cos(z_ii)
			   else
			     phi = acos(1._r8/(tan(zen)*tan(z_ii)))
			     ath = cos(zen)*cos(z_ii)*(1._r8 + 2._r8/PI*(tan(phi)-phi))
			   end if
               hth = hth + ath * 1._r8*PI180
             end do
             hth = hth /(10._r8*PI180)	!mean value for each 10 degree interval: very important!
             Gfunc = Gfunc + hth * lad(p,iv,il)
          end do
		  gfunc_sky(p,iv,iz) = max(mpe, Gfunc)
        end do
       end do

      end if

    end do   ! end of pft loop

    end associate
   end subroutine Gfunc_all


   !-----------------------------------------------------------------------
   subroutine Beta_Dist (mean, std, dens)
     !
     ! !DESCRIPTION: (Y.Fan 01.12.2014)
       ! this program uses the beta distribution to compute the probability density
       ! distribution for a known mean leaf inclination angle
       ! starting from the top of the canopy to bottom where MEAN between 0 and PI/2
       ! after Goel and Strebel (1984) and Wang et al. (2007) and some errors modified

     ! !ARGUMENTS:
     real(r8), intent(in)  :: mean   ! mean leaf inclination angle (LAD) per canopy layer
     real(r8), intent(in)  :: std    ! standard deviation of LAD per canopy layer
     real(r8), intent(out) :: dens(1:9) ! leaf angle distribution probability density per canopy layer
     !
     ! !LOCAL VARIABLES:
     real, PARAMETER :: PI = 3.1415926_r8
     real, PARAMETER :: PI180 = 0.017453292_r8   ! PI divided by 180, radians per degree
     real(r8), parameter :: mpe = 1.e-06_r8 ! prevents overflow for division by zero
     integer  :: il
     real(r8) :: ang, mean_t, std_t, var_t, cons
     real(r8) :: nu, mu, numu, beta, nu1, mu1

    !-----------------------------------------------------------------------
	    ! 10<mean<80 for 9 angle classes, to ensure integral of dens(:) is close to 1
        mean_t = min(max(2._r8*mean/PI, 10._r8/90),80._r8/90)
		! 5<std<10 degree for 9 angle classes to ensure integral of dens(:) is close to 1
        std_t = min(max(2._r8*std/PI, 5._r8/90), 10._r8/90)
        var_t = mean_t*mean_t + std_t*std_t

        nu = mean_t*(mean_t - var_t)/(std_t*std_t)
        mu = nu*(1._r8/mean_t - 1._r8)
        numu = nu + mu

       !The gamma function is defined for all complex numbers except the negative integers and zero.
       !The result of gamma function is nonzero everywhere along the real line
        beta = (gamma(nu) * gamma(mu)) / max(gamma(numu), mpe)
        mu1 = mu - 1._r8
        nu1 = nu - 1._r8

        !Compute probability distribution for 9 angle classes
        !between 5 and 85 degrees, with increments of 10 degrees
        cons = 1._r8 / 9
        do il =1, 9
           ang = (10._r8 * il - 5._r8)*PI180
           dens(il) = cons * ((1._r8 - 2._r8*ang/PI)**mu1) * &
		          ((2._r8*ang/PI)**nu1) / max(beta, mpe)
        end do

   end subroutine Beta_Dist

   !-----------------------------------------------------------------------
!   subroutine G_FUNC (bounds, filter_vegsol, num_vegsol, iv, iz, bang)
     !
     ! !DESCRIPTION: (Y.Fan 01.12.2014)
       !This subroutine computes the G Function according to
       !the algorithms of Lemeur (1973, Agric. Meteorol. 12: 229-247).
       !Computes G for each sky sector, as needed to compute the transmission of diffuse light
       !G varies with canopy height to account for vertical variations in leaf angles

	   !G_Func for each incident zenith angle
	   !for direct beam, bang= solar elevation angle
	   !for diffuse, bang= divide sky as 5 to 85 degree sectors

       !This following subroutine computes the G Function according to Lemeur (1973)
       ! Midpoints for azimuth intervals (PI/8), 16 azimuth sections
!        do i = 1, 17
!          azi(i) = (2*i - 3)*PI/16._r8
!        end do
!
!        do i = 1, 16
!          del_azi(i) = PI/8._r8
!          del_sin(i) = sin(azi(i+1)) - sin(azi(i))
!          a_dens(i) = 1._r8/16._r8   !leaf azimuth density function (assume symmetric distribution here)
!        end do

!		PPP = 0.0
!		do ia = 1, 9 !leaf inclination a angle: 5 to 85 degree
!		aang = ((ia - 1) * 10._r8 + 5._r8) * PI180
!		x = cos(aang)*sin(bang)
!		y = sin(aang)*cos(bang)
!		if (aang <= bang) then
!		  do i = 1, 16 !leaf azimuth classes
!			PGF(i) = x * del_azi(i) + y * del_sin(i)
!		  end do
!		else
!		  TT0 = 2._r8 * atan(sqrt(max(0._r8,(1._r8 + x/y)/(1._r8 - x/y))))
!		  TT1 = 2._r8 * PI - TT0
!		  do i = 1, 16
!			if (azi(i+1) <= TT0) then
!			   PGF(i) = x * del_azi(i) + y * del_sin(i)
!			   continue
!			else
!			   if (azi(i+1) <= TT1) then
!				 if (azi(i) >= TT0) then
!				   PGF(i) = -x * del_azi(i) - y * del_sin(i)
!				   continue
!				 else
!				   PGF(i) = x * (TT0-azi(i)) + y * (sin(TT0) - sin(azi(i))) - &
!						  (x *(azi(i+1)-TT0) +  y *(sin(azi(i+1)) - sin(TT0)))
!				   continue
!				 end if
!			   else
!				 if (azi(i) >= TT1) then
!				   PGF(i) = x * del_azi(i) + y * del_sin(i)
!				   continue
!				 else
!				   PGF(i) = (x *(azi(i+1)-TT1) +  y *(sin(azi(i+1)) - sin(TT1))) - &
!						  (x * (TT1-azi(i)) + y * (sin(TT1) - sin(azi(i))))
!				   continue
!				 end if
!			   end if
!			end if
!		  end do
!		end if
!		PP = 0._r8
!		do i = 1, 16
!		  PP = PP + (PGF(i)*adens(i))
!		end do
!		PPP = PPP + (PP*bdens(ia)*9._r8/PI)
!		end do
!		Gfunc(iv,iz) = PPP

!   end subroutine G_FUNC

end module SurfaceAlbedoMod
