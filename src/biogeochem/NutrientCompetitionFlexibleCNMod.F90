module NutrientCompetitionFlexibleCNMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to do soil nutrient competition dynamics
  !
  ! FIXME(bja, 2015-08) This module was copied from
  ! NutrientCompetitionCLM45default then flexible cn modifications
  ! were added for the clm50 nitrogen science changes (r120). There is
  ! a significant amount of duplicate code between the two
  ! modules. They need to be reexamined and the common code pulled out
  ! into a common base class.
  !
  ! created by Jinyun Tang, Sep 8, 2014
  ! modified by Mariana Vertenstein, Nov 15, 2014
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : patch
  use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
  use NutrientCompetitionMethodMod, only : params_inst
  use clm_varctl          , only : iulog
  use pftconMod           , only : mxnp  !max number of phytomers for multilayer structure
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: nutrient_competition_FlexibleCN_type
  !
  type, extends(nutrient_competition_method_type) :: nutrient_competition_FlexibleCN_type
     private
     real(r8), pointer :: actual_leafcn(:)                    ! leaf CN ratio used by flexible CN
     real(r8), pointer :: actual_storage_leafcn(:)            ! storage leaf CN ratio used by flexible CN
     real(r8), pointer :: afruitn(:,:),aleafn(:,:)            ! fruit and leaf allocation for each phytomer
   contains
     ! public methocs
     procedure, public :: Init                                ! Initialization
     procedure, public :: calc_plant_nutrient_competition     ! calculate nutrient yield rate from competition
     procedure, public :: calc_plant_nutrient_demand          ! calculate plant nutrient demand
     !
     ! private methods
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: calc_plant_cn_alloc
     procedure, private :: calc_plant_nitrogen_demand
  end type nutrient_competition_FlexibleCN_type
  !
  interface nutrient_competition_FlexibleCN_type
     ! initialize a new nutrient_competition_FlexibleCN_type object
     module procedure constructor
  end interface nutrient_competition_FlexibleCN_type
  !

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  type(nutrient_competition_FlexibleCN_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type nutrient_competition_FlexibleCN_type.
    ! For now, this is simply a place-holder.
  end function constructor

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize the class
    !
    class(nutrient_competition_FlexibleCN_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for the class data
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)

    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type) :: this
    type(bounds_type), intent(in) :: bounds

    allocate(this%actual_leafcn(bounds%begp:bounds%endp))         ; this%actual_leafcn(:)         = nan
    allocate(this%actual_storage_leafcn(bounds%begp:bounds%endp)) ; this%actual_storage_leafcn(:) = nan
    !!!!new allocation parameters for multilayer phytomer structure (Y.Fan)
    if ( mxnp > 0 )then
       allocate(this%afruitn(bounds%begp:bounds%endp,1:mxnp))
       allocate(this%aleafn(bounds%begp:bounds%endp,1:mxnp))
       this%afruitn(bounds%begp:bounds%endp,1:mxnp) = nan
       this%aleafn(bounds%begp:bounds%endp,1:mxnp) = nan
    end if
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Send data to history file
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    use clm_varcon     , only : spval
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%actual_leafcn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFCN', units='gC/gN', &
         avgflag='A', long_name='Leaf CN ratio used for flexible CN', &
         ptr_patch=this%actual_leafcn )
    this%actual_storage_leafcn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFCN_STORAGE', units='gC/gN', &
         avgflag='A', long_name='Storage Leaf CN ratio used for flexible CN', &
         ptr_patch=this%actual_storage_leafcn, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine calc_plant_nutrient_competition (this, &
       bounds, num_soilp, filter_soilp, &
       cnveg_state_inst, crop_inst, canopystate_inst, cnveg_carbonstate_inst, &
       cnveg_carbonflux_inst, &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst, &
       aroot, arepr, fpg_col)
    !
    ! !USES:
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    real(r8), intent(in)    :: aroot   (bounds%begp:)
    real(r8), intent(in)    :: arepr   (bounds%begp:)
    real(r8), intent(in)    :: fpg_col (bounds%begc:)

    call this%calc_plant_cn_alloc(bounds, num_soilp, filter_soilp,   &
         cnveg_state_inst, crop_inst, canopystate_inst, &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
         c14_cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
         soilbiogeochem_nitrogenstate_inst, &
         aroot=aroot(bounds%begp:bounds%endp),                               &
         arepr=arepr(bounds%begp:bounds%endp),                               &
         fpg_col=fpg_col(bounds%begc:bounds%endc))

  end subroutine calc_plant_nutrient_competition

!-----------------------------------------------------------------------
  subroutine calc_plant_cn_alloc(this, bounds, num_soilp, filter_soilp,   &
       cnveg_state_inst, crop_inst, canopystate_inst, &
       cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
       c14_cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst, &
       aroot, arepr, fpg_col)
    !
    ! !USES:
    use pftconMod             , only : pftcon, npcropmin
    use clm_varctl            , only : use_c13, use_c14, carbon_resp_opt
    use clm_varctl            , only : downreg_opt
    use clm_varctl            , only : CN_residual_opt
    use clm_varctl            , only : CN_partition_opt
    use clm_time_manager       , only : get_step_size
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
    use CNSharedParamsMod     , only : use_fun
    use CNPrecisionControlMod , only : n_min
    use clm_varcon            , only : spval
    
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)  , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    real(r8)                        , intent(in)    :: aroot(bounds%begp:)
    real(r8)                        , intent(in)    :: arepr(bounds%begp:)
    real(r8)                        , intent(in)    :: fpg_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p            ! indices
    integer  :: fp                 ! lake filter patch index
    real(r8) :: f1,f2,f3,f4,g1,g2  ! allocation parameters
    real(r8) :: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8) :: fcur               ! fraction of current psn displayed as growth
    real(r8) :: gresp_storage      ! temporary variable for growth resp to storage
    real(r8) :: nlc                ! temporary variable for total new leaf carbon allocation
    real(r8) :: f5                 ! grain allocation parameter
    real(r8) :: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8) :: dt                 ! model time step
    real(r8):: fsmn(bounds%begp:bounds%endp)  ! A emperate variable for adjusting FUN uptakes
	
    real(r8):: frootcn_storage_actual       										
    real(r8):: frootcn_actual               	
    real(r8):: livestemcn_storage_actual    	
    real(r8):: livestemcn_actual            	
    real(r8):: livecrootcn_storage_actual   			
    real(r8):: livecrootcn_actual           	
    real(r8):: leafcn_max                   	
    real(r8):: frootcn_max                  	
    real(r8):: livewdcn_max  
    real(r8):: frac_resp    
    real(r8) :: npool_to_leafn_demand                           (bounds%begp:bounds%endp)
    real(r8) :: npool_to_leafn_supply                           (bounds%begp:bounds%endp)
    real(r8) :: npool_to_leafn_storage_demand                   (bounds%begp:bounds%endp)
    real(r8) :: npool_to_leafn_storage_supply                   (bounds%begp:bounds%endp)
    real(r8) :: npool_to_frootn_demand                          (bounds%begp:bounds%endp)
    real(r8) :: npool_to_frootn_supply                          (bounds%begp:bounds%endp)
    real(r8) :: npool_to_frootn_storage_demand                  (bounds%begp:bounds%endp)
    real(r8) :: npool_to_frootn_storage_supply                  (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livestemn_demand                       (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livestemn_supply                       (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livestemn_storage_demand               (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livestemn_storage_supply               (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livecrootn_demand                      (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livecrootn_supply                      (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livecrootn_storage_demand              (bounds%begp:bounds%endp)
    real(r8) :: npool_to_livecrootn_storage_supply              (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadstemn_demand                       (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadstemn_supply                       (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadstemn_storage_demand               (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadstemn_storage_supply               (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadcrootn_demand                      (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadcrootn_supply                      (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadcrootn_storage_demand              (bounds%begp:bounds%endp)
    real(r8) :: npool_to_deadcrootn_storage_supply              (bounds%begp:bounds%endp)
    real(r8) :: npool_to_grainn_demand                          (bounds%begp:bounds%endp)
    real(r8) :: npool_to_grainn_supply                          (bounds%begp:bounds%endp)
    real(r8) :: npool_to_grainn_storage_demand                  (bounds%begp:bounds%endp)
    real(r8) :: npool_to_grainn_storage_supply                  (bounds%begp:bounds%endp)
    real(r8) :: total_plant_Ndemand                             (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_leafn                        (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_leafn_storage        (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_frootn               (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_frootn_storage       (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_livestemn            (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_livestemn_storage    (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_deadstemn            (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_deadstemn_storage    (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_livecrootn           (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_livecrootn_storage   (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_deadcrootn           (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_deadcrootn_storage   (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_grainn               (bounds%begp:bounds%endp)
    real(r8) :: frNdemand_npool_to_grainn_storage       (bounds%begp:bounds%endp)

    ! -----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(aroot)   == (/bounds%endp/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(arepr)   == (/bounds%endp/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(fpg_col) == (/bounds%endc/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%actual_storage_leafcn) >= (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%actual_storage_leafcn) <= (/bounds%begp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%aleafn) >= (/bounds%endp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%aleafn) <= (/bounds%begp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%afruitn) >= (/bounds%endp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%afruitn) <= (/bounds%begp,mxnp/)), errMsg(sourcefile, __LINE__))
	
    associate(                                                                                       &
         fpg                          => fpg_col                                                   , & ! Input:  [real(r8) (:)   ]  fraction of potential gpp (no units)

         ivt                          => patch%itype                                                 , & ! Input:  [integer  (:) ]  patch vegetation type

         woody                        => pftcon%woody                                              , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf                   => pftcon%froot_leaf                                         , & ! Input:  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem                   => pftcon%croot_stem                                         , & ! Input:  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf                    => pftcon%stem_leaf                                          , & ! Input:  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd                      => pftcon%flivewd                                            , & ! Input:  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn                       => pftcon%leafcn                                             , & ! Input:  leaf C:N (gC/gN)
         frootcn                      => pftcon%frootcn                                            , & ! Input:  fine root C:N (gC/gN)
         livewdcn                     => pftcon%livewdcn                                           , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn                     => pftcon%deadwdcn                                           , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         fcur2                        => pftcon%fcur                                               , & ! Input:  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
         graincn                      => pftcon%graincn                                            , & ! Input:  grain C:N (gC/gN)
         grperc                       => pftcon%grperc                                             , & ! Input:  growth respiration parameter
         grpnow                       => pftcon%grpnow                                             , & ! Input:  growth respiration parameter
         evergreen                    => pftcon%evergreen                                          , & ! Input:  binary flag for evergreen leaf habit (0 or 1)
         phytomer                     => pftcon%phytomer                                   , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
         !leafcnr                   =>    pftcon%leafcnr                                            , & ! Input:  [real(r8) (:)]  range of departure from default leaf C:N ratio, used to determine max/min C:N (gC/gN)
         !frootcnr                  =>    pftcon%frootcnr                                           , & ! Input:  [real(r8) (:)]  range of departure from default fineroot C:N ratio, used to determine max/min C:N (gC/gN)
         !livewdcnr                 =>    pftcon%livewdcnr                                          , & ! Input:  [real(r8) (:)]  range of departure from default livestem C:N ratio, used to determine max/min C:N (gC/gN)
         !graincnr                  =>    pftcon%graincnr                                           , & ! Input:  [real(r8) (:)]  range of departure from default grain C:N ratio, used to determine max/min C:N (gC/gN)
         huileafnp                    => crop_inst%huileafnp_patch                , & ! Input:  [real(r8) (:,:)]  hui needed for initiation of successive phytomers
         huilfexpnp                   => crop_inst%huilfexpnp_patch               , & ! Input:  [real(r8) (:,:)]  hui needed for leaf expansion of successive phytomers
         huilfmatnp                   => crop_inst%huilfmatnp_patch               , & ! Input:  [real(r8) (:,:)]  hui needed for leaf maturity of successive phytomers
         pleafc                    => cnveg_carbonstate_inst%pleafc_patch                   , & ! InOut:  [real(r8) (:,:)]  (gC/m2) phytomer leaf C
         pleafn                    => cnveg_nitrogenstate_inst%pleafn_patch                 , & ! InOut:  [real(r8) (:,:)]  (gN/m2) phytomer leaf N
         cpool_to_pleafc           => cnveg_carbonflux_inst%cpool_to_pleafc_patch           , & ! InOut:  [real(r8) (:,:)]  allocation to phytomer leaf C (gC/m2/s)
         cpool_to_pleafc_storage   => cnveg_carbonflux_inst%cpool_to_pleafc_storage_patch   , & ! Input:  [real(r8) (:,:)]
         cpool_to_pgrainc          => cnveg_carbonflux_inst%cpool_to_pgrainc_patch          , & ! InOut:  [real(r8) (:,:)]  allocation to phytomer grain C (gC/m2/s)
         npool_to_pleafn           => cnveg_nitrogenflux_inst%npool_to_pleafn_patch         , & ! InOut:  [real(r8) (:,:)]  allocation to phytomer leaf N (gN/m2/s)
         npool_to_pleafn_storage   => cnveg_nitrogenflux_inst%npool_to_pleafn_storage_patch  , & ! Input:  [real(r8) (:,:)]
         npool_to_pgrainn          => cnveg_nitrogenflux_inst%npool_to_pgrainn_patch         , & ! InOut:  [real(r8) (:,:)]  allocation to phytomer grain N (gN/m2/s)

         croplive                     => crop_inst%croplive_patch                                  , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         peaklai                      => cnveg_state_inst%peaklai_patch                            , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         aleaf                        => cnveg_state_inst%aleaf_patch                              , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnveg_state_inst%astem_patch                              , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         c_allometry                  => cnveg_state_inst%c_allometry_patch                        , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry                  => cnveg_state_inst%n_allometry_patch                        , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         downreg                      => cnveg_state_inst%downreg_patch                            , & ! Output: [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)

         annsum_npp                   => cnveg_carbonflux_inst%annsum_npp_patch                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         gpp                          => cnveg_carbonflux_inst%gpp_before_downreg_patch            , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                       => cnveg_carbonflux_inst%availc_patch                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         excess_cflux                 => cnveg_carbonflux_inst%excess_cflux_patch                  , & ! Output: [real(r8) (:)   ]  C flux not allocated due to downregulation (gC/m2/s)
         plant_calloc                 => cnveg_carbonflux_inst%plant_calloc_patch                  , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)
         npp_growth                  => cnveg_carbonflux_inst%npp_growth_patch              , & ! Output:  [real(r8) (:) ] C for growth in FUN. g/m2/s
         cpool_to_resp                => cnveg_carbonflux_inst%cpool_to_resp_patch                 , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_resp          => cnveg_carbonflux_inst%cpool_to_leafc_resp_patch              , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_storage_resp  => cnveg_carbonflux_inst%cpool_to_leafc_storage_resp_patch      , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_resp         => cnveg_carbonflux_inst%cpool_to_frootc_resp_patch             , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_storage_resp => cnveg_carbonflux_inst%cpool_to_frootc_storage_resp_patch     , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_resp     => cnveg_carbonflux_inst%cpool_to_livecrootc_resp_patch         , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_storage_resp  => cnveg_carbonflux_inst%cpool_to_livecrootc_storage_resp_patch , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_resp      => cnveg_carbonflux_inst%cpool_to_livestemc_resp_patch          , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_storage_resp => cnveg_carbonflux_inst%cpool_to_livestemc_storage_resp_patch  , & ! Output: [real(r8) (:)   ]
         psnsun_to_cpool              => cnveg_carbonflux_inst%psnsun_to_cpool_patch               , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool            => cnveg_carbonflux_inst%psnshade_to_cpool_patch             , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc               => cnveg_carbonflux_inst%cpool_to_leafc_patch                , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_storage       => cnveg_carbonflux_inst%cpool_to_leafc_storage_patch        , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc              => cnveg_carbonflux_inst%cpool_to_frootc_patch               , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_storage      => cnveg_carbonflux_inst%cpool_to_frootc_storage_patch       , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc           => cnveg_carbonflux_inst%cpool_to_livestemc_patch            , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_storage   => cnveg_carbonflux_inst%cpool_to_livestemc_storage_patch    , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc           => cnveg_carbonflux_inst%cpool_to_deadstemc_patch            , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc_storage   => cnveg_carbonflux_inst%cpool_to_deadstemc_storage_patch    , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc          => cnveg_carbonflux_inst%cpool_to_livecrootc_patch           , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_storage  => cnveg_carbonflux_inst%cpool_to_livecrootc_storage_patch   , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc          => cnveg_carbonflux_inst%cpool_to_deadcrootc_patch           , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc_storage  => cnveg_carbonflux_inst%cpool_to_deadcrootc_storage_patch   , & ! Output: [real(r8) (:)   ]
         cpool_to_gresp_storage       => cnveg_carbonflux_inst%cpool_to_gresp_storage_patch        , & ! Output: [real(r8) (:)   ]  allocation to growth respiration storage (gC/m2/s)
         cpool_to_grainc              => cnveg_carbonflux_inst%cpool_to_grainc_patch               , & ! Output: [real(r8) (:)   ]  allocation to grain C (gC/m2/s)
         cpool_to_grainc_storage      => cnveg_carbonflux_inst%cpool_to_grainc_storage_patch       , & ! Output: [real(r8) (:)   ]  allocation to grain C storage (gC/m2/s)
         
         laisun                       => canopystate_inst%laisun_patch  , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                       => canopystate_inst%laisha_patch  , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index
         smin_no3_vr                  => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col         , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3
         leafn                        => cnveg_nitrogenstate_inst%leafn_patch                      , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         leafn_storage                => cnveg_nitrogenstate_inst%leafn_storage_patch              , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         npool                        => cnveg_nitrogenstate_inst%npool_patch                      , & ! Input:  [real(r8) (:)   ]  (gN/m2) temporary plant N pool
         plant_ndemand                => cnveg_nitrogenflux_inst%plant_ndemand_patch               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_nalloc                 => cnveg_nitrogenflux_inst%plant_nalloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)
         npool_to_grainn              => cnveg_nitrogenflux_inst%npool_to_grainn_patch             , & ! Output: [real(r8) (:)   ]  allocation to grain N (gN/m2/s)
         npool_to_grainn_storage      => cnveg_nitrogenflux_inst%npool_to_grainn_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to grain N storage (gN/m2/s)
         retransn_to_npool            => cnveg_nitrogenflux_inst%retransn_to_npool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         sminn_to_npool               => cnveg_nitrogenflux_inst%sminn_to_npool_patch              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
         npool_to_leafn               => cnveg_nitrogenflux_inst%npool_to_leafn_patch              , & ! Output: [real(r8) (:)   ]  allocation to leaf N (gN/m2/s)
         npool_to_leafn_storage       => cnveg_nitrogenflux_inst%npool_to_leafn_storage_patch      , & ! Output: [real(r8) (:)   ]  allocation to leaf N storage (gN/m2/s)
         npool_to_frootn              => cnveg_nitrogenflux_inst%npool_to_frootn_patch             , & ! Output: [real(r8) (:)   ]  allocation to fine root N (gN/m2/s)
         npool_to_frootn_storage      => cnveg_nitrogenflux_inst%npool_to_frootn_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to fine root N storage (gN/m2/s)
         npool_to_livestemn           => cnveg_nitrogenflux_inst%npool_to_livestemn_patch          , & ! Output: [real(r8) (:)   ]
         npool_to_livestemn_storage   => cnveg_nitrogenflux_inst%npool_to_livestemn_storage_patch  , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn           => cnveg_nitrogenflux_inst%npool_to_deadstemn_patch          , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn_storage   => cnveg_nitrogenflux_inst%npool_to_deadstemn_storage_patch  , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn          => cnveg_nitrogenflux_inst%npool_to_livecrootn_patch         , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn_storage  => cnveg_nitrogenflux_inst%npool_to_livecrootn_storage_patch , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn          => cnveg_nitrogenflux_inst%npool_to_deadcrootn_patch         , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn_storage  => cnveg_nitrogenflux_inst%npool_to_deadcrootn_storage_patch , & ! Output: [real(r8) (:)   ]
         Npassive                     => cnveg_nitrogenflux_inst%Npassive_patch                    , & ! Output:  [real(r8) (:) ]  Passive N uptake (gN/m2/s)
         Nfix                         => cnveg_nitrogenflux_inst%Nfix_patch                        , & ! Output:  [real(r8) (:) ]  Symbiotic BNF (gN/m2/s)
         Nactive                      => cnveg_nitrogenflux_inst%Nactive_patch                     , & ! Output:  [real(r8) (:) ] Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc                      => cnveg_nitrogenflux_inst%Nnonmyc_patch                     , & ! Output:  [real(r8) (:) ] Non-mycorrhizal N uptake (gN/m2/s)
         Nam                          => cnveg_nitrogenflux_inst%Nam_patch                         , & ! Output:  [real(r8) (:) ]  AM uptake (gN/m2/s)
         Necm                         => cnveg_nitrogenflux_inst%Necm_patch                        , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
         sminn_to_plant_fun           => cnveg_nitrogenflux_inst%sminn_to_plant_fun_patch            & ! Output:  [real(r8) (:) ]  Total soil N uptake of FUN (gN/m2/s)

         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! patch loop to distribute the available N between the competing patches
      ! on the basis of relative demand, and allocate C and N to new growth and storage

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         ! set some local allocation variables
         f1 = froot_leaf(ivt(p))
         f2 = croot_stem(ivt(p))

         ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
         ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
         ! There was an error in this formula in previous version, where the coefficient
         ! was 0.004 instead of 0.0025.
         ! This variable allocation is only for trees. Shrubs have a constant
         ! allocation as specified in the pft-physiology file.  The value is also used
         ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
         if (stem_leaf(ivt(p)) == -1._r8) then
            f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
         else
            f3 = stem_leaf(ivt(p))
         end if

         f4   = flivewd(ivt(p))
         g1   = grperc(ivt(p))
         g2   = grpnow(ivt(p))
         cnl  = leafcn(ivt(p))
         cnfr = frootcn(ivt(p))
         cnlw = livewdcn(ivt(p))
         cndw = deadwdcn(ivt(p))
         fcur = fcur2(ivt(p))

         if (.not. downreg_opt) then
            if (evergreen(ivt(p)) == 1._r8) then
               fcur = 0.0_r8
            end if
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            if (croplive(p)) then
               f1 = aroot(p) / aleaf(p)
               f3 = astem(p) / aleaf(p)
               f5 = arepr(p) / aleaf(p)
               g1 = 0.25_r8
            else
               f1 = 0._r8
               f3 = 0._r8
               f5 = 0._r8
               g1 = 0.25_r8
            end if
         end if

         ! increase fcur linearly with ndays_active, until fcur reaches 1.0 at
         ! ndays_active = days/year.  This prevents the continued storage of C and N.
         ! turning off this correction (PET, 12/11/03), instead using bgtr in
         ! phenology algorithm.
         
         
         if(use_fun)then ! if we are using FUN, we get the N available from there.
           sminn_to_npool(p) = sminn_to_plant_fun(p)
         else ! no FUN. :( we get N available from the FPG calculation in soilbiogeochemistry competition. 
           sminn_to_npool(p) = plant_ndemand(p) * fpg(c)        
         endif
       
         plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)
         
         if(.not.use_fun)then
	         if (downreg_opt) then
	            ! calculate the associated carbon allocation, and the excess
	            ! carbon flux that must be accounted for through downregulation
	            plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
	            excess_cflux(p) = availc(p) - plant_calloc(p)

	            ! reduce gpp fluxes due to N limitation
	            if (gpp(p) > 0.0_r8) then
	               downreg(p) = excess_cflux(p)/gpp(p)

	               psnsun_to_cpool(p)   = psnsun_to_cpool(p)  *(1._r8 - downreg(p))
	               psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1._r8 - downreg(p))

	               if ( use_c13 ) then
	                  c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
	                       c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
	                  c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
	                       c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
	               endif
	               if ( use_c14 ) then
	                  c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
	                       c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
	                  c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
	                       c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
	               endif
	            end if
	         end if
         end if

         
         if(use_fun)then
            plant_calloc(p) = npp_growth(p) 
         else
            if (.not. downreg_opt) then
               plant_calloc(p) = availc(p)
            end if         
         end if
        
         ! calculate the amount of new leaf C dictated by these allocation
         ! decisions, and calculate the daily fluxes of C and N to current
         ! growth and storage pools

         ! fcur is the proportion of this day's growth that is displayed now,
         ! the remainder going into storage for display next year through the
         ! transfer pools

         nlc = plant_calloc(p) / c_allometry(p)
         cpool_to_leafc(p)          = nlc * fcur
         cpool_to_leafc_storage(p)  = nlc * (1._r8 - fcur)
         cpool_to_frootc(p)         = nlc * f1 * fcur
         cpool_to_frootc_storage(p) = nlc * f1 * (1._r8 - fcur)
         if (woody(ivt(p)) == 1._r8) then
            cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_grainc(p)             = nlc * f5 * fcur
            cpool_to_grainc_storage(p)     = nlc * f5 * (1._r8 -fcur)
         end if
	 !phytomer-based allocation (Y.Fan)
	 !fcur is used differently than other PFTs
       if (phytomer(ivt(p)) > 0) then 
         cpool_to_leafc(p)              = nlc * fcur
         cpool_to_leafc_storage(p)      = nlc * (1._r8 - fcur) !only allow leaf storage growth
         cpool_to_frootc(p)             = nlc * f1
	 cpool_to_frootc_storage(p)     = 0.0_r8
         cpool_to_livestemc(p)          = nlc * f3 * f4
	 cpool_to_livestemc_storage(p)  = 0.0_r8
         cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4)
	 cpool_to_deadstemc_storage(p)  = 0.0_r8
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4
	 cpool_to_livecrootc_storage(p) = 0.0_r8
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4)
	 cpool_to_deadcrootc_storage(p) = 0.0_r8
         cpool_to_grainc(p)             = nlc * f5
	 cpool_to_grainc_storage(p)     = 0.0_r8

         where (hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:))
            cpool_to_pleafc_storage(p,:) = cpool_to_leafc_storage(p) * this%aleafn(p,:)
         elsewhere (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:))
            cpool_to_pleafc(p,:) = cpool_to_leafc(p) * this%aleafn(p,:)
         endwhere
         cpool_to_pgrainc(p,:) = cpool_to_grainc(p) * this%afruitn(p,:)
       end if	
	   
         if (downreg_opt) then
            ! corresponding N fluxes
            npool_to_leafn(p)          = (nlc / cnl) * fcur
            npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur)
            npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
            npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
            if (woody(ivt(p)) == 1._r8) then
               npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops
               cng = graincn(ivt(p))
               npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_grainn(p)             = (nlc * f5 / cng) * fcur
               npool_to_grainn_storage(p)     = (nlc * f5 / cng) * (1._r8 -fcur)
            end if
			
	    if (phytomer(ivt(p)) > 0) then ! phytomer-based allocation (Y.Fan)
		 cng = graincn(ivt(p))
		 !only use fcur for storage growth of spear leaf of oil palm
		 !erase other storage terms to zero
		 npool_to_leafn(p)          = (nlc / cnl) * fcur
		 npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur) 
		 npool_to_frootn(p)         = (nlc * f1 / cnfr)
		 npool_to_frootn_storage(p) = 0.0_r8
		 
		 npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw)
		 npool_to_livestemn_storage(p)  = 0.0_r8
		 npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw)
		 npool_to_deadstemn_storage(p)  = 0.0_r8
		 npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw)
		 npool_to_livecrootn_storage(p) = 0.0_r8
		 npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw)
		 npool_to_deadcrootn_storage(p) = 0.0_r8
		 npool_to_grainn(p)             = (nlc * f5 / cng)
     		 npool_to_grainn_storage(p)     = 0.0_r8
		 
		 where (hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:))
			npool_to_pleafn_storage(p,:) = npool_to_leafn_storage(p) * this%aleafn(p,:)
		 elsewhere (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:))
			npool_to_pleafn(p,:) = npool_to_leafn(p) * this%aleafn(p,:)
		 endwhere
		 npool_to_pgrainn(p,:) = npool_to_grainn(p) * this%afruitn(p,:)
		end if	
			
         end if

         if (downreg_opt .eqv. .false. .AND. CN_partition_opt == 0) then

            ! N transfer depends on supply and demand
            npool_to_frootn_demand(p)         = (nlc * f1 / cnfr) * fcur
            npool_to_frootn_supply(p)         = npool(p)/dt * fcur
            npool_to_frootn(p) = max(min(npool_to_frootn_supply(p),npool_to_frootn_demand(p)),0.0_r8)

            npool_to_frootn_storage_demand(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
            npool_to_frootn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur)
            npool_to_frootn_storage(p) = max(min(npool_to_frootn_storage_supply(p),npool_to_frootn_storage_demand(p)),0.0_r8)

            npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
            npool_to_leafn_supply(p)   = npool(p)/dt * fcur - npool_to_frootn(p)
            npool_to_leafn(p)   = max(min(npool_to_leafn_supply(p),npool_to_leafn_demand(p)),0.0_r8)

            npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
            npool_to_leafn_storage_supply(p)  = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p)
            npool_to_leafn_storage(p)  = max(min(npool_to_leafn_storage_supply(p),npool_to_leafn_storage_demand(p)),0.0_r8)

            if (CN_residual_opt == 1) then
               npool_to_leafn(p)   = max(npool_to_leafn_supply(p),0.0_r8)
               npool_to_leafn_storage(p)  = max(npool_to_leafn_storage_supply(p),0.0_r8)
            end if

            if (woody(ivt(p)) == 1._r8) then
               npool_to_livestemn_demand(p) = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)
               npool_to_livestemn(p) = max(min(npool_to_livestemn_supply(p),npool_to_livestemn_demand(p)),0.0_r8)

               npool_to_livestemn_storage_demand(p) = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_livestemn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p)
               npool_to_livestemn_storage(p) = max(min(npool_to_livestemn_storage_supply(p), &
                    npool_to_livestemn_storage_demand(p)),0.0_r8)

               npool_to_livecrootn_demand(p) = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p)
               npool_to_livecrootn(p) = max(min(npool_to_livecrootn_supply(p),npool_to_livecrootn_demand(p)),0.0_r8)

               npool_to_livecrootn_storage_demand(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_livecrootn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p)
               npool_to_livecrootn_storage(p) = max(min(npool_to_livecrootn_storage_supply(p), &
                    npool_to_livecrootn_storage_demand(p)),0.0_r8)

               npool_to_deadstemn_demand(p) = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p)
               npool_to_deadstemn(p) = max(min(npool_to_deadstemn_supply(p),npool_to_deadstemn_demand(p)),0.0_r8)

               npool_to_deadstemn_storage_demand(p) = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_deadstemn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)
               npool_to_deadstemn_storage(p) = max(min(npool_to_deadstemn_storage_supply(p), &
                    npool_to_deadstemn_storage_demand(p)),0.0_r8)

               npool_to_deadcrootn_demand(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p)
               npool_to_deadcrootn(p) = max(min(npool_to_deadcrootn_supply(p),npool_to_deadcrootn_demand(p)),0.0_r8)

               npool_to_deadcrootn_storage_demand(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_deadcrootn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  - npool_to_deadstemn_storage(p)
               npool_to_deadcrootn_storage(p) = max(min(npool_to_deadcrootn_storage_supply(p), &
                    npool_to_deadcrootn_storage_demand(p)),0.0_r8)

               npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
               npool_to_leafn_supply(p)   = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p) - npool_to_deadcrootn(p)
               npool_to_leafn(p)   = max(min(npool_to_leafn_supply(p),npool_to_leafn_demand(p)),0.0_r8)

               npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
               npool_to_leafn_storage_supply(p)  = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  - npool_to_deadstemn_storage(p) - &
                    npool_to_deadcrootn_storage(p)
               npool_to_leafn_storage(p)  = max(min(npool_to_leafn_storage_supply(p),&
                    npool_to_leafn_storage_demand(p)),0.0_r8)

               if (CN_residual_opt == 1) then
                  npool_to_leafn(p)   = max(npool_to_leafn_supply(p),0.0_r8)
                  npool_to_leafn_storage(p)  = max(npool_to_leafn_storage_supply(p),0.0_r8)
               end if

            end if

            if (ivt(p) >= npcropmin) then     ! skip 2 generic crops
               cng = graincn(ivt(p))
               npool_to_livestemn_demand(p)          = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)
               npool_to_livestemn(p) = max(min(npool_to_livestemn_supply(p),npool_to_livestemn_demand(p)),0.0_r8)

               npool_to_livestemn_storage_demand(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_livestemn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p)
               npool_to_livestemn_storage(p) = max(min(npool_to_livestemn_storage_supply(p), &
                    npool_to_livestemn_storage_demand(p)),0.0_r8)

               npool_to_livecrootn_demand(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p)
               npool_to_livecrootn(p) = max(min(npool_to_livecrootn_supply(p),npool_to_livecrootn_demand(p)),0.0_r8)

               npool_to_livecrootn_storage_demand(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_livecrootn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p)
               npool_to_livecrootn_storage(p) = max(min(npool_to_livecrootn_storage_supply(p), &
                    npool_to_livecrootn_storage_demand(p)),0.0_r8)

               npool_to_deadstemn_demand(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p)
               npool_to_deadstemn(p) = max(min(npool_to_deadstemn_supply(p), npool_to_deadstemn_demand(p)), 0.0_r8)

               npool_to_deadstemn_storage_demand(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_deadstemn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)
               npool_to_deadstemn_storage(p) = max(min(npool_to_deadstemn_storage_supply(p), &
                    npool_to_deadstemn_storage_demand(p)),0.0_r8)

               npool_to_deadcrootn_demand(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p)
               npool_to_deadcrootn(p) = max(min(npool_to_deadcrootn_supply(p), npool_to_deadcrootn_demand(p)), 0.0_r8)

               npool_to_deadcrootn_storage_demand(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_deadcrootn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  - npool_to_deadstemn_storage(p)
               npool_to_deadcrootn_storage(p) = max(min(npool_to_deadcrootn_storage_supply(p), &
                    npool_to_deadcrootn_storage_demand(p)),0.0_r8)

               npool_to_grainn_demand(p)             = (nlc * f5 / cng) * fcur
               npool_to_grainn_supply(p) = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p) - npool_to_deadcrootn(p)
               npool_to_grainn(p) = max(min(npool_to_grainn_supply(p), npool_to_grainn_demand(p)), 0.0_r8)

               npool_to_grainn_storage_demand(p)     = (nlc * f5 / cng) * (1._r8 -fcur)
               npool_to_grainn_storage_supply(p) = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  - npool_to_deadstemn_storage(p) - &
                    npool_to_deadcrootn_storage(p)
               npool_to_grainn_storage(p) = max(min(npool_to_grainn_storage_supply(p), npool_to_grainn_storage_demand(p)), &
                    0.0_r8)

               npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
               npool_to_leafn_supply(p)   = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p) - npool_to_deadcrootn(p) - npool_to_grainn(p)
               npool_to_leafn(p)   = max(min(npool_to_leafn_supply(p), npool_to_leafn_demand(p)), 0.0_r8)

               npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
               npool_to_leafn_storage_supply(p)  = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  &
                    - npool_to_deadstemn_storage(p) - npool_to_deadcrootn_storage(p) - npool_to_grainn_storage(p)
               npool_to_leafn_storage(p)  = max(min(npool_to_leafn_storage_supply(p), npool_to_leafn_storage_demand(p)), &
                    0.0_r8)

               if (CN_residual_opt == 1) then
                  npool_to_leafn(p)   = max(npool_to_leafn_supply(p),0.0_r8)
                  npool_to_leafn_storage(p)  = max(npool_to_leafn_storage_supply(p),0.0_r8)
               end if

            end if

	    if (phytomer(ivt(p)) > 0) then ! phytomer-based allocation (Y.Fan)
	       cng = graincn(ivt(p))
	       !calculate allocation for froot, live/deadstem, live/deadcroot, grain, and then leaf at the end
               npool_to_frootn_demand(p)         = (nlc * f1 / cnfr) 
               npool_to_frootn_supply(p)         = npool(p)/dt
               npool_to_frootn(p) = max(min(npool_to_frootn_supply(p),npool_to_frootn_demand(p)),0.0_r8)
               npool_to_frootn_storage(p) = 0.0_r8  !only phytomer leaf has non-zero storage term 
	       npool_to_frootn_storage_demand(p) = 0.0_r8
	       npool_to_frootn_storage_supply(p) = 0.0_r8
   
               npool_to_livestemn_demand(p)          = (nlc * f3 * f4 / cnlw) 
               npool_to_livestemn_supply(p) = npool(p)/dt - npool_to_frootn(p)
               npool_to_livestemn(p) = max(min(npool_to_livestemn_supply(p),npool_to_livestemn_demand(p)),0.0_r8)
               npool_to_livestemn_storage(p) = 0.0_r8
	       npool_to_livestemn_storage_demand(p)  = 0.0_r8
	       npool_to_livestemn_storage_supply(p) = 0.0_r8

               npool_to_livecrootn_demand(p)         = (nlc * f2 * f3 * f4 / cnlw) 
               npool_to_livecrootn_supply(p) = npool(p)/dt - npool_to_frootn(p)  -  npool_to_livestemn(p)
               npool_to_livecrootn(p) = max(min(npool_to_livecrootn_supply(p),npool_to_livecrootn_demand(p)),0.0_r8)
               npool_to_livecrootn_storage(p) = 0.0_r8
               npool_to_livecrootn_storage_demand(p) = 0.0_r8
	       npool_to_livecrootn_storage_supply(p) = 0.0_r8
	   
               npool_to_deadstemn_demand(p)          = (nlc * f3 * (1._r8 - f4) / cndw)
               npool_to_deadstemn_supply(p) = npool(p)/dt - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p)
               npool_to_deadstemn(p) = max(min(npool_to_deadstemn_supply(p), npool_to_deadstemn_demand(p)), 0.0_r8)
               npool_to_deadstemn_storage(p) = 0.0_r8
	       npool_to_deadstemn_storage_demand(p)  = 0.0_r8
	       npool_to_deadstemn_storage_supply(p)  = 0.0_r8

               npool_to_deadcrootn_demand(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) 
               npool_to_deadcrootn_supply(p) = npool(p)/dt - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p)
               npool_to_deadcrootn(p) = max(min(npool_to_deadcrootn_supply(p), npool_to_deadcrootn_demand(p)), 0.0_r8)
               npool_to_deadcrootn_storage(p) = 0.0_r8
               npool_to_deadcrootn_storage_demand(p) = 0.0_r8
	       npool_to_deadcrootn_storage_supply(p) = 0.0_r8
			   
               npool_to_grainn_demand(p)             = (nlc * f5 / cng)
               npool_to_grainn_supply(p) = npool(p)/dt - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p) - npool_to_deadcrootn(p)
               npool_to_grainn(p) = max(min(npool_to_grainn_supply(p), npool_to_grainn_demand(p)), 0.0_r8)
               npool_to_grainn_storage(p) = 0.0_r8
               npool_to_grainn_storage_demand(p)     = 0.0_r8
	       npool_to_grainn_storage_supply(p)     = 0.0_r8
	   
	       !only use fcur for storage growth of spear leaf of oil palm 
               npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
               npool_to_leafn_supply(p)   = npool(p)/dt * fcur - npool_to_frootn(p)  -  npool_to_livestemn(p) - &
                    npool_to_livecrootn(p) -  npool_to_deadstemn(p) - npool_to_deadcrootn(p) - npool_to_grainn(p)
               npool_to_leafn(p)   = max(min(npool_to_leafn_supply(p), npool_to_leafn_demand(p)), 0.0_r8)

               npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
               npool_to_leafn_storage_supply(p)  = npool(p)/dt * (1._r8 - fcur) - npool_to_frootn_storage(p) - &
                    npool_to_livestemn_storage(p) - npool_to_livecrootn_storage(p)  &
                    - npool_to_deadstemn_storage(p) - npool_to_deadcrootn_storage(p) - npool_to_grainn_storage(p)
               npool_to_leafn_storage(p)  = max(min(npool_to_leafn_storage_supply(p), npool_to_leafn_storage_demand(p)), &
                    0.0_r8)

               if (CN_residual_opt == 1) then
                  npool_to_leafn(p)   = max(npool_to_leafn_supply(p),0.0_r8)
                  npool_to_leafn_storage(p)  = max(npool_to_leafn_storage_supply(p),0.0_r8)
               end if

              !storage allocation to spear leaf phytomers
	       where (hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:))
		    npool_to_pleafn_storage(p,:) = npool_to_leafn_storage(p) * this%aleafn(p,:)
	       elsewhere (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:))
		    npool_to_pleafn(p,:) = npool_to_leafn(p) * this%aleafn(p,:)
	       endwhere
	       npool_to_pgrainn(p,:) = npool_to_grainn(p) * this%afruitn(p,:)
	    end if	
			
         end if
         

         ! Calculate the amount of carbon that needs to go into growth
         ! respiration storage to satisfy all of the storage growth demands.
         ! Allows for the fraction of growth respiration that is released at the
         ! time of fixation, versus the remaining fraction that is stored for
         ! release at the time of display. Note that all the growth respiration
         ! fluxes that get released on a given timestep are calculated in growth_resp(),
         ! but that the storage of C for growth resp during display of transferred
         ! growth is assigned here.

         gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p)
         if (woody(ivt(p)) == 1._r8) then
            gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
            gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)

            gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
            gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
         end if
         if (ivt(p) >= npcropmin) then     ! skip 2 generic crops
            gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
            gresp_storage = gresp_storage + cpool_to_grainc_storage(p)
         end if

         !oil palm is woody crop, above two clauses will double count livestemc_storage
         !have to rewrite the total storage of all vegetation parts for oil palm
         if (phytomer(ivt(p)) > 0) then
            gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p) &
                                    + cpool_to_livestemc_storage(p) &
                                    + cpool_to_deadstemc_storage(p) &
                                    + cpool_to_livecrootc_storage(p) &
                                    + cpool_to_deadcrootc_storage(p) &
                                    + cpool_to_grainc_storage(p) 
         end if

         cpool_to_gresp_storage(p) = gresp_storage * g1 * (1._r8 - g2)


         ! computing 1.) fractional N demand and 2.) N allocation after uptake for different plant parts
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (downreg_opt .eqv. .false. .AND. CN_partition_opt == 1) then

            ! computing nitrogen demand for different pools based on carbon allocated and CN ratio
            npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
            npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
            npool_to_frootn_demand(p)         = (nlc * f1 / cnfr) * fcur
            npool_to_frootn_storage_demand(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
            if (woody(ivt(p)) == 1._r8) then

               npool_to_livestemn_demand(p)          = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_storage_demand(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadstemn_demand(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_storage_demand(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_livecrootn_demand(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_storage_demand(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadcrootn_demand(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_storage_demand(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops

               cng = graincn(ivt(p))
               npool_to_livestemn_demand(p)          = (nlc * f3 * f4 / cnlw) * fcur
               npool_to_livestemn_storage_demand(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadstemn_demand(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadstemn_storage_demand(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_livecrootn_demand(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
               npool_to_livecrootn_storage_demand(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
               npool_to_deadcrootn_demand(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
               npool_to_deadcrootn_storage_demand(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
               npool_to_grainn_demand(p)             = (nlc * f5 / cng) * fcur
               npool_to_grainn_storage_demand(p)     = (nlc * f5 / cng) * (1._r8 -fcur)
            end if
	    if (phytomer(ivt(p)) > 0) then ! phytomer-based allocation (Y.Fan)
	       cng = graincn(ivt(p))
               npool_to_leafn_demand(p)          = (nlc / cnl) * fcur
               npool_to_leafn_storage_demand(p)  = (nlc / cnl) * (1._r8 - fcur)
               npool_to_frootn_demand(p)         = (nlc * f1 / cnfr) 
               npool_to_frootn_storage_demand(p) = 0.0_r8 !only phytomer leaf has non-zero storage term 
		   
               npool_to_livestemn_demand(p)          = (nlc * f3 * f4 / cnlw) 
               npool_to_livestemn_storage_demand(p) = 0.0_r8			   
               npool_to_livecrootn_demand(p)         = (nlc * f2 * f3 * f4 / cnlw) 
               npool_to_livecrootn_storage_demand(p) = 0.0_r8			   
               npool_to_deadstemn_demand(p)          = (nlc * f3 * (1._r8 - f4) / cndw)
               npool_to_deadstemn_storage_demand(p) = 0.0_r8			   
               npool_to_deadcrootn_demand(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) 
               npool_to_deadcrootn_storage_demand(p) = 0.0_r8			   
               npool_to_grainn_demand(p)             = (nlc * f5 / cng)
               npool_to_grainn_storage_demand(p) = 0.0_r8			   
	    end if

            ! computing 1.) fractional N demand for different plant parts
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            total_plant_Ndemand(p) = npool_to_leafn_demand(p) + npool_to_leafn_storage_demand(p) + &
                 npool_to_frootn_demand(p) + npool_to_frootn_storage_demand(p)

            if (woody(ivt(p)) == 1._r8) then

               total_plant_Ndemand(p) = npool_to_leafn_demand(p) + npool_to_leafn_storage_demand(p) + &
                    npool_to_frootn_demand(p) + npool_to_frootn_storage_demand(p) + &
                    npool_to_livestemn_demand(p) + npool_to_livestemn_storage_demand(p) + npool_to_deadstemn_demand(p) + &
                    npool_to_deadstemn_storage_demand(p)  + &
                    npool_to_livecrootn_demand(p) + npool_to_livecrootn_storage_demand(p) + npool_to_deadcrootn_demand(p) + &
                    npool_to_deadcrootn_storage_demand(p)

            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops

               total_plant_Ndemand(p) = npool_to_leafn_demand(p) + npool_to_leafn_storage_demand(p) + &
                    npool_to_frootn_demand(p) + npool_to_frootn_storage_demand(p) + &
                    npool_to_livestemn_demand(p) + npool_to_livestemn_storage_demand(p) + npool_to_deadstemn_demand(p) + &
                    npool_to_deadstemn_storage_demand(p)  + &
                    npool_to_livecrootn_demand(p) + npool_to_livecrootn_storage_demand(p) + npool_to_deadcrootn_demand(p) + &
                    npool_to_deadcrootn_storage_demand(p) + &
                 npool_to_grainn_demand(p) + npool_to_grainn_storage_demand(p)

            end if
            !!!!!!!!!!!!!!!!!!!!
            !if (phytomer(ivt(p)) > 0) then ! oil palm is woody crop, total demand is calculated like above crops
            !below crop clauses include all vegetative parts: leaf, froot, live/deadstem, live/deadcroot, grain
            !no need to rewrite these items for oil palm
            !!!!!!!!!!!!!!!!!!!

            if (total_plant_Ndemand(p) == 0.0_r8) then    ! removing division by zero

                frNdemand_npool_to_leafn(p) = 0.0_r8
                frNdemand_npool_to_leafn_storage(p) = 0.0_r8
                frNdemand_npool_to_frootn(p) = 0.0_r8
                frNdemand_npool_to_frootn_storage(p) = 0.0_r8
                if (woody(ivt(p)) == 1._r8) then

                   frNdemand_npool_to_livestemn(p) = 0.0_r8
                   frNdemand_npool_to_livestemn_storage(p)  = 0.0_r8
                   frNdemand_npool_to_deadstemn(p) = 0.0_r8
                   frNdemand_npool_to_deadstemn_storage(p) = 0.0_r8
                   frNdemand_npool_to_livecrootn(p) = 0.0_r8
                   frNdemand_npool_to_livecrootn_storage(p) = 0.0_r8
                   frNdemand_npool_to_deadcrootn(p) = 0.0_r8
                   frNdemand_npool_to_deadcrootn_storage(p) = 0.0_r8
                end if
                if (ivt(p) >= npcropmin) then ! skip 2 generic crops

                   frNdemand_npool_to_livestemn(p) = 0.0_r8
                   frNdemand_npool_to_livestemn_storage(p) = 0.0_r8
                   frNdemand_npool_to_deadstemn(p) = 0.0_r8
                   frNdemand_npool_to_deadstemn_storage(p) = 0.0_r8
                   frNdemand_npool_to_livecrootn(p) = 0.0_r8
                   frNdemand_npool_to_livecrootn_storage(p) = 0.0_r8
                   frNdemand_npool_to_deadcrootn(p) = 0.0_r8
                   frNdemand_npool_to_deadcrootn_storage(p) = 0.0_r8
                   frNdemand_npool_to_grainn(p) = 0.0_r8
                   frNdemand_npool_to_grainn_storage(p) = 0.0_r8
                end if

            else

               frNdemand_npool_to_leafn(p) = npool_to_leafn_demand(p) / total_plant_Ndemand(p)
               frNdemand_npool_to_leafn_storage(p) = npool_to_leafn_storage_demand(p) / total_plant_Ndemand(p)
               frNdemand_npool_to_frootn(p) = npool_to_frootn_demand(p) / total_plant_Ndemand(p)
               frNdemand_npool_to_frootn_storage(p) = npool_to_frootn_storage_demand(p) / total_plant_Ndemand(p)
               if (woody(ivt(p)) == 1._r8) then

                  frNdemand_npool_to_livestemn(p) = npool_to_livestemn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livestemn_storage(p)  = npool_to_livestemn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadstemn(p) = npool_to_deadstemn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadstemn_storage(p) = npool_to_deadstemn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livecrootn(p) = npool_to_livecrootn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livecrootn_storage(p) = npool_to_livecrootn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadcrootn(p) = npool_to_deadcrootn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadcrootn_storage(p) = npool_to_deadcrootn_storage_demand(p) / total_plant_Ndemand(p)
               end if
               if (ivt(p) >= npcropmin) then ! skip 2 generic crops

                  frNdemand_npool_to_livestemn(p) = npool_to_livestemn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livestemn_storage(p) = npool_to_livestemn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadstemn(p) = npool_to_deadstemn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadstemn_storage(p) = npool_to_deadstemn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livecrootn(p) = npool_to_livecrootn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_livecrootn_storage(p) = npool_to_livecrootn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadcrootn(p) = npool_to_deadcrootn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_deadcrootn_storage(p) = npool_to_deadcrootn_storage_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_grainn(p) = npool_to_grainn_demand(p) / total_plant_Ndemand(p)
                  frNdemand_npool_to_grainn_storage(p) = npool_to_grainn_storage_demand(p) / total_plant_Ndemand(p)
               end if

            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! computing N allocation for different plant parts
            ! allocating allocation to different plant parts in proportion to the fractional demand
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            npool_to_leafn(p) = frNdemand_npool_to_leafn(p) * npool(p) / dt
            npool_to_leafn_storage(p)  = frNdemand_npool_to_leafn_storage(p) * npool(p) / dt
            npool_to_frootn(p) = frNdemand_npool_to_frootn(p) * npool(p) / dt
            npool_to_frootn_storage(p) = frNdemand_npool_to_frootn_storage(p) * npool(p) / dt
            if (woody(ivt(p)) == 1._r8) then
               npool_to_livestemn(p) = frNdemand_npool_to_livestemn(p) * npool(p) / dt
               npool_to_livestemn_storage(p) = frNdemand_npool_to_livestemn_storage(p) * npool(p) / dt
               npool_to_deadstemn(p) = frNdemand_npool_to_deadstemn(p) * npool(p) / dt
               npool_to_deadstemn_storage(p) = frNdemand_npool_to_deadstemn_storage(p) * npool(p) / dt
               npool_to_livecrootn(p) = frNdemand_npool_to_livecrootn(p) * npool(p) / dt
               npool_to_livecrootn_storage(p) = frNdemand_npool_to_livecrootn_storage(p) * npool(p) / dt
               npool_to_deadcrootn(p) = frNdemand_npool_to_deadcrootn(p) * npool(p) / dt
               npool_to_deadcrootn_storage(p) = frNdemand_npool_to_deadcrootn_storage(p) * npool(p) / dt
            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops
               npool_to_livestemn(p) = frNdemand_npool_to_livestemn(p) * npool(p) / dt
               npool_to_livestemn_storage(p) = frNdemand_npool_to_livestemn_storage(p) * npool(p) / dt
               npool_to_deadstemn(p) = frNdemand_npool_to_deadstemn(p) * npool(p) / dt
               npool_to_deadstemn_storage(p) = frNdemand_npool_to_deadstemn_storage(p) * npool(p) / dt
               npool_to_livecrootn(p) = frNdemand_npool_to_livecrootn(p) * npool(p) / dt
               npool_to_livecrootn_storage(p) = frNdemand_npool_to_livecrootn_storage(p) * npool(p) / dt
               npool_to_deadcrootn(p) = frNdemand_npool_to_deadcrootn(p) * npool(p) / dt
               npool_to_deadcrootn_storage(p) = frNdemand_npool_to_deadcrootn_storage(p) * npool(p) / dt
               npool_to_grainn(p) = frNdemand_npool_to_grainn(p) * npool(p) / dt
               npool_to_grainn_storage(p) = frNdemand_npool_to_grainn_storage(p) * npool(p) / dt
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

            cpool_to_resp(p) = 0.0_r8
            cpool_to_leafc_resp(p) = 0.0_r8
            cpool_to_leafc_storage_resp(p) = 0.0_r8
            cpool_to_frootc_resp(p) = 0.0_r8
            cpool_to_frootc_storage_resp(p) = 0.0_r8
            cpool_to_livecrootc_resp(p) = 0.0_r8
            cpool_to_livecrootc_storage_resp(p) = 0.0_r8
            cpool_to_livestemc_resp(p) = 0.0_r8
            cpool_to_livestemc_storage_resp(p) = 0.0_r8

            if ( laisun(p)+laisha(p) > 0.0_r8 ) then
               if (cnveg_nitrogenstate_inst%leafn_storage_patch(p) == 0.0_r8 ) then
                  ! to avoid division by zero, and also to make actual_leafncn(p) a very large number if leafn(p) is zero 
                  this%actual_storage_leafcn(p) = spval
               else                        	
                  ! leaf CN ratio 
                  this%actual_storage_leafcn(p) = cnveg_carbonstate_inst%leafc_storage_patch(p)  &
                  / cnveg_nitrogenstate_inst%leafn_storage_patch(p)
               end if
            end if
            
            if (carbon_resp_opt == 1 .AND. laisun(p)+laisha(p) > 0.0_r8) then   
               ! computing carbon to nitrogen ratio of different plant parts 

          
               if (cnveg_nitrogenstate_inst%frootn_storage_patch(p) == 0.0_r8) then
                  ! to avoid division by zero, and also to make frootcn_actual(p) a very large number if frootc(p) is zero 
                  frootcn_actual = cnveg_carbonstate_inst%frootc_storage_patch(p) / n_min
               else
                  ! fine root CN ratio 
                  frootcn_actual = cnveg_carbonstate_inst%frootc_storage_patch(p) / cnveg_nitrogenstate_inst%frootn_storage_patch(p)
               end if
      
               if (woody(ivt(p)) == 1._r8) then  

                  if (cnveg_nitrogenstate_inst%livestemn_storage_patch(p) == 0.0_r8) then
                     ! to avoid division by zero, and also to make livestemcn_actual(p) a very large number if livestemc(p) is zero 
                     livestemcn_actual = cnveg_carbonstate_inst%livestemc_storage_patch(p) / n_min
                  else                        										
                     ! live stem CN ratio 
                     livestemcn_actual =  cnveg_carbonstate_inst%livestemc_storage_patch(p) / &
                          cnveg_nitrogenstate_inst%livestemn_storage_patch(p)
                  end if 														 														
              
                  if (cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) == 0.0_r8) then
                     ! to avoid division by zero, and also to make livecrootcn_actual(p) a very large number if livecrootc(p) is zero 
                     livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / n_min
                  else
                     ! live coarse root CN ratio
                     livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livecrootn_storage_patch(p)
                  end if              
               end if   
          
               if (ivt(p) >= npcropmin) then ! skip 2 generic crops    
														              
                  if (cnveg_nitrogenstate_inst%livestemn_storage_patch(p) == 0.0_r8) then
                     ! to avoid division by zero, and also to make livestemcn_actual(p) a very large number if livestemc(p) is zero 
                     livestemcn_actual = cnveg_carbonstate_inst%livestemc_storage_patch(p) / n_min
                  else                        										
                     ! live stem CN ratio 
                     livestemcn_actual =  cnveg_carbonstate_inst%livestemc_storage_patch(p) / &
                          cnveg_nitrogenstate_inst%livestemn_storage_patch(p)
                  end if
              
                  if (cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) == 0.0_r8) then
                     ! to avoid division by zero, and also to make livecrootcn_actual(p) a very large number if livecrootc(p) is zero 
                     livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / n_min
                  else
                     ! live coarse root CN ratio
                     livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livecrootn_storage_patch(p)
                  end if          
               end if
                          
               leafcn_max = leafcn(ivt(p)) + 15.0_r8
               frootcn_max = frootcn(ivt(p)) + 15.0_r8
            
               ! Note that for high CN ratio stress the plant part does not retranslocate nitrogen as the plant part will need the N 
               ! if high leaf CN ratio (i.e., high leaf C compared to N) then turnover extra C           
               if (this%actual_storage_leafcn(p) > leafcn_max) then   
               
                  frac_resp =  (this%actual_storage_leafcn(p) - leafcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))
               
                  cpool_to_leafc_resp(p)          = frac_resp * cpool_to_leafc(p) 
                  cpool_to_leafc_storage_resp(p)  = frac_resp * cpool_to_leafc_storage(p)
                  
                  !cpool_to_leafc(p) = cpool_to_leafc(p) - cpool_to_leafc_resp(p) 
                  !cpool_to_leafc_storage(p) = cpool_to_leafc_storage(p) - cpool_to_leafc_storage_resp(p)
                  
               end if    

               ! if high fine root CN ratio (i.e., high fine root C compared to N) then turnover extra C 
               if (frootcn_actual > frootcn_max) then     
               
                  frac_resp =  (frootcn_actual - frootcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))
                
                  cpool_to_frootc_resp(p)         = frac_resp * cpool_to_frootc(p) 
                  cpool_to_frootc_storage_resp(p) = frac_resp * cpool_to_frootc_storage(p)
                  
                  !cpool_to_frootc(p) = cpool_to_frootc(p) - cpool_to_frootc_resp(p)    
                  !cpool_to_frootc_storage(p) = cpool_to_frootc_storage(p) - cpool_to_frootc_storage_resp(p)
                 
               end if    
					
		       if (woody(ivt(p)) == 1._r8) then  
		  
		          livewdcn_max = livewdcn(ivt(p)) + 15.0_r8

                  ! if high coarse root CN ratio (i.e., high coarse root C compared to N) then turnover extra C 
                  if (livecrootcn_actual > livewdcn_max) then   
                  
                     frac_resp =  (livecrootcn_actual - livewdcn_max) / 10.0_r8
                     frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))  
                  
                     cpool_to_livecrootc_resp(p)         = frac_resp * cpool_to_livecrootc(p) 
                     cpool_to_livecrootc_storage_resp(p) = frac_resp * cpool_to_livecrootc_storage(p)
                     
                     !cpool_to_livecrootc(p) = cpool_to_livecrootc(p) - cpool_to_livecrootc_resp(p)     
                     !cpool_to_livecrootc_storage(p) = cpool_to_livecrootc_storage(p) - cpool_to_livecrootc_storage_resp(p) 
                  
                  end if    
          
                  ! if high stem CN ratio (i.e., high stem C compared to N) then turnover extra C 
                  if (livestemcn_actual > livewdcn_max) then  
                  
                     frac_resp =  (livestemcn_actual - livewdcn_max) / 10.0_r8
                     frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))     
                   
                     cpool_to_livestemc_resp(p)          = frac_resp * cpool_to_livestemc(p) 
                     cpool_to_livestemc_storage_resp(p)  = frac_resp * cpool_to_livestemc_storage(p)
                     
                     !cpool_to_livestemc(p) = cpool_to_livestemc(p) - cpool_to_livestemc_resp(p)    
                     !cpool_to_livestemc_storage(p) = cpool_to_livestemc_storage(p) - cpool_to_livestemc_storage_resp(p)
                     
                  end if    
              
               end if  
          
               if (ivt(p) >= npcropmin) then ! skip 2 generic crops    

                  livewdcn_max = livewdcn(ivt(p)) + 15.0_r8

                  ! if high coarse root CN ratio (i.e., high coarse root C compared to N) then turnover extra C 
                  if (livecrootcn_actual > livewdcn_max) then     
                  
                     frac_resp =  (livecrootcn_actual - livewdcn_max) / 10.0_r8
                     frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))  
                  
                     cpool_to_livecrootc_resp(p)         = frac_resp * cpool_to_livecrootc(p) 
                     cpool_to_livecrootc_storage_resp(p) = frac_resp * cpool_to_livecrootc_storage(p)
                     
                     !cpool_to_livecrootc(p) = cpool_to_livecrootc(p) - cpool_to_livecrootc_resp(p)     
                     !cpool_to_livecrootc_storage(p) = cpool_to_livecrootc_storage(p) - cpool_to_livecrootc_storage_resp(p) 
                     
                  end if    
          
                  ! if high stem CN ratio (i.e., high stem C compared to N) then turnover extra C 
                  if (livestemcn_actual > livewdcn_max) then       

                     frac_resp =  (livestemcn_actual - livewdcn_max) / 10.0_r8
                     frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))     
                   
                     cpool_to_livestemc_resp(p)          = frac_resp * cpool_to_livestemc(p) 
                     cpool_to_livestemc_storage_resp(p)  = frac_resp * cpool_to_livestemc_storage(p)
                     
                     !cpool_to_livestemc(p) = cpool_to_livestemc(p) - cpool_to_livestemc_resp(p)    
                     !cpool_to_livestemc_storage(p) = cpool_to_livestemc_storage(p) - cpool_to_livestemc_storage_resp(p)

                  end if    
              
               end if  
               
               cpool_to_resp(p) = cpool_to_leafc_resp(p) + cpool_to_leafc_storage_resp(p)  + cpool_to_frootc_resp(p) + &
                    cpool_to_frootc_storage_resp(p) + cpool_to_livecrootc_resp(p) + cpool_to_livecrootc_storage_resp(p) + &
                    cpool_to_livestemc_resp(p) + cpool_to_livestemc_storage_resp(p)

            end if   ! end of if (carbon_resp_opt == 1 .AND. laisun(p)+laisha(p) > 0.0_r8) then   

            !if (cnveg_nitrogenstate_inst%leafn_storage_patch(p) < n_min .or. laisun(p)+laisha(p) <= 0.0_r8) then
               !! to make output on history missing value
               !this%actual_storage_leafcn(p) = spval
            !end if
									
         end if    ! end of if (downreg_opt .eqv. .false. .AND. CN_partition_opt == 1) then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end do ! end patch loop

    end associate

  end subroutine calc_plant_cn_alloc

! -----------------------------------------------------------------------
  subroutine calc_plant_nutrient_demand(this, bounds,  num_soilp, filter_soilp,&
       photosyns_inst, crop_inst, canopystate_inst,                            &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,        &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,                   &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst, &
       aroot, arepr)
    !
    ! !USES:
    use CanopyStateType        , only : canopystate_type
    use PhotosynthesisMod      , only : photosyns_type
    use CropType               , only : crop_type
    use CNVegStateType         , only : cnveg_state_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
    use SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use EnergyFluxType         , only : energyflux_type     !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)            , intent(in)    :: photosyns_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)   , intent(in) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(energyflux_type)           , intent(in)    :: energyflux_inst
    real(r8)                        , intent(out)   :: aroot(bounds%begp:)
    real(r8)                        , intent(out)   :: arepr(bounds%begp:)
    !-----------------------------------------------------------------------

    call this%calc_plant_nitrogen_demand(bounds,  num_soilp, filter_soilp, &
       photosyns_inst, crop_inst, canopystate_inst,                        &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,    &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,               &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                  &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst, &
       aroot=aroot(bounds%begp:bounds%endp),                               &
       arepr=arepr(bounds%begp:bounds%endp))

  end subroutine calc_plant_nutrient_demand

  !-----------------------------------------------------------------------
  subroutine calc_plant_nitrogen_demand(this, bounds,  num_soilp, filter_soilp, &
       photosyns_inst, crop_inst, canopystate_inst,                             &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,         &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,                    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst, &
       aroot, arepr)
    !
    ! !USES:
    use pftconMod              , only : npcropmin, pftcon
    use pftconMod              , only : ntmp_soybean, nirrig_tmp_soybean
    use pftconMod              , only : ntrp_soybean, nirrig_trp_soybean
    use pftconMod              , only : noilpalm, nirrig_oilpalm, mxnp
! mxmat, a_par, b_par 
    use clm_varcon             , only : secspday, dzsoi_decomp
    use clm_varctl             , only : use_c13, use_c14
    use clm_varctl             , only : nscalar_opt, plant_ndemand_opt, substrate_term_opt, temp_scalar_opt
    use clm_varpar             , only : nlevdecomp
    use clm_time_manager       , only : get_step_size, get_days_per_year, get_curr_calday
    use CanopyStateType        , only : canopystate_type
    use PhotosynthesisMod      , only : photosyns_type
    use CropType               , only : crop_type
    use CNVegStateType         , only : cnveg_state_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
    use SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use EnergyFluxType         , only : energyflux_type     !
    use CNSharedParamsMod      , only : use_fun
    use CNPrecisionControlMod  , only : n_min
    use clm_varcon             , only : spval
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)            , intent(in)    :: photosyns_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)   , intent(in) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    type(energyflux_type)           , intent(in)    :: energyflux_inst
    real(r8)                        , intent(out)   :: aroot(bounds%begp:)
    real(r8)                        , intent(out)   :: arepr(bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c, p, j                                 ! indices
    integer  :: fp                                         ! lake filter patch index
    real(r8) :: mr                                         ! maintenance respiration (gC/m2/s)
    real(r8) :: f1, f2, f3, f4, g1, g2                     ! allocation parameters
    real(r8) :: cnl, cnfr, cnlw, cndw                      ! C:N ratios for leaf, fine root, and wood
    real(r8) :: curmr, curmr_ratio                         ! xsmrpool temporary variables
    real(r8) :: f5                                         ! grain allocation                   parameter
    real(r8) :: cng                                        ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8) :: fleaf                                      ! fraction allocated to leaf
    real(r8) :: t1                                         ! temporary variable
    real(r8) :: dt                                         ! model time step
    real(r8) :: dayscrecover                               ! number of days to recover negative cpool
    real(r8) :: f_N              (bounds%begp:bounds%endp)
    real(r8) :: Kmin
    real(r8) :: leafcn_max
    real(r8) :: leafcn_min
    real(r8) :: nscalar
    real(r8) :: sminn_total
    real(r8) :: substrate_term
    real(r8) :: temp_scalar
    real(r8) :: Vmax_N
    real(r8) :: allocation_leaf  (bounds%begp:bounds%endp)
    real(r8) :: allocation_stem  (bounds%begp:bounds%endp)
    real(r8) :: allocation_froot (bounds%begp:bounds%endp)
    !new oil palm-related variables
    integer jday                    ! julian day of the year
    !real(r8) astem0                 !initial allocation to stem
    real(r8) pi,n                   !indices
    real(r8) dayspyr                ! days per year
    real(r8) gddperday
    real(r8) :: psum_exp, psum_bud   !total leaf C/N sink size for expanded/unexpanded phytomers
    real(r8) :: rleafn (bounds%begp:bounds%endp,1:mxnp)   !fruit and leaf allocation ratio for each phytomer
    real(r8) :: rfruitn(bounds%begp:bounds%endp,1:mxnp)
	
    ! -----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(aroot) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(arepr) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%actual_leafcn) >= (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%actual_leafcn) <= (/bounds%begp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%aleafn) >= (/bounds%endp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%aleafn) <= (/bounds%begp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(this%afruitn) >= (/bounds%endp,mxnp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((lbound(this%afruitn) <= (/bounds%begp,mxnp/)), errMsg(sourcefile, __LINE__))

    associate(                                                                        &
         ivt                   => patch%itype                                        ,  & ! Input:  [integer  (:) ]  patch vegetation type

         woody                 => pftcon%woody                                     ,  & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf            => pftcon%froot_leaf                                ,  & ! Input:  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem            => pftcon%croot_stem                                ,  & ! Input:  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf             => pftcon%stem_leaf                                 ,  & ! Input:  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd               => pftcon%flivewd                                   ,  & ! Input:  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn                => pftcon%leafcn                                    ,  & ! Input:  leaf C:N (gC/gN)
         frootcn               => pftcon%frootcn                                    , & ! Input:  fine root C:N (gC/gN)
         livewdcn              => pftcon%livewdcn                                   , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn              => pftcon%deadwdcn                                   , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         graincn               => pftcon%graincn                                    , & ! Input:  grain C:N (gC/gN)
         fleafcn               => pftcon%fleafcn                                    , & ! Input:  leaf c:n during organ fill
         ffrootcn              => pftcon%ffrootcn                                   , & ! Input:  froot c:n during organ fill
         fstemcn               => pftcon%fstemcn                                    , & ! Input:  stem c:n during organ fill
         bfact                 => pftcon%bfact                                      , & ! Input:  parameter used below
         aleaff                => pftcon%aleaff                                     , & ! Input:  parameter used below
         arootf                => pftcon%arootf                                     , & ! Input:  parameter used below
         astemf                => pftcon%astemf                                     , & ! Input:  parameter used below
         arooti                => pftcon%arooti                                     , & ! Input:  parameter used below
         fleafi                => pftcon%fleafi                                     , & ! Input:  parameter used below
         allconsl              => pftcon%allconsl                                   , & ! Input:  parameter used below
         allconss              => pftcon%allconss                                   , & ! Input:  parameter used below
         grperc                => pftcon%grperc                                     , & ! Input:  parameter used below
         grpnow                => pftcon%grpnow                                     , & ! Input:  parameter used below
         declfact              => pftcon%declfact                                   , & ! Input:
         season_decid          => pftcon%season_decid                               , & ! Input:  binary flag for seasonal-deciduous leaf habit (0 or 1)
         stress_decid          => pftcon%stress_decid                               , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)
         psnsun                => photosyns_inst%psnsun_patch                       , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         psnsha                => photosyns_inst%psnsha_patch                       , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun            => photosyns_inst%c13_psnsun_patch                   , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsha            => photosyns_inst%c13_psnsha_patch                   , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsun            => photosyns_inst%c14_psnsun_patch                   , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsha            => photosyns_inst%c14_psnsha_patch                   , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

         laisun                => canopystate_inst%laisun_patch                     , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                => canopystate_inst%laisha_patch                     , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index

         hui                   => crop_inst%gddplant_patch                          , & ! Input:  [real(r8) (:)   ]  =gdd since planting (gddplant)
         leafout               => crop_inst%gddtsoi_patch                           , & ! Input:  [real(r8) (:)   ]  =gdd from top soil layer temperature
         croplive              => crop_inst%croplive_patch                          , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         gddmaturity           => cnveg_state_inst%gddmaturity_patch                , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest
         huileaf               => cnveg_state_inst%huileaf_patch                    , & ! Input:  [real(r8) (:)   ]  heat unit index needed from planting to leaf emergence
         huigrain              => cnveg_state_inst%huigrain_patch                   , & ! Input:  [real(r8) (:)   ]  same to reach vegetative maturity
         peaklai               => cnveg_state_inst%peaklai_patch                    , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         aleafi                => cnveg_state_inst%aleafi_patch                     , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         astemi                => cnveg_state_inst%astemi_patch                     , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         aleaf                 => cnveg_state_inst%aleaf_patch                      , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                 => cnveg_state_inst%astem_patch                      , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         grain_flag            => cnveg_state_inst%grain_flag_patch                 , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not
         c_allometry           => cnveg_state_inst%c_allometry_patch                , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry           => cnveg_state_inst%n_allometry_patch                , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         tempsum_potential_gpp => cnveg_state_inst%tempsum_potential_gpp_patch      , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP
         tempmax_retransn      => cnveg_state_inst%tempmax_retransn_patch           , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)
         annsum_potential_gpp  => cnveg_state_inst%annsum_potential_gpp_patch       , & ! Output: [real(r8) (:)   ]  annual sum of potential GPP
         annmax_retransn       => cnveg_state_inst%annmax_retransn_patch            , & ! Output: [real(r8) (:)   ]  annual max of retranslocated N pool

         perennial             => pftcon%perennial                                  , & ! Input:  [integer (:)]   binary flag for perennial crop phenology (1=perennial, 0=not perennial) (added by Y.Fan)
         phytomer              => pftcon%phytomer                                   , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
         leaf_long             => pftcon%leaf_long                                  , & ! Input:  [real(r8) (:)]  leaf longevity (yrs)
         slatop                => pftcon%slatop                                     , & ! Input:  [real(r8) (:)]  specific leaf area at top of canopy, projected area basis [m^2/gC]
         mxmat                 =>    pftcon%mxmat                                   , & ! Input:
         a_par                 =>    pftcon%a_par                                   , & ! Input:
         b_par                 =>    pftcon%b_par                                   , & ! Input:
         gddmaturity2          => crop_inst%gddmaturity2_patch               , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest since previous harvest (Y.Fan)
         huigrain2             => crop_inst%huigrain2_patch                  , & ! Input:  [real(r8) (:)]  gdd needed from last harvest to start of next grainfill (Y.Fan)
         idpp                  => crop_inst%idpp_patch                       , & ! Input:  [integer (:)]  days past planting (Y.Fan)
         idpp2                 => crop_inst%idpp2_patch                      , & ! InOut:  [integer (:)]  Saved idpp from phase2 before grainfill starts

         gdd15                 => temperature_inst%gdd15_patch               , & ! Input:  [real(r8) (:)]  growing deg. days base 15 deg C (ddays) (Y.Fan)
        !gdd1520               => temperature_inst%gdd1520_patch             , & ! Input:  [real(r8) (:)]  20 yr mean of gdd15
         aleaf0                => crop_inst%aleaf0_patch                     , & ! Input:  [real(r8) (:)]  initial leaf allocation coefficient
         livep                 => crop_inst%livep_patch                      , & ! Input:  [logical(r8) (:,:)]  Flag, true if this phytomer is alive
         plaipeak                  => crop_inst%plaipeak_patch                 , & ! Input:  [integer (:,:)]   Flag, 1: max allowed lai per phytomer; 0: not at max
         huileafnp                 => crop_inst%huileafnp_patch                , & ! Input:  [real(r8) (:,:)]  hui needed for initiation of successive phytomers
         huilfexpnp                => crop_inst%huilfexpnp_patch               , & ! Input:  [real(r8) (:,:)]  hui needed for leaf expansion of successive phytomers
         huilfmatnp                => crop_inst%huilfmatnp_patch               , & ! Input:  [real(r8) (:,:)]  hui needed for leaf maturity of successive phytomers
         huigrnnp                  => crop_inst%huigrnnp_patch                 , & ! Input:  [real(r8) (:,:)]  hui needed for start of grainfill of successive phytomers
         grnmatnp                  => crop_inst%grnmatnp_patch                 , & ! Input:  [real(r8) (:,:)]  hui needed for grain maturity of successive phytomers
         pgrainc                   => cnveg_carbonstate_inst%pgrainc_patch     , & ! InOut:  [real(r8) (:,:)]  (gC/m2) phytomer grain C
         pgrainn                   => cnveg_nitrogenstate_inst%pgrainn_patch   , & ! InOut:  [real(r8) (:,:)]  (gN/m2) phytomer grain N
         tlai                      => cnveg_state_inst%tlai_patch              , & ! Input:  [real(r8) (:)] one-sided leaf area index, no burying by snow
         plai                      => crop_inst%plai_patch              , & ! Input:  [real(r8) (:,:)]  one-sided leaf area index of each phytomer

         xsmrpool              => cnveg_carbonstate_inst%xsmrpool_patch             , & ! Input:  [real(r8) (:)   ]  (gC/m2) temporary photosynthate C pool
         leafc                 => cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:)   ]
         frootc                => cnveg_carbonstate_inst%frootc_patch               , & ! Input:  [real(r8) (:)   ]
         livestemc             => cnveg_carbonstate_inst%livestemc_patch            , & ! Input:  [real(r8) (:)   ]
         livecrootc            => cnveg_carbonstate_inst%livecrootc_patch           , & ! Input:  [real(r8) (:)   ]
         retransn              => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N

         monsum_npp            => cnveg_carbonflux_inst%monsum_npp_patch            , & ! InOut:  [real(r8) (:)   ]  monthly sum NPP (gC/m2/yr)
         annsum_npp            => cnveg_carbonflux_inst%annsum_npp_patch            , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         leaf_mr               => cnveg_carbonflux_inst%leaf_mr_patch               , & ! Input:  [real(r8) (:)   ]
         froot_mr              => cnveg_carbonflux_inst%froot_mr_patch              , & ! Input:  [real(r8) (:)   ]
         livestem_mr           => cnveg_carbonflux_inst%livestem_mr_patch           , & ! Input:  [real(r8) (:)   ]
         livecroot_mr          => cnveg_carbonflux_inst%livecroot_mr_patch          , & ! Input:  [real(r8) (:)   ]
         grain_mr              => cnveg_carbonflux_inst%grain_mr_patch              , & ! Input:  [real(r8) (:)   ]
         gpp                   => cnveg_carbonflux_inst%gpp_before_downreg_patch    , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                => cnveg_carbonflux_inst%availc_patch                , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         xsmrpool_recover      => cnveg_carbonflux_inst%xsmrpool_recover_patch      , & ! Output: [real(r8) (:)   ]  C flux assigned to recovery of negative cpool (gC/m2/s)
         psnsun_to_cpool       => cnveg_carbonflux_inst%psnsun_to_cpool_patch       , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool     => cnveg_carbonflux_inst%psnshade_to_cpool_patch     , & ! Output: [real(r8) (:)   ]
         leaf_curmr            => cnveg_carbonflux_inst%leaf_curmr_patch            , & ! Output: [real(r8) (:)   ]
         froot_curmr           => cnveg_carbonflux_inst%froot_curmr_patch           , & ! Output: [real(r8) (:)   ]
         livestem_curmr        => cnveg_carbonflux_inst%livestem_curmr_patch        , & ! Output: [real(r8) (:)   ]
         livecroot_curmr       => cnveg_carbonflux_inst%livecroot_curmr_patch       , & ! Output: [real(r8) (:)   ]
         grain_curmr           => cnveg_carbonflux_inst%grain_curmr_patch           , & ! Output: [real(r8) (:)   ]
         leaf_xsmr             => cnveg_carbonflux_inst%leaf_xsmr_patch             , & ! Output: [real(r8) (:)   ]
         froot_xsmr            => cnveg_carbonflux_inst%froot_xsmr_patch            , & ! Output: [real(r8) (:)   ]
         livestem_xsmr         => cnveg_carbonflux_inst%livestem_xsmr_patch         , & ! Output: [real(r8) (:)   ]
         livecroot_xsmr        => cnveg_carbonflux_inst%livecroot_xsmr_patch        , & ! Output: [real(r8) (:)   ]
         grain_xsmr            => cnveg_carbonflux_inst%grain_xsmr_patch            , & ! Output: [real(r8) (:)   ]
         cpool_to_xsmrpool     => cnveg_carbonflux_inst%cpool_to_xsmrpool_patch     , & ! Output: [real(r8) (:)   ]

         leafn                 => cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         plant_ndemand         => cnveg_nitrogenflux_inst%plant_ndemand_patch       , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         avail_retransn        => cnveg_nitrogenflux_inst%avail_retransn_patch      , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         retransn_to_npool     => cnveg_nitrogenflux_inst%retransn_to_npool_patch   , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         sminn_to_npool        => cnveg_nitrogenflux_inst%sminn_to_npool_patch      , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
         leafn_to_retransn     => cnveg_nitrogenflux_inst%leafn_to_retransn_patch   , & ! Output: [real(r8) (:)   ]
         frootn_to_retransn    => cnveg_nitrogenflux_inst%frootn_to_retransn_patch  , & ! Output: [real(r8) (:)   ]
         livestemn_to_retransn => cnveg_nitrogenflux_inst%livestemn_to_retransn_patch,& ! Output: [real(r8) (:)   ]
         livestemn             => cnveg_nitrogenstate_inst%livestemn_patch          , & ! Input:  [real(r8) (:)   ]  (gN/m2) livestem N
         frootn                => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:)   ]  (gN/m2) fine root N
         sminn_vr              => soilbiogeochem_nitrogenstate_inst%sminn_vr_col    , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
         btran                 => energyflux_inst%btran_patch                       , & ! Input: [real(r8) (:)    ]  transpiration wetness factor (0 to 1)
         t_scalar              => soilbiogeochem_carbonflux_inst%t_scalar_col        & ! Input:  [real(r8) (:,:) ]  soil temperature scalar for decomp

         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      dayspyr = get_days_per_year()
      jday    = get_curr_calday()
      gddperday = max(10._r8, gdd15(p)/jday)

      ! set number of days to recover negative cpool
      dayscrecover = params_inst%dayscrecover     ! loop over patches to assess the total plant N demand
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! get the time step total gross photosynthesis
         ! this is coming from the canopy fluxes code, and is the
         ! gpp that is used to control stomatal conductance.
         ! For the nitrogen downregulation code, this is assumed
         ! to be the potential gpp, and the actual gpp will be
         ! reduced due to N limitation.

         ! Convert psn from umol/m2/s -> gC/m2/s

         ! The input psn (psnsun and psnsha) are expressed per unit LAI
         ! in the sunlit and shaded canopy, respectively. These need to be
         ! scaled by laisun and laisha to get the total gpp for allocation

         ! Note that no associate statement is used for the isotope carbon fluxes below
         ! since they are not always allocated AND nag compiler will complain if you try to
         ! to have an associate statement with unallocated memory

         psnsun_to_cpool(p)   = psnsun(p) * laisun(p) * 12.011e-6_r8
         psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6_r8

         if ( use_c13 ) then
            c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
            c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
         endif

         if ( use_c14 ) then
            c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
            c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
         endif

         gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

         ! get the time step total maintenance respiration
         ! These fluxes should already be in gC/m2/s

         mr = leaf_mr(p) + froot_mr(p)
         if (woody(ivt(p)) == 1.0_r8) then
            mr = mr + livestem_mr(p) + livecroot_mr(p)
	    !include all mr terms for woody crop types (Y.Fan)
	    if (ivt(p) >= npcropmin .and. croplive(p)) mr = mr + grain_mr(p)
         else if (ivt(p) >= npcropmin) then
            if (croplive(p)) mr = mr + livestem_mr(p) + grain_mr(p)
         end if
         !ensure no negative mr. negative mr will cause availc > 0 even when gpp=0
         !this is added in case some plant c/n pools have small negative values after crop die or at crop rotation (Y.Fan 2016)
         mr = max(mr, 0._r8)

         ! carbon flux available for allocation
         availc(p) = gpp(p) - mr

         ! new code added for isotope calculations, 7/1/05, PET
         ! If mr > gpp, then some mr comes from gpp, the rest comes from
         ! cpool (xsmr)
         if (mr > 0._r8 .and. availc(p) < 0._r8) then
            curmr = gpp(p)
            curmr_ratio = curmr / mr
         else
            curmr_ratio = 1._r8
         end if
         leaf_curmr(p)      = leaf_mr(p) * curmr_ratio
         leaf_xsmr(p)       = leaf_mr(p) - leaf_curmr(p)
         froot_curmr(p)     = froot_mr(p) * curmr_ratio
         froot_xsmr(p)      = froot_mr(p) - froot_curmr(p)
         livestem_curmr(p)  = livestem_mr(p) * curmr_ratio
         livestem_xsmr(p)   = livestem_mr(p) - livestem_curmr(p)
         livecroot_curmr(p) = livecroot_mr(p) * curmr_ratio
         livecroot_xsmr(p)  = livecroot_mr(p) - livecroot_curmr(p)
         grain_curmr(p)     = grain_mr(p) * curmr_ratio
         grain_xsmr(p)      = grain_mr(p) - grain_curmr(p)

         ! no allocation when available c is negative
         availc(p) = max(availc(p),0.0_r8)

         ! test for an xsmrpool deficit
         if (xsmrpool(p) < 0.0_r8) then
            ! Running a deficit in the xsmrpool, so the first priority is to let
            ! some availc from this timestep accumulate in xsmrpool.
            ! Determine rate of recovery for xsmrpool deficit

            xsmrpool_recover(p) = -xsmrpool(p)/(dayscrecover*secspday)
            if (xsmrpool_recover(p) < availc(p)) then
               ! available carbon reduced by amount for xsmrpool recovery
               availc(p) = availc(p) - xsmrpool_recover(p)
            else
               ! all of the available carbon goes to xsmrpool recovery
               xsmrpool_recover(p) = availc(p)
               availc(p) = 0.0_r8
            end if
            cpool_to_xsmrpool(p) = xsmrpool_recover(p)
         end if

         f1 = froot_leaf(ivt(p))
         f2 = croot_stem(ivt(p))

         ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
         ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
         ! This variable allocation is only for trees. Shrubs have a constant
         ! allocation as specified in the pft-physiology file.  The value is also used
         ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

         if (stem_leaf(ivt(p)) == -1._r8) then
            f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
         else
            f3 = stem_leaf(ivt(p))
         end if

         f4   = flivewd(ivt(p))
         g1   = grperc(ivt(p))
         g2   = grpnow(ivt(p))
         cnl  = leafcn(ivt(p))
         cnfr = frootcn(ivt(p))
         cnlw = livewdcn(ivt(p))
         cndw = deadwdcn(ivt(p))

         ! calculate f1 to f5 for prog crops following AgroIBIS subr phenocrop

         f5 = 0._r8 ! continued intializations from above

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops

            if (croplive(p)) then
			
	     !!for multilayer phytmoer structure, calculate f1, f3, f5 according to subroutine PalmPhenology
             if (phytomer(ivt(p)) > 0) then
	        !Sub-PFT (phytomer) leaf allocation
	        ! ==================

		!calculate leaf sink size and update leaf alloc ratio for each phytomer
		where (livep(p,:) == .false. .or. plaipeak(p,:) == 1) !adjust lai for each phytomer
			rleafn(p,:) = 0._r8
		elsewhere (hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:)) !pre-expansion growth
			rleafn(p,:) = 1._r8	!use flat rate
		elsewhere (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:)) !post-expansion
			!   rleafn(p,:) = (1._r8 - (hui(p)-huilfexpnp(p,:))/ &   !simple linear decline
			!                 (huilfmatnp(p,:)-huilfexpnp(p,:)))
			rleafn(p,:) = 1._r8	!use flat rate
		elsewhere
			rleafn(p,:) = 0._r8
		endwhere

		!calculate bunch sink size and update fruit allocation for each phytomer
		!assume bunch development increases until end of fruit-fill (summit of oil synthesis)
		!where (pgrainc(p,:) >= grnmx(ivt(p)))   !fruit may also stop growth before phenological
		!    rfruitn(p,:) = 0._r8                 !maturity when a max bunch C size is reached
		where (hui(p) >= huigrnnp(p,:) .and. hui(p) < grnmatnp(p,:) &
				.and. plai(p,:) > 0._r8)
			rfruitn(p,:) = (hui(p)- huigrnnp(p,:))/(grnmatnp(p,:)- huigrnnp(p,:))   !simple linear
			!   rfruitn(p,:) = (1- rsex(p,:))*(1- rabort(p,:)) !use flat rate, only rescale by rsex/rabort
		elsewhere
			rfruitn(p,:) = 0._r8
		endwhere
		where (huigrnnp(p,:) < huigrain(p)) rfruitn(p,:) = 0._r8 !the initial bunches before fruiting threshold are mostly male (Corley and Tinker 2003)


		  !PFT level allocation
		  ! ==================
		  !Although root allocation is not constant across soil environments,
		  !here assume root:leaf ratio is more or less constant, do not consider coarse root for oil palm

		  !oil palm stem growth is also closely related to leaf/phytomer growth

		  !above-ground vegetative growth is more sink-limited (more or less constant at given developmental stage, leaf number and size)
		  !while yield is source-limited (e.g. depend on excess assimilates available)
		  if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p)) then
			 arepr(p) = 0._r8

			 if (sum(rleafn(p,:)) == 0._r8) then !when plai all reach laimx too early, all allocation goes to root
				aleaf(p) = 1.e-5_r8
				astem(p) = 0._r8
				aroot(p) = 1._r8 - aleaf(p) - astem(p)
			 else
				aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -  &
							  (arooti(ivt(p)) - arootf(ivt(p))) *  &
						   min(1._r8, real(idpp(p))/mxmat(ivt(p)))))
				aleaf(p) = max(1.e-5_r8, (1._r8 - aroot(p))*fleafi(ivt(p)))
				astem(p) = 1._r8 - aleaf(p) - aroot(p)
			 end if
			 grain_flag(p) = 0._r8

			 !astemi(p) = astem(p) ! save for use by equations after shift
			 aleafi(p) = aleaf(p) ! to reproductive phenology stage
			 idpp2(p)  = idpp(p)
		  else if (hui(p) >= huigrain(p)) then
			 grain_flag(p) = 1._r8

			 !assume allocation rates apprach the final values and stablize after 1/2 of max age,
			 aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
						   (arooti(ivt(p)) - arootf(ivt(p))) * &
							min(1._r8, real(idpp(p))/mxmat(ivt(p)))))

			 aleaf(p) = max(1.e-5_r8, min(1._r8, aleafi(p) - &
						   (aleafi(p) - aleaff(ivt(p))) * &
							max(0._r8, min(1._r8, real(idpp(p) - idpp2(p))/ &
						   (declfact(ivt(p))*mxmat(ivt(p))- idpp2(p))))**allconsl(ivt(p)) ))

			 astem(p) = 1._r8 - aleaf(p) - aroot(p)

			!palm is perennial evergreen and yields every month, using last year NPP will delay the yield dynamics. 23.01.2015
			
			!use monthly sum of NPP so that fruit allocation dynamics is not delayed 20.04.2015
			 arepr(p) = 2._r8/(1.0_r8 + exp(-b_par(ivt(p))*(monsum_npp(p) - 100._r8))) - a_par(ivt(p))

			!when plai all reach laimx too early, allocation goes to root and fruit
			 if (sum(rleafn(p,:)) == 0._r8) then
				!arepr(p) = 1.e5_r8  !the same effect as setting aleaf(p) = 1.e-5_r8
				aleaf(p) = 1.e-5_r8
				astem(p) = 0._r8
			 end if
			 if (sum(rfruitn(p,:)) == 0._r8) then
				arepr(p) = 0._r8
			 end if

		  else !pre-emergence
			 aleaf(p) = 1.e-5_r8
			 astem(p) = 0._r8
			 aroot(p) = 0._r8
			 arepr(p) = 0._r8

		  end if


		  !Phytomer level allocation
		  ! =======================

		  !Sub-PFT level allocation to leaf/grain of phytomers
		  !calculate the total sink size of expended / bud phytomers
		  psum_exp = sum(rleafn(p,:), mask=(hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:)))
		  psum_bud = sum(rleafn(p,:), mask=(hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:)))

		  where (hui(p) >= huileafnp(p,:) .and. hui(p) < huilfexpnp(p,:)) !pre-expansion growth
			 this%aleafn(p,:) = rleafn(p,:) / max(1.e-5_r8, psum_bud)
		  elsewhere (hui(p) >= huilfexpnp(p,:) .and. hui(p) < huilfmatnp(p,:)) !post-expansion
			 this%aleafn(p,:) = rleafn(p,:) / max(1.e-5_r8, psum_exp)
		  endwhere
		  this%afruitn(p,:) = rfruitn(p,:) / max(1.e-5_r8, sum(rfruitn(p,:)))

             else   ! .not phytomer structure


               ! same phases appear in subroutine CropPhenology

               ! Phase 1 completed:
               ! ==================
               ! if hui is less than the number of gdd needed for filling of grain
               ! leaf emergence also has to have taken place for lai changes to occur
               ! and carbon assimilation
               ! Next phase: leaf emergence to start of leaf decline

               if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p)) then

                  ! allocation rules for crops based on maturity and linear decrease
                  ! of amount allocated to roots over course of the growing season

                  if (peaklai(p) == 1) then ! lai at maximum allowed
                     arepr(p) = 0._r8
                     aleaf(p) = 1.e-5_r8
                     astem(p) = 0._r8
		     !for perennial crops
		     if (perennial(ivt(p)) == 1) then
                        aleaf(p) = aleaff(ivt(p))
                        astem(p) = astemf(ivt(p))
                     end if
                     aroot(p) = 1._r8 - arepr(p) - aleaf(p) - astem(p)
                  else
                     arepr(p) = 0._r8
                     aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -   &
                          (arooti(ivt(p)) - arootf(ivt(p))) *  &
                          min(1._r8, hui(p)/gddmaturity(p))))
                     fleaf = fleafi(ivt(p)) * (exp(-bfact(ivt(p))) -         &
                          exp(-bfact(ivt(p))*hui(p)/huigrain(p))) / &
                          (exp(-bfact(ivt(p)))-1) ! fraction alloc to leaf (from J Norman alloc curve)
                     aleaf(p) = max(1.e-5_r8, (1._r8 - aroot(p)) * fleaf)
					 
		    if (perennial(ivt(p)) == 1 ) then
                     aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
                                   (arooti(ivt(p)) - arootf(ivt(p))) * &
                                   min(1._r8, real(idpp(p))/real(mxmat(ivt(p))))))
                     aleaf0(p)= max(1.e-5_r8, (1._r8 - arooti(ivt(p)))*fleafi(ivt(p)))
                     aleaf(p) = max(1.e-5_r8, min(1._r8, aleaf0(p) - &
                                   (aleaf0(p) - aleaff(ivt(p))) * &
                                   min(1._r8, real(idpp(p))/real(mxmat(ivt(p))))))
		    end if
					 
                     astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
                  end if

                  ! AgroIBIS included here an immediate adjustment to aleaf & astem if the
                  ! predicted lai from the above allocation coefficients exceeded laimx.
                  ! We have decided to live with lais slightly higher than laimx by
                  ! enforcing the cap in the following tstep through the peaklai logic above.

                  astemi(p) = astem(p) ! save for use by equations after shift
                  aleafi(p) = aleaf(p) ! to reproductive phenology stage begins
                  grain_flag(p) = 0._r8 ! setting to 0 while in phase 2

                  ! Phase 2 completed:
                  ! ==================
                  ! shift allocation either when enough gdd are accumulated or maximum number
                  ! of days has elapsed since planting
               
	       !for perennial crops add a continuous phenological cycle at annual time step (Y.Fan)
	       else if (hui(p) >= huigrain(p) .and. perennial(ivt(p)) == 1) then
		    if (hui(p) >= huigrain2(p) .and. hui(p) < gddmaturity2(p)) then
			  if (peaklai(p) == 1) then !
				 if (perennial(ivt(p)) == 1) then
					aleaf(p) = aleaff(ivt(p))
					astem(p) = astemf(ivt(p))
				 else
					aleaf(p) = 1.e-5_r8
					astem(p) = 0._r8
				 end if
				 aroot(p) = arootf(ivt(p))
			  else
				!root allocation continue decrease to the base level arootf until the end of life
				 aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
							   (arooti(ivt(p)) - arootf(ivt(p))) * &
							   min(1._r8, real(idpp(p))/real(mxmat(ivt(p))))))
				!leaf /stem allocation decreases to the base level during grainfill
				 if (aleafi(p) > aleaff(ivt(p))) then
					aleaf(p) = max(1.e-5_r8, max(aleaff(ivt(p)), aleaf(p) * &
							 (1._r8 - min((hui(p)- huigrain2(p))/           &
							 ((gddmaturity2(p)*declfact(ivt(p)))-           &
							 huigrain2(p)),1._r8)**allconsl(ivt(p)) )))
				 end if
				 if (astemi(p) > astemf(ivt(p))) then
					astem(p) = max(0._r8, max(astemf(ivt(p)), astem(p) * &
							 (1._r8 - min((hui(p)- huigrain2(p))/        &
							 ((gddmaturity2(p)*declfact(ivt(p)))-        &
							 huigrain2(p)),1._r8)**allconss(ivt(p)) )))
				 end if
			  end if
			  arepr(p) = 1._r8 - aroot(p) - astem(p) - aleaf(p)
			  !retranslocation starts at the begining of grnfill for one step like corn and cereals
			  if (grain_flag(p) == 0._r8) then
				 t1 = 1 / dt
				 leafn_to_retransn(p) = t1 * ((leafc(p) / leafcn(ivt(p))) - (leafc(p) / fleafcn(ivt(p))))
				 livestemn_to_retransn(p) = t1 * ((livestemc(p) / livewdcn(ivt(p))) - (livestemc(p) / fstemcn(ivt(p))))
				 if (ffrootcn(ivt(p)) > 0._r8) then
					frootn_to_retransn(p) = t1 * ((frootc(p) / frootcn(ivt(p))) - (frootc(p) / ffrootcn(ivt(p))))
				 else
					frootn_to_retransn(p) = 0._r8
				 end if
				 grain_flag(p) = 1._r8
			  end if

		    else !outside of grainfill period, root and leaf alloc continue decline through ageing
			  arepr(p) = 0._r8
			  grain_flag(p) = 0._r8 !set to 0 before the next grainfill
			  if (peaklai(p) == 1) then ! lai at maximum allowed
				 if (perennial(ivt(p)) == 1) then
					aleaf(p) = aleaff(ivt(p))
					astem(p) = astemf(ivt(p))
				 else
					aleaf(p) = 1.e-5_r8
					astem(p) = 0._r8
				 end if
				 aroot(p) = 1._r8 - arepr(p) - aleaf(p) - astem(p)
			  else
				 aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
							   (arooti(ivt(p)) - arootf(ivt(p))) * &
							   min(1._r8, real(idpp(p))/real(mxmat(ivt(p))))))
				 aleaf(p) = max(1.e-5_r8, min(1._r8, aleaf0(p) - &
							   (aleaf0(p) - aleaff(ivt(p))) * &
							   min(1._r8, real(idpp(p))/real(mxmat(ivt(p))))))
				 astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
			  end if
		    end if
		   
               else if (hui(p) >= huigrain(p)) then
                  aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
                       (arooti(ivt(p)) - arootf(ivt(p))) * min(1._r8, hui(p)/gddmaturity(p))))
                  if (astemi(p) > astemf(ivt(p))) then
                     astem(p) = max(0._r8, max(astemf(ivt(p)), astem(p) * &
                          (1._r8 - min((hui(p)-                 &
                          huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                          huigrain(p)),1._r8)**allconss(ivt(p)) )))
                  end if

                  ! If crops have hit peaklai, then set leaf allocation to small value
                  if (peaklai(p) == 1) then 
                     aleaf(p) = 1.e-5_r8
                  else if (aleafi(p) > aleaff(ivt(p))) then
                     aleaf(p) = max(1.e-5_r8, max(aleaff(ivt(p)), aleaf(p) * &
                          (1._r8 - min((hui(p)-                    &
                          huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                          huigrain(p)),1._r8)**allconsl(ivt(p)) )))
                  end if

                  !Beth's retranslocation of leafn, stemn, rootn to organ
                  !Filter excess plant N to retransn pool for organ N
                  !Only do one time then hold grain_flag till onset next season

                  ! slevis: Will astem ever = astemf exactly?
                  ! Beth's response: ...looks like astem can equal astemf under the right circumstances.
                  !It might be worth a rewrite to capture what I was trying to do, but the retranslocation for
                  !corn and wheat begins at the beginning of the grain fill stage, but for soybean I was holding it
                  !until after the leaf and stem decline were complete. Looking at how astem is calculated, once the
                  !stem decline is near complete, astem should (usually) be set to astemf. The reason for holding off
                  !on soybean is that the retranslocation scheme begins at the beginning of the grain phase, when the
                  !leaf and stem are still growing, but declining. Since carbon is still getting allocated and now
                  !there is more nitrogen available, the nitrogen can be diverted from grain. For corn and wheat
                  !the impact was probably enough to boost productivity, but for soybean the nitrogen was better off
                  !fulfilling the grain fill. It seems that if the peak lai is reached for soybean though that this
                  !would be bypassed altogether, not the intended outcome. I checked several of my output files and
                  !they all seemed to be going through the retranslocation loop for soybean - good news.

                  if (astem(p) == astemf(ivt(p)) .or. &
                       (ivt(p) /= ntmp_soybean .and. ivt(p) /= nirrig_tmp_soybean .and.&
                        ivt(p) /= ntrp_soybean .and. ivt(p) /= nirrig_trp_soybean)) then
                     if (grain_flag(p) == 0._r8) then
                        t1 = 1 / dt
                        leafn_to_retransn(p) = t1 * max(leafn(p)- (leafc(p) / fleafcn(ivt(p))),0._r8)
                        livestemn_to_retransn(p) = t1 * max(livestemn(p) - (livestemc(p) / fstemcn(ivt(p))),0._r8)
                        frootn_to_retransn(p) = 0._r8
                        if (ffrootcn(ivt(p)) > 0._r8) then
                           frootn_to_retransn(p) = t1 * max(frootn(p) - (frootc(p) / ffrootcn(ivt(p))),0._r8)
                        end if
                        grain_flag(p) = 1._r8
                     end if
                  end if

                  arepr(p) = 1._r8 - aroot(p) - astem(p) - aleaf(p)
		 !NOTE: B. Drewniak's N restranslocation during organ development might not check for N balance,
		 !Because in the phenology model, a high N content (leafcn) is used again for leaf residue at offset (offset means harvest for annual crops)
		 !fleafcn should be used for vegetative residues at harvest, if Beth Drewniak's N restranslocation is used. (Y.Fan 2015)

               else                   ! pre emergence
                  aleaf(p) = 1.e-5_r8 ! allocation coefficients should be irrelevant
                  astem(p) = 0._r8    ! because crops have no live carbon pools;
                  aroot(p) = 0._r8    ! this applies to this "else" and to the "else"
                  arepr(p) = 0._r8    ! a few lines down
               end if


	     end if !oil palm .or. normal crops

               f1 = aroot(p) / aleaf(p)
               f3 = astem(p) / aleaf(p)
               f5 = arepr(p) / aleaf(p)
               g1 = 0.25_r8

            else   ! .not croplive
               f1 = 0._r8
               f3 = 0._r8
               f5 = 0._r8
               g1 = 0.25_r8
            end if
         end if

         ! based on available C, use constant allometric relationships to
         ! determine N requirements
         if(use_fun)then ! In FUN, growth respiration is not part of the allometry calculation. 
	         if (woody(ivt(p)) == 1.0_r8) then
	            c_allometry(p) = (1._r8)*(1._r8+f1+f3*(1._r8+f2))
	            n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
	                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
	         else if (ivt(p) >= npcropmin) then ! skip generic crops
	            cng = graincn(ivt(p))
	            c_allometry(p) = (1._r8)*(1._r8+f1+f5+f3*(1._r8+f2))
	            n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
	                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
	         else
	            c_allometry(p) = 1._r8+f1
	            n_allometry(p) = 1._r8/cnl + f1/cnfr
	         end if           
        else !no FUN.
	         if (woody(ivt(p)) == 1.0_r8) then
	            c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
	            n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
	                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
	         else if (ivt(p) >= npcropmin) then ! skip generic crops
	            cng = graincn(ivt(p))
	            c_allometry(p) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
	            n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
	                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
	         else
	            c_allometry(p) = 1._r8+g1+f1+f1*g1
	            n_allometry(p) = 1._r8/cnl + f1/cnfr
	         end if
         end if !FUN

         ! when we have "if (leafn(p) == 0.0_r8)" below then we
         ! have floating overflow (out of floating point range)
         ! error in "actual_leafcn(p) = leafc(p) / leafn(p)"
         if (leafn(p) < n_min ) then
            ! to avoid division by zero, and to set leafcn to missing value for history files
            this%actual_leafcn(p) = spval
         else                        	
            ! leaf CN ratio 
            this%actual_leafcn(p) = leafc(p)  / leafn(p)
         end if
            

         if (nscalar_opt) then

            leafcn_min = leafcn(ivt(p)) - 10.0_r8
            leafcn_max = leafcn(ivt(p)) + 10.0_r8

            this%actual_leafcn(p) = max( this%actual_leafcn(p), leafcn_min-0.0001_r8 )
            this%actual_leafcn(p) = min( this%actual_leafcn(p), leafcn_max )

            nscalar = (this%actual_leafcn(p) - leafcn_min ) / (leafcn_max - leafcn_min)  ! Nitrogen scaler factor
            nscalar = min( max(0.0_r8, nscalar), 1.0_r8 )
         else ! if (nscalar_opt == .false.) then
            nscalar = 1.0_r8
         end if


         if (substrate_term_opt) then
            c = patch%column(p)
            sminn_total = 0.0_r8
            do j = 1, nlevdecomp
               sminn_total = sminn_total + sminn_vr(c,j) * dzsoi_decomp(j)
            end do
            Kmin = 1.0_r8
            substrate_term = sminn_total / (sminn_total + Kmin)
         else ! if (substrate_term_opt == .false) then
            substrate_term = 1.0_r8
         end if

         if (.not. temp_scalar_opt) then
            temp_scalar = 1.0_r8
         else !(temp_scalar_opt == .true.) then
            c = patch%column(p)
            temp_scalar=t_scalar(c,1)
            temp_scalar = min( max(0.0_r8, temp_scalar), 1.0_r8 )
         end if
         
         if(use_fun)then ! in FUN, plant_ndemand is just used as a maximum draw on soil N pools. 
             plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))
         else !FUN
	         if (plant_ndemand_opt == 0) then
	            plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))
	         else if (plant_ndemand_opt == 1) then
	            plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p)) * substrate_term
	         else if (plant_ndemand_opt == 2) then      ! N uptake happens at day time only

	            if (gpp(p) > 0.0_r8) then
	               Vmax_N = 2.7E-8_r8
	               plant_ndemand(p) = Vmax_N * frootc(p) * substrate_term * temp_scalar  * nscalar
	            else
	               plant_ndemand(p) = 0.0_r8
	            end if
	         else if (plant_ndemand_opt == 3) then      ! N uptake happens at day and night time

	            if (laisun(p)+laisha(p) > 0.0_r8) then
	               Vmax_N = 2.7E-8_r8
	               plant_ndemand(p) =  Vmax_N  * frootc(p) * substrate_term * temp_scalar * nscalar
	            else
	               plant_ndemand(p) = 0.0_r8
	            end if

	            if (this%actual_leafcn(p) < leafcn_min )then
	               plant_ndemand(p) = 0.0_r8
	            end if

	         end if
         end if  !FUN

         !if (leafn(p) < n_min ) then
            !! to set leafcn to missing value for history files
            !this%actual_leafcn(p) = spval
         !end if

         ! retranslocated N deployment depends on seasonal cycle of potential GPP
         ! (requires one year run to accumulate demand)

         tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

         ! Adding the following line to carry max retransn info to CN Annual Update
         tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))

         ! Beth's code: crops pull from retransn pool only during grain fill;
         !              retransn pool has N from leaves, stems, and roots for
         !              retranslocation

         !palm pulls from retransn pool continuously after fruiting starts (Y.Fan)
         if (phytomer(ivt(p)) > 0) then !pulls gradually from retransn pool accroding to gpp like trees
            avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
         else if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
         !if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
            avail_retransn(p) = plant_ndemand(p)
         else if (ivt(p) < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
            avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
         else
            avail_retransn(p) = 0.0_r8
         end if

         ! make sure available retrans N doesn't exceed storage
         avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)

         ! modify plant N demand according to the availability of
         ! retranslocated N
         ! take from retransn pool at most the flux required to meet
         ! plant ndemand

         if (plant_ndemand(p) > avail_retransn(p)) then
            retransn_to_npool(p) = avail_retransn(p)
         else
            retransn_to_npool(p) = plant_ndemand(p)
         end if

         if ( .not. use_fun ) then
            plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)
         end if

      end do ! end patch loop

    end associate

  end subroutine calc_plant_nitrogen_demand

end module NutrientCompetitionFlexibleCNMod
