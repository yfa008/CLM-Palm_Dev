module CNVegStructUpdateMod

  !-----------------------------------------------------------------------
  ! Module for vegetation structure updates (LAI, SAI, htop, hbot)
  !
  ! !USES:
  use shr_kind_mod         , only: r8 => shr_kind_r8
  use shr_const_mod        , only : SHR_CONST_PI
  use clm_varctl           , only : iulog, use_cndv
  use CNDVType             , only : dgv_ecophyscon    
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use FrictionVelocityMod  , only : frictionvel_type
  use CNDVType             , only : dgvs_type
  use CNVegStateType       , only : cnveg_state_type
  use CropType             , only : crop_type
  use CNVegCarbonStateType , only : cnveg_carbonstate_type
  use CanopyStateType      , only : canopystate_type
  use PatchType            , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNVegStructUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNVegStructUpdate(num_soilp, filter_soilp, &
       waterdiagnosticbulk_inst, frictionvel_inst, dgvs_inst, cnveg_state_inst, crop_inst, &
       cnveg_carbonstate_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, use C state variables and epc to diagnose
    ! vegetation structure (LAI, SAI, height)
    !
    ! !USES:
    use pftconMod        , only : noveg, nc3crop, nc3irrig, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub
    use pftconMod        , only : npcropmin 
    use pftconMod        , only : ntmp_corn, nirrig_tmp_corn
    use pftconMod        , only : ntrp_corn, nirrig_trp_corn
    use pftconMod        , only : nsugarcane, nirrig_sugarcane
    use pftconMod        , only : pftcon
    use pftconMod        , only : noilpalm, nirrig_oilpalm !(Y.Fan)
    use pftconMod        , only : ncocoa, nirrig_cocoa !(Ashehad added)
    use clm_varpar       , only : nlevcan ! N canopy layers: nlevcan =1 big leaf; nlevcan >1 multilayer canopy (Y.Fan)
    use clm_varctl       , only : spinup_state
    use clm_time_manager , only : get_rad_step_size, get_days_per_year ! added get_days (Y.Fan)
    !
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilp       ! number of column soil points in patch filter
    integer                      , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(waterdiagnosticbulk_type)        , intent(in)    :: waterdiagnosticbulk_inst
    type(frictionvel_type)       , intent(in)    :: frictionvel_inst
    type(dgvs_type)              , intent(in)    :: dgvs_inst
    type(cnveg_state_type)       , intent(inout) :: cnveg_state_inst
    type(crop_type)              , intent(in)    :: crop_inst
    type(cnveg_carbonstate_type) , intent(in)    :: cnveg_carbonstate_inst
    type(canopystate_type)       , intent(inout) :: canopystate_inst
    !
    ! !REVISION HISTORY:
    ! 10/28/03: Created by Peter Thornton
    ! 2/29/08, David Lawrence: revised snow burial fraction for short vegetation
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g      ! indices
    integer  :: fp         ! lake filter indices
    real(r8) :: taper      ! ratio of height:radius_breast_height (tree allometry)
    real(r8) :: stocking   ! #stems / ha (stocking density)
    real(r8) :: ol         ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb         ! fraction of canopy layer covered by snow
    real(r8) :: tlai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_min   ! PATCH derived minimum tsai
    real(r8) :: tsai_alpha ! monthly decay rate of tsai
    real(r8) :: dt         ! radiation time step (sec)
    integer :: n,np1,np2,i      !indices
    real(r8):: plaimx       ! max lai per phytomer rescaled with planting density

    !real(r8) :: sla         ! specific leaf area at different canopy depths, temporary
    real(r8), parameter :: dtsmonth = 2592000._r8 ! number of seconds in a 30 day month (60x60x24x30)
    !-----------------------------------------------------------------------
    ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835
    !
    ! tsai(p) = max( tsai_alpha(ivt(p))*tsai_old + max(tlai_old-tlai(p),0_r8), tsai_min(ivt(p)) )
    ! notes:
    ! * RHS tsai & tlai are from previous timestep
    ! * should create tsai_alpha(ivt(p)) & tsai_min(ivt(p)) in pftconMod.F90 - slevis
    ! * all non-crop patches use same values:
    !   crop    tsai_alpha,tsai_min = 0.0,0.1
    !   noncrop tsai_alpha,tsai_min = 0.5,1.0  (includes bare soil and urban)
    !-------------------------------------------------------------------------------
    
    associate(                                                            & 
         ivt                =>  patch%itype                               , & ! Input:  [integer  (:) ] patch vegetation type                                

         woody              =>  pftcon%woody                            , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         slatop             =>  pftcon%slatop                           , & ! Input:  specific leaf area at top of canopy, projected area basis [m^2/gC]
         dsladlai           =>  pftcon%dsladlai                         , & ! Input:  dSLA/dLAI, projected area basis [m^2/gC]           
         z0mr               =>  pftcon%z0mr                             , & ! Input:  ratio of momentum roughness length to canopy top height (-)
         displar            =>  pftcon%displar                          , & ! Input:  ratio of displacement height to canopy top height (-)
         dwood              =>  pftcon%dwood                            , & ! Input:  density of wood (gC/m^3)                          
         ztopmx             =>  pftcon%ztopmx                           , & ! Input:
         laimx              =>  pftcon%laimx                            , & ! Input:
         
         allom2             =>  dgv_ecophyscon%allom2                   , & ! Input:  [real(r8) (:) ] ecophys const                                     
         allom3             =>  dgv_ecophyscon%allom3                   , & ! Input:  [real(r8) (:) ] ecophys const                                     

         nind               =>  dgvs_inst%nind_patch                    , & ! Input:  [real(r8) (:) ] number of individuals (#/m**2)                    
         fpcgrid            =>  dgvs_inst%fpcgrid_patch                 , & ! Input:  [real(r8) (:) ] fractional area of patch (pft area/nat veg area)    

         snow_depth         =>  waterdiagnosticbulk_inst%snow_depth_col          , & ! Input:  [real(r8) (:) ] snow height (m)                                   

         forc_hgt_u_patch   =>  frictionvel_inst%forc_hgt_u_patch       , & ! Input:  [real(r8) (:) ] observational height of wind at patch-level [m]     

        perennial           =>  pftcon%perennial         , & ! Input:  [integer (:)]  binary flag for perennial crop phenology (1=perennial, 0=not perennial) (added by Y.Fan)
        phytomer            =>  pftcon%phytomer          , & ! Input:  [integer (:)]   total number of phytomers in life time (if >0 use phytomer phenology)
      !  hui                 =>  pps%gddplant             , & ! Input:  [real(r8) (:)]  =gdd since planting (gddplant)
      !  huigrnnp            =>  pps%huigrnnp             , & ! Input:  [real(r8) (:,:)]  hui needed for start of grainfill of successive phytomers
      !  huilfmatnp          =>  pps%huilfmatnp           , & ! Input: [real(r8) (:,:)]  hui needed for leaf maturity of successive phytomers
      !  leaf_long           =>  pftcon%leaf_long         , & ! Input:  [real(r8) (:)]  leaf longevity (yrs)
        lfdays              =>  crop_inst%lfdays_patch               , & ! Input:  [integer (:,:)]  days past leaf emergence for each phytomer (Y.Fan)
        np                  =>  crop_inst%np_patch                   , & ! Input:  [integer (:)]   total number of phytomers having appeared so far
        rankp               =>  crop_inst%rankp_patch                , & ! Input:  [integer (:,:)]  rank of phytomers 1=youngest and 0=dead
        !livep               =>  pps%livep                , & ! Input:  [real(r8) (:,:)] Flag, true if this phytomer is alive
        psla                =>  crop_inst%psla_patch                 , & ! Input:  [real(r8) (:,:)] specific leaf area of each phytomer (added by Y.Fan)
        plai                =>  crop_inst%plai_patch                 , & ! Input:  [real(r8) (:,:)] one-sided leaf area index of each phytomer
        pleafc              =>  cnveg_carbonstate_inst%pleafc_patch               , & ! Input:  [real(r8) (:,:)] (gC/m2) phytomer leaf C
        plaipeak            =>  crop_inst%plaipeak_patch         , & ! Output: [integer  (:) ] Flag, 1: max allowed lai; 0: not at max (Y.Fan)

         leafc              =>  cnveg_carbonstate_inst%leafc_patch      , & ! Input:  [real(r8) (:) ] (gC/m2) leaf C
         deadstemc          =>  cnveg_carbonstate_inst%deadstemc_patch  , & ! Input:  [real(r8) (:) ] (gC/m2) dead stem C

         farea_burned       =>  cnveg_state_inst%farea_burned_col       , & ! Input:  [real(r8) (:) ] F. Li and S. Levis                                 
         htmx               =>  cnveg_state_inst%htmx_patch             , & ! Output: [real(r8) (:) ] max hgt attained by a crop during yr (m)          
         peaklai            =>  cnveg_state_inst%peaklai_patch          , & ! Output: [integer  (:) ] 1: max allowed lai; 0: not at max                  

         harvdate           =>  crop_inst%harvdate_patch                , & ! Input:  [integer  (:) ] harvest date                                       

         ! *** Key Output from CN***
         tlai               =>  canopystate_inst%tlai_patch             , & ! Output: [real(r8) (:) ] one-sided leaf area index, no burying by snow      
         tsai               =>  canopystate_inst%tsai_patch             , & ! Output: [real(r8) (:) ] one-sided stem area index, no burying by snow      
         htop               =>  canopystate_inst%htop_patch             , & ! Output: [real(r8) (:) ] canopy top (m)                                     
         hbot               =>  canopystate_inst%hbot_patch             , & ! Output: [real(r8) (:) ] canopy bottom (m)                                  
         elai               =>  canopystate_inst%elai_patch             , & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow    
         esai               =>  canopystate_inst%esai_patch             , & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow    
         frac_veg_nosno_alb =>  canopystate_inst%frac_veg_nosno_alb_patch & ! Output: [integer  (:) ] frac of vegetation not covered by snow [-]         
         )

      dt = real( get_rad_step_size(), r8 )

      ! constant allometric parameters
      taper = 200._r8
      stocking = 1000._r8

      ! convert from stems/ha -> stems/m^2
      stocking = stocking / 10000._r8

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         if (ivt(p) /= noveg) then

            tlai_old = tlai(p) ! n-1 value
            tsai_old = tsai(p) ! n-1 value

            ! update the leaf area index based on leafC and SLA
            ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923. 
            if (dsladlai(ivt(p)) > 0._r8) then
               tlai(p) = (slatop(ivt(p))*(exp(leafc(p)*dsladlai(ivt(p))) - 1._r8))/dsladlai(ivt(p))
            else
               tlai(p) = slatop(ivt(p)) * leafc(p)
            end if
            tlai(p) = max(0._r8, tlai(p))

            ! update the stem area index and height based on LAI, stem mass, and veg type.
            ! With the exception of htop for woody vegetation, this follows the DGVM logic.

            ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835 (see notes)
            ! Assumes doalb time step .eq. CLM time step, SAI min and monthly decay factor
            ! alpha are set by PFT, and alpha is scaled to CLM time step by multiplying by
            ! dt and dividing by dtsmonth (seconds in average 30 day month)
            ! tsai_min scaled by 0.5 to match MODIS satellite derived values
            if (ivt(p) == nc3crop .or. ivt(p) == nc3irrig) then ! generic crops

               tsai_alpha = 1.0_r8-1.0_r8*dt/dtsmonth
               tsai_min = 0.1_r8
            else
               tsai_alpha = 1.0_r8-0.5_r8*dt/dtsmonth
               tsai_min = 1.0_r8
            end if
            tsai_min = tsai_min * 0.5_r8
            tsai(p) = max(tsai_alpha*tsai_old+max(tlai_old-tlai(p),0._r8),tsai_min)
          ! new calculations for phytomer-based canopy structure
          ! for oil palm the stocking (planting density) and wood density values can be adjusted here (Y.Fan)
          if (phytomer(ivt(p)) > 0) then
             dwood(ivt(p)) = 1.0e5_r8   !wood density (gC/m3); cn:2.5e5_r8; lpj:2.0e5
             stocking = 150._r8         !planting density (stems/ha); 150 is standard stocking for oil palm plantation
             taper = 20._r8             !stem height 0-5m, diameter 20 to 75 cm

             !laimx is per phytomer for oil palm, not per PFT
             plaimx = laimx(ivt(p)) *min(stocking/150._r8, 2._r8) !rescale plaimx when planting is too sparse/too dense

             !dsladlai is set to zero, effectively psla(p,:) = slatop(ivt(p)); do not use SLA scaling for oil palm!
             !Oil palm SLA actually decrease slightly, instead of increase, from canopy top to bottom
             !SLA scaling was only intented for trees to give a canopy gradient in foliage nitrogen in the CLM4
             !CLM4.5 use Kn instead of dsladlai to scale foliage nitrogen profile

             np1 = sum(minloc(rankp(p,:), mask=(rankp(p,:) > 0))) !the top youngest expanded phytomer
             np2 = sum(maxloc(rankp(p,:), mask=(rankp(p,:) > 0))) !the bottom expanded phytomer, sequence np2 < np1
             psla(p,np2:np1) = (/ ((slatop(ivt(p)) + (np1-i)*dsladlai(ivt(p))), i=np2, np1, 1) /)
!             !or below function
!             do n = np1, np2, -1
!                if (hui(p) < huilfmatnp(p,n) .and. plai(p,n) < laimx(ivt(p))) &
!                   psla(p,n) = slatop(ivt(p)) + dsladlai(ivt(p))*sum(plai(p,n:np1)) !SLA stops increase when C allocation stops
!                !else if (rankp(p,n) == 0) then
!                !   exit
!                !end if
!             end do

             !only calculate LAI for expanded phytomers
             where (rankp(p,:) > 0)
                plai(p,:) = psla(p,:) * pleafc(p,:)
             elsewhere
                plai(p,:) = 0._r8
             endwhere
             where (plai(p,:) >= plaimx)
                plaipeak(p,:) = 1 ! used in CNAllocation
             elsewhere
                plaipeak(p,:) = 0
             endwhere
             tlai(p) = sum(plai(p,:))

             ! canopy integration of TLAI
!             if (nlevcan == 1) then !for one-layer sun/shade big-leaf canopy
!                tlai(p) = (slatop(ivt(p))*(exp(sum(pleafc(p,:))*dsladlai(ivt(p))) - 1._r8))/dsladlai(ivt(p))
!             else			       !for multilayer canopy, each phytomer is considered a layer
!                tlai(p) = sum(plai(p,:))
!             end if

	    !Stem area index depends on planting density
             tsai_min = 0.005_r8 !the same sai threshold used in STATICEcosysDynMod
             !stocking = stocking / 10000._r8 !convert from stems/ha -> stems/m^2
             !tsai(p) = 0.05_r8 * tlai(p) * stocking
             !the above sai might be too small. don't rescale sai again by stocking. LAI is already limited by stocking (Y.Fan 2015.11)
             !sai should consider the vertical thickness (height) of stem as well as the rachis of phytomers. LAI and SAI are affected by stocking to the same degree.
             tsai(p) = max(0.1_r8 * tlai(p), tsai_min) !0.1 is the similar ratio as trees

             !canopy top and bottom height follow trees
             stocking = stocking / 10000._r8 !convert from stems/ha -> stems/m^2
             !correct height calculation if doing accelerated spinup
             if (spinup_state == 2) then
                 htop(p) = ((3._r8 * deadstemc(p) * 10._r8 * taper * taper)/ &
                    (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
             else
                 htop(p) = ((3._r8 * deadstemc(p) * taper * taper)/ &
                     (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
             end if
             htop(p) = max(0.05_r8, min(htop(p),(forc_hgt_u_patch(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8))
             htop(p) = min(ztopmx(ivt(p)),htop(p))
             hbot(p) = max(0._r8, htop(p)-2._r8)  !assume palm canopy thichness == 2m
            ! "stubble" after harvest
            ! if (harvdate(p) < 999 .and. tlai(p) == 0._r8) then
            !    tsai(p) = 0.25_r8*(1._r8-farea_burned(c)*0.90_r8)
            ! end if

          else if (woody(ivt(p)) == 1._r8) then


            !if (woody(ivt(p)) == 1._r8) then

               ! trees and shrubs

               ! if shrubs have a squat taper 
               if (ivt(p) >= nbrdlf_evr_shrub .and. ivt(p) <= nbrdlf_dcd_brl_shrub) then
                  taper = 10._r8
                  ! otherwise have a tall taper
               else
                  taper = 200._r8
               end if

               ! trees and shrubs for now have a very simple allometry, with hard-wired
               ! stem taper (height:radius) and hard-wired stocking density (#individuals/area)
               if (use_cndv) then

                  if (fpcgrid(p) > 0._r8 .and. nind(p) > 0._r8) then

                     stocking = nind(p)/fpcgrid(p) !#ind/m2 nat veg area -> #ind/m2 patch area
                     htop(p) = allom2(ivt(p)) * ( (24._r8 * deadstemc(p) / &
                          (SHR_CONST_PI * stocking * dwood(ivt(p)) * taper))**(1._r8/3._r8) )**allom3(ivt(p)) ! lpj's htop w/ cn's stemdiam

                  else
                     htop(p) = 0._r8
                  end if

               else
                  !correct height calculation if doing accelerated spinup
                  if (spinup_state == 2) then
                    htop(p) = ((3._r8 * deadstemc(p) * 10._r8 * taper * taper)/ &
                         (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
                  else
                    htop(p) = ((3._r8 * deadstemc(p) * taper * taper)/ &
                         (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
                  end if

               endif

               ! Peter Thornton, 5/3/2004
               ! Adding test to keep htop from getting too close to forcing height for windspeed
               ! Also added for grass, below, although it is not likely to ever be an issue.
               htop(p) = min(htop(p),(forc_hgt_u_patch(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
               ! This becomes an issue when fire mortality is pushing deadstemc
               ! to 0.0.
               htop(p) = max(htop(p), 0.01_r8)

               hbot(p) = max(0._r8, min(3._r8, htop(p)-1._r8))

            else if (ivt(p) >= npcropmin) then ! prognostic crops

               if (tlai(p) >= laimx(ivt(p))) peaklai(p) = 1 ! used in CNAllocation
              if (perennial(ivt(p)) == 1) then
                if (tlai(p) >= laimx(ivt(p))) then
                   peaklai(p) = 1
                else  !allow leaf regrowth for perennial crops if tlai drop below laimx (Y.Fan)
                   peaklai(p) = 0
                end if
              end if

               if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
                   ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
                   ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane) then
                  tsai(p) = 0.1_r8 * tlai(p)
               else
                  tsai(p) = 0.2_r8 * tlai(p)
               end if

               ! "stubble" after harvest
               if (harvdate(p) < 999 .and. tlai(p) == 0._r8) then
                  tsai(p) = 0.25_r8*(1._r8-farea_burned(c)*0.90_r8)    !changed by F. Li and S. Levis
                  htmx(p) = 0._r8
                  peaklai(p) = 0
               end if
               !if (harvdate(p) < 999 .and. tlai(p) > 0._r8) write(iulog,*) 'CNVegStructUpdate: tlai>0 after harvest!' ! remove after initial debugging?

               ! canopy top and bottom heights
               htop(p) = ztopmx(ivt(p)) * (min(tlai(p)/(laimx(ivt(p))-1._r8),1._r8))**2
               htmx(p) = max(htmx(p), htop(p))
               htop(p) = max(0.05_r8, max(htmx(p),htop(p)))
               hbot(p) = 0.02_r8

            else ! generic crops and ...

               ! grasses

               ! height for grasses depends only on LAI
               htop(p) = max(0.25_r8, tlai(p) * 0.25_r8)

               htop(p) = min(htop(p),(forc_hgt_u_patch(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
               htop(p) = max(htop(p), 0.01_r8)

               hbot(p) = max(0.0_r8, min(0.05_r8, htop(p)-0.20_r8))
            end if

         else

            tlai(p) = 0._r8
            tsai(p) = 0._r8
            htop(p) = 0._r8
            hbot(p) = 0._r8

         end if

         ! adjust lai and sai for burying by snow. 
         ! snow burial fraction for short vegetation (e.g. grasses) as in
         ! Wang and Zeng, 2007.
       !!This condition is not reasonable for many tall crops (Y.Fan)
       !not necessary as htop > 0.25 and hbot>0.05 for grasses have been defined above; 
       !the first condition works for most PFTs
         if (ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
            ol = min( max(snow_depth(c)-hbot(p), 0._r8), htop(p)-hbot(p))
            fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
         else if (ivt(p) >= npcropmin) then ! prognostic crops (including oil palm)
            ol = min( max(snow_depth(c)-hbot(p), 0._r8), htop(p)-hbot(p))
            fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
         else
            fb = 1._r8 - max(min(snow_depth(c),0.2_r8),0._r8)/0.2_r8   ! 0.2m is assumed
            !depth of snow required for complete burial of grasses
         endif

         elai(p) = max(tlai(p)*fb, 0.0_r8)
         esai(p) = max(tsai(p)*fb, 0.0_r8)

         ! Fraction of vegetation free of snow
         if ((elai(p) + esai(p)) > 0._r8) then
            frac_veg_nosno_alb(p) = 1
         else
            frac_veg_nosno_alb(p) = 0
         end if

      end do

    end associate 

 end subroutine CNVegStructUpdate

end module CNVegStructUpdateMod
