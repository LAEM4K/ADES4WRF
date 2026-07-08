! ****************************************************************************
! *  Source: A new semi-empirical dust emission scheme for East Asia in
!            WRF-Chem constrained by in situ tower observations
! *  Author: Lee, Jae-Hyeong (jaehyeonglee@kongju.ac.kr), Lee, Sang-Hyun (sanghyun@kongju.ac.kr)
! *  This scheme is enabled at dust_opt = 6
! ****************************************************************************

MODULE NEW_DUST

  USE module_data_gocart_dust

  IMPLICIT NONE  

  INTRINSIC max, min

CONTAINS
  SUBROUTINE new_dust_driver(ktau,dt,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,erod_dri,dustin,snowc,tsk,zs,   &
         ivgtyp,isltyp,vegfra,lai_vegmask,xland,xlat,xlong,gsw,dx,g,emis_dust,      &
         ust,znt,clay_wrf,sand_wrf,clay_nga,sand_nga,                     &
         tot_dust,tot_edust,                                               &
         new_source,                                              &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  USE module_data_sorgam, ONLY: factnuma,factnumc,soilfac

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: julday, ktau,                            &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                              moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                           chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,           &
         INTENT(INOUT ) ::                                                 &
         emis_dust
   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,     &
         INTENT(IN   ) ::                            smois
   REAL, DIMENSION( config_flags%num_soil_layers ) ,                       &
         INTENT(IN   ) ::                            zs
   REAL, DIMENSION( ims:ime , jms:jme, ndcls )             ,               &
         INTENT(IN   ) ::                            erod,erod_dri
   REAL, DIMENSION( ims:ime , jms:jme, 5 )                 ,               &
         INTENT(INOUT) ::                            dustin
   REAL, DIMENSION( ims:ime , jms:jme )                    ,               &
         INTENT(IN   ) ::                                                  &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                     vegfra,               &
                                                     lai_vegmask,          &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong,                &
                                                     ust,                  &
                                                     znt,                  &
                                                     clay_wrf,             &
                                                     sand_wrf,             &
                                                     clay_nga,             &
                                                     sand_nga,             &
                                                     snowc,                &
                                                     tsk                    
  REAL, DIMENSION( ims:ime, 12, jms:jme),                                  &
        INTENT(IN   ) ::                             new_source            
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                          &
         INTENT(IN   ) ::                                                  &
                                                     alt,                  &
                                                     t_phy,                &
                                                     dz8w,p8w,             &
                                                     u_phy,v_phy,rho_phy
  REAL, DIMENSION( ims:ime , jms:jme )                    ,                &
         INTENT(  OUT) ::                            tot_edust
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                          &
         INTENT(  OUT) ::                            tot_dust
  REAL, INTENT(IN   ) :: dt,dx,g

  ! Local variables

  INTEGER :: nmx,smx,i,j,k,imx,jmx,lmx,lhave
  INTEGER,DIMENSION (1,1) :: ilwi
  REAL*8, DIMENSION (1,1) :: erodtot
  REAL*8, DIMENSION (1,1) :: vegmask
  REAL*8, DIMENSION (1,1) :: gravsm
  REAL*8, DIMENSION (1,1) :: rfreeze,tskin,snowin,roughnessin,nsrcin
  REAL, DIMENSION( ims:ime , jms:jme ) :: clay,sand
  REAL*8, DIMENSION (1,1) :: drylimit
  REAL*8, DIMENSION (5)   :: tc,bems
  REAL*8, DIMENSION (1,1) :: airden,airmas,ustar
  REAL*8, DIMENSION (1) :: dxy
  REAL*8, DIMENSION (3) :: massfrac
  REAL*8                :: volsm
  REAL :: conver,converi
  REAL :: ustart
  REAL*8 :: zwant
  REAL*8 :: cf,rough_cor_in, lambda     
  REAL*8, PARAMETER  :: c_lambda = 0.35 

  conver=1.e-9
  converi=1.e9

  ! Number of dust bins

  imx=1
  jmx=1
  lmx=1
  nmx=ndust
  smx=ngsalt

  k=kts
  DO j=jts,jte
  DO i=its,ite



    ! Don't do dust over water!!!
    IF (xland(i,j) .lt. 1.5 .AND. (new_source(i, config_flags%start_month, j)) .GT. 0.) THEN
      ilwi(1,1)=1

      ! Total concentration at lowest model level. This is still hardcoded
      ! for 5 bins.             
      IF(config_flags%chem_opt == CB05_SORG_AQ_KPP .or. &
         config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP ) then
        tc(1)=0.0
        tc(2)=0.0
        tc(3)=0.0
        tc(4)=0.0
        tc(5)=0.0
      ELSE
        tc(1)=chem(i,kts,j,p_dust_1)*conver
        tc(2)=chem(i,kts,j,p_dust_2)*conver
        tc(3)=chem(i,kts,j,p_dust_3)*conver
        tc(4)=chem(i,kts,j,p_dust_4)*conver
        tc(5)=chem(i,kts,j,p_dust_5)*conver
      END IF

      !new dust: soil freezing
      tskin(1,1) = tsk(i,j)-273.15
     IF (tskin(1,1) .LE. -5.) THEN
       rfreeze(1,1) = 1-(tskin(1,1)/(-100.))
     ELSE IF (tskin(1,1) .GT. -5.) THEN
       rfreeze(1,1) = 1.
     END IF

     ! ! Air mass and density at lowest model level.
      airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g
      airden(1,1)=rho_phy(i,kts,j)
      ustar(1,1)=ust(i,j)
      dxy(1)=dx*dx
      snowin(1,1) = snowc(i,j)


      clay(i,j) =clay_wrf(i,j)
      sand(i,j) =sand_wrf(i,j)

      ! Mass fractions of clay, silt, and sand.

      massfrac(1)=clay(i,j)
      massfrac(2)=1-(clay(i,j)+sand(i,j))
      massfrac(3)=sand(i,j)

      ! Do not allow areas with bedrock, lava, or land-ice to loft.
      !==AFWA==
      IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
          isltyp(i,j) .eq. 18) then
        ilwi(1,1)=0
      ENDIF

      ! Volumetric soil moisture can be tuned here with a dust_smtune
      ! set in the namelist.

      volsm=max(smois(i,1,j),0.)

      ! Calculate gravimetric soil moisture.

      gravsm(1,1)=100*volsm/((1.-porosity(isltyp(i,j)))*(2.65*(1-clay(i,j))+2.50*clay(i,j)))
      drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
 

      cf=vegfra(i,j)
      cf = cf/100.d0
      IF (cf .GT. 0.95) THEN
          cf = 0.95
      ENDIF

      !new dust: vegetation-roughness drag (Raupach, 1992)
       lambda = - c_lambda*dlog( 1.d0 - cf )
       call r_c(lambda, rough_cor_in)
       roughnessin(1,1)= rough_cor_in

      ! Call dust emission routine.
      call source_dust(imx, jmx, lmx, nmx, smx, dt, tc, ustar, massfrac, &
                       ilwi, dxy, gravsm, volsm, airden, airmas,rfreeze, snowin, roughnessin,&
                       bems, ustart, g, drylimit)

      IF(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
        dustin(i,j,1:5)=tc(1:5)*converi
      ELSE IF(config_flags%chem_opt == CB05_SORG_AQ_KPP .or.  &
              config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP) then
!KW            chem(i,kts,j,p_p25i)=chem(i,kts,j,p_p25i) &
!                         +.25*(tc(1)+.286*tc(2))*converi
!            chem(i,kts,j,p_p25i)=max(chem(i,kts,j,p_p25i),1.e-16)
!            chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) &
!                         +.75*(tc(1)+.286*tc(2))*converi
!            chem(i,kts,j,p_p25j)=max(chem(i,kts,j,p_p25j),1.e-16)
!            chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) &
!                         +(.714*tc(2)+tc(3))*converi
!            chem(i,kts,j,p_soila)=max(chem(i,kts,j,p_soila),1.e-16)
         chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) + 0.03*sum(tc(1:5))*converi
         chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) + 0.97*1.02*sum(tc(1:5))*converi

         chem(i,kts,j,p_ac0) =  chem(i,kts,j,p_ac0) +  0.03*sum(tc(1:5))*converi*factnuma*soilfac
         chem(i,kts,j,p_corn) =  chem(i,kts,j,p_corn) + 0.97*1.02*sum(tc(1:5))*converi*factnumc*soilfac
      ELSE
        chem(i,kts,j,p_dust_1)=tc(1)*converi
        chem(i,kts,j,p_dust_2)=tc(2)*converi
        chem(i,kts,j,p_dust_3)=tc(3)*converi
        chem(i,kts,j,p_dust_4)=tc(4)*converi
        chem(i,kts,j,p_dust_5)=tc(5)*converi
      ENDIF


      ! For output diagnostics (g m^-2 s^-1)

      emis_dust(i,1,j,p_edust1)=bems(1)
      emis_dust(i,1,j,p_edust2)=bems(2)
      emis_dust(i,1,j,p_edust3)=bems(3)
      emis_dust(i,1,j,p_edust4)=bems(4)
      emis_dust(i,1,j,p_edust5)=bems(5)

      ! Diagnostic total emitted dust (g m^-2 s^-1)

      tot_edust(i,j)=(bems(1)+bems(2)+bems(3)+bems(4)+bems(5))

    ENDIF

    ! Cumulative dust concentration (ug m^-3) and visibility (m) diagnostics
    ! Note visibility is capped at 20 km.  Simple visibility algorithm based
    ! on mean particle diameter for each dust bin - perfect spheres - perfect
    ! blackbodies.

    DO k=kts,kte
      tot_dust(i,k,j)=(chem(i,k,j,p_dust_1)+chem(i,k,j,p_dust_2)+       &
                       chem(i,k,j,p_dust_3)+chem(i,k,j,p_dust_4)+       &
                       chem(i,k,j,p_dust_5))*rho_phy(i,k,j)
    ENDDO

  ENDDO
  ENDDO

  END SUBROUTINE new_dust_driver

  SUBROUTINE source_dust(imx, jmx, lmx, nmx, smx, dt1, tc, ustar, massfrac,&
                         ilwi, dxy, gravsm, volsm, airden, airmas, rfreeze,snowin, roughnessin,&
                         bems, ustart, g0, drylimit)

  ! ****************************************************************************
  ! *  Source: A new semi-empirical dust emission scheme for East Asia in
  !            WRF-Chem constrained by in situ tower observations
  ! *  Author: Lee, Jae-Hyeong (jaehyeonglee@kongju.ac.kr), Lee, Sang-Hyun (sanghyun@kongju.ac.kr)
  ! ****************************************************************************

  INTEGER, INTENT(IN)   :: nmx,imx,jmx,lmx,smx
  INTEGER, INTENT(IN)   :: ilwi(imx,jmx)
  REAL*8, INTENT(IN)    :: snowin(imx,jmx)      !new dust
  REAL*8, INTENT(IN)    :: roughnessin(imx,jmx) !new dust
  REAL*8, INTENT(IN)    :: rfreeze(imx,jmx)     !new dust
  REAL*8, INTENT(IN)    :: ustar(imx,jmx)
  REAL*8, INTENT(IN)    :: gravsm(imx,jmx)
  REAL*8, INTENT(IN)    :: drylimit(imx,jmx) 
  REAL*8, INTENT(IN)    :: dxy(jmx)
  REAL*8, INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8, INTENT(OUT)   :: bems(imx,jmx,nmx) 
  REAL, INTENT(IN)    :: g0,dt1
  REAL, INTENT(OUT)   :: ustart
  REAL*8, INTENT(IN)    :: volsm

  REAL*8    :: den(smx), diam(smx)
  REAL*8    :: dvol(nmx), dlndp(nmx)
! REAL*8    :: distr_dust(nmx)
  REAL*8    :: dsurface(smx), ds_rel(smx)
  REAL*8    :: massfrac(3)
  REAL*8    :: u_ts0, u_ts, dsrc, srce, dmass, dvol_tot
  REAL*8    :: emit, emit_vol
  REAL      :: rhoa, g
  REAL*8    :: salt, stotal
  INTEGER   :: i, j, m, n, s

! Global tuning constant, alpha.  Sandblasting mass efficiency, beta.
! Beta maxes out for clay fractions above 0.2 => betamax.

  REAL, PARAMETER :: betamax=5.25E-4
  REAL, PARAMETER :: snow_gamma=0.027 !new dust (Jae-Hyeong Lee and Sang-Hyun Lee,KNU)
  REAL, PARAMETER :: von_karman=0.4   !new dust (Jae-Hyeong Lee and Sang-Hyun Lee,KNU)
  REAL*8 :: beta

! Rather than relying on simple log-linear regression between dust emission flux
! and friction velocity, the new scheme introduces a refined dust flux
! formulation by directly optimizing it against field measurements.
 REAL, PARAMETER :: cmb  = 0.18    ! new dust (Jae-Hyeong Lee and Sang-Hyun Lee,KNU)
 REAL, PARAMETER :: uidx = 1.72    ! new dust (Jae-Hyeong Lee and Sang-Hyun Lee,KNU)
! REAL, PARAMETER :: cmb=1.0    ! Marticorena et al,1997
! REAL, PARAMETER :: cmb=2.61   ! White,1979

! Parameters used in Kok distribution function. Advise not to play with 
! these without the expressed written consent of someone who knows what
! they're doing. (See Kok, 2010 PNAS for details on this scheme).

! REAL, PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
! REAL, PARAMETER :: gsd_dust=3.0     ! geometric standard deviation
! REAL, PARAMETER :: lambda=12.0D-6   ! crack propogation length (m)
! REAL, PARAMETER :: cv=12.62D-6      ! normalization constant

! Calculate saltation surface area distribution from sand, silt, and clay
! mass fractions and saltation bin fraction. This will later become a 
! modifier to the total saltation flux.  The reasoning here is that the 
! size and availability of saltators affects saltation efficiency. Based
! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

  DO n=1,smx
    dmass=massfrac(spoint(n))*frac_salt(n)
    dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))  
  ENDDO
  
! The following equation yields relative surface area fraction.  It will only
! work if you are representing the "full range" of all three soil classes.
! For this reason alone, we have incorporated particle sizes that encompass
! the clay class, to account for its relative area over the basal
! surface, even though these smaller bins would be unlikely to play any large
! role in the actual saltation process.

  stotal=SUM(dsurface(:))
  DO n=1,smx
    ds_rel(n)=dsurface(n)/stotal
  ENDDO

! Calculate total dust emission due to saltation of sand sized particles.
! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
! u_ts0 for vegetation roughness drag to get threshold friction velocity (u_ts). Then
! Additional consideration is given to the critical friction speed effect due to snow cover. Then
! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally, 
! calculate total dust emission (tot_emit)
! The new scheme completely eliminates the conventional dependence on static source strength functions (S) 

  g = g0*1.0E2                          ! (cm s^-2)
  emit=0.0

  DO n = 1, smx                         ! Loop over saltation bins
    den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
    diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
    DO i = 1,imx
      DO j = 1,jmx
        rhoa = airden(i,j,1)*1.0D-3     ! (g cm^-3)
 
        ! Threshold friction velocity as a function of the dust density and
        ! diameter from Bagnold (1941) (m s^-1).

        u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 

        ! Friction velocity threshold correction function based on physical
        ! properties related to moisture tension. Soil moisture greater than
        ! dry limit serves to increase threshold friction velocity (making
        ! it more difficult to loft dust). When soil moisture has not reached
        ! dry limit, treat as dry (no correction to threshold friction
        ! velocity). GAC

        ! Calculate threshold friction velocity. If volumetric (gravimetric)
        ! water content is less than the drylimit, then the threshold friction
        ! velocity (u_ts) will be equal to the dry threshold friction velocity
        ! (u_ts0). EDH

        IF (gravsm(i,j) > drylimit(i,j)) THEN
          u_ts = MAX(0.0D+0,roughnessin(1,1)*u_ts0*(sqrt(1.0+1.21*(gravsm(i,j)-drylimit(i,j))**0.68)))
        ELSE
          u_ts = u_ts0*roughnessin(1,1)  !new dust - roughness
        END IF 

        !new dust, the snow cover effect ∆u_(*,t) is incorporated following
        !Kurosaki and Mikami (2004), where a 10% snow cover fraction increases
        !the threshold by approximately 0.015 m s-1.
        IF (snowin(i,j) .GT. 0.) THEN
          u_ts = u_ts + (von_karman*snow_gamma * (snowin(i,j)*100.))/LOG(10/0.01)
        ELSE
          u_ts = u_ts
        END IF
        

        ! Bin 7 threshold friction velocity for diagnostic dust lofting
        ! potential product

        IF (n .eq. 7) THEN  ! Saltation bin 7 is small sand
          ustart = u_ts
        ENDIF
 
        ! Saltation flux (kg m^-1 s^-1) from MB95
        ! ds_rel is the relative surface area distribution
        ! new dust. Based on the independent non-linear best-fit curve derived
        ! from the filtered monitoring tower data, we determined a
        ! proportionality constant (α_new) of 0.18 and a friction velocity
        ! exponent (β_new) of 1.72.

        IF (ustar(i,j) .gt. u_ts .and. ilwi(i,j) == 1) THEN
          salt = cmb*ds_rel(n)*(airden(i,j,1)/g0)*(ustar(i,j)**uidx)* &
                 (1. + u_ts/ustar(i,j))*(1. - (u_ts**2)/(ustar(i,j)**2))
        ELSE 
          salt = 0.D0
        ENDIF
 
        ! Calculate total vertical mass flux (note beta has units of m^-1)
        ! Beta acts to tone down dust in areas with so few dust-sized particles
        ! that the lofting efficiency decreases.  Otherwise, super sandy zones
        ! would be huge dust producers, which is generally not the case.
        ! Equation derived from wind-tunnel experiments (see MB95).
 
        beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
        IF (beta .gt. betamax) THEN
          beta=betamax
        ENDIF
        emit=emit+salt*beta         ! (kg m^-2 s^-1)
      END DO
    END DO
  END DO                                 ! End do over saltation bins
   
  emit =  rfreeze(1,1) * emit   !new dust soil feezing effect

  ! Now distribute total vertical emission into dust bins and update
  ! concentration.

  DO n=1,nmx  ! Loop over all dust bins
    DO i=1,imx
      DO j=1,jmx

        ! Calculate total mass emitted

        dsrc = emit*distr_dust(n)*dxy(j)*dt1  ! (kg)
        IF (dsrc < 0.0) dsrc = 0.0

        ! Update dust mixing ratio at first model level.

        tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1) ! (kg/kg)
        bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
      END DO
    END DO
  END DO

  END SUBROUTINE source_dust
!*****************************************************************************
  SUBROUTINE r_c (x,r)
!
!   Y. Shao 17-07-92
!   CORRECTION FUNCTION FOR UST(D) BASED ON Raupach et al. (1992)
!   x = frontal area index
!
!   R_C = (1 - sig m x)^{1/2} (1 + m beta x)^{1/2}
!         Note I deife R_C = u_{*tR}/u_{*tS}
!         While Raupach et al. defined 
!                      R_C = u_{*tS}/u_{*tR} and their R function is
!   R_C = (1 - sig m x)^{-1/2} (1 + m beta x)^{-1/2}
!                      
!   sig   : basal to frontal area; about 1
!   m     : parameter less than 1; about 0.5
!   beta  : a ratio of drag coef.; about 90. ; changed to 200 (recommendation by
!   Y. Shao based on data by Gillies et al.)
!
      real(8) :: xc
      real(8), intent(in) :: x
      real(8), intent(out) :: r
      real(8), parameter :: sig=1., m=0.5, beta=200.
!
      xc = 1./(sig*m)
      IF (x.ge.xc) THEN
        r   = 999.           ! Full covered surface
      ELSE
        r   = dsqrt(1.-sig*m*x)*dsqrt(1.+m*beta*x)
      ENDIF
!
 END SUBROUTINE r_c



END MODULE NEW_DUST
