!> @file photolysis_model_mod.f90
!------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2016 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
! Description:
! ------------
!> photolysis models and interfaces (Adapted from photolysis_model_mod.f90)
!> @todo Alles!
!------------------------------------------------------------------------------!
 MODULE kchem_photolysis
 
    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, pt, q, ql, zu, zw

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  time_since_reference_point

    USE pegrid,             ONLY: myid, threads_per_task

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb_s_inner, nzb, nzt

    USE radiation_model_mod,                                                   &
        ONLY:  calc_zenith, day, day_init, lat, lambda, lon,                   &
               time_utc, time_utc_init, zenith

    USE control_parameters,                                                    &
        ONLY:  initializing_actions           
!       ONLY:  cloud_droplets, cloud_physics, g, initializing_actions,         &
!              large_scale_forcing, lsf_surf, phi, pt_surface, rho_surface,    &
!              surface_pressure, time_since_reference_point

    USE kchem_kpp,                                                             &
        ONLY: nphot, phot_names, phot 

    USE kchem_driver,                                                          &
        ONLY: phot_frequen

!   REAL(wp), DIMENSION(0:0) ::  zenith,         & !< cosine of solar zenith angle
!                                sun_dir_lat,    & !< solar directional vector in latitudes
!                                sun_dir_lon       !< solar directional vector in longitudes
    USE kinds

#if defined ( __netcdf )
    USE NETCDF
#endif


    IMPLICIT NONE

    CHARACTER(10) :: photolysis_scheme = 'constant' ! photolysis_scheme should be namelist input later 
!                                        'constant', 
!                                        'simple' (Vila-Guerau de Arellano et al., 2011, J. Geophys. Res., 116, D07304, doi:10.1029/2010JD014857)  STILL  NOT IMPLEMENTED
!                                        'fastj'  (Wild et al., 2000, J. Atmos. Chem., 37, 245-282) STILL NOT IMPLEMENTED

!


!   LOGICAL ::  unscheduled_photolysis_calls = .TRUE., & !< flag parameter indicating whether additional calls of the photolysis code are allowed
!               constant_albedo = .FALSE.,            & !< flag parameter indicating whether the albedo may change depending on zenith
!               force_photolysis_call = .FALSE.,       & !< flag parameter for unscheduled photolysis calls
!               photolysis = .FALSE.,                  & !< flag parameter indicating whether the photolysis model is used
!               sun_up    = .TRUE.,                   & !< flag parameter indicating whether the sun is up or down
!               photolysis = .TRUE.,                 & !< flag parameter indicing whether shortwave photolysis shall be calculated
!               sun_direction = .FALSE.                 !< flag parameter indicing whether solar direction shall be calculated



    REAL(wp) :: time_photolysis = 0.0_wp,         & !< time since last call of photolysis code
                dt_photolysis = 0.0_wp,           & !< hotolysis model timestep
                skip_time_do_photolysis = 0.0_wp    !< Radiation model is not called before this time


!
!-- Variables and parameters used in Fast-J only
! ....

!   INTERFACE photolysis_check_parameters
!      MODULE PROCEDURE photolysis_check_parameters
!   END INTERFACE photolysis_check_parameters
  
    INTERFACE photolysis_constant
       MODULE PROCEDURE photolysis_constant
    END INTERFACE photolysis_constant
 
!   INTERFACE photolysis_simple
!      MODULE PROCEDURE photolysis_simple
!   END INTERFACE photolysis_simple
!
!   INTERFACE photolysis_fastj
!      MODULE PROCEDURE photolysis_fastj
!   END INTERFACE photolysis_fastj
!
    INTERFACE photolysis_control
       MODULE PROCEDURE photolysis_control
    END INTERFACE photolysis_control

!   INTERFACE photolysis_header
!      MODULE PROCEDURE photolysis_header
!   END INTERFACE photolysis_header 
! 
    INTERFACE photolysis_init
       MODULE PROCEDURE photolysis_init
    END INTERFACE photolysis_init

!   INTERFACE photolysis_parin
!      MODULE PROCEDURE photolysis_parin
!   END INTERFACE photolysis_parin
    
!   INTERFACE photolysis_read_restart_data
!      MODULE PROCEDURE photolysis_read_restart_data
!   END INTERFACE photolysis_read_restart_data

!   INTERFACE photolysis_last_actions
!      MODULE PROCEDURE photolysis_last_actions
!   END INTERFACE photolysis_last_actions

    SAVE

    PRIVATE

!
!-- Public functions / NEEDS SORTING
!   PUBLIC photolysis_init
!          photolysis_check_parameters, photolysis_control,                      &
!          photolysis_header, photolysis_init, photolysis_parin !,                 &
!          photolysis_define_netcdf_grid, photolysis_last_actions,               &
!          photolysis_read_restart_data, photolysis_data_output_mask
 !   PUBLIC  photolysis_control
    
!

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the photolysis schemes
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_control
 
 
       IMPLICIT NONE


       SELECT CASE ( TRIM( photolysis_scheme ) )

          CASE ( 'constant' )
             CALL photolysis_constant
          
!         CASE ( 'simple' )
!            CALL photolysis_simple
       
!         CASE ( 'fastj' )
!            CALL photolysis_fastj

          CASE DEFAULT

       END SELECT


    END SUBROUTINE photolysis_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the photolysis model
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_init
    
       IMPLICIT NONE

!
!--    Fix net photolysis in case of radiation_scheme = 'constant'
       IF ( photolysis_scheme == 'constant' )  THEN
!
!--    Calculate orbital constants
       ELSE
!         decl_1 = SIN(23.45_wp * pi / 180.0_wp)
!         decl_2 = 2.0_wp * pi / 365.0_wp
!         decl_3 = decl_2 * 81.0_wp
!         lat    = phi * pi / 180.0_wp
!         lon    = lambda * pi / 180.0_wp
       ENDIF

!      IF ( photolysis_scheme == 'simple'  .OR.                              &
!           photolysis_scheme == 'constant')  THEN
!
!      ELSE
!      ENDIF
!
!         ALLOCATE ( alpha(nysg:nyng,nxlg:nxrg) )
!
!--    Perform user actions if required
!      CALL user_init_photolysis

!
!--    Calculate radiative fluxes at model start
       IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

          SELECT CASE ( photolysis_scheme )
!            CASE ( 'fastj' )
!               CALL photolysis_fastj
!            CASE ( 'simple' )
!               CALL photolysis_simple
             CASE ( 'constant' )
                CALL photolysis_constant
             CASE DEFAULT
          END SELECT

       ENDIF

       RETURN

    END SUBROUTINE photolysis_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_constant


       IMPLICIT NONE

       INTEGER(iwp) :: i, j, k   !< loop indices
       INTEGER(iwp) :: iphot,i_avail     !< loop indix for photolysis reaction
!      REAL(wp)     :: exn,   &  !< Exner functions at surface
!                      pt1       !< potential temperature at first grid level
       REAL(wp)     :: cosz     ! cosine of Zenith angle - zenith angle should be input later


       REAL(wp),DIMENSION(:,:,:), ALLOCATABLE                   :: phot_const

       CHARACTER(LEN=10), PARAMETER, DIMENSION(25) :: names =  (/                      &
                     'J_O31D    ','J_O33P    ','J_NO2     ','J_HNO3    ','J_RCHO    ', &
                     'J         ','J         ','J         ','J         ','J         ', &
                     'J         ','J         ','J         ','J         ','J         ', &
                     'J         ','J         ','J         ','J         ','J         ', &
                     'J         ','J         ','J         ','J         ','J         ' /)
! Photolysis frequency at zenith angle 0
       REAL(wp), PARAMETER, DIMENSION(25) :: phot0 =  (/                               &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-3_wp, 0.000E00_wp, 3.734E-05_wp, &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-3_wp, 0.000E00_wp, 3.734E-05_wp, &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-3_wp, 0.000E00_wp, 3.734E-05_wp, &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-3_wp, 0.000E00_wp, 3.734E-05_wp, &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-3_wp, 0.000E00_wp, 3.734E-05_wp /)

       ALLOCATE ( phot_const(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!
       cosz = 0.7_wp    ! 45 deg  Should come from Namelist

!--    Prescribe fixed photolysis frequencies  [1/s]
       DO i_avail = 1,25   
          phot_const(:,:,:) = phot0(i_avail) * cosz 

       DO iphot = 1,nphot
        IF ( names(i_avail) == trim(PHOT_NAMES(iphot) ) ) then
         if(myid == 0)  write(06,*) iphot, i_avail,PHOT_NAMES(iphot),names(i_avail)
          DO i = nxlg, nxrg
             DO j = nysg, nyng

                phot_frequen(iphot)%freq(nzb_s_inner(j,i)+1:nzt,j,i) =    &
                       phot_const(nzb_s_inner(j,i)+1:nzt,j,i)

 !              phot_frequen(iphot)%name = PHOT_NAMES(iphot)

             ENDDO
          ENDDO
         if(myid == 0)  write(06,*) phot_frequen(iphot)%freq(1,5,5)
       ENDIF
       ENDDO
       ENDDO

    END SUBROUTINE photolysis_constant

!
 END MODULE kchem_photolysis
