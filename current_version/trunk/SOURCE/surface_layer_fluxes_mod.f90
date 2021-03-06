!> @file surface_layer_fluxes_mod.f90
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
! Copyright 1997-2017 Leibniz Universitaet Hannover
!
!------------------------------------------------------------------------------!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: surface_layer_fluxes_mod.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2382 2017-09-01 12:20:53Z basit
! replaced 'chemistry' with 'air_chemistry' 
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! air_chemistry 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Account for urban surface model in computation of vertical kinematic heatflux
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! Minor bug, declaration of look-up index as INTEGER
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
! 
! 1929 2016-06-09 16:25:25Z suehring
! Bugfix: avoid segmentation fault in case one grid point is horizontally 
! completely surrounded by topography
! 
! 1920 2016-05-30 10:50:15Z suehring
! Avoid segmentation fault (see change in 1915) by different initialization of 
! us instead of adding a very small number in the denominator
!
! 1915 2016-05-27 11:05:02Z suehring
! Bugfix: avoid segmentation fault in case of most_method = 'circular' at first
! timestep
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1788 2016-03-10 11:01:04Z maronga
! Added parameter z0q which replaces z0h in the similarity functions for
! humidity.
! Syntax layout improved.
! 
! 1757 2016-02-22 15:49:32Z maronga
! Minor fixes.
!
! 1749 2016-02-09 12:19:56Z raasch
! further OpenACC adjustments
!
! 1747 2016-02-08 12:25:53Z raasch
! adjustments for OpenACC usage
!
! 1709 2015-11-04 14:47:01Z maronga
! Bugfix: division by zero could occur when calculating rib at low wind speeds
! Bugfix: calculation of uv_total for neutral = .T., initial value for ol for
! neutral = .T.
! 
! 1705 2015-11-02 14:28:56Z maronga
! Typo removed
!
! 1697 2015-10-28 17:14:10Z raasch
! FORTRAN and OpenMP errors removed
!
! 1696 2015-10-27 10:03:34Z maronga
! Modularized and completely re-written version of prandtl_fluxes.f90. In the
! course of the re-writing two additional methods have been implemented. See
! updated description.
!
! 1551 2015-03-03 14:18:16Z maronga
! Removed land surface model part. The surface fluxes are now always calculated
! within prandtl_fluxes, based on the given surface temperature/humidity (which 
! is either provided by the land surface model, by large scale forcing data, or
! directly prescribed by the user.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Adapted for land surface model
! 
! 1494 2014-11-21 17:14:03Z maronga
! Bugfixes: qs is now calculated before calculation of Rif. calculation of
! buoyancy flux in Rif corrected (added missing humidity term), allow use of 
! topography for coupled runs (not tested)
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! Bugfix: calculation of turbulent fluxes of rain water content (qrsws) and rain 
! drop concentration (nrsws) added
! 
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1257 2013-11-08 15:18:40Z raasch
! openACC "kernels do" replaced by "kernels loop", "loop independent" added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! OpenACC statements added
!
! 978 2012-08-09 08:28:32Z fricke
! roughness length for scalar quantities z0h added
!
! Revision 1.1  1998/01/23 10:06:06  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diagnostic computation of vertical fluxes in the constant flux layer from the
!> values of the variables at grid point k=1. Three different methods are
!> available:
!> 1) the "old" version (most_method = 'circular') which is fast, but inaccurate
!> 2) a Newton iteration method (most_method = 'newton'), which is accurate, but
!>    slower
!> 3) a method using a lookup table which is fast and accurate. Note, however,
!>    that this method cannot be used in case of roughness heterogeneity
!>
!> @todo (re)move large_scale_forcing actions
!> @todo check/optimize OpenMP and OpenACC directives
!------------------------------------------------------------------------------!
 MODULE surface_layer_fluxes_mod

    USE arrays_3d,                                                             &
        ONLY:  e, kh, nr, nrs, nrsws, ol, pt, q, ql, qr, qrs, qrsws, qs, qsws, &
               rs, rss, rssws, s, shf, ss, ssws, ts, u, us, usws, v, vpt, vsws, zu, zw, z0,    &  !bK added rs, rss, rssws
               z0h, z0q, drho_air_zw, rho_air_zw

    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, pt_d_t

    USE constants,                                                             &
        ONLY:  pi

    USE cpulog

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, cloud_physics, constant_chemistryflux,           &
               constant_heatflux, constant_scalarflux,                         &     !bK added air_chemistry, constant_chemistryflux
               constant_waterflux, coupling_mode, g, humidity, ibc_e_b,        &
               ibc_pt_b, initializing_actions, kappa,                          &
               intermediate_timestep_count,                                    &
               intermediate_timestep_count_max, large_scale_forcing, lsf_surf, &
               message_string, microphysics_seifert, most_method, neutral,     &
               passive_scalar, pt_surface, q_surface, run_coupled,             &
               surface_pressure, simulated_time, terminate_run,                &
               urban_surface,                                                  &
               zeta_max, zeta_min 

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb_s_inner,        &
               nzb_u_inner, nzb_v_inner

    USE kinds

    USE pegrid

    USE land_surface_model_mod,                                                &
        ONLY:  land_surface, skip_time_do_lsm

        

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< loop index x direction
    INTEGER(iwp) ::  j              !< loop index y direction
    INTEGER(iwp) ::  k              !< loop index z direction
    INTEGER(iwp) ::  l_bnd  = 7500  !< Lookup table index of the last time step

    INTEGER(iwp), PARAMETER     :: num_steps = 15000  !< number of steps in the lookup table

    LOGICAL      ::  coupled_run  !< Flag for coupled atmosphere-ocean runs

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pt1,      & !< Potential temperature at first grid level (required for cloud_physics = .T.)
                                             qv1,      & !< Specific humidity at first grid level (required for cloud_physics = .T.)
                                             uv_total    !< Total velocity at first grid level

    REAL(wp), DIMENSION(0:num_steps-1) :: rib_tab,  & !< Lookup table bulk Richardson number
                                          ol_tab      !< Lookup table values of L

    REAL(wp)     ::  e_s,               & !< Saturation water vapor pressure
                     ol_max = 1.0E6_wp, & !< Maximum Obukhov length
                     rib_max,           & !< Maximum Richardson number in lookup table
                     rib_min,           & !< Minimum Richardson number in lookup table
                     z_mo                 !< Height of the constant flux layer where MOST is assumed


    SAVE

    PRIVATE

    PUBLIC init_surface_layer_fluxes, pt1, qv1, surface_layer_fluxes, uv_total

    INTERFACE init_surface_layer_fluxes
       MODULE PROCEDURE init_surface_layer_fluxes
    END INTERFACE init_surface_layer_fluxes

    INTERFACE surface_layer_fluxes
       MODULE PROCEDURE surface_layer_fluxes
    END INTERFACE surface_layer_fluxes


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Main routine to compute the surface fluxes
!------------------------------------------------------------------------------!
    SUBROUTINE surface_layer_fluxes

       IMPLICIT NONE

!
!--    In case cloud physics is used, it is required to derive potential
!--    temperature and specific humidity at first grid level from the fields pt 
!--    and q
       IF ( cloud_physics )  THEN
          CALL calc_pt_q
       ENDIF

!
!--    First, calculate the new Obukhov length, then new friction velocity,
!--    followed by the new scaling parameters (th*, q*, etc.), and the new
!--    surface fluxes if required. The old routine ("circular") requires a 
!--    different order of calls as the scaling parameters from the previous time
!--    steps are used to calculate the Obukhov length

!
!--    Depending on setting of most_method use the "old" routine
       IF ( most_method == 'circular' )  THEN

          CALL calc_scaling_parameters

          CALL calc_uv_total

          IF ( .NOT. neutral )  THEN
             CALL calc_ol
          ENDIF

          CALL calc_us

          CALL calc_surface_fluxes

!
!--    Use either Newton iteration or a lookup table for the bulk Richardson
!--    number to calculate the Obukhov length 
       ELSEIF ( most_method == 'newton'  .OR.  most_method == 'lookup' )  THEN

          CALL calc_uv_total

          IF ( .NOT. neutral )  THEN
             CALL calc_ol
          ENDIF

          CALL calc_us

          CALL calc_scaling_parameters

          CALL calc_surface_fluxes

       ENDIF

    END SUBROUTINE surface_layer_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing actions for the surface layer routine. Basically, this involves
!> the preparation of a lookup table for the the bulk Richardson number vs 
!> Obukhov length L when using the lookup table method.
!------------------------------------------------------------------------------!
    SUBROUTINE init_surface_layer_fluxes

       IMPLICIT NONE

       INTEGER(iwp) :: l,          & !< Index for loop to create lookup table
                       num_steps_n   !< Number of non-stretched zeta steps

       LOGICAL :: terminate_run_l = .FALSE.    !< Flag to terminate run (global)

       REAL(wp), PARAMETER ::  zeta_stretch = -10.0_wp !< Start of stretching in the free convection limit
                               
       REAL(wp), DIMENSION(:), ALLOCATABLE :: zeta_tmp


       REAL(wp) :: zeta_step,            & !< Increment of zeta
                   regr      = 1.01_wp,  & !< Stretching factor of zeta_step in the free convection limit
                   regr_old  = 1.0E9_wp, & !< Stretching factor of last iteration step
                   z0h_min   = 0.0_wp,   & !< Minimum value of z0h to create table
                   z0_min    = 0.0_wp      !< Minimum value of z0 to create table
!
!--    When cloud physics is used, arrays for storing potential temperature and
!--    specific humidity at first grid level are required
       IF ( cloud_physics )  THEN
          ALLOCATE ( pt1(nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( qv1(nysg:nyng,nxlg:nxrg) )
       ENDIF

!
!--    Allocate field for storing the horizontal velocity
       ALLOCATE ( uv_total(nysg:nyng,nxlg:nxrg) )


!
!--    In case of runs with neutral statification, set Obukhov length to a
!--    large value
       IF ( neutral ) ol = 1.0E10_wp

       IF ( most_method == 'lookup' )  THEN

!
!--       Check for roughness heterogeneity. In that case terminate run and
!--       inform user
          IF ( MINVAL( z0h ) /= MAXVAL( z0h )  .OR.                            &
               MINVAL( z0  ) /= MAXVAL( z0  ) )  THEN
             terminate_run_l = .TRUE.
          ENDIF

#if defined( __parallel )
!
!--       Make a logical OR for all processes. Force termiation of model if result
!--       is TRUE
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( terminate_run_l, terminate_run, 1, MPI_LOGICAL,  &
                              MPI_LOR, comm2d, ierr )
#else
          terminate_run = terminate_run_l
#endif

          IF ( terminate_run )  THEN
             message_string = 'most_method = "lookup" cannot be used in ' //   &
                              'combination with a prescribed roughness '  //   &
                              'heterogeneity'
             CALL message( 'surface_layer_fluxes', 'PA0417', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE(  zeta_tmp(0:num_steps-1) )

!
!--       Use the lowest possible value for z_mo
          k    = MINVAL(nzb_s_inner)
          z_mo = zu(k+1) - zw(k)

!
!--       Calculate z/L range from zeta_stretch to zeta_max using 90% of the
!--       available steps (num_steps). The calculation is done with negative
!--       values of zeta in order to simplify the stretching in the free
!--       convection limit for the remaining 10% of steps.
          zeta_tmp(0) = - zeta_max
          num_steps_n = ( num_steps * 9 / 10 ) - 1
          zeta_step   = (zeta_max - zeta_stretch) / REAL(num_steps_n)

          DO l = 1, num_steps_n
             zeta_tmp(l) = zeta_tmp(l-1) + zeta_step
          ENDDO

!
!--       Calculate stretching factor for the free convection range
          DO  WHILE ( ABS( (regr-regr_old) / regr_old ) > 1.0E-10_wp )
             regr_old = regr
             regr = ( 1.0_wp - ( -zeta_min / zeta_step ) * ( 1.0_wp - regr )   &
                    )**( 10.0_wp / REAL(num_steps) )
          ENDDO

!
!--       Calculate z/L range from zeta_min to zeta_stretch
          DO l = num_steps_n+1, num_steps-1
             zeta_tmp(l) = zeta_tmp(l-1) + zeta_step
             zeta_step = zeta_step * regr
          ENDDO

!
!--       Save roughness lengths to temporary variables
          z0h_min = z0h(nys,nxl)
          z0_min  = z0(nys,nxl)
          
!
!--       Calculate lookup table for the Richardson number versus Obukhov length
!--       The Richardson number (rib) is defined depending on the choice of 
!--       boundary conditions for temperature
          IF ( ibc_pt_b == 1 )  THEN
             DO l = 0, num_steps-1
                ol_tab(l)  = - z_mo / zeta_tmp(num_steps-1-l)
                rib_tab(l) = z_mo / ol_tab(l)  / ( LOG( z_mo / z0_min )        &
                                                - psi_m( z_mo / ol_tab(l) )    &
                                                + psi_m( z0_min / ol_tab(l) )  &
                                                  )**3
             ENDDO  
          ELSE
             DO l = 0, num_steps-1
                ol_tab(l)  = - z_mo / zeta_tmp(num_steps-1-l)
                rib_tab(l) = z_mo / ol_tab(l)  * ( LOG( z_mo / z0h_min )       &
                                              - psi_h( z_mo / ol_tab(l) )      &
                                              + psi_h( z0h_min / ol_tab(l) )   &
                                            )                                  &
                                          / ( LOG( z_mo / z0_min )             &
                                              - psi_m( z_mo / ol_tab(l) )      &
                                              + psi_m( z0_min / ol_tab(l) )    &
                                            )**2
             ENDDO
          ENDIF

!
!--       Determine minimum values of rib in the lookup table. Set upper limit 
!--       to critical Richardson number (0.25)
          rib_min  = MINVAL(rib_tab)
          rib_max  = 0.25 !MAXVAL(rib_tab)

          DEALLOCATE( zeta_tmp )
       ENDIF

    END SUBROUTINE init_surface_layer_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface). This is required by all methods
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uv_total

       IMPLICIT NONE


       !$OMP PARALLEL DO PRIVATE( k )
       !$acc kernels loop present( nzb_s_inner, u, uv_total, v ) private( j, k )
       DO  i = nxl, nxr
          DO  j = nys, nyn

             k   = nzb_s_inner(j,i)
             uv_total(j,i) = SQRT( ( 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1)      &
                                         - u(k,j,i)   - u(k,j,i+1) ) )**2 +    &
                              ( 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i)           &
                                         - v(k,j,i)   - v(k,j+1,i) ) )**2 )

!
!--          For too small values of the local wind, MOST does not work. A
!--          threshold value is thus set if required
!            uv_total(j,i) = MAX(0.01_wp,uv_total(j,i))

          ENDDO
       ENDDO

!
!--    Values of uv_total need to be exchanged at the ghost boundaries
       !$acc update host( uv_total )
       CALL exchange_horiz_2d( uv_total )
       !$acc update device( uv_total )

    END SUBROUTINE calc_uv_total


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L) and Richardson flux number (z/L)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_ol

       IMPLICIT NONE

       INTEGER(iwp) :: iter,  & !< Newton iteration step
                       l        !< look index

       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) :: rib !< Bulk Richardson number

       REAL(wp)     :: f,      & !< Function for Newton iteration: f = Ri - [...]/[...]^2 = 0
                       f_d_ol, & !< Derivative of f
                       ol_l,   & !< Lower bound of L for Newton iteration
                       ol_m,   & !< Previous value of L for Newton iteration
                       ol_old, & !< Previous time step value of L
                       ol_u      !< Upper bound of L for Newton iteration

       IF ( TRIM( most_method ) /= 'circular' )  THEN
      
          !$acc data present( nzb_s_inner, pt, q, qsws, rib, shf, uv_total, vpt, zu, zw )

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !$acc kernels loop private( j, k, z_mo )
          DO  i = nxl, nxr
             DO  j = nys, nyn

                k   = nzb_s_inner(j,i)
                z_mo = zu(k+1) - zw(k)

!
!--             Evaluate bulk Richardson number (calculation depends on
!--             definition based on setting of boundary conditions
                IF ( ibc_pt_b /= 1 )  THEN
                   IF ( humidity )  THEN
                      rib(j,i) = g * z_mo * ( vpt(k+1,j,i) - vpt(k,j,i) )      &
                           / ( uv_total(j,i)**2 * vpt(k+1,j,i) + 1.0E-20_wp )
                   ELSE
                      rib(j,i) = g * z_mo * ( pt(k+1,j,i) - pt(k,j,i) )        &
                           / ( uv_total(j,i)**2 * pt(k+1,j,i)  + 1.0E-20_wp )
                   ENDIF      
                ELSE
!
!--                When using Neumann boundary conditions, the buoyancy flux
!--                is required but cannot be calculated at the surface, as pt
!--                and q are not known at the surface. Hence the values at 
!--                first grid level are used to estimate the buoyancy flux
                   IF ( humidity )  THEN
                      rib(j,i) = - g * z_mo * ( ( 1.0_wp + 0.61_wp             &
                                 * q(k+1,j,i) ) * shf(j,i) + 0.61_wp           &
                                 * pt(k+1,j,i) * qsws(j,i) ) * drho_air_zw(k)  &
                                 / ( uv_total(j,i)**3 * vpt(k+1,j,i) * kappa**2&
                                 + 1.0E-20_wp)
                   ELSE
                      rib(j,i) = - g * z_mo * shf(j,i) * drho_air_zw(k)        &
                           / ( uv_total(j,i)**3 * pt(k+1,j,i) * kappa**2       &
                           + 1.0E-20_wp )
                   ENDIF
                ENDIF  
     
             ENDDO
          ENDDO 
          !$acc end data

       ENDIF

!
!--    Calculate the Obukhov length either using a Newton iteration
!--    method, via a lookup table, or using the old circular way
       IF ( TRIM( most_method ) == 'newton' )  THEN

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !# WARNING: does not work on GPU so far because of DO-loop with
          !#          undetermined iterations
          !!!!!!$acc kernels loop
          DO  i = nxl, nxr
             DO  j = nys, nyn

                k   = nzb_s_inner(j,i)
                z_mo = zu(k+1) - zw(k)

!
!--             Store current value in case the Newton iteration fails
                ol_old = ol(j,i)

!
!--             Ensure that the bulk Richardson number and the Obukhov 
!--             lengtH have the same sign
                IF ( rib(j,i) * ol(j,i) < 0.0_wp  .OR.                         &
                     ABS( ol(j,i) ) == ol_max )  THEN
                   IF ( rib(j,i) > 0.0_wp ) ol(j,i) =  0.01_wp
                   IF ( rib(j,i) < 0.0_wp ) ol(j,i) = -0.01_wp
                ENDIF
!
!--             Iteration to find Obukhov length
                iter = 0
                DO
                   iter = iter + 1
!
!--                In case of divergence, use the value of the previous time step
                   IF ( iter > 1000 )  THEN
                      ol(j,i) = ol_old
                      EXIT
                   ENDIF

                   ol_m = ol(j,i)
                   ol_l = ol_m - 0.001_wp * ol_m
                   ol_u = ol_m + 0.001_wp * ol_m


                   IF ( ibc_pt_b /= 1 )  THEN
!
!--                   Calculate f = Ri - [...]/[...]^2 = 0
                      f = rib(j,i) - ( z_mo / ol_m ) * ( LOG( z_mo / z0h(j,i) )&
                                                    - psi_h( z_mo / ol_m )     &
                                                    + psi_h( z0h(j,i) / ol_m ) &
                                                   )                           &
                                                 / ( LOG( z_mo / z0(j,i) )     &
                                                    - psi_m( z_mo / ol_m )     &
                                                    + psi_m( z0(j,i) / ol_m )  &
                                                    )**2

!
!--                    Calculate df/dL
                       f_d_ol = ( - ( z_mo / ol_u ) * ( LOG( z_mo / z0h(j,i) ) &
                                                   - psi_h( z_mo / ol_u )      &
                                                   + psi_h( z0h(j,i) / ol_u )  &
                                                 )                             &
                                               / ( LOG( z_mo / z0(j,i) )       &
                                                   - psi_m( z_mo / ol_u )      &
                                                   + psi_m( z0(j,i) / ol_u )   &
                                                 )**2                          &
                              + ( z_mo / ol_l ) * ( LOG( z_mo / z0h(j,i) )     &
                                                   - psi_h( z_mo / ol_l )      &
                                                   + psi_h( z0h(j,i) / ol_l )  &
                                                 )                             &
                                               / ( LOG( z_mo / z0(j,i) )       &
                                                   - psi_m( z_mo / ol_l )      &
                                                   + psi_m( z0(j,i) / ol_l )   &
                                                 )**2                          &
                                ) / ( ol_u - ol_l )
                   ELSE
!
!--                   Calculate f = Ri - 1 /[...]^3 = 0
                      f = rib(j,i) - ( z_mo / ol_m ) / ( LOG( z_mo / z0(j,i) )&
                                                    - psi_m( z_mo / ol_m )    &
                                                    + psi_m( z0(j,i) / ol_m ) &
                                                       )**3

!
!--                   Calculate df/dL
                      f_d_ol = ( - ( z_mo / ol_u ) / ( LOG( z_mo / z0(j,i) )  &
                                                   - psi_m( z_mo / ol_u )     &
                                                   + psi_m( z0(j,i) / ol_u )  &
                                                 )**3                         &
                              + ( z_mo / ol_l ) / ( LOG( z_mo / z0(j,i) )     &
                                                   - psi_m( z_mo / ol_l )     &
                                                   + psi_m( z0(j,i) / ol_l )  &
                                                 )**3                         &
                                     ) / ( ol_u - ol_l )
                   ENDIF
!
!--                Calculate new L
                   ol(j,i) = ol_m - f / f_d_ol

!
!--                Ensure that the bulk Richardson number and the Obukhov 
!--                length have the same sign and ensure convergence.
                   IF ( ol(j,i) * ol_m < 0.0_wp )  ol(j,i) = ol_m * 0.5_wp

!
!--                If unrealistic value occurs, set L to the maximum
!--                value that is allowed
                   IF ( ABS( ol(j,i) ) > ol_max )  THEN
                      ol(j,i) = ol_max
                      EXIT
                   ENDIF
!
!--                Check for convergence
                   IF ( ABS( ( ol(j,i) - ol_m ) / ol(j,i) ) < 1.0E-4_wp )  THEN
                      EXIT
                   ELSE
                      CYCLE
                   ENDIF

                ENDDO
                       
             ENDDO
          ENDDO

       ELSEIF ( TRIM( most_method ) == 'lookup' )  THEN

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !# WARNING: does not work on GPU so far because of DO  WHILE construct
          !!!!!!$acc kernels loop
          DO  i = nxl, nxr
             DO  j = nys, nyn

!
!--             If the bulk Richardson number is outside the range of the lookup
!--             table, set it to the exceeding threshold value
                IF ( rib(j,i) < rib_min )  rib(j,i) = rib_min
                IF ( rib(j,i) > rib_max )  rib(j,i) = rib_max

!
!--             Find the correct index bounds for linear interpolation. As the
!--             Richardson number will not differ very much from time step to
!--             time step , use the index from the last step and search in the 
!--             correct direction
                l = l_bnd
                IF ( rib_tab(l) - rib(j,i) > 0.0_wp )  THEN
                   DO  WHILE ( rib_tab(l-1) - rib(j,i) > 0.0_wp  .AND.  l > 0 )
                      l = l-1
                   ENDDO
                ELSE
                   DO  WHILE ( rib_tab(l) - rib(j,i) < 0.0_wp                  &
                              .AND.  l < num_steps-1 )
                      l = l+1
                   ENDDO
                ENDIF
                l_bnd = l

!
!--             Linear interpolation to find the correct value of z/L
                ol(j,i) = ( ol_tab(l-1) + ( ol_tab(l) - ol_tab(l-1) )          &
                            / (  rib_tab(l) - rib_tab(l-1) )                   &
                            * ( rib(j,i) - rib_tab(l-1) ) )

             ENDDO
          ENDDO

       ELSEIF ( TRIM( most_method ) == 'circular' )  THEN

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !$acc kernels loop present( nzb_s_inner, ol, pt, pt1, q, ql, qs, qv1, ts, us, vpt, zu, zw ) private( j, k, z_mo )
          DO  i = nxl, nxr
             DO  j = nys, nyn

                k   = nzb_s_inner(j,i)
                z_mo = zu(k+1) - zw(k)

                IF ( .NOT. humidity )  THEN
                   ol(j,i) =  ( pt(k+1,j,i) *  us(j,i)**2 ) / ( kappa * g      &
                              * ts(j,i) + 1E-30_wp )
                ELSEIF ( cloud_physics )  THEN

                   ol(j,i) =  ( vpt(k+1,j,i) * us(j,i)**2 ) / ( kappa * g      &
                              * ( ts(j,i) + 0.61_wp * pt1(j,i) * qs(j,i)       &
                              + 0.61_wp * qv1(j,i) * ts(j,i) - ts(j,i)         &
                              * ql(k+1,j,i) ) + 1E-30_wp )
                ELSE
                   ol(j,i) =  ( vpt(k+1,j,i) *  us(j,i)**2 ) / ( kappa * g     &
                              * ( ts(j,i) + 0.61_wp * pt(k+1,j,i) * qs(j,i)    &
                                  + 0.61_wp * q(k+1,j,i) * ts(j,i) ) + 1E-30_wp )
                ENDIF
!
!--             Limit the value range of the Obukhov length.
!--             This is necessary for very small velocities (u,v --> 0), because
!--             the absolute value of ol can then become very small, which in
!--             consequence would result in very large shear stresses and very
!--             small momentum fluxes (both are generally unrealistic).
                IF ( ( z_mo / ( ol(j,i) + 1E-30_wp ) ) < zeta_min )            &
                   ol(j,i) = z_mo / zeta_min
                IF ( ( z_mo / ( ol(j,i) + 1E-30_wp ) ) > zeta_max )            &
                   ol(j,i) = z_mo / zeta_max
             ENDDO
          ENDDO

       ENDIF

!
!--    Values of ol at ghost point locations are needed for the evaluation
!--    of usws and vsws.
       !$acc update host( ol )
       CALL exchange_horiz_2d( ol )
       !$acc update device( ol )

    END SUBROUTINE calc_ol

!
!-- Calculate friction velocity u*
    SUBROUTINE calc_us

       IMPLICIT NONE

       !$OMP PARALLEL DO PRIVATE( k, z_mo )
       !$acc kernels loop present( nzb_s_inner, ol, us, uv_total, zu, zw, z0 ) private( j, k, z_mo )
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng

             k   = nzb_s_inner(j,i)+1
             z_mo = zu(k+1) - zw(k)

!
!--          Compute u* at the scalars' grid points
             us(j,i) = kappa * uv_total(j,i) / ( LOG( z_mo / z0(j,i) )         &
                                          - psi_m( z_mo / ol(j,i) )            &
                                          + psi_m( z0(j,i) / ol(j,i) ) )
          ENDDO
       ENDDO

    END SUBROUTINE calc_us

!
!-- Calculate potential temperature and specific humidity at first grid level
    SUBROUTINE calc_pt_q

       IMPLICIT NONE

       !$acc kernels loop present( nzb_s_inner, pt, pt1, pt_d_t, q, ql, qv1 ) private( j, k )
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             k   = nzb_s_inner(j,i)+1
             pt1(j,i) = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
             qv1(j,i) = q(k,j,i) - ql(k,j,i)
          ENDDO
       ENDDO

    END SUBROUTINE calc_pt_q

!
!-- Calculate the other MOST scaling parameters theta*, q*, (qr*, nr*)
    SUBROUTINE calc_scaling_parameters

       IMPLICIT NONE

!
!--    Data information for accelerators
       !$acc data present( e, nrsws, nzb_u_inner, nzb_v_inner, nzb_s_inner, pt )  &
       !$acc      present( q, qs, qsws, qrsws, shf, ts, u, us, usws, v )     &
       !$acc      present( vpt, vsws, zu, zw, z0, z0h )
! 
!--    Compute theta*
       IF ( constant_heatflux )  THEN

!
!--       For a given heat flux in the surface layer:
          !$OMP PARALLEL DO
          !$acc kernels loop private( j, k )
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                k   = nzb_s_inner(j,i)
                ts(j,i) = -shf(j,i) * drho_air_zw(k) / ( us(j,i) + 1E-30_wp )
!
!--             ts must be limited, because otherwise overflow may occur in case
!--             of us=0 when computing ol further below
                IF ( ts(j,i) < -1.05E5_wp )  ts(j,i) = -1.0E5_wp
                IF ( ts(j,i) >   1.0E5_wp )  ts(j,i) =  1.0E5_wp
             ENDDO
          ENDDO

       ELSE
!
!--       For a given surface temperature:
          IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
             !$OMP PARALLEL DO
             !$acc kernels loop private( j, k )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   k = nzb_s_inner(j,i)
                   pt(k,j,i) = pt_surface
                ENDDO
             ENDDO
          ENDIF

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !$acc kernels loop present( nzb_s_inner, ol, pt, pt1, ts, zu, zw, z0h ) private( j, k, z_mo )
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng

                k   = nzb_s_inner(j,i)
                z_mo = zu(k+1) - zw(k)

                IF ( cloud_physics )  THEN
                   ts(j,i) = kappa * ( pt1(j,i) - pt(k,j,i) )                  &
                                     / ( LOG( z_mo / z0h(j,i) )                &
                                         - psi_h( z_mo / ol(j,i) )             &
                                         + psi_h( z0h(j,i) / ol(j,i) ) )
                ELSE
                   ts(j,i) = kappa * ( pt(k+1,j,i) - pt(k,j,i) )               &
                                     / ( LOG( z_mo / z0h(j,i) )                &
                                         - psi_h( z_mo / ol(j,i) )             &
                                         + psi_h( z0h(j,i) / ol(j,i) ) )
                ENDIF

             ENDDO
          ENDDO
       ENDIF

!
!--    If required compute q*
       IF ( humidity )  THEN
          IF ( constant_waterflux )  THEN
!
!--          For a given water flux in the surface layer
             !$OMP PARALLEL DO
             !$acc kernels loop private( j )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   k   = nzb_s_inner(j,i)
                   qs(j,i) = -qsws(j,i) * drho_air_zw(k) / ( us(j,i) + 1E-30_wp )
                ENDDO
             ENDDO

          ELSE
             coupled_run = ( coupling_mode == 'atmosphere_to_ocean'  .AND.     &
                             run_coupled )

             IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
                !$OMP PARALLEL DO
                !$acc kernels loop private( j, k )
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      k = nzb_s_inner(j,i)
                      q(k,j,i) = q_surface
                   ENDDO
                ENDDO
             ENDIF

             !$OMP PARALLEL DO PRIVATE( e_s, k, z_mo )
             !$acc kernels loop independent present( nzb_s_inner, ol, pt, q, qs, qv1, zu, zw, z0q ) private( e_s, j, k, z_mo )
             DO  i = nxlg, nxrg
                !$acc loop independent
                DO  j = nysg, nyng

                   k   = nzb_s_inner(j,i)
                   z_mo = zu(k+1) - zw(k)

!
!--                Assume saturation for atmosphere coupled to ocean (but not
!--                in case of precursor runs)
                   IF ( coupled_run )  THEN
                      e_s = 6.1_wp * &
                              EXP( 0.07_wp * ( MIN(pt(k,j,i),pt(k+1,j,i))      &
                                               - 273.15_wp ) )
                      q(k,j,i) = 0.622_wp * e_s / ( surface_pressure - e_s )
                   ENDIF

                   IF ( cloud_physics )  THEN
                      qs(j,i) = kappa * ( qv1(j,i) - q(k,j,i) )                &
                                        / ( LOG( z_mo / z0q(j,i) )             &
                                            - psi_h( z_mo / ol(j,i) )          &
                                            + psi_h( z0q(j,i) / ol(j,i) ) )

                   ELSE
                      qs(j,i) = kappa * ( q(k+1,j,i) - q(k,j,i) )              &
                                        / ( LOG( z_mo / z0q(j,i) )             &
                                            - psi_h( z_mo / ol(j,i) )          &
                                            + psi_h( z0q(j,i) / ol(j,i) ) )
                   ENDIF

                ENDDO
             ENDDO
          ENDIF
       ENDIF
       
!
!--    If required compute s*
       IF ( passive_scalar )  THEN
          IF ( constant_scalarflux )  THEN
!
!--          For a given water flux in the surface layer
             !$OMP PARALLEL DO
             !$acc kernels loop private( j )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ss(j,i) = -ssws(j,i) / ( us(j,i) + 1E-30_wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( air_chemistry )  THEN                               !bK added this nested IF block
          IF ( constant_chemistryflux )  THEN
!
!--          For a given water flux in the surface layer
             !$OMP PARALLEL DO
             !$acc kernels loop private( j )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   rss(j,i) = -rssws(j,i) / ( us(j,i) + 1E-30_wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF






!
!--    If required compute qr* and nr*
       IF ( cloud_physics  .AND.  microphysics_seifert )   &
       THEN

          !$OMP PARALLEL DO PRIVATE( k, z_mo )
          !$acc kernels loop independent present( nr, nrs, nzb_s_inner, ol, qr, qrs, zu, zw, z0q ) private( j, k, z_mo )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng

                k   = nzb_s_inner(j,i)
                z_mo = zu(k+1) - zw(k)

                qrs(j,i) = kappa * ( qr(k+1,j,i) - qr(k,j,i) )              &
                                 / ( LOG( z_mo / z0q(j,i) )                 &
                                     - psi_h( z_mo / ol(j,i) )              &
                                     + psi_h( z0q(j,i) / ol(j,i) ) )

                nrs(j,i) = kappa * ( nr(k+1,j,i) - nr(k,j,i) )              &
                                 / ( LOG( z_mo / z0q(j,i) )                 &
                                     - psi_h( z_mo / ol(j,i) )              &
                                     + psi_h( z0q(j,i) / ol(j,i) ) )
             ENDDO
          ENDDO

       ENDIF
       !$acc end data

    END SUBROUTINE calc_scaling_parameters



!
!-- Calculate surface fluxes usws, vsws, shf, qsws, (qrsws, nrsws)
    SUBROUTINE calc_surface_fluxes

       IMPLICIT NONE

       REAL(wp) :: ol_mid !< Grid-interpolated L

!
!--    Compute u'w' for the total model domain.
!--    First compute the corresponding component of u* and square it.
       !$OMP PARALLEL DO PRIVATE( k, ol_mid, z_mo )
       !$acc kernels loop present( nzb_u_inner, ol, u, us, usws, zu, zw, z0 ) private( j, k, z_mo )
       DO  i = nxl, nxr
          DO  j = nys, nyn

             k   = nzb_u_inner(j,i)
             z_mo = zu(k+1) - zw(k)
!
!--          Compute bulk Obukhov length for this point
             ol_mid = 0.5_wp * ( ol(j,i-1) + ol(j,i) )

             IF ( ol_mid == 0.0_wp )  THEN
                ol_mid = MIN(ol(j,i-1), ol(j,i))
             ENDIF

             usws(j,i) = kappa * ( u(k+1,j,i) - u(k,j,i) )                     &
                                 / ( LOG( z_mo / z0(j,i) )                     &
                                     - psi_m( z_mo / ol_mid )                  &
                                     + psi_m( z0(j,i) / ol_mid ) )

             usws(j,i) = -usws(j,i) * 0.5_wp * ( us(j,i-1) + us(j,i) )         &
                                    * rho_air_zw(k)
          ENDDO
       ENDDO

!
!--    Compute v'w' for the total model domain.
!--    First compute the corresponding component of u* and square it.
       !$OMP PARALLEL DO PRIVATE( k, ol_mid, z_mo )
       !$acc kernels loop present( nzb_v_inner, ol, v, us, vsws, zu, zw, z0 ) private( j, k, ol_mid, z_mo )
       DO  i = nxl, nxr
          DO  j = nys, nyn

             k   = nzb_v_inner(j,i)
             z_mo = zu(k+1) - zw(k)
!
!--          Compute bulk Obukhov length for this point
             ol_mid = 0.5_wp * ( ol(j-1,i) + ol(j,i) )

             IF ( ol_mid == 0.0_wp )  THEN
                ol_mid = MIN(ol(j-1,i), ol(j-1,i))
             ENDIF

             vsws(j,i) = kappa * ( v(k+1,j,i) - v(k,j,i) )                     &
                                 / ( LOG( z_mo / z0(j,i) )                     &
                                     - psi_m( z_mo / ol_mid )                  &
                                     + psi_m( z0(j,i) / ol_mid ) )

             vsws(j,i) = -vsws(j,i) * 0.5_wp * ( us(j,i-1) + us(j,i) )         &
                                    * rho_air_zw(k)

          ENDDO
       ENDDO

!
!--    Exchange the boundaries for the momentum fluxes (is this still required?)
       !$acc update host( usws, vsws )
       CALL exchange_horiz_2d( usws )
       CALL exchange_horiz_2d( vsws )
       !$acc update device( usws, vsws )

!
!--    Compute the vertical kinematic heat flux
       IF (  .NOT.  constant_heatflux  .AND.  ( simulated_time <=            &
            skip_time_do_lsm  .OR.  .NOT.  land_surface )  .AND.             &
            .NOT.  urban_surface )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( shf, ts, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                k   = nzb_s_inner(j,i)
                shf(j,i) = -ts(j,i) * us(j,i) * rho_air_zw(k)
             ENDDO
          ENDDO

       ENDIF

!
!--    Compute the vertical water flux
       IF (  .NOT.  constant_waterflux  .AND.  humidity  .AND.                 &
             ( simulated_time <= skip_time_do_lsm                              &
            .OR.  .NOT.  land_surface ) )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( qs, qsws, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                k   = nzb_s_inner(j,i)
                qsws(j,i) = -qs(j,i) * us(j,i) * rho_air_zw(k)
             ENDDO
          ENDDO

       ENDIF
       
!
!--    Compute the vertical scalar flux
       IF (  .NOT.  constant_scalarflux  .AND.  passive_scalar  .AND.          &
             ( simulated_time <= skip_time_do_lsm                              &
            .OR.  .NOT.  land_surface ) )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( qs, qsws, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                ssws(j,i) = -ss(j,i) * us(j,i)
             ENDDO
          ENDDO

       ENDIF       
!
!--    Compute the vertical chemistry flux
       IF (  .NOT.  constant_chemistryflux  .AND.  air_chemistry  .AND.            &        !bK adde this IF block 
             ( simulated_time <= skip_time_do_lsm                              &
            .OR.  .NOT.  land_surface ) )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( qs, qsws, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                rssws(j,i) = -rss(j,i) * us(j,i)
             ENDDO
          ENDDO

       ENDIF       

!
!--    Compute (turbulent) fluxes of rain water content and rain drop conc.
       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( nrs, nrsws, qrs, qrsws, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                qrsws(j,i) = -qrs(j,i) * us(j,i)
                nrsws(j,i) = -nrs(j,i) * us(j,i)
             ENDDO
          ENDDO
       ENDIF

!
!--    Bottom boundary condition for the TKE
       IF ( ibc_e_b == 2 )  THEN
          !$OMP PARALLEL DO
          !$acc kernels loop independent present( e, nzb_s_inner, us )
          DO  i = nxlg, nxrg
             !$acc loop independent
             DO  j = nysg, nyng
                k = nzb_s_inner(j,i)
                e(k+1,j,i) = ( us(j,i) / 0.1_wp )**2
!
!--             As a test: cm = 0.4
!               e(k+1,j,i) = ( us(j,i) / 0.4_wp )**2
                e(k,j,i)   = e(k+1,j,i)
             ENDDO
          ENDDO
       ENDIF

    END SUBROUTINE calc_surface_fluxes


!
!-- Integrated stability function for momentum
    FUNCTION psi_m( zeta ) 
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_m !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( SQRT( 1.0_wp  - 16.0_wp * zeta ) )
          psi_m = pi * 0.5_wp - 2.0_wp * ATAN( x ) + LOG( ( 1.0_wp + x )**2    &
                  * ( 1.0_wp + x**2 ) * 0.125_wp )
       ELSE

          psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta         &
                   - bc_d_d
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_m = - 5.0_wp * zeta

       ENDIF

    END FUNCTION psi_m


!
!-- Integrated stability function for heat and moisture
    FUNCTION psi_h( zeta ) 
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_h !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( 1.0_wp  - 16.0_wp * zeta )
          psi_h = 2.0_wp * LOG( (1.0_wp + x ) / 2.0_wp )
       ELSE
          psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_wp          &
                  + 0.66666666666_wp * a * zeta )**1.5_wp - bc_d_d             &
                  + 1.0_wp
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_h = - 5.0_wp * zeta
       ENDIF

    END FUNCTION psi_h

 END MODULE surface_layer_fluxes_mod
