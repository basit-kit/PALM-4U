!> @file init_1d_model.f90
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
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: init_1d_model.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1960 2016-07-12 16:34:24Z suehring
! Remove passive_scalar from IF-statements, as 1D-scalar profile is effectively
! not used.
! Formatting adjustment
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_flush replaced by FORTRAN statement
!
! 1709 2015-11-04 14:47:01Z maronga
! Set initial time step to 10 s to avoid instability of the 1d model for small
! grid spacings
!
! 1697 2015-10-28 17:14:10Z raasch
! small E- and F-FORMAT changes to avoid informative compiler messages about
! insufficient field width
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer. rif is replaced by ol and zeta.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
! 
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! adjustment of mixing length to the Prandtl mixing length at first grid point
! above ground removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! roughness length for scalar quantities z0h1d added
!
! Revision 1.1  1998/03/09 16:22:10  raasch
! Initial revision
!
!
! Description:
! ------------
!> 1D-model to initialize the 3D-arrays.
!> The temperature profile is set as steady and a corresponding steady solution 
!> of the wind profile is being computed.
!> All subroutines required can be found within this file.
!>
!> @todo harmonize code with new surface_layer_fluxes module
!> @bug 1D model crashes when using small grid spacings in the order of 1 m
!------------------------------------------------------------------------------!
 SUBROUTINE init_1d_model
 

    USE arrays_3d,                                                             &
        ONLY:  l_grid, ug, u_init, vg, v_init, zu
    
    USE indices,                                                               &
        ONLY:  nzb, nzt
    
    USE kinds
    
    USE model_1d,                                                              &
        ONLY:  e1d, e1d_p, kh1d, km1d, l1d, l_black, qs1d, rif1d,              &
               simulated_time_1d, te_e, te_em, te_u, te_um, te_v, te_vm, ts1d, &
               u1d, u1d_p, us1d, usws1d, v1d, v1d_p, vsws1d, z01d, z0h1d
    
    USE control_parameters,                                                    &
        ONLY:  constant_diffusion, constant_flux_layer, f, humidity, kappa,    &
               km_constant, mixing_length_1d, prandtl_number,                  &
               roughness_length, simulated_time_chr, z0h_factor

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string  !<
    
    INTEGER(iwp) ::  k  !< 
    
    REAL(wp) ::  lambda !<

!
!-- Allocate required 1D-arrays
    ALLOCATE( e1d(nzb:nzt+1),    e1d_p(nzb:nzt+1),                             &
              kh1d(nzb:nzt+1),   km1d(nzb:nzt+1),                              &
              l_black(nzb:nzt+1), l1d(nzb:nzt+1),                              &
              rif1d(nzb:nzt+1),   te_e(nzb:nzt+1),                             &
              te_em(nzb:nzt+1),  te_u(nzb:nzt+1),    te_um(nzb:nzt+1),         &
              te_v(nzb:nzt+1),   te_vm(nzb:nzt+1),    u1d(nzb:nzt+1),          &
              u1d_p(nzb:nzt+1),  v1d(nzb:nzt+1),                               &
              v1d_p(nzb:nzt+1) )

!
!-- Initialize arrays
    IF ( constant_diffusion )  THEN
       km1d = km_constant
       kh1d = km_constant / prandtl_number
    ELSE
       e1d = 0.0_wp; e1d_p = 0.0_wp
       kh1d = 0.0_wp; km1d = 0.0_wp
       rif1d = 0.0_wp
!
!--    Compute the mixing length
       l_black(nzb) = 0.0_wp

       IF ( TRIM( mixing_length_1d ) == 'blackadar' )  THEN
!
!--       Blackadar mixing length
          IF ( f /= 0.0_wp )  THEN
             lambda = 2.7E-4_wp * SQRT( ug(nzt+1)**2 + vg(nzt+1)**2 ) /        &
                               ABS( f ) + 1E-10_wp
          ELSE
             lambda = 30.0_wp
          ENDIF

          DO  k = nzb+1, nzt+1
             l_black(k) = kappa * zu(k) / ( 1.0_wp + kappa * zu(k) / lambda )
          ENDDO

       ELSEIF ( TRIM( mixing_length_1d ) == 'as_in_3d_model' )  THEN
!
!--       Use the same mixing length as in 3D model
          l_black(1:nzt) = l_grid
          l_black(nzt+1) = l_black(nzt)

       ENDIF
    ENDIF
    l1d   = l_black
    u1d   = u_init
    u1d_p = u_init
    v1d   = v_init
    v1d_p = v_init

!
!-- Set initial horizontal velocities at the lowest grid levels to a very small
!-- value in order to avoid too small time steps caused by the diffusion limit
!-- in the initial phase of a run (at k=1, dz/2 occurs in the limiting formula!)
    u1d(0:1)   = 0.1_wp
    u1d_p(0:1) = 0.1_wp
    v1d(0:1)   = 0.1_wp
    v1d_p(0:1) = 0.1_wp

!
!-- For u*, theta* and the momentum fluxes plausible values are set
    IF ( constant_flux_layer )  THEN
       us1d = 0.1_wp   ! without initial friction the flow would not change
    ELSE
       e1d(nzb+1)  = 1.0_wp
       km1d(nzb+1) = 1.0_wp
       us1d = 0.0_wp
    ENDIF
    ts1d = 0.0_wp
    usws1d = 0.0_wp
    vsws1d = 0.0_wp
    z01d  = roughness_length
    z0h1d = z0h_factor * z01d 
    IF ( humidity )  qs1d = 0.0_wp

!
!-- Tendencies must be preset in order to avoid runtime errors within the
!-- first Runge-Kutta step
    te_em = 0.0_wp
    te_um = 0.0_wp
    te_vm = 0.0_wp

!
!-- Set start time in hh:mm:ss - format
    simulated_time_chr = time_to_string( simulated_time_1d )

!
!-- Integrate the 1D-model equations using the leap-frog scheme
    CALL time_integration_1d


 END SUBROUTINE init_1d_model



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Leap-frog time differencing scheme for the 1D-model.
!------------------------------------------------------------------------------!
 
 SUBROUTINE time_integration_1d


    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzw, l_grid, pt_init, q_init, ug, vg, zu
        
    USE control_parameters,                                                    &
        ONLY:  constant_diffusion, constant_flux_layer, dissipation_1d,        &
               humidity, intermediate_timestep_count,                          &
               intermediate_timestep_count_max, f, g, ibc_e_b, kappa,          &  
               mixing_length_1d,                                               &
               simulated_time_chr, timestep_scheme, tsc, zeta_max, zeta_min
                
    USE indices,                                                               &
        ONLY:  nzb, nzb_diff, nzt
        
    USE kinds
    
    USE model_1d,                                                              &
        ONLY:  current_timestep_number_1d, damp_level_ind_1d, dt_1d,           &
               dt_pr_1d, dt_run_control_1d, e1d, e1d_p, end_time_1d,           &
               kh1d, km1d, l1d, l_black, qs1d, rif1d, simulated_time_1d,       &
               stop_dt_1d, te_e, te_em, te_u, te_um, te_v, te_vm, time_pr_1d,  &
               ts1d, time_run_control_1d, u1d, u1d_p, us1d, usws1d, v1d,       &
               v1d_p, vsws1d, z01d, z0h1d
        
    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string  !<
    
    INTEGER(iwp) ::  k  !<
    
    REAL(wp) ::  a            !<
    REAL(wp) ::  b            !<
    REAL(wp) ::  dissipation  !<
    REAL(wp) ::  dpt_dz       !<
    REAL(wp) ::  flux         !<
    REAL(wp) ::  kmzm         !<
    REAL(wp) ::  kmzp         !<
    REAL(wp) ::  l_stable     !<
    REAL(wp) ::  pt_0         !<
    REAL(wp) ::  uv_total     !<

!
!-- Determine the time step at the start of a 1D-simulation and
!-- determine and printout quantities used for run control
    dt_1d = 10.0_wp
    CALL run_control_1d

!
!-- Start of time loop
    DO  WHILE ( simulated_time_1d < end_time_1d  .AND.  .NOT. stop_dt_1d )

!
!--    Depending on the timestep scheme, carry out one or more intermediate
!--    timesteps

       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

          CALL timestep_scheme_steering

!
!--       Compute all tendency terms. If a Prandtl-layer is simulated, k starts 
!--       at nzb+2.
          DO  k = nzb_diff, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
!
!--          u-component
             te_u(k) =  f * ( v1d(k) - vg(k) ) + ( &
                              kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1) &
                            - kmzm * ( u1d(k) - u1d(k-1) ) * ddzu(k)   &
                                                 ) * ddzw(k)
!
!--          v-component
             te_v(k) = -f * ( u1d(k) - ug(k) ) + (                     &
                              kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1) &
                            - kmzm * ( v1d(k) - v1d(k-1) ) * ddzu(k)   &
                                                 ) * ddzw(k)
          ENDDO
          IF ( .NOT. constant_diffusion )  THEN
             DO  k = nzb_diff, nzt
!
!--             TKE
                kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
                kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                  &
                            0.61_wp * pt_init(k) *                             &
                            ( q_init(k+1) - q_init(k-1) ) ) * dd2zu(k)
                ENDIF

                IF ( dissipation_1d == 'detering' )  THEN
!
!--                According to Detering, c_e=0.064
                   dissipation = 0.064_wp * e1d(k) * SQRT( e1d(k) ) / l1d(k)
                ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                   dissipation = ( 0.19_wp + 0.74_wp * l1d(k) / l_grid(k) )    &
                                 * e1d(k) * SQRT( e1d(k) ) / l1d(k)
                ENDIF

                te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2&
                                    + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2&
                                    )                                          &
                                    - g / pt_0 * kh1d(k) * flux                &
                                    +            (                             &
                                     kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)  &
                                   - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)    &
                                                 ) * ddzw(k)                   &
                                   - dissipation
             ENDDO
          ENDIF

!
!--       Tendency terms at the top of the Prandtl-layer.
!--       Finite differences of the momentum fluxes are computed using half the 
!--       normal grid length (2.0*ddzw(k)) for the sake of enhanced accuracy 
          IF ( constant_flux_layer )  THEN

             k = nzb+1
             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
             IF ( .NOT. humidity )  THEN
                pt_0 = pt_init(k)
                flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
             ELSE
                pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                     &
                         0.61_wp * pt_init(k) * ( q_init(k+1) - q_init(k-1) )  &
                       ) * dd2zu(k)
             ENDIF

             IF ( dissipation_1d == 'detering' )  THEN
!
!--             According to Detering, c_e=0.064
                dissipation = 0.064_wp * e1d(k) * SQRT( e1d(k) ) / l1d(k)
             ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                dissipation = ( 0.19_wp + 0.74_wp * l1d(k) / l_grid(k) )       &
                              * e1d(k) * SQRT( e1d(k) ) / l1d(k)
             ENDIF

!
!--          u-component
             te_u(k) = f * ( v1d(k) - vg(k) ) + (                              &
                       kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1) + usws1d       &
                                                ) * 2.0_wp * ddzw(k)
!
!--          v-component
             te_v(k) = -f * ( u1d(k) - ug(k) ) + (                             &
                       kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1) + vsws1d       &
                                                 ) * 2.0_wp * ddzw(k)
!
!--          TKE
             te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2   &
                                 + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2   &
                                 )                                             &
                                 - g / pt_0 * kh1d(k) * flux                   &
                                 +           (                                 &
                                  kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)     &
                                - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)       &
                                              ) * ddzw(k)                      &
                                - dissipation
          ENDIF

!
!--       Prognostic equations for all 1D variables
          DO  k = nzb+1, nzt

             u1d_p(k) = u1d(k) + dt_1d * ( tsc(2) * te_u(k) + &
                                           tsc(3) * te_um(k) )
             v1d_p(k) = v1d(k) + dt_1d * ( tsc(2) * te_v(k) + &
                                           tsc(3) * te_vm(k) )

          ENDDO
          IF ( .NOT. constant_diffusion )  THEN
             DO  k = nzb+1, nzt

                e1d_p(k) = e1d(k) + dt_1d * ( tsc(2) * te_e(k) + &
                                              tsc(3) * te_em(k) )

             ENDDO
!
!--          Eliminate negative TKE values, which can result from the
!--          integration due to numerical inaccuracies. In such cases the TKE
!--          value is reduced to 10 percent of its old value.
             WHERE ( e1d_p < 0.0_wp )  e1d_p = 0.1_wp * e1d
          ENDIF

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' ) THEN
             IF ( intermediate_timestep_count == 1 )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = te_u(k)
                   te_vm(k) = te_v(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO k = nzb+1, nzt
                      te_em(k) = te_e(k)
                   ENDDO
                ENDIF

             ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = -9.5625_wp * te_u(k) + 5.3125_wp * te_um(k)
                   te_vm(k) = -9.5625_wp * te_v(k) + 5.3125_wp * te_vm(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO k = nzb+1, nzt
                      te_em(k) = -9.5625_wp * te_e(k) + 5.3125_wp * te_em(k)
                   ENDDO
                ENDIF

             ENDIF
          ENDIF


!
!--       Boundary conditions for the prognostic variables.
!--       At the top boundary (nzt+1) u,v and e keep their initial values 
!--       (ug(nzt+1), vg(nzt+1), 0), at the bottom boundary the mirror 
!--       boundary condition applies to u and v.
!--       The boundary condition for e is set further below ( (u*/cm)**2 ).
         ! u1d_p(nzb) = -u1d_p(nzb+1)
         ! v1d_p(nzb) = -v1d_p(nzb+1)

          u1d_p(nzb) = 0.0_wp
          v1d_p(nzb) = 0.0_wp

!
!--       Swap the time levels in preparation for the next time step.
          u1d  = u1d_p
          v1d  = v1d_p
          IF ( .NOT. constant_diffusion )  THEN
             e1d  = e1d_p
          ENDIF

!
!--       Compute diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          First compute the vertical fluxes in the Prandtl-layer
             IF ( constant_flux_layer )  THEN
!
!--             Compute theta* using Rif numbers of the previous time step
                IF ( rif1d(1) >= 0.0_wp )  THEN
!
!--                Stable stratification
                   ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /          &
                          ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * rif1d(nzb+1) * &
                                          ( zu(nzb+1) - z0h1d ) / zu(nzb+1)    &
                          )
                ELSE
!
!--                Unstable stratification
                   a = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) )
                   b = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /                 &
                       zu(nzb+1) * z0h1d )
!
!--                In the borderline case the formula for stable stratification 
!--                must be applied, because otherwise a zero division would
!--                occur in the argument of the logarithm.
                   IF ( a == 0.0_wp  .OR.  b == 0.0_wp )  THEN
                      ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /       &
                             ( LOG( zu(nzb+1) / z0h1d ) +                      &
                               5.0_wp * rif1d(nzb+1) *                         &
                               ( zu(nzb+1) - z0h1d ) / zu(nzb+1)               &
                             )
                   ELSE
                      ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /       &
                             LOG( (a-1.0_wp) / (a+1.0_wp) *                    &
                                  (b+1.0_wp) / (b-1.0_wp) )
                   ENDIF
                ENDIF

             ENDIF    ! constant_flux_layer

!
!--          Compute the Richardson-flux numbers,
!--          first at the top of the Prandtl-layer using u* of the previous
!--          time step (+1E-30, if u* = 0), then in the remaining area. There
!--          the rif-numbers of the previous time step are used.

             IF ( constant_flux_layer )  THEN
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(nzb+1)
                   flux = ts1d
                ELSE
                   pt_0 = pt_init(nzb+1) * ( 1.0_wp + 0.61_wp * q_init(nzb+1) )
                   flux = ts1d + 0.61_wp * pt_init(k) * qs1d
                ENDIF
                rif1d(nzb+1) = zu(nzb+1) * kappa * g * flux / &
                               ( pt_0 * ( us1d**2 + 1E-30_wp ) )
             ENDIF

             DO  k = nzb_diff, nzt
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) )                    &
                            + 0.61_wp * pt_init(k)                             &
                            * ( q_init(k+1) - q_init(k-1) )                    &
                          ) * dd2zu(k)
                ENDIF
                IF ( rif1d(k) >= 0.0_wp )  THEN
                   rif1d(k) = g / pt_0 * flux /                                &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2     &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2     &
                               + 1E-30_wp                                      &
                              )
                ELSE
                   rif1d(k) = g / pt_0 * flux /                                &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2     &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2     &
                               + 1E-30_wp                                      &
                              ) * ( 1.0_wp - 16.0_wp * rif1d(k) )**0.25_wp
                ENDIF
             ENDDO
!
!--          Richardson-numbers must remain restricted to a realistic value
!--          range. It is exceeded excessively for very small velocities
!--          (u,v --> 0).
             WHERE ( rif1d < zeta_min )  rif1d = zeta_min
             WHERE ( rif1d > zeta_max )  rif1d = zeta_max

!
!--          Compute u* from the absolute velocity value
             IF ( constant_flux_layer )  THEN
                uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )

                IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
!
!--                Stable stratification
                   us1d = kappa * uv_total / (                                 &
                             LOG( zu(nzb+1) / z01d ) + 5.0_wp * rif1d(nzb+1) * &
                                              ( zu(nzb+1) - z01d ) / zu(nzb+1) &
                                             )
                ELSE
!
!--                Unstable stratification
                   a = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) ) )
                   b = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /  &
                                                     zu(nzb+1) * z01d ) )
!
!--                In the borderline case the formula for stable stratification 
!--                must be applied, because otherwise a zero division would
!--                occur in the argument of the logarithm.
                   IF ( a == 1.0_wp  .OR.  b == 1.0_wp )  THEN
                      us1d = kappa * uv_total / (                              &
                             LOG( zu(nzb+1) / z01d ) +                         &
                             5.0_wp * rif1d(nzb+1) * ( zu(nzb+1) - z01d ) /    &
                                                  zu(nzb+1) )
                   ELSE
                      us1d = kappa * uv_total / (                              &
                                 LOG( (1.0_wp+b) / (1.0_wp-b) * (1.0_wp-a) /   &
                                      (1.0_wp+a) ) +                           &
                                 2.0_wp * ( ATAN( b ) - ATAN( a ) )            &
                                                )
                   ENDIF
                ENDIF

!
!--             Compute the momentum fluxes for the diffusion terms
                usws1d  = - u1d(nzb+1) / uv_total * us1d**2
                vsws1d  = - v1d(nzb+1) / uv_total * us1d**2

!
!--             Boundary condition for the turbulent kinetic energy at the top
!--             of the Prandtl-layer. c_m = 0.4 according to Detering.
!--             Additional Neumann condition de/dz = 0 at nzb is set to ensure
!--             compatibility with the 3D model.
                IF ( ibc_e_b == 2 )  THEN
                   e1d(nzb+1) = ( us1d / 0.1_wp )**2
!                  e1d(nzb+1) = ( us1d / 0.4_wp )**2  !not used so far, see also
                                                      !prandtl_fluxes
                ENDIF
                e1d(nzb) = e1d(nzb+1)

                IF ( humidity ) THEN
!
!--                Compute q* 
                   IF ( rif1d(1) >= 0.0_wp )  THEN
!
!--                   Stable stratification
                      qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /         &
                          ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * rif1d(nzb+1) * &
                                          ( zu(nzb+1) - z0h1d ) / zu(nzb+1)    &
                          )
                   ELSE
!
!--                   Unstable stratification
                      a = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) )
                      b = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /              &
                                         zu(nzb+1) * z0h1d )
!
!--                   In the borderline case the formula for stable stratification 
!--                   must be applied, because otherwise a zero division would
!--                   occur in the argument of the logarithm.
                      IF ( a == 1.0_wp  .OR.  b == 1.0_wp )  THEN
                         qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /      &
                                ( LOG( zu(nzb+1) / z0h1d ) +                   &
                                  5.0_wp * rif1d(nzb+1) *                      &
                                  ( zu(nzb+1) - z0h1d ) / zu(nzb+1)            &
                                )
                      ELSE
                         qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /      &
                                LOG( (a-1.0_wp) / (a+1.0_wp) *                 &
                                     (b+1.0_wp) / (b-1.0_wp) )
                      ENDIF
                   ENDIF                
                ELSE
                   qs1d = 0.0_wp
                ENDIF             

             ENDIF   !  constant_flux_layer

!
!--          Compute the diabatic mixing length
             IF ( mixing_length_1d == 'blackadar' )  THEN
                DO  k = nzb+1, nzt
                   IF ( rif1d(k) >= 0.0_wp )  THEN
                      l1d(k) = l_black(k) / ( 1.0_wp + 5.0_wp * rif1d(k) )
                   ELSE
                      l1d(k) = l_black(k) *                                    &
                               ( 1.0_wp - 16.0_wp * rif1d(k) )**0.25_wp
                   ENDIF
                   l1d(k) = l_black(k)
                ENDDO

             ELSEIF ( mixing_length_1d == 'as_in_3d_model' )  THEN
                DO  k = nzb+1, nzt
                   dpt_dz = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                   IF ( dpt_dz > 0.0_wp )  THEN
                      l_stable = 0.76_wp * SQRT( e1d(k) ) /                    &
                                     SQRT( g / pt_init(k) * dpt_dz ) + 1E-5_wp
                   ELSE
                      l_stable = l_grid(k)
                   ENDIF
                   l1d(k) = MIN( l_grid(k), l_stable )
                ENDDO
             ENDIF

!
!--          Compute the diffusion coefficients for momentum via the
!--          corresponding Prandtl-layer relationship and according to
!--          Prandtl-Kolmogorov, respectively. The unstable stratification is
!--          computed via the adiabatic mixing length, for the unstability has
!--          already been taken account of via the TKE (cf. also Diss.).
             IF ( constant_flux_layer )  THEN
                IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) /                    &
                                 ( 1.0_wp + 5.0_wp * rif1d(nzb+1) )
                ELSE
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) *                    &
                                 ( 1.0_wp - 16.0_wp * rif1d(nzb+1) )**0.25_wp
                ENDIF
             ENDIF
             DO  k = nzb_diff, nzt
!                km1d(k) = 0.4 * SQRT( e1d(k) ) !changed: adjustment to 3D-model
                km1d(k) = 0.1_wp * SQRT( e1d(k) )
                IF ( rif1d(k) >= 0.0_wp )  THEN
                   km1d(k) = km1d(k) * l1d(k)
                ELSE
                   km1d(k) = km1d(k) * l_black(k)
                ENDIF
             ENDDO

!
!--          Add damping layer
             DO  k = damp_level_ind_1d+1, nzt+1
                km1d(k) = 1.1_wp * km1d(k-1)
                km1d(k) = MIN( km1d(k), 10.0_wp )
             ENDDO

!
!--          Compute the diffusion coefficient for heat via the relationship
!--          kh = phim / phih * km
             DO  k = nzb+1, nzt
                IF ( rif1d(k) >= 0.0_wp )  THEN
                   kh1d(k) = km1d(k)
                ELSE
                   kh1d(k) = km1d(k) * ( 1.0_wp - 16.0_wp * rif1d(k) )**0.25_wp
                ENDIF
             ENDDO

          ENDIF   ! .NOT. constant_diffusion

       ENDDO   ! intermediate step loop

!
!--    Increment simulated time and output times
       current_timestep_number_1d = current_timestep_number_1d + 1
       simulated_time_1d          = simulated_time_1d + dt_1d
       simulated_time_chr         = time_to_string( simulated_time_1d )
       time_pr_1d                 = time_pr_1d          + dt_1d
       time_run_control_1d        = time_run_control_1d + dt_1d

!
!--    Determine and print out quantities for run control
       IF ( time_run_control_1d >= dt_run_control_1d )  THEN
          CALL run_control_1d
          time_run_control_1d = time_run_control_1d - dt_run_control_1d
       ENDIF

!
!--    Profile output on file
       IF ( time_pr_1d >= dt_pr_1d )  THEN
          CALL print_1d_model
          time_pr_1d = time_pr_1d - dt_pr_1d
       ENDIF

!
!--    Determine size of next time step
       CALL timestep_1d

    ENDDO   ! time loop


 END SUBROUTINE time_integration_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute and print out quantities for run control of the 1D model.
!------------------------------------------------------------------------------!
 
 SUBROUTINE run_control_1d


    USE constants,                                                             &
        ONLY:  pi
        
    USE indices,                                                               &
        ONLY:  nzb, nzt
        
    USE kinds
    
    USE model_1d,                                                              &
        ONLY:  current_timestep_number_1d, dt_1d, run_control_header_1d, u1d,  &
               us1d, v1d
    
    USE pegrid
    
    USE control_parameters,                                                    &
        ONLY:  simulated_time_chr

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !<
    
    REAL(wp) ::  alpha 
    REAL(wp) ::  energy 
    REAL(wp) ::  umax
    REAL(wp) ::  uv_total 
    REAL(wp) ::  vmax

!
!-- Output
    IF ( myid == 0 )  THEN
!
!--    If necessary, write header
       IF ( .NOT. run_control_header_1d )  THEN
          CALL check_open( 15 )
          WRITE ( 15, 100 )
          run_control_header_1d = .TRUE.
       ENDIF

!
!--    Compute control quantities
!--    grid level nzp is excluded due to mirror boundary condition
       umax = 0.0_wp; vmax = 0.0_wp; energy = 0.0_wp
       DO  k = nzb+1, nzt+1
          umax = MAX( ABS( umax ), ABS( u1d(k) ) )
          vmax = MAX( ABS( vmax ), ABS( v1d(k) ) )
          energy = energy + 0.5_wp * ( u1d(k)**2 + v1d(k)**2 )
       ENDDO
       energy = energy / REAL( nzt - nzb + 1, KIND=wp )

       uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )
       IF ( ABS( v1d(nzb+1) ) < 1.0E-5_wp )  THEN
          alpha = ACOS( SIGN( 1.0_wp , u1d(nzb+1) ) )
       ELSE
          alpha = ACOS( u1d(nzb+1) / uv_total )
          IF ( v1d(nzb+1) <= 0.0_wp )  alpha = 2.0_wp * pi - alpha
       ENDIF
       alpha = alpha / ( 2.0_wp * pi ) * 360.0_wp

       WRITE ( 15, 101 )  current_timestep_number_1d, simulated_time_chr, &
                          dt_1d, umax, vmax, us1d, alpha, energy
!
!--    Write buffer contents to disc immediately
       FLUSH( 15 )

    ENDIF

!
!-- formats
100 FORMAT (///'1D-Zeitschrittkontrollausgaben:'/ &
              &'------------------------------'// &
           &'ITER.  HH:MM:SS    DT      UMAX   VMAX    U*   ALPHA   ENERG.'/ &
           &'-------------------------------------------------------------')
101 FORMAT (I5,2X,A9,1X,F6.2,2X,F6.2,1X,F6.2,1X,F6.3,2X,F5.1,2X,F7.2)


 END SUBROUTINE run_control_1d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the time step w.r.t. the diffusion criterion
!------------------------------------------------------------------------------!
 
 SUBROUTINE timestep_1d


    USE arrays_3d,                                                             &
        ONLY:  dzu, zu
        
    USE indices,                                                               &
        ONLY:  nzb, nzt
    
    USE kinds
    
    USE model_1d,                                                              &
        ONLY:  dt_1d, dt_max_1d, km1d, old_dt_1d, stop_dt_1d
    
    USE pegrid
    
    USE control_parameters,                                                    &
        ONLY:  message_string

    IMPLICIT NONE

    INTEGER(iwp) ::  k !<
    
    REAL(wp) ::  div      !<
    REAL(wp) ::  dt_diff  !<
    REAL(wp) ::  fac      !<
    REAL(wp) ::  value    !<


!
!-- Compute the currently feasible time step according to the diffusion
!-- criterion. At nzb+1 the half grid length is used.
    fac = 0.35_wp
    dt_diff = dt_max_1d
    DO  k = nzb+2, nzt
       value   = fac * dzu(k) * dzu(k) / ( km1d(k) + 1E-20_wp )
       dt_diff = MIN( value, dt_diff )
    ENDDO
    value   = fac * zu(nzb+1) * zu(nzb+1) / ( km1d(nzb+1) + 1E-20_wp )
    dt_1d = MIN( value, dt_diff )

!
!-- Set flag when the time step becomes too small
    IF ( dt_1d < ( 0.00001_wp * dt_max_1d ) )  THEN
       stop_dt_1d = .TRUE.

       WRITE( message_string, * ) 'timestep has exceeded the lower limit &', &
                                  'dt_1d = ',dt_1d,' s   simulation stopped!'
       CALL message( 'timestep_1d', 'PA0192', 1, 2, 0, 6, 0 )
       
    ENDIF

!
!-- A more or less simple new time step value is obtained taking only the
!-- first two significant digits
    div = 1000.0_wp
    DO  WHILE ( dt_1d < div )
       div = div / 10.0_wp
    ENDDO
    dt_1d = NINT( dt_1d * 100.0_wp / div ) * div / 100.0_wp

    old_dt_1d = dt_1d


 END SUBROUTINE timestep_1d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> List output of profiles from the 1D-model
!------------------------------------------------------------------------------!
 
 SUBROUTINE print_1d_model


    USE arrays_3d,                                                             &
        ONLY:  pt_init, zu
        
    USE indices,                                                               &
        ONLY:  nzb, nzt
        
    USE kinds
    
    USE model_1d,                                                              &
        ONLY:  e1d, kh1d, km1d, l1d, rif1d, u1d, v1d
    
    USE pegrid
    
    USE control_parameters,                                                    &
        ONLY:  run_description_header, simulated_time_chr

    IMPLICIT NONE


    INTEGER(iwp) ::  k  !<


    IF ( myid == 0 )  THEN
!
!--    Open list output file for profiles from the 1D-model
       CALL check_open( 17 )

!
!--    Write Header
       WRITE ( 17, 100 )  TRIM( run_description_header ), &
                          TRIM( simulated_time_chr )
       WRITE ( 17, 101 )

!
!--    Write the values
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )
       DO  k = nzt+1, nzb, -1
          WRITE ( 17, 103)  k, zu(k), u1d(k), v1d(k), pt_init(k), e1d(k), &
                            rif1d(k), km1d(k), kh1d(k), l1d(k), zu(k), k
       ENDDO
       WRITE ( 17, 101 )
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )

!
!--    Write buffer contents to disc immediately
       FLUSH( 17 )

    ENDIF

!
!-- Formats
100 FORMAT (//1X,A/1X,10('-')/' 1d-model profiles'/ &
            ' Time: ',A)
101 FORMAT (1X,79('-'))
102 FORMAT ('   k     zu      u      v     pt      e    rif    Km    Kh     ', &
            'l      zu      k')
103 FORMAT (1X,I4,1X,F7.1,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F5.2,1X,F5.2, &
            1X,F5.2,1X,F6.2,1X,F7.1,2X,I4)


 END SUBROUTINE print_1d_model
