!> @file production_e.f90
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
! $Id: production_e.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1374 2014-04-25 12:55:07Z raasch
! nzb_s_outer removed from acc-present-list
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1342 2014-03-26 17:04:47Z kanani
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
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1179 2013-06-14 05:57:58Z raasch
! use_reference renamed use_single_reference_value
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1007 2012-09-19 14:30:36Z franke
! Bugfix: calculation of buoyancy production has to consider the liquid water
! mixing ratio in case of cloud droplets
!
! 940 2012-07-09 14:31:00Z raasch
! TKE production by buoyancy can be switched off in case of runs with pure
! neutral stratification
!
! Revision 1.1  1997/09/19 07:45:35  raasch
! Initial revision
!
!
! Description:
! ------------
!> Production terms (shear + buoyancy) of the TKE.
!> @warning The case with constant_flux_layer = F and use_surface_fluxes = T is
!>          not considered well!
!------------------------------------------------------------------------------!
 MODULE production_e_mod
 

    USE wall_fluxes_mod,                                                       &
        ONLY:  wall_fluxes_e, wall_fluxes_e_acc

    USE kinds

    PRIVATE
    PUBLIC production_e, production_e_acc, production_e_init

    LOGICAL, SAVE ::  first_call = .TRUE.  !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  u_0  !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  v_0  !<

    INTERFACE production_e
       MODULE PROCEDURE production_e
       MODULE PROCEDURE production_e_ij
    END INTERFACE production_e
    
    INTERFACE production_e_acc
       MODULE PROCEDURE production_e_acc
    END INTERFACE production_e_acc

    INTERFACE production_e_init
       MODULE PROCEDURE production_e_init
    END INTERFACE production_e_init
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE production_e

       USE arrays_3d,                                                          &
           ONLY:  ddzw, dd2zu, kh, km, pt, q, ql, qsws, qswst, rho_ocean, shf,       &
                  tend, tswst, u, v, vpt, w

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, l_d_r, pt_d_t, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  cloud_droplets, cloud_physics, constant_flux_layer, g,       &
                  humidity, kappa, neutral, ocean, pt_reference,               &
                  rho_reference, use_single_reference_value,                   &
                  use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                     &
           ONLY:  ddx, dx, ddy, dy, wall_e_x, wall_e_y

       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzb_diff_s_inner,                   &
                   nzb_diff_s_outer, nzb_s_inner, nzt, nzt_diff

       IMPLICIT NONE

       INTEGER(iwp) ::  i           !<
       INTEGER(iwp) ::  j           !<
       INTEGER(iwp) ::  k           !<

       REAL(wp)     ::  def         !<
       REAL(wp)     ::  dudx        !<
       REAL(wp)     ::  dudy        !<
       REAL(wp)     ::  dudz        !<
       REAL(wp)     ::  dvdx        !<
       REAL(wp)     ::  dvdy        !<
       REAL(wp)     ::  dvdz        !<
       REAL(wp)     ::  dwdx        !<
       REAL(wp)     ::  dwdy        !<
       REAL(wp)     ::  dwdz        !<
       REAL(wp)     ::  k1          !<
       REAL(wp)     ::  k2          !<
       REAL(wp)     ::  km_neutral  !<
       REAL(wp)     ::  theta       !<
       REAL(wp)     ::  temp        !<

!       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  usvs, vsus, wsus, wsvs
       REAL(wp), DIMENSION(nzb:nzt+1) ::  usvs  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  vsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  wsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  wsvs  !<

!
!--    First calculate horizontal momentum flux u'v', w'v', v'u', w'u' at
!--    vertical walls, if neccessary
!--    So far, results are slightly different from the ij-Version.
!--    Therefore, ij-Version is called further below within the ij-loops.
!       IF ( topography /= 'flat' )  THEN
!          CALL wall_fluxes_e( usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, wall_e_y )
!          CALL wall_fluxes_e( wsvs, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, wall_e_y )
!          CALL wall_fluxes_e( vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, wall_e_x )
!          CALL wall_fluxes_e( wsus, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, wall_e_x )
!       ENDIF


       DO  i = nxl, nxr

!
!--       Calculate TKE production by shear
          DO  j = nys, nyn
             DO  k = nzb_diff_s_outer(j,i), nzt

                dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                    u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                    u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                    v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                    v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                    w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                    w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                      dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                      dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                IF ( def < 0.0_wp )  def = 0.0_wp

                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

             ENDDO
          ENDDO

          IF ( constant_flux_layer )  THEN

!
!--          Position beneath wall
!--          (2) - Will allways be executed.
!--          'bottom and wall: use u_0,v_0 and wall functions'
             DO  j = nys, nyn

                IF ( ( wall_e_x(j,i) /= 0.0_wp ) .OR. ( wall_e_y(j,i) /= 0.0_wp ) ) &
                THEN

                   k = nzb_diff_s_inner(j,i) - 1
                   dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
                   dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                     u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)
                   dvdy = ( v(k,j+1,i) - v(k,j,i) ) * ddy
                   dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                     v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)
                   dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

                   IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--                   Inconsistency removed: as the thermal stratification is 
!--                   not taken into account for the evaluation of the wall 
!--                   fluxes at vertical walls, the eddy viscosity km must not 
!--                   be used for the evaluation of the velocity gradients dudy 
!--                   and dwdy
!--                   Note: The validity of the new method has not yet been 
!--                         shown, as so far no suitable data for a validation 
!--                         has been available
                      CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                          usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
                      CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                          wsvs, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp )
                      km_neutral = kappa * ( usvs(k)**2 + wsvs(k)**2 )**0.25_wp * &
                                   0.5_wp * dy 
                      IF ( km_neutral > 0.0_wp )  THEN
                         dudy = - wall_e_y(j,i) * usvs(k) / km_neutral
                         dwdy = - wall_e_y(j,i) * wsvs(k) / km_neutral
                      ELSE
                         dudy = 0.0_wp
                         dwdy = 0.0_wp
                      ENDIF
                   ELSE
                      dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                         u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                      dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                         w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                   ENDIF

                   IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--                   Inconsistency removed: as the thermal stratification is 
!--                   not taken into account for the evaluation of the wall 
!--                   fluxes at vertical walls, the eddy viscosity km must not 
!--                   be used for the evaluation of the velocity gradients dvdx 
!--                   and dwdx
!--                   Note: The validity of the new method has not yet been 
!--                         shown, as so far no suitable data for a validation 
!--                         has been available
                      CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                          vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp )
                      CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                          wsus, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp )
                      km_neutral = kappa * ( vsus(k)**2 + wsus(k)**2 )**0.25_wp * &
                                   0.5_wp * dx
                      IF ( km_neutral > 0.0_wp )  THEN
                         dvdx = - wall_e_x(j,i) * vsus(k) / km_neutral
                         dwdx = - wall_e_x(j,i) * wsus(k) / km_neutral
                      ELSE
                         dvdx = 0.0_wp
                         dwdx = 0.0_wp
                      ENDIF
                   ELSE
                      dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                         v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                      dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                         w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                   ENDIF

                   def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                         dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                         dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                   IF ( def < 0.0_wp )  def = 0.0_wp

                   tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def


!
!--                (3) - will be executed only, if there is at least one level
!--                between (2) and (4), i.e. the topography must have a
!--                minimum height of 2 dz. Wall fluxes for this case have
!--                already been calculated for (2).
!--                'wall only: use wall functions'

                   DO  k = nzb_diff_s_inner(j,i), nzb_diff_s_outer(j,i)-2

                      dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
                      dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                        u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
                      dvdy =          ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                      dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                        v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
                      dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

                      IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--                      Inconsistency removed: as the thermal stratification 
!--                      is not taken into account for the evaluation of the 
!--                      wall fluxes at vertical walls, the eddy viscosity km 
!--                      must not be used for the evaluation of the velocity 
!--                      gradients dudy and dwdy
!--                      Note: The validity of the new method has not yet 
!--                            been shown, as so far no suitable data for a 
!--                            validation has been available
                         km_neutral = kappa * ( usvs(k)**2 + & 
                                                wsvs(k)**2 )**0.25_wp * 0.5_wp * dy
                         IF ( km_neutral > 0.0_wp )  THEN
                            dudy = - wall_e_y(j,i) * usvs(k) / km_neutral
                            dwdy = - wall_e_y(j,i) * wsvs(k) / km_neutral
                         ELSE
                            dudy = 0.0_wp
                            dwdy = 0.0_wp
                         ENDIF
                      ELSE
                         dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                            u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                         dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                            w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                      ENDIF

                      IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--                      Inconsistency removed: as the thermal stratification 
!--                      is not taken into account for the evaluation of the 
!--                      wall fluxes at vertical walls, the eddy viscosity km 
!--                      must not be used for the evaluation of the velocity 
!--                      gradients dvdx and dwdx
!--                      Note: The validity of the new method has not yet 
!--                            been shown, as so far no suitable data for a 
!--                            validation has been available
                         km_neutral = kappa * ( vsus(k)**2 + & 
                                                wsus(k)**2 )**0.25_wp * 0.5_wp * dx
                         IF ( km_neutral > 0.0_wp )  THEN
                            dvdx = - wall_e_x(j,i) * vsus(k) / km_neutral
                            dwdx = - wall_e_x(j,i) * wsus(k) / km_neutral
                         ELSE
                            dvdx = 0.0_wp
                            dwdx = 0.0_wp
                         ENDIF
                      ELSE
                         dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                            v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                         dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                            w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                      ENDIF

                      def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                           dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 +  &
                           dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                      IF ( def < 0.0_wp )  def = 0.0_wp

                      tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                   ENDDO

                ENDIF

             ENDDO

!
!--          (4) - will allways be executed.
!--          'special case: free atmosphere' (as for case (0))
             DO  j = nys, nyn

                IF ( ( wall_e_x(j,i) /= 0.0_wp ) .OR. ( wall_e_y(j,i) /= 0.0_wp ) ) &
                THEN

                   k = nzb_diff_s_outer(j,i)-1

                   dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                   dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                       u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                   dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                       u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                   dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                       v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                   dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                   dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                       v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                   dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                       w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                   dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                       w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                   dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                   def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                         dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                         dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                   IF ( def < 0.0_wp )  def = 0.0_wp

                   tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                ENDIF

             ENDDO

!
!--          Position without adjacent wall
!--          (1) - will allways be executed. 
!--          'bottom only: use u_0,v_0'
             DO  j = nys, nyn

                IF ( ( wall_e_x(j,i) == 0.0_wp ) .AND. ( wall_e_y(j,i) == 0.0_wp ) ) &
                THEN

                   k = nzb_diff_s_inner(j,i)-1

                   dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                   dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                       u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                   dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                       u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)

                   dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                       v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                   dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                   dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                       v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)

                   dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                       w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                   dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                       w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                   dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                   def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                         dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                         dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                   IF ( def < 0.0_wp )  def = 0.0_wp

                   tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                ENDIF

             ENDDO

          ELSEIF ( use_surface_fluxes )  THEN

             DO  j = nys, nyn

                k = nzb_diff_s_outer(j,i)-1

                dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                    u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                   u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                    v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                    v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                    w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                    w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                      dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                      dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                IF ( def < 0.0_wp )  def = 0.0_wp

                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

             ENDDO

          ENDIF

!
!--       If required, calculate TKE production by buoyancy
          IF ( .NOT. neutral )  THEN

             IF ( .NOT. humidity )  THEN

                IF ( use_single_reference_value )  THEN

                   IF ( ocean )  THEN
!
!--                   So far in the ocean no special treatment of density flux
!--                   in the bottom and top surface layer
                      DO  j = nys, nyn
                         DO  k = nzb_s_inner(j,i)+1, nzt
                            tend(k,j,i) = tend(k,j,i) +                     &
                                          kh(k,j,i) * g / rho_reference *   &
                                          ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * &
                                          dd2zu(k)
                         ENDDO
                      ENDDO

                   ELSE

                      DO  j = nys, nyn
                         DO  k = nzb_diff_s_inner(j,i), nzt_diff
                            tend(k,j,i) = tend(k,j,i) -                   &
                                          kh(k,j,i) * g / pt_reference *  &
                                          ( pt(k+1,j,i) - pt(k-1,j,i) ) * &
                                          dd2zu(k)
                         ENDDO

                         IF ( use_surface_fluxes )  THEN
                            k = nzb_diff_s_inner(j,i)-1
                            tend(k,j,i) = tend(k,j,i) + g / pt_reference * &
                                                        shf(j,i)
                         ENDIF

                         IF ( use_top_fluxes )  THEN
                            k = nzt
                            tend(k,j,i) = tend(k,j,i) + g / pt_reference * &
                                                        tswst(j,i)
                         ENDIF
                      ENDDO

                   ENDIF

                ELSE

                   IF ( ocean )  THEN
!
!--                   So far in the ocean no special treatment of density flux
!--                   in the bottom and top surface layer
                      DO  j = nys, nyn
                         DO  k = nzb_s_inner(j,i)+1, nzt
                            tend(k,j,i) = tend(k,j,i) +                     &
                                          kh(k,j,i) * g / rho_ocean(k,j,i) *      &
                                          ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * &
                                          dd2zu(k)
                         ENDDO
                      ENDDO

                   ELSE

                      DO  j = nys, nyn
                         DO  k = nzb_diff_s_inner(j,i), nzt_diff
                            tend(k,j,i) = tend(k,j,i) -                   &
                                          kh(k,j,i) * g / pt(k,j,i) *     &
                                          ( pt(k+1,j,i) - pt(k-1,j,i) ) * &
                                          dd2zu(k)
                         ENDDO

                         IF ( use_surface_fluxes )  THEN
                            k = nzb_diff_s_inner(j,i)-1
                            tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * &
                                                        shf(j,i)
                         ENDIF

                         IF ( use_top_fluxes )  THEN
                            k = nzt
                            tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * &
                                                        tswst(j,i)
                         ENDIF
                      ENDDO

                   ENDIF

                ENDIF

             ELSE

                DO  j = nys, nyn

                   DO  k = nzb_diff_s_inner(j,i), nzt_diff

                      IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                         tend(k,j,i) = tend(k,j,i) - kh(k,j,i) *               &
                                         g / vpt(k,j,i) *                      &
                                         ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) +  &
                                           k2 * ( q(k+1,j,i) - q(k-1,j,i) )    &
                                         ) * dd2zu(k)
                      ELSE IF ( cloud_physics )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                            temp  = theta * t_d_pt(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                          ( q(k,j,i) - ql(k,j,i) ) *          &
                                 ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /        &
                                 ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *          &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                         tend(k,j,i) = tend(k,j,i) - kh(k,j,i) *               &
                                         g / vpt(k,j,i) *                      &
                                         ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) +  &
                                           k2 * ( q(k+1,j,i) - q(k-1,j,i) )    &
                                         ) * dd2zu(k)
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                         tend(k,j,i) = tend(k,j,i) -                          & 
                                       kh(k,j,i) * g / vpt(k,j,i) *           &
                                       ( k1 * ( pt(k+1,j,i)- pt(k-1,j,i) ) +  &
                                         k2 * ( q(k+1,j,i) -  q(k-1,j,i) ) -  &
                                         pt(k,j,i) * ( ql(k+1,j,i) -          &
                                         ql(k-1,j,i) ) ) * dd2zu(k)
                      ENDIF

                   ENDDO

                ENDDO

                IF ( use_surface_fluxes )  THEN

                   DO  j = nys, nyn

                      k = nzb_diff_s_inner(j,i)-1

                      IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE IF ( cloud_physics )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                            temp  = theta * t_d_pt(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                          ( q(k,j,i) - ql(k,j,i) ) *           &
                                 ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /      &
                                 ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *        &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ENDIF

                      tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
                                            ( k1* shf(j,i) + k2 * qsws(j,i) )
                   ENDDO

                ENDIF

                IF ( use_top_fluxes )  THEN

                   DO  j = nys, nyn

                      k = nzt

                      IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE IF ( cloud_physics )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                            temp  = theta * t_d_pt(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                          ( q(k,j,i) - ql(k,j,i) ) *           &
                                 ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /      &
                                 ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *        &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ENDIF

                      tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
                                            ( k1* tswst(j,i) + k2 * qswst(j,i) )
                   ENDDO

                ENDIF

             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE production_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE production_e_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzw, dd2zu, kh, km, pt, q, ql, qsws, qswst, rho_ocean, shf,       &
                  tend, tswst, u, v, vpt, w

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, l_d_r, pt_d_t, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  cloud_droplets, cloud_physics, constant_flux_layer, g,       &
                  humidity, kappa, neutral, ocean, pt_reference,               &
                  rho_reference, topography, use_single_reference_value,       &
                  use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                     &
           ONLY:  ddx, dx, ddy, dy, wall_e_x, wall_e_y

       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nys, nyn, nzb,  &
                  nzb_diff_s_inner, nzb_diff_s_outer, nzb_s_inner, nzt,        &
                  nzt_diff

       IMPLICIT NONE

       INTEGER(iwp) ::  i           !<
       INTEGER(iwp) ::  j           !<
       INTEGER(iwp) ::  k           !<

       REAL(wp)     ::  def         !<
       REAL(wp)     ::  dudx        !<
       REAL(wp)     ::  dudy        !<
       REAL(wp)     ::  dudz        !<
       REAL(wp)     ::  dvdx        !<
       REAL(wp)     ::  dvdy        !<
       REAL(wp)     ::  dvdz        !<
       REAL(wp)     ::  dwdx        !<
       REAL(wp)     ::  dwdy        !<
       REAL(wp)     ::  dwdz        !<
       REAL(wp)     ::  k1          !<
       REAL(wp)     ::  k2          !<
       REAL(wp)     ::  km_neutral  !<
       REAL(wp)     ::  theta       !<
       REAL(wp)     ::  temp        !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  usvs  !<
       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  vsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  wsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  wsvs  !<
       !$acc declare create ( usvs, vsus, wsus, wsvs )

!
!--    First calculate horizontal momentum flux u'v', w'v', v'u', w'u' at
!--    vertical walls, if neccessary
!--    CAUTION: results are slightly different from the ij-version!!
!--    ij-version should be called further below within the ij-loops!!
       IF ( topography /= 'flat' )  THEN
          CALL wall_fluxes_e_acc( usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, wall_e_y )
          CALL wall_fluxes_e_acc( wsvs, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, wall_e_y )
          CALL wall_fluxes_e_acc( vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, wall_e_x )
          CALL wall_fluxes_e_acc( wsus, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, wall_e_x )
       ENDIF


!
!--    Calculate TKE production by shear
       !$acc kernels present( ddzw, dd2zu, kh, km, nzb_diff_s_inner, nzb_diff_s_outer ) &
       !$acc         present( nzb_s_inner, pt, q, ql, qsws, qswst, rho_ocean )                &
       !$acc         present( shf, tend, tswst, u, v, vpt, w, wall_e_x, wall_e_y )      &
       !$acc         copyin( u_0, v_0 )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = 1, nzt

                IF ( k >= nzb_diff_s_outer(j,i) )  THEN

                   dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                   dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                       u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                   dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                       u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                   dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                       v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                   dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                   dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                       v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                   dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                       w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                   dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                       w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                   dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                   def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                         dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                         dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                   IF ( def < 0.0_wp )  def = 0.0_wp

                   tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       IF ( constant_flux_layer )  THEN

!
!--       Position beneath wall
!--       (2) - Will allways be executed.
!--       'bottom and wall: use u_0,v_0 and wall functions'
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = 1, nzt

                   IF ( ( wall_e_x(j,i) /= 0.0_wp ).OR.( wall_e_y(j,i) /= 0.0_wp ) ) &
                   THEN

                      IF ( k == nzb_diff_s_inner(j,i) - 1 )  THEN
                         dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
                         dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                           u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)
                         dvdy = ( v(k,j+1,i) - v(k,j,i) ) * ddy
                         dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                           v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)
                         dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

                         IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--                         Inconsistency removed: as the thermal stratification is
!--                         not taken into account for the evaluation of the wall
!--                         fluxes at vertical walls, the eddy viscosity km must not
!--                         be used for the evaluation of the velocity gradients dudy
!--                         and dwdy
!--                         Note: The validity of the new method has not yet been
!--                               shown, as so far no suitable data for a validation
!--                               has been available
!                            CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
!                                                usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
!                            CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
!                                                wsvs, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp )
                            km_neutral = kappa *                                    &
                                        ( usvs(k,j,i)**2 + wsvs(k,j,i)**2 )**0.25_wp * &
                                         0.5_wp * dy
                            IF ( km_neutral > 0.0_wp )  THEN
                               dudy = - wall_e_y(j,i) * usvs(k,j,i) / km_neutral
                               dwdy = - wall_e_y(j,i) * wsvs(k,j,i) / km_neutral
                            ELSE
                               dudy = 0.0_wp
                               dwdy = 0.0_wp
                            ENDIF
                         ELSE
                            dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                               u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                            dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                               w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                         ENDIF

                         IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--                         Inconsistency removed: as the thermal stratification is
!--                         not taken into account for the evaluation of the wall
!--                         fluxes at vertical walls, the eddy viscosity km must not
!--                         be used for the evaluation of the velocity gradients dvdx
!--                         and dwdx
!--                         Note: The validity of the new method has not yet been
!--                               shown, as so far no suitable data for a validation
!--                               has been available
!                            CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
!                                                vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp )
!                            CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
!                                                wsus, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp )
                            km_neutral = kappa *                                     &
                                         ( vsus(k,j,i)**2 + wsus(k,j,i)**2 )**0.25_wp * &
                                         0.5_wp * dx
                            IF ( km_neutral > 0.0_wp )  THEN
                               dvdx = - wall_e_x(j,i) * vsus(k,j,i) / km_neutral
                               dwdx = - wall_e_x(j,i) * wsus(k,j,i) / km_neutral
                            ELSE
                               dvdx = 0.0_wp
                               dwdx = 0.0_wp
                            ENDIF
                         ELSE
                            dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                               v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                            dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                               w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                         ENDIF

                         def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                               dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                               dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                         IF ( def < 0.0_wp )  def = 0.0_wp

                         tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                      ENDIF
!
!--                   (3) - will be executed only, if there is at least one level
!--                   between (2) and (4), i.e. the topography must have a
!--                   minimum height of 2 dz. Wall fluxes for this case have
!--                   already been calculated for (2).
!--                   'wall only: use wall functions'

                      IF ( k >= nzb_diff_s_inner(j,i)  .AND.  &
                           k <= nzb_diff_s_outer(j,i)-2 )  THEN

                         dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
                         dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                           u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
                         dvdy =          ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                         dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                           v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
                         dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

                         IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--                         Inconsistency removed: as the thermal stratification
!--                         is not taken into account for the evaluation of the
!--                         wall fluxes at vertical walls, the eddy viscosity km
!--                         must not be used for the evaluation of the velocity
!--                         gradients dudy and dwdy
!--                         Note: The validity of the new method has not yet
!--                               been shown, as so far no suitable data for a
!--                               validation has been available
                            km_neutral = kappa * ( usvs(k,j,i)**2 + &
                                                   wsvs(k,j,i)**2 )**0.25_wp * 0.5_wp * dy
                            IF ( km_neutral > 0.0_wp )  THEN
                               dudy = - wall_e_y(j,i) * usvs(k,j,i) / km_neutral
                               dwdy = - wall_e_y(j,i) * wsvs(k,j,i) / km_neutral
                            ELSE
                               dudy = 0.0_wp
                               dwdy = 0.0_wp
                            ENDIF
                         ELSE
                            dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                               u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                            dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                               w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                         ENDIF

                         IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--                         Inconsistency removed: as the thermal stratification
!--                         is not taken into account for the evaluation of the
!--                         wall fluxes at vertical walls, the eddy viscosity km
!--                         must not be used for the evaluation of the velocity
!--                         gradients dvdx and dwdx
!--                         Note: The validity of the new method has not yet
!--                               been shown, as so far no suitable data for a
!--                               validation has been available
                            km_neutral = kappa * ( vsus(k,j,i)**2 + &
                                                   wsus(k,j,i)**2 )**0.25_wp * 0.5_wp * dx
                            IF ( km_neutral > 0.0_wp )  THEN
                               dvdx = - wall_e_x(j,i) * vsus(k,j,i) / km_neutral
                               dwdx = - wall_e_x(j,i) * wsus(k,j,i) / km_neutral
                            ELSE
                               dvdx = 0.0_wp
                               dwdx = 0.0_wp
                            ENDIF
                         ELSE
                            dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                               v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                            dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                               w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                         ENDIF

                         def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                              dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 +  &
                              dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                         IF ( def < 0.0_wp )  def = 0.0_wp

                         tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                      ENDIF

!
!--                   (4) - will allways be executed.
!--                   'special case: free atmosphere' (as for case (0))
                      IF ( k == nzb_diff_s_outer(j,i)-1 )  THEN

                         dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                         dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                             u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                         dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                             u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                         dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                             v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                         dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                         dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                             v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                         dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                             w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                         dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                             w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                         dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                         def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                               dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                               dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                         IF ( def < 0.0_wp )  def = 0.0_wp

                         tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                      ENDIF

                   ENDIF

                ENDDO
             ENDDO
          ENDDO

!
!--       Position without adjacent wall
!--       (1) - will allways be executed.
!--       'bottom only: use u_0,v_0'
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = 1, nzt

                   IF ( ( wall_e_x(j,i) == 0.0_wp ) .AND. ( wall_e_y(j,i) == 0.0_wp ) ) &
                   THEN

                      IF ( k == nzb_diff_s_inner(j,i)-1 )  THEN

                         dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                         dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                             u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                         dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                             u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)

                         dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                             v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                         dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                         dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                             v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)

                         dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                             w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                         dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                             w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                         dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                         def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                               dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                               dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                         IF ( def < 0.0_wp )  def = 0.0_wp

                         tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                      ENDIF

                   ENDIF

                ENDDO
             ENDDO
          ENDDO

       ELSEIF ( use_surface_fluxes )  THEN

          DO  i = i_left, i_right
             DO  j = j_south, j_north
                 DO  k = 1, nzt

                   IF ( k == nzb_diff_s_outer(j,i)-1 )  THEN

                      dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                      dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                          u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                      dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                          u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

                      dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                          v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                      dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                      dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                          v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

                      dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                          w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                      dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                          w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                      dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

                      def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                            dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                            dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                      IF ( def < 0.0_wp )  def = 0.0_wp

                      tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

                   ENDIF

                ENDDO
             ENDDO
          ENDDO

       ENDIF

!
!--    If required, calculate TKE production by buoyancy
       IF ( .NOT. neutral )  THEN

          IF ( .NOT. humidity )  THEN

             IF ( use_single_reference_value )  THEN

                IF ( ocean )  THEN
!
!--                So far in the ocean no special treatment of density flux
!--                in the bottom and top surface layer
                   DO  i = i_left, i_right
                      DO  j = j_south, j_north
                         DO  k = 1, nzt
                            IF ( k > nzb_s_inner(j,i) )  THEN
                               tend(k,j,i) = tend(k,j,i) +                     &
                                             kh(k,j,i) * g / rho_reference *   &
                                             ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * &
                                             dd2zu(k)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  i = i_left, i_right
                      DO  j = j_south, j_north
                         DO  k = 1, nzt_diff
                            IF ( k >= nzb_diff_s_inner(j,i) )  THEN
                               tend(k,j,i) = tend(k,j,i) -                   &
                                             kh(k,j,i) * g / pt_reference *  &
                                             ( pt(k+1,j,i) - pt(k-1,j,i) ) * &
                                             dd2zu(k)
                            ENDIF

                            IF ( k == nzb_diff_s_inner(j,i)-1  .AND.  &
                                 use_surface_fluxes )  THEN
                               tend(k,j,i) = tend(k,j,i) + g / pt_reference * &
                                                           shf(j,i)
                            ENDIF

                            IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                               tend(k,j,i) = tend(k,j,i) + g / pt_reference * &
                                                           tswst(j,i)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF

             ELSE

                IF ( ocean )  THEN
!
!--                So far in the ocean no special treatment of density flux
!--                in the bottom and top surface layer
                   DO  i = i_left, i_right
                      DO  j = j_south, j_north
                         DO  k = 1, nzt
                            IF ( k > nzb_s_inner(j,i) )  THEN
                               tend(k,j,i) = tend(k,j,i) +                     &
                                             kh(k,j,i) * g / rho_ocean(k,j,i) *      &
                                             ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * &
                                             dd2zu(k)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  i = i_left, i_right
                      DO  j = j_south, j_north
                         DO  k = 1, nzt_diff
                            IF( k >= nzb_diff_s_inner(j,i) )  THEN
                               tend(k,j,i) = tend(k,j,i) -                   &
                                             kh(k,j,i) * g / pt(k,j,i) *     &
                                             ( pt(k+1,j,i) - pt(k-1,j,i) ) * &
                                             dd2zu(k)
                            ENDIF

                            IF (  k == nzb_diff_s_inner(j,i)-1  .AND.  &
                                  use_surface_fluxes )  THEN
                               tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * &
                                                           shf(j,i)
                            ENDIF

                            IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                               tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * &
                                                           tswst(j,i)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF

             ENDIF

          ELSE
!
!++          This part gives the PGI compiler problems in the previous loop
!++          even without any acc statements????
!             STOP '+++ production_e problems with acc-directives'
!             !acc loop
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   !acc loop vector
!                   DO  k = 1, nzt_diff
!
!                      IF ( k >= nzb_diff_s_inner(j,i) )  THEN
!
!                         IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
!                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                            k2 = 0.61_wp * pt(k,j,i)
!                            tend(k,j,i) = tend(k,j,i) - kh(k,j,i) *               &
!                                            g / vpt(k,j,i) *                      &
!                                            ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) +  &
!                                              k2 * ( q(k+1,j,i) - q(k-1,j,i) )    &
!                                            ) * dd2zu(k)
!                         ELSE IF ( cloud_physics )  THEN
!                            IF ( ql(k,j,i) == 0.0_wp )  THEN
!                               k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                               k2 = 0.61_wp * pt(k,j,i)
!                            ELSE
!                               theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
!                               temp  = theta * t_d_pt(k)
!                               k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                 &
!                                             ( q(k,j,i) - ql(k,j,i) ) *          &
!                                    ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /        &
!                                    ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *          &
!                                    ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
!                               k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
!                            ENDIF
!                            tend(k,j,i) = tend(k,j,i) - kh(k,j,i) *               &
!                                            g / vpt(k,j,i) *                      &
!                                            ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) +  &
!                                              k2 * ( q(k+1,j,i) - q(k-1,j,i) )    &
!                                            ) * dd2zu(k)
!                         ELSE IF ( cloud_droplets )  THEN
!                            k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
!                            k2 = 0.61_wp * pt(k,j,i)
!                            tend(k,j,i) = tend(k,j,i) -                          &
!                                          kh(k,j,i) * g / vpt(k,j,i) *           &
!                                          ( k1 * ( pt(k+1,j,i)- pt(k-1,j,i) ) +  &
!                                            k2 * ( q(k+1,j,i) -  q(k-1,j,i) ) -  &
!                                            pt(k,j,i) * ( ql(k+1,j,i) -          &
!                                            ql(k-1,j,i) ) ) * dd2zu(k)
!                         ENDIF
!
!                      ENDIF
!
!                   ENDDO
!                ENDDO
!             ENDDO
!

!!++          Next two loops are probably very inefficiently parallellized
!!++          and will require better optimization
!             IF ( use_surface_fluxes )  THEN
!
!                !acc loop
!                DO  i = nxl, nxr
!                   DO  j = nys, nyn
!                      !acc loop vector
!                      DO  k = 1, nzt_diff
!
!                         IF ( k == nzb_diff_s_inner(j,i)-1 )  THEN
!
!                            IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
!                               k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                               k2 = 0.61_wp * pt(k,j,i)
!                            ELSE IF ( cloud_physics )  THEN
!                               IF ( ql(k,j,i) == 0.0_wp )  THEN
!                                  k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                                  k2 = 0.61_wp * pt(k,j,i)
!                               ELSE
!                                  theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
!                                  temp  = theta * t_d_pt(k)
!                                  k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *        &
!                                                ( q(k,j,i) - ql(k,j,i) ) *    &
!                                       ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /&
!                                       ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp * &
!                                       ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
!                                  k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
!                               ENDIF
!                            ELSE IF ( cloud_droplets )  THEN
!                               k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
!                               k2 = 0.61_wp * pt(k,j,i)
!                            ENDIF
!
!                            tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
!                                                  ( k1* shf(j,i) + k2 * qsws(j,i) )
!                         ENDIF
!
!                      ENDDO
!                   ENDDO
!                ENDDO
!
!             ENDIF
!
!             IF ( use_top_fluxes )  THEN
!
!                !acc loop
!                DO  i = nxl, nxr
!                   DO  j = nys, nyn
!                      !acc loop vector
!                      DO  k = 1, nzt
!                         IF ( k == nzt )  THEN
!
!                            IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets ) THEN
!                               k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                               k2 = 0.61_wp * pt(k,j,i)
!                            ELSE IF ( cloud_physics )  THEN
!                               IF ( ql(k,j,i) == 0.0_wp )  THEN
!                                  k1 = 1.0_wp + 0.61_wp * q(k,j,i)
!                                  k2 = 0.61_wp * pt(k,j,i)
!                               ELSE
!                                  theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
!                                  temp  = theta * t_d_pt(k)
!                                  k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *        &
!                                                ( q(k,j,i) - ql(k,j,i) ) *    &
!                                       ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /&
!                                       ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp * &
!                                       ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
!                                  k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
!                               ENDIF
!                            ELSE IF ( cloud_droplets )  THEN
!                               k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
!                               k2 = 0.61_wp * pt(k,j,i)
!                            ENDIF
!
!                            tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
!                                                  ( k1* tswst(j,i) + k2 * qswst(j,i) )
!
!                         ENDIF
!
!                      ENDDO
!                   ENDDO
!                ENDDO
!
!             ENDIF

          ENDIF

       ENDIF
       !$acc end kernels

    END SUBROUTINE production_e_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE production_e_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzw, dd2zu, kh, km, pt, q, ql, qsws, qswst, rho_ocean, shf,       &
                  tend, tswst, u, v, vpt, w

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, l_d_r, pt_d_t, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  cloud_droplets, cloud_physics, constant_flux_layer, g,       &
                  humidity, kappa, neutral, ocean, pt_reference,               &
                  rho_reference, use_single_reference_value,                   &
                  use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                     &
           ONLY:  ddx, dx, ddy, dy, wall_e_x, wall_e_y

       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzb_diff_s_inner,                   &
                  nzb_diff_s_outer, nzb_s_inner, nzt, nzt_diff

       IMPLICIT NONE

       INTEGER(iwp) ::  i           !<
       INTEGER(iwp) ::  j           !<
       INTEGER(iwp) ::  k           !<

       REAL(wp)     ::  def         !<
       REAL(wp)     ::  dudx        !<
       REAL(wp)     ::  dudy        !<
       REAL(wp)     ::  dudz        !<
       REAL(wp)     ::  dvdx        !<
       REAL(wp)     ::  dvdy        !<
       REAL(wp)     ::  dvdz        !<
       REAL(wp)     ::  dwdx        !<
       REAL(wp)     ::  dwdy        !<
       REAL(wp)     ::  dwdz        !<
       REAL(wp)     ::  k1          !<
       REAL(wp)     ::  k2          !<
       REAL(wp)     ::  km_neutral  !<
       REAL(wp)     ::  theta       !<
       REAL(wp)     ::  temp        !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  usvs  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  vsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  wsus  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  wsvs  !<

!
!--    Calculate TKE production by shear
       DO  k = nzb_diff_s_outer(j,i), nzt

          dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
          dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                              u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
          dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                              u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

          dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                              v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
          dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
          dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                              v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

          dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                              w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
          dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                              w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
          dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

          def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 )                       &
                + dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + dvdz**2    &
                + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

          IF ( def < 0.0_wp )  def = 0.0_wp

          tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

       ENDDO

       IF ( constant_flux_layer )  THEN

          IF ( ( wall_e_x(j,i) /= 0.0_wp ) .OR. ( wall_e_y(j,i) /= 0.0_wp ) )  THEN

!
!--          Position beneath wall
!--          (2) - Will allways be executed.
!--          'bottom and wall: use u_0,v_0 and wall functions'
             k = nzb_diff_s_inner(j,i)-1

             dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
             dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                               u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)
             dvdy =          ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                               v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)
             dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

             IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--             Inconsistency removed: as the thermal stratification 
!--             is not taken into account for the evaluation of the 
!--             wall fluxes at vertical walls, the eddy viscosity km 
!--             must not be used for the evaluation of the velocity 
!--             gradients dudy and dwdy
!--             Note: The validity of the new method has not yet 
!--                   been shown, as so far no suitable data for a 
!--                   validation has been available
                CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                    usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
                CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                    wsvs, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp )
                km_neutral = kappa * ( usvs(k)**2 + wsvs(k)**2 )**0.25_wp * &
                             0.5_wp * dy
                IF ( km_neutral > 0.0_wp )  THEN
                   dudy = - wall_e_y(j,i) * usvs(k) / km_neutral
                   dwdy = - wall_e_y(j,i) * wsvs(k) / km_neutral
                ELSE
                   dudy = 0.0_wp
                   dwdy = 0.0_wp
                ENDIF
             ELSE
                dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                   u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                   w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             ENDIF

             IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--             Inconsistency removed: as the thermal stratification 
!--             is not taken into account for the evaluation of the 
!--             wall fluxes at vertical walls, the eddy viscosity km 
!--             must not be used for the evaluation of the velocity 
!--             gradients dvdx and dwdx
!--             Note: The validity of the new method has not yet 
!--                   been shown, as so far no suitable data for a 
!--                   validation has been available
                CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                    vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp )
                CALL wall_fluxes_e( i, j, k, nzb_diff_s_outer(j,i)-2, &
                                    wsus, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp )
                km_neutral = kappa * ( vsus(k)**2 + wsus(k)**2 )**0.25_wp * & 
                             0.5_wp * dx
                IF ( km_neutral > 0.0_wp )  THEN
                   dvdx = - wall_e_x(j,i) * vsus(k) / km_neutral
                   dwdx = - wall_e_x(j,i) * wsus(k) / km_neutral
                ELSE
                   dvdx = 0.0_wp
                   dwdx = 0.0_wp
                ENDIF
             ELSE
                dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                   v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                   w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             ENDIF

             def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                   dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                   dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

             IF ( def < 0.0_wp )  def = 0.0_wp

             tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

!
!--          (3) - will be executed only, if there is at least one level
!--          between (2) and (4), i.e. the topography must have a
!--          minimum height of 2 dz. Wall fluxes for this case have
!--          already been calculated for (2).
!--          'wall only: use wall functions'
             DO  k = nzb_diff_s_inner(j,i), nzb_diff_s_outer(j,i)-2

                dudx = ( u(k,j,i+1) - u(k,j,i) ) * ddx
                dudz = 0.5_wp * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                  u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
                dvdy =          ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                dvdz = 0.5_wp * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                  v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
                dwdz = ( w(k,j,i) - w(k-1,j,i) ) * ddzw(k)

                IF ( wall_e_y(j,i) /= 0.0_wp )  THEN
!
!--                Inconsistency removed: as the thermal stratification 
!--                is not taken into account for the evaluation of the 
!--                wall fluxes at vertical walls, the eddy viscosity km 
!--                must not be used for the evaluation of the velocity 
!--                gradients dudy and dwdy
!--                Note: The validity of the new method has not yet 
!--                      been shown, as so far no suitable data for a 
!--                      validation has been available
                   km_neutral = kappa * ( usvs(k)**2 + & 
                                          wsvs(k)**2 )**0.25_wp * 0.5_wp * dy
                   IF ( km_neutral > 0.0_wp )  THEN
                      dudy = - wall_e_y(j,i) * usvs(k) / km_neutral
                      dwdy = - wall_e_y(j,i) * wsvs(k) / km_neutral
                   ELSE
                      dudy = 0.0_wp
                      dwdy = 0.0_wp
                   ENDIF
                ELSE
                   dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                      u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                   dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                      w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                ENDIF

                IF ( wall_e_x(j,i) /= 0.0_wp )  THEN
!
!--                Inconsistency removed: as the thermal stratification 
!--                is not taken into account for the evaluation of the 
!--                wall fluxes at vertical walls, the eddy viscosity km 
!--                must not be used for the evaluation of the velocity 
!--                gradients dvdx and dwdx
!--                Note: The validity of the new method has not yet 
!--                      been shown, as so far no suitable data for a 
!--                      validation has been available
                   km_neutral = kappa * ( vsus(k)**2 + & 
                                          wsus(k)**2 )**0.25_wp * 0.5_wp * dx
                   IF ( km_neutral > 0.0_wp )  THEN
                      dvdx = - wall_e_x(j,i) * vsus(k) / km_neutral
                      dwdx = - wall_e_x(j,i) * wsus(k) / km_neutral
                   ELSE
                      dvdx = 0.0_wp
                      dwdx = 0.0_wp
                   ENDIF
                ELSE
                   dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                      v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                   dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                      w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                ENDIF

                def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                      dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                      dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

                IF ( def < 0.0_wp )  def = 0.0_wp

                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

             ENDDO

!
!--          (4) - will allways be executed.
!--          'special case: free atmosphere' (as for case (0))
             k = nzb_diff_s_outer(j,i)-1

             dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                 u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                 u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

             dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                 v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                 v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

             dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                 w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                 w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

             def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +        &
                   dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                   dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

             IF ( def < 0.0_wp )  def = 0.0_wp

             tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

          ELSE

!
!--          Position without adjacent wall
!--          (1) - will allways be executed. 
!--          'bottom only: use u_0,v_0'
             k = nzb_diff_s_inner(j,i)-1

             dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                 u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                 u_0(j,i)   - u_0(j,i+1)   ) * dd2zu(k)

             dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                 v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                 v_0(j,i)   - v_0(j+1,i)   ) * dd2zu(k)

             dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                 w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                 w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

             def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 )                       &
                   + dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + dvdz**2 &
                   + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

             IF ( def < 0.0_wp )  def = 0.0_wp

             tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

          ENDIF

       ELSEIF ( use_surface_fluxes )  THEN

          k = nzb_diff_s_outer(j,i)-1

          dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
          dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                              u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
          dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                              u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

          dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                              v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
          dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
          dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                              v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

          dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                              w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
          dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                              w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
          dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

          def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                dudy**2 + dvdx**2 + dwdx**2 + dwdy**2 + dudz**2 + &
                dvdz**2 + 2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )

          IF ( def < 0.0_wp )  def = 0.0_wp

          tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

       ENDIF

!
!--    If required, calculate TKE production by buoyancy
       IF ( .NOT. neutral )  THEN

          IF ( .NOT. humidity )  THEN

             IF ( use_single_reference_value )  THEN

                IF ( ocean )  THEN
!
!--                So far in the ocean no special treatment of density flux in
!--                the bottom and top surface layer
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) +                   &
                                    kh(k,j,i) * g / rho_reference * &
                                    ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * dd2zu(k)
                   ENDDO

                ELSE

                   DO  k = nzb_diff_s_inner(j,i), nzt_diff
                      tend(k,j,i) = tend(k,j,i) -                  &
                                    kh(k,j,i) * g / pt_reference * &
                                    ( pt(k+1,j,i) - pt(k-1,j,i) ) * dd2zu(k)
                   ENDDO

                   IF ( use_surface_fluxes )  THEN
                      k = nzb_diff_s_inner(j,i)-1
                      tend(k,j,i) = tend(k,j,i) + g / pt_reference * shf(j,i)
                   ENDIF

                   IF ( use_top_fluxes )  THEN
                      k = nzt
                      tend(k,j,i) = tend(k,j,i) + g / pt_reference * tswst(j,i)
                   ENDIF

                ENDIF

             ELSE

                IF ( ocean )  THEN
!
!--                So far in the ocean no special treatment of density flux in
!--                the bottom and top surface layer
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) +                &
                                    kh(k,j,i) * g / rho_ocean(k,j,i) * &
                                    ( rho_ocean(k+1,j,i) - rho_ocean(k-1,j,i) ) * dd2zu(k)
                   ENDDO

                ELSE

                   DO  k = nzb_diff_s_inner(j,i), nzt_diff
                      tend(k,j,i) = tend(k,j,i) -               &
                                    kh(k,j,i) * g / pt(k,j,i) * &
                                    ( pt(k+1,j,i) - pt(k-1,j,i) ) * dd2zu(k)
                   ENDDO

                   IF ( use_surface_fluxes )  THEN
                      k = nzb_diff_s_inner(j,i)-1
                      tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * shf(j,i)
                   ENDIF

                   IF ( use_top_fluxes )  THEN
                      k = nzt
                      tend(k,j,i) = tend(k,j,i) + g / pt(k,j,i) * tswst(j,i)
                   ENDIF

                ENDIF

             ENDIF

          ELSE

             DO  k = nzb_diff_s_inner(j,i), nzt_diff

                IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                   tend(k,j,i) = tend(k,j,i) - kh(k,j,i) * g / vpt(k,j,i) *   &
                                         ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) + &
                                           k2 * ( q(k+1,j,i) - q(k-1,j,i) )   &
                                         ) * dd2zu(k)
                ELSE IF ( cloud_physics )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                      temp  = theta * t_d_pt(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                 &
                                    ( q(k,j,i) - ql(k,j,i) ) *          &
                           ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /        &
                           ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *          &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                   tend(k,j,i) = tend(k,j,i) - kh(k,j,i) * g / vpt(k,j,i) *   &
                                         ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) + &
                                           k2 * ( q(k+1,j,i) - q(k-1,j,i) )   &
                                         ) * dd2zu(k)
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                   tend(k,j,i) = tend(k,j,i) - kh(k,j,i) * g / vpt(k,j,i) *  &
                                     ( k1 * ( pt(k+1,j,i)-pt(k-1,j,i) ) +    &
                                       k2 * ( q(k+1,j,i) - q(k-1,j,i) ) -    &
                                       pt(k,j,i) * ( ql(k+1,j,i) -           &
                                                     ql(k-1,j,i) ) ) * dd2zu(k)
                ENDIF
             ENDDO

             IF ( use_surface_fluxes )  THEN
                k = nzb_diff_s_inner(j,i)-1

                IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE IF ( cloud_physics )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                      temp  = theta * t_d_pt(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                 &
                                    ( q(k,j,i) - ql(k,j,i) ) *          &
                           ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /        &
                           ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *          &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ENDIF

                tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
                                            ( k1* shf(j,i) + k2 * qsws(j,i) )
             ENDIF

             IF ( use_top_fluxes )  THEN
                k = nzt

                IF ( .NOT. cloud_physics .AND. .NOT. cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE IF ( cloud_physics )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + pt_d_t(k) * l_d_cp * ql(k,j,i)
                      temp  = theta * t_d_pt(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                 &
                                    ( q(k,j,i) - ql(k,j,i) ) *          &
                           ( 1.0_wp + 0.622_wp * l_d_r / temp ) ) /        &
                           ( 1.0_wp + 0.622_wp * l_d_r * l_d_cp *          &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( l_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ENDIF

                tend(k,j,i) = tend(k,j,i) + g / vpt(k,j,i) * &
                                            ( k1* tswst(j,i) + k2 * qswst(j,i) )
             ENDIF

          ENDIF

       ENDIF

    END SUBROUTINE production_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE production_e_init

       USE arrays_3d,                                                          &
           ONLY:  kh, km, u, us, usws, v, vsws, zu

       USE control_parameters,                                                 &
           ONLY:  constant_flux_layer, kappa

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb_u_inner,     &
                  nzb_v_inner

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !<
       INTEGER(iwp) ::  j   !<
       INTEGER(iwp) ::  ku  !<
       INTEGER(iwp) ::  kv  !<

       IF ( constant_flux_layer )  THEN

          IF ( first_call )  THEN
             ALLOCATE( u_0(nysg:nyng,nxlg:nxrg), v_0(nysg:nyng,nxlg:nxrg) )
             u_0 = 0.0_wp   ! just to avoid access of uninitialized memory
             v_0 = 0.0_wp   ! within exchange_horiz_2d
             first_call = .FALSE.
          ENDIF

!
!--       Calculate a virtual velocity at the surface in a way that the
!--       vertical velocity gradient at k = 1 (u(k+1)-u_0) matches the
!--       Prandtl law (-w'u'/km). This gradient is used in the TKE shear
!--       production term at k=1 (see production_e_ij).
!--       The velocity gradient has to be limited in case of too small km
!--       (otherwise the timestep may be significantly reduced by large
!--       surface winds).
!--       Upper bounds are nxr+1 and nyn+1 because otherwise these values are
!--       not available in case of non-cyclic boundary conditions.
!--       WARNING: the exact analytical solution would require the determination
!--                of the eddy diffusivity by km = u* * kappa * zp / phi_m.
          !$OMP PARALLEL DO PRIVATE( ku, kv )
          DO  i = nxl, nxr+1
             DO  j = nys, nyn+1

                ku = nzb_u_inner(j,i)+1
                kv = nzb_v_inner(j,i)+1

                u_0(j,i) = u(ku+1,j,i) + usws(j,i) * ( zu(ku+1) - zu(ku-1) ) / &
                                 ( 0.5_wp * ( km(ku,j,i) + km(ku,j,i-1) ) +    &
                                   1.0E-20_wp )
!                                  ( us(j,i) * kappa * zu(1) )
                v_0(j,i) = v(kv+1,j,i) + vsws(j,i) * ( zu(kv+1) - zu(kv-1) ) / &
                                 ( 0.5_wp * ( km(kv,j,i) + km(kv,j-1,i) ) +    &
                                   1.0E-20_wp )
!                                  ( us(j,i) * kappa * zu(1) )

                IF ( ABS( u(ku+1,j,i) - u_0(j,i) )  > &
                     ABS( u(ku+1,j,i) - u(ku-1,j,i) ) )  u_0(j,i) = u(ku-1,j,i)
                IF ( ABS( v(kv+1,j,i) - v_0(j,i) )  > &
                     ABS( v(kv+1,j,i) - v(kv-1,j,i) ) )  v_0(j,i) = v(kv-1,j,i)

             ENDDO
          ENDDO

          CALL exchange_horiz_2d( u_0 )
          CALL exchange_horiz_2d( v_0 )

       ENDIF

    END SUBROUTINE production_e_init

 END MODULE production_e_mod
