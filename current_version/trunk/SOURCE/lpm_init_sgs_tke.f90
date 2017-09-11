!> @file lpm_init_sgs_tke.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: lpm_init_sgs_tke.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1929 2016-06-09 16:25:25Z suehring
! sgs_wfu_par, sgs_wfv_par and sgs_wfw_par are replaced by sgs_wf_par
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
! Kind definition added to all floating point numbers.
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Calculates quantities required for considering the SGS velocity fluctuations
!> in the particle transport by a stochastic approach. The respective
!> quantities are: SGS-TKE gradients and horizontally averaged profiles of the
!> SGS TKE and the resolved-scale velocity variances. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_init_sgs_tke
 

    USE arrays_3d,                                                             &
        ONLY:  de_dx, de_dy, de_dz, diss, e, u, v, w, zu

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nbgp, ngp_2dh_outer, nxl, nxr, nyn, nys, nzb, nzb_s_inner,      &
               nzb_s_outer, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  sgs_wf_part

    USE pegrid

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, sums, sums_l

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  k      !<


!
!-- TKE gradient along x and y
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1

             IF ( k <= nzb_s_inner(j,i-1)  .AND.  k > nzb_s_inner(j,i)  .AND.  &
                  k  > nzb_s_inner(j,i+1) )                                    &
             THEN
                de_dx(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i+1) - e(k,j,i) ) * ddx
             ELSEIF ( k  > nzb_s_inner(j,i-1)  .AND.  k > nzb_s_inner(j,i)     &
                      .AND.  k <= nzb_s_inner(j,i+1) )                         &
             THEN
                de_dx(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i) - e(k,j,i-1) ) * ddx
             ELSEIF ( k < nzb_s_inner(j,i)  .AND.  k < nzb_s_inner(j,i+1) )    &
             THEN
                de_dx(k,j,i) = 0.0_wp
             ELSEIF ( k < nzb_s_inner(j,i-1)  .AND.  k < nzb_s_inner(j,i) )    &
             THEN
                de_dx(k,j,i) = 0.0_wp
             ELSE
                de_dx(k,j,i) = sgs_wf_part * ( e(k,j,i+1) - e(k,j,i-1) ) * ddx
             ENDIF

             IF ( k <= nzb_s_inner(j-1,i)  .AND.  k > nzb_s_inner(j,i)  .AND.  &
                  k  > nzb_s_inner(j+1,i) )                                    &
             THEN
                de_dy(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j+1,i) - e(k,j,i) ) * ddy
             ELSEIF ( k  > nzb_s_inner(j-1,i)  .AND.  k  > nzb_s_inner(j,i)    &
                      .AND.  k <= nzb_s_inner(j+1,i) )                         &
             THEN
                de_dy(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i) - e(k,j-1,i) ) * ddy
             ELSEIF ( k < nzb_s_inner(j,i)  .AND.  k < nzb_s_inner(j+1,i) )    &
             THEN
                de_dy(k,j,i) = 0.0_wp
             ELSEIF ( k < nzb_s_inner(j-1,i)  .AND.  k < nzb_s_inner(j,i) )    &
             THEN
                de_dy(k,j,i) = 0.0_wp
             ELSE
                de_dy(k,j,i) = sgs_wf_part * ( e(k,j+1,i) - e(k,j-1,i) ) * ddy
             ENDIF

          ENDDO
       ENDDO
    ENDDO

!
!-- TKE gradient along z, including bottom and top boundary conditions
    DO  i = nxl, nxr
       DO  j = nys, nyn

          DO  k = nzb_s_inner(j,i)+2, nzt-1
             de_dz(k,j,i)  = 2.0_wp * sgs_wf_part *                            &
                             ( e(k+1,j,i) - e(k-1,j,i) ) / ( zu(k+1)-zu(k-1) )
          ENDDO

          k = nzb_s_inner(j,i)
          de_dz(nzb:k,j,i) = 0.0_wp
          de_dz(k+1,j,i)   = 2.0_wp * sgs_wf_part *                            &
                             ( e(k+2,j,i) - e(k+1,j,i) ) / ( zu(k+2) - zu(k+1) )
          de_dz(nzt,j,i)   = 0.0_wp
          de_dz(nzt+1,j,i) = 0.0_wp
       ENDDO
    ENDDO


!
!-- Lateral boundary conditions
    CALL exchange_horiz( de_dx, nbgp )
    CALL exchange_horiz( de_dy, nbgp )
    CALL exchange_horiz( de_dz, nbgp )
    CALL exchange_horiz( diss, nbgp  )


!
!-- Calculate the horizontally averaged profiles of SGS TKE and resolved
!-- velocity variances (they may have been already calculated in routine
!-- flow_statistics).
    IF ( .NOT. flow_statistics_called )  THEN

!
!--    First calculate horizontally averaged profiles of the horizontal
!--    velocities.
       sums_l(:,1,0) = 0.0_wp
       sums_l(:,2,0) = 0.0_wp

       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb_s_outer(j,i), nzt+1
                sums_l(k,1,0)  = sums_l(k,1,0)  + u(k,j,i)
                sums_l(k,2,0)  = sums_l(k,2,0)  + v(k,j,i)
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,1,0), sums(nzb,1), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,2,0), sums(nzb,2), nzt+2-nzb, &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       sums(:,1) = sums_l(:,1,0)
       sums(:,2) = sums_l(:,2,0)
#endif

!
!--    Final values are obtained by division by the total number of grid
!--    points used for the summation.
       hom(:,1,1,0) = sums(:,1) / ngp_2dh_outer(:,0)   ! u
       hom(:,1,2,0) = sums(:,2) / ngp_2dh_outer(:,0)   ! v

!
!--    Now calculate the profiles of SGS TKE and the resolved-scale
!--    velocity variances
       sums_l(:,8,0)  = 0.0_wp
       sums_l(:,30,0) = 0.0_wp
       sums_l(:,31,0) = 0.0_wp
       sums_l(:,32,0) = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_outer(j,i), nzt+1
                sums_l(k,8,0)  = sums_l(k,8,0)  + e(k,j,i)
                sums_l(k,30,0) = sums_l(k,30,0) + ( u(k,j,i) - hom(k,1,1,0) )**2
                sums_l(k,31,0) = sums_l(k,31,0) + ( v(k,j,i) - hom(k,1,2,0) )**2
                sums_l(k,32,0) = sums_l(k,32,0) + w(k,j,i)**2
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,8,0), sums(nzb,8), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,30,0), sums(nzb,30), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,31,0), sums(nzb,31), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,32,0), sums(nzb,32), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )

#else
       sums(:,8)  = sums_l(:,8,0)
       sums(:,30) = sums_l(:,30,0)
       sums(:,31) = sums_l(:,31,0)
       sums(:,32) = sums_l(:,32,0)
#endif

!
!--    Final values are obtained by division by the total number of grid
!--    points used for the summation.
       hom(:,1,8,0)  = sums(:,8)  / ngp_2dh_outer(:,0)   ! e
       hom(:,1,30,0) = sums(:,30) / ngp_2dh_outer(:,0)   ! u*2
       hom(:,1,31,0) = sums(:,31) / ngp_2dh_outer(:,0)   ! v*2 
       hom(:,1,32,0) = sums(:,32) / ngp_2dh_outer(:,0)   ! w*2

    ENDIF

 END SUBROUTINE lpm_init_sgs_tke
