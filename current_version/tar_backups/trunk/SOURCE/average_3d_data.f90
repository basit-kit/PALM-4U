!> @file average_3d_data.f90
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
! 
! Former revisions:
! -----------------
! $Id: average_3d_data.f90 2032 2016-10-21 15:13:51Z knoop $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! introduced control parameter varnamelength for LEN of trimvar.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added support for new urban surface model (temporary modifications of 
! SELECT CASE ( ) necessary, see variable trimvar),
! added comments in variable declaration section
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1972 2016-07-26 07:52:02Z maronga
! Output of land surface quantities is now done directly in the respective module
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length and radiative heating rates for RRTMG. 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
! 
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for land surface and radiation model parameters.
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
! 1318 2014-03-17 13:35:16Z raasch
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! +qc
!
! 1053 2012-11-13 17:11:03Z hoffmann
! averaging of nr, qr added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h_av
!
! Revision 1.1  2006/02/23 09:48:58  raasch
! Initial revision
!
!
! Description:
! ------------
!> Time-averaging of 3d-data-arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE average_3d_data
 

    USE averaging

    USE control_parameters,                                                    &
        ONLY:  average_count_3d, doav, doav_n, urban_surface, varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE land_surface_model_mod,                                                &
        ONLY:  land_surface, lsm_3d_data_averaging

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_3d_data_averaging

    USE urban_surface_mod,                                                     &
        ONLY:  usm_average_3d_data

#ifdef KPP_CHEM
    USE kchem_driver,                                                          &
        ONLY: chem_species, chem_species_av, kchem_integrate, NSPEC, NVAR,     &   !bK NVAR, SPC_NAMES added bk pe1
              use_kpp_chemistry, SPC_NAMES, kchem_3d_data_averaging
    USE arrays_3d,                                                             &   !bK added moduel arrays_3d  
        ONLY: rssws, rsswst, rs_p, rs, trs_m
#endif


    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< running index
    INTEGER(iwp) ::  ii !< running index
    INTEGER(iwp) ::  j  !< running index
    INTEGER(iwp) ::  k  !< running index

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(35),'average_3d_data','start')

!
!-- Check, if averaging is necessary
    IF ( average_count_3d <= 1 )  RETURN

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n

!
!--    Temporary solution to account for data output within the new urban 
!--    surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
       trimvar = TRIM( doav(ii) )
       IF ( urban_surface  .AND.  trimvar(1:4) == 'usm_' )  THEN
          trimvar = 'usm_output'
       ENDIF

        print*,' fm #21 avg_3d_data, trimvar conents are: ',trimvar     !bK debug
!   
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )

          CASE ( 'e' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      e_av(k,j,i) = e_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qsws*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_av(j,i) = qsws_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'lpt' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      lpt_av(k,j,i) = lpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'lwp*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lwp_av(j,i) = lwp_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'nr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      nr_av(k,j,i) = nr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

         CASE ( 'ol*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ol_av(j,i) = ol_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'p' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      p_av(k,j,i) = p_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'pc' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      pc_av(k,j,i) = pc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'pr' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      pr_av(k,j,i) = pr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'prr*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   precipitation_rate_av(j,i) = precipitation_rate_av(j,i) /   &
                                                REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'pt' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      pt_av(k,j,i) = pt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'q' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      q_av(k,j,i) = q_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qc' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      qc_av(k,j,i) = qc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'ql' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      ql_av(k,j,i) = ql_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'ql_c' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      ql_c_av(k,j,i) = ql_c_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'ql_v' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      ql_v_av(k,j,i) = ql_v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'ql_vp' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      ql_vp_av(k,j,i) = ql_vp_av(k,j,i) /                      &
                                        REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      qr_av(k,j,i) = qr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qv' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      qv_av(k,j,i) = qv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rho_ocean' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 's' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      s_av(k,j,i) = s_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

!          CASE ( 'rs' )                                                             !bK added this CASE block
!           CASE('k_NO','k_NO2','kc_O3','kc_RCHO')
!            print*,'fm #21, ln343, val of i is  ',i
!             DO  i = nxlg, nxrg
!                DO  j = nysg, nyng
!                   DO  k = nzb, nzt+1
!                      rs_av(k,j,i) = rs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
!                   ENDDO
!                ENDDO
!             ENDDO

          CASE ( 'sa' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      sa_av(k,j,i) = sa_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

         CASE ( 'shf*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   shf_av(j,i) = shf_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'ssws*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ssws_av(j,i) = ssws_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 't*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ts_av(j,i) = ts_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'u' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      u_av(k,j,i) = u_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'u*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   us_av(j,i) = us_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'v' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      v_av(k,j,i) = v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'vpt' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      vpt_av(k,j,i) = vpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'w' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      w_av(k,j,i) = w_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'z0*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   z0_av(j,i) = z0_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'z0h*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   z0h_av(j,i) = z0h_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO
!             
!--       Block of urban surface model outputs   
          CASE ( 'usm_output' )
             CALL usm_average_3d_data( 'average', doav(ii) )

          CASE DEFAULT
!
!--          Land surface quantity
             IF ( land_surface )  THEN
                CALL lsm_3d_data_averaging( 'average', doav(ii) )
             ENDIF

!
!--          Radiation quantity
             IF ( radiation )  THEN
                CALL radiation_3d_data_averaging( 'average', doav(ii) )
             ENDIF

!
!--         kchem chemisry
#ifdef KPP_CHEM
            IF ( use_kpp_chemistry ) THEN
                 CALL kchem_3d_data_averaging('average',doav(ii) )
            ENDIF
#endif
!
!--          User-defined quantity
             CALL user_3d_data_averaging( 'average', doav(ii) )

       END SELECT

    ENDDO

!
!-- Reset the counter
    average_count_3d = 0.0

    CALL cpu_log( log_point(35), 'average_3d_data', 'stop' )


 END SUBROUTINE average_3d_data
