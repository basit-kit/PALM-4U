!> @file read_3d_binary.f90
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
! $Id: read_3d_binary.f90 2032 2016-10-21 15:13:51Z knoop $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho_av to rho_ocean_av
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Bugfix: read of land surface data only when module is switched on.
! Radiation parts are now done in the respective module.
! Binary version increased to 4.5.
! 
! 1972 2016-07-26 07:52:02Z maronga
! Land surface parts are now done in the respective module
! 
! 1849 2016-04-08 11:33:18Z hoffmann 
! prr, precipitation_amount moved to arrays_3d
!
! 1833 2016-04-07 14:23:03Z raasch
! statistics module replaced by spectra module
!
! 1808 2016-04-05 19:44:00Z raasch
! test output removed
!
! 1788 2016-03-10 11:01:04Z maronga
! Added z0q and z0q_av
! 
! 1757 2016-02-22 15:49:32Z maronga
! Bugfix in allocation of radiative heating rates
! 
! 1709 2015-11-04 14:47:01Z maronga
! Added rad_lw_out_change_0, increased binary_version
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of radiative heating rates and Obukhov length. Removed output of
! rif.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for binary input of land surface and radiation model data. In the
! course of this work, new temporary arrays tmp_3d_soil1, tmp_3d_soil2 were
! required.
! 
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
! 
! 1400 2014-05-09 14:03:54Z knoop
! reading of arrays for random_generator_parallel added
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
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr) 
! of the two-moment cloud physics scheme:
! +prr, prr_av, *, *_av, *s, *sws, *swst
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h, z0h_av
!
! Revision 1.1  2004/04/30 12:47:27  raasch
! Initial revision
!
!
! Description:
! ------------
!> Binary input of variables and arrays from restart file
!------------------------------------------------------------------------------!
 SUBROUTINE read_3d_binary
 

    USE arrays_3d,                                                             &
        ONLY:  e, kh, km, ol, p, pt, q, ql, qc, nr, nrs, nrsws, nrswst,        &
               prr, precipitation_amount, qr,                                  &
               qrs, qrsws, qrswst, qs, qsws, qswst, rs, rss, rssws, rsswst, s, sa, saswsb, saswst,     &  !bK added rs, rss, rssws, rsswst
               ss, ssws, sswst, rif_wall, shf, ss, ts, tswst, u, u_m_l, u_m_n, &
               u_m_r, u_m_s, us, usws, uswst, v, v_m_l, v_m_n, v_m_r, v_m_s,   &
               vpt, vsws, vswst, w, w_m_l, w_m_n, w_m_r, w_m_s, z0, z0h, z0q

    USE averaging

    USE control_parameters,                                                    &
        ONLY:  iran, message_string, outflow_l, outflow_n, outflow_r, outflow_s, &
               chemistry 

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nx_on_file, ny, nys, nysg, nyn, &
               nyng, ny_on_file, nzb, nzt

    USE kinds

    USE land_surface_model_mod,                                                &
        ONLY:  land_surface, lsm_read_restart_data

    USE particle_attributes,                                                   &
        ONLY:  iran_part

    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY: radiation, radiation_read_restart_data

    USE random_function_mod,                                                   &
        ONLY:  random_iv, random_iy

    USE random_generator_parallel,                                             &
        ONLY:  id_random_array, seq_random_array

    USE spectra_mod,                                                           &
        ONLY:  spectrum_x, spectrum_y

    USE kchem_driver,                                                          &
        ONLY: chem_species, NSPEC, use_kpp_chemistry, kchem_read_restart_data


    IMPLICIT NONE

    CHARACTER (LEN=7)  ::  myid_char_save
    CHARACTER (LEN=10) ::  binary_version
    CHARACTER (LEN=10) ::  version_on_file
    CHARACTER (LEN=20) ::  field_chr

    INTEGER(iwp) ::  files_to_be_opened  !<
    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  myid_on_file        !<
    INTEGER(iwp) ::  numprocs_on_file    !<
    INTEGER(iwp) ::  nxlc                !<
    INTEGER(iwp) ::  nxlf                !<
    INTEGER(iwp) ::  nxlpr               !<
    INTEGER(iwp) ::  nxl_on_file         !<
    INTEGER(iwp) ::  nxrc                !<
    INTEGER(iwp) ::  nxrf                !<
    INTEGER(iwp) ::  nxrpr               !<
    INTEGER(iwp) ::  nxr_on_file         !<
    INTEGER(iwp) ::  nync                !<
    INTEGER(iwp) ::  nynf                !<
    INTEGER(iwp) ::  nynpr               !<
    INTEGER(iwp) ::  nyn_on_file         !<
    INTEGER(iwp) ::  nysc                !<
    INTEGER(iwp) ::  nysf                !<
    INTEGER(iwp) ::  nyspr               !<
    INTEGER(iwp) ::  nys_on_file         !<
    INTEGER(iwp) ::  nzb_on_file         !<
    INTEGER(iwp) ::  nzt_on_file         !<
    INTEGER(iwp) ::  offset_x            !<
    INTEGER(iwp) ::  offset_y            !<
    INTEGER(iwp) ::  shift_x             !<
    INTEGER(iwp) ::  shift_y             !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  file_list       !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  overlap_count   !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxlfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxrfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nynfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nysfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_xa  !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_ya  !<

#ifdef KPP_CHEM
    INTEGER(iwp) ::  kc_n                   !<   
#endif 

    REAL(wp) ::  rdummy

    REAL(wp), DIMENSION(:,:), ALLOCATABLE     ::  tmp_2d      !< temporary array for storing 2D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3d      !< temporary array for storing 3D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwul   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwun   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwur   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwus   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvl   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvn   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvr   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvs   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwl   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwn   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwr   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwws   !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tmp_4d     !<


!
!-- Read data from previous model run.
    CALL cpu_log( log_point_s(14), 'read_3d_binary', 'start' )

!
!-- Check which of the restart files contain data needed for the subdomain
!-- of this PE
    files_to_be_opened = 0

    DO  i = 1, numprocs_previous_run
!
!--    Store array bounds of the previous run ("pr") in temporary scalars
       nxlpr = hor_index_bounds_previous_run(1,i-1)
       nxrpr = hor_index_bounds_previous_run(2,i-1)
       nyspr = hor_index_bounds_previous_run(3,i-1)
       nynpr = hor_index_bounds_previous_run(4,i-1)

!
!--    Determine the offsets. They may be non-zero in case that the total domain
!--    on file is smaller than the current total domain.
       offset_x = ( nxl / ( nx_on_file + 1 ) ) * ( nx_on_file + 1 )
       offset_y = ( nys / ( ny_on_file + 1 ) ) * ( ny_on_file + 1 )

!
!--    Start with this offset and then check, if the subdomain on file
!--    matches another time(s) in the current subdomain by shifting it
!--    for nx_on_file+1, ny_on_file+1 respectively
    
       shift_y = 0
       j       = 0
       DO WHILE (  nyspr+shift_y <= nyn-offset_y )
          
          IF ( nynpr+shift_y >= nys-offset_y ) THEN 

             shift_x = 0
             DO WHILE ( nxlpr+shift_x <= nxr-offset_x )
                
                IF ( nxrpr+shift_x >= nxl-offset_x ) THEN
                   j = j +1
                   IF ( j > 1000 )  THEN
!
!--                   Array bound exceeded
                      message_string = 'data from subdomain of previous' // &
                                       ' run mapped more than 1000 times'
                      CALL message( 'read_3d_binary', 'PA0284', 2, 2, -1,   &
                                       6, 1 )
                   ENDIF

                   IF ( j == 1 )  THEN
                      files_to_be_opened = files_to_be_opened + 1
                      file_list(files_to_be_opened) = i-1
                   ENDIF
                      
                   offset_xa(files_to_be_opened,j) = offset_x + shift_x
                   offset_ya(files_to_be_opened,j) = offset_y + shift_y
!
!--                Index bounds of overlapping data
                   nxlfa(files_to_be_opened,j) = MAX( nxl-offset_x-shift_x, nxlpr )
                   nxrfa(files_to_be_opened,j) = MIN( nxr-offset_x-shift_x, nxrpr )
                   nysfa(files_to_be_opened,j) = MAX( nys-offset_y-shift_y, nyspr )
                   nynfa(files_to_be_opened,j) = MIN( nyn-offset_y-shift_y, nynpr )

                ENDIF

                shift_x = shift_x + ( nx_on_file + 1 )
             ENDDO
       
          ENDIF
             
          shift_y = shift_y + ( ny_on_file + 1 )             
       ENDDO
          
       IF ( j > 0 )  overlap_count(files_to_be_opened) = j
          
    ENDDO
   
!
!-- Save the id-string of the current process, since myid_char may now be used
!-- to open files created by PEs with other id.
    myid_char_save = myid_char

    IF ( files_to_be_opened /= 1  .OR.  numprocs /= numprocs_previous_run ) &
    THEN
       WRITE( message_string, * ) 'number of PEs or virtual PE-grid changed ', &
                        'in restart run&  PE', myid, ' will read from files ', &
                         file_list(1:files_to_be_opened)
       CALL message( 'read_3d_binary', 'PA0285', 0, 0, 0, 6, 0 )
    ENDIF

!
!-- Read data from all restart files determined above
    DO  i = 1, files_to_be_opened

        print*,'41@331, files to be opened .. i .. ',i, 'and files_to_be_opened is .. ',files_to_be_opened     !bK
       j = file_list(i)
!
!--    Set the filename (underscore followed by four digit processor id)
       WRITE (myid_char,'(''_'',I6.6)')  j

!
!--    Open the restart file. If this file has been created by PE0 (_000000),
!--    the global variables at the beginning of the file have to be skipped
!--    first.
       CALL check_open( 13 )
       IF ( j == 0 )  CALL skip_var_list

!
!--    First compare the version numbers
       READ ( 13 )  version_on_file
       binary_version = '4.5'
       IF ( TRIM( version_on_file ) /= TRIM( binary_version ) )  THEN
          WRITE( message_string, * ) 'version mismatch concerning data ',      &
                      'from prior run',                                        &
                      '&version on file    = "', TRIM( version_on_file ), '"', &
                      '&version in program = "', TRIM( binary_version ), '"'
          CALL message( 'read_3d_binary', 'PA0286', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Read number of processors, processor-id, and array ranges.
!--    Compare the array ranges with those stored in the index bound array.
       READ ( 13 )  numprocs_on_file, myid_on_file, nxl_on_file, nxr_on_file, &
                    nys_on_file, nyn_on_file, nzb_on_file, nzt_on_file

       IF ( nxl_on_file /= hor_index_bounds_previous_run(1,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nxl on ',  &
                            'restart file "', myid_char, '"',             &
                            '&nxl = ', nxl_on_file, ' but it should be',  &
                            '&= ', hor_index_bounds_previous_run(1,j),    &
                            '&from the index bound information array'
          CALL message( 'read_3d_binary', 'PA0287', 2, 2, -1, 6, 1 )
       ENDIF

       IF ( nxr_on_file /= hor_index_bounds_previous_run(2,j) )  THEN
           WRITE( message_string, * ) 'problem with index bound nxr on ',   &
                               'restart file "', myid_char, '"'  ,          &
                               '&nxr = ', nxr_on_file, ' but it should be', &
                               '&= ', hor_index_bounds_previous_run(2,j),   &
                               '&from the index bound information array' 
          CALL message( 'read_3d_binary', 'PA0288', 2, 2, -1, 6, 1 )

       ENDIF

       IF ( nys_on_file /= hor_index_bounds_previous_run(3,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nys on ',      &
                                 'restart file "', myid_char, '"',            &
                                 '&nys = ', nys_on_file, ' but it should be', &
                                 '&= ', hor_index_bounds_previous_run(3,j),   &
                                     '&from the index bound information array'
          CALL message( 'read_3d_binary', 'PA0289', 2, 2, -1, 6, 1 ) 
       ENDIF

       IF ( nyn_on_file /= hor_index_bounds_previous_run(4,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nyn on ',    &
                               'restart file "', myid_char, '"',            &
                               '&nyn = ', nyn_on_file, ' but it should be', &
                               '&= ', hor_index_bounds_previous_run(4,j),   &
                               '&from the index bound information array'
          CALL message( 'read_3d_binary', 'PA0290', 2, 2, -1, 6, 1 ) 
       ENDIF

       IF ( nzb_on_file /= nzb )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     '&from prior run on PE ', myid,           &
                                     '&nzb on file = ', nzb_on_file,           &
                                     '&nzb         = ', nzb
          CALL message( 'read_3d_binary', 'PA0291', 1, 2, 0, 6, 0 )  
       ENDIF

       IF ( nzt_on_file /= nzt )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     '&from prior run on PE ', myid,           &
                                     '&nzt on file = ', nzt_on_file,           &
                                     '&nzt         = ', nzt
          CALL message( 'read_3d_binary', 'PA0292', 1, 2, 0, 6, 0 )  
       ENDIF

!
!--    Allocate temporary arrays sized as the arrays on the restart file
       ALLOCATE( tmp_2d(nys_on_file-nbgp:nyn_on_file+nbgp,                     &
                        nxl_on_file-nbgp:nxr_on_file+nbgp),                    &
                 tmp_3d(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,           &
                        nxl_on_file-nbgp:nxr_on_file+nbgp) )

!
!--    Read arrays
!--    ATTENTION: If the following read commands have been altered, the 
!--    ---------- version number of the variable binary_version must be altered,
!--               too. Furthermore, the output list of arrays in write_3d_binary
!--               must also be altered accordingly.
       READ ( 13 )  field_chr
       DO  WHILE ( TRIM( field_chr ) /= '*** end ***' )

!
!--       Map data on file as often as needed (data are read only for k=1)
          DO  k = 1, overlap_count(i)
            if(myid==0) print*,'fm #41@434 . field_chr ',field_chr,'.. k .. ',k,'..files_to_be_opened .. ',files_to_be_opened 

! 
!--          Get the index range of the subdomain on file which overlap with the
!--          current subdomain
             nxlf = nxlfa(i,k)
             nxlc = nxlfa(i,k) + offset_xa(i,k)
             nxrf = nxrfa(i,k)
             nxrc = nxrfa(i,k) + offset_xa(i,k)
             nysf = nysfa(i,k)
             nysc = nysfa(i,k) + offset_ya(i,k)
             nynf = nynfa(i,k)
             nync = nynfa(i,k) + offset_ya(i,k)


             SELECT CASE ( TRIM( field_chr ) )

                CASE ( 'e' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'e_av' )
                   IF ( .NOT. ALLOCATED( e_av ) )  THEN
                      ALLOCATE( e_av(nzb:nzt+1,nys-nbgp:nyn+nbgp,nxl-nbgp:nxr+nbgp) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                            tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'iran' ) ! matching random numbers is still unresolved
                                ! issue
                   IF ( k == 1 )  READ ( 13 )  iran, iran_part

                CASE ( 'kh' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   kh(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'km' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   km(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'lpt_av' )
                   IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                      ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   lpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'lwp_av' )
                   IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                      ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   lwp_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nr' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nr_av' )
                   IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                      ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nrs' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   nrs(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nrsws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   nrsws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nrswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   nrswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                CASE ( 'ol' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ol(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                 tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p_av' )
                   IF ( .NOT. ALLOCATED( p_av ) )  THEN
                      ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pc_av' )
                   IF ( .NOT. ALLOCATED( pc_av ) )  THEN
                      ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pr_av' )
                   IF ( .NOT. ALLOCATED( pr_av ) )  THEN
                      ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'prr' )
                   IF ( .NOT. ALLOCATED( prr ) )  THEN
                      ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   prr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'prr_av' )
                   IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                      ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   prr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'precipitation_amount' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   precipitation_amount(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'precipitation_rate_a' )
                   IF ( .NOT. ALLOCATED( precipitation_rate_av ) )  THEN
                      ALLOCATE( precipitation_rate_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   precipitation_rate_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt_av' )
                   IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                      ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'q' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   q(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'q_av' )
                   IF ( .NOT. ALLOCATED( q_av ) )  THEN
                      ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   q_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                     tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qc' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                       tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qc_av' )
                   IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                      ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                       tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                       tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_av' )
                   IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                      ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                       tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qr' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qr_av' )
                   IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                      ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qrs' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qrs(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qrsws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qrsws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qrswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qrswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_c_av' )
                   IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
                      ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_c_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                        tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_v_av' )
                   IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
                      ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                        tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_vp_av' )
                   IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
                      ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_vp_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                        tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qs' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qs(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws_av' )
                   IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                      ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qv_av' )
                   IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                      ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qv_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'random_iv' )  ! still unresolved issue
                   IF ( k == 1 )  READ ( 13 )  random_iv
                   IF ( k == 1 )  READ ( 13 )  random_iy
                   
                CASE ( 'seq_random_array' )  ! still unresolved issue
                   IF ( k == 1 )  READ ( 13 )  id_random_array
                   IF ( k == 1 )  READ ( 13 )  seq_random_array

                CASE ( 'rho_ocean_av' )
                   IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                      ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rho_ocean_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rif_wall' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_4d(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp, &
                                       nxl_on_file-nbgp:nxr_on_file+nbgp,1:4) )
                      READ ( 13 )  tmp_4d
                   ENDIF
                   rif_wall(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp,:) = &
                            tmp_4d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp,:)

                CASE ( 's' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   s(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                CASE ( 's_av' )
                   IF ( .NOT. ALLOCATED( s_av ) )  THEN
                      ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   s_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                 tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

!#ifdef KPP_CHEM
!                CASE ( field_chr(1:3)== 'kc_O3' )                                        !bK added this CASE block
!                   IF ( k == 1 )  READ ( 13 )  tmp_3d
!                    print*,'fm #41@755 read_3d_binar, reading rs data'                        !bK debug
!                   rs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
!                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!                 CASE ( 'kc_NO2' )                                                               !bK added this CASE block
!                   IF ( k == 1 )  READ ( 13 )  tmp_3d
!                    print*, 'fm #41@761 read_3d_binary NO2 is read in rs ..'                        !bK debug
!                    rs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
!                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!                CASE ( 'rs_av' )                                                            !bK added this CASE block
!                   IF ( .NOT. ALLOCATED( rs_av ) )  THEN
!                      ALLOCATE( rs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg))
!                   ENDIF
!                   IF ( k == 1 )  READ ( 13 )  tmp_3d
!                        write(9,*) 'read_3d_binary, reading rs_av'                          !bK debug
!                   rs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
!                                 tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!#endif

                CASE ( 'sa' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sa_av' )
                   IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                      ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'saswsb' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   saswsb(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'saswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   saswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_av' )
                   IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                      ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'spectrum_x' )
                   IF ( k == 1 )  THEN
                      IF ( nx_on_file /= nx )  THEN
                         message_string = 'read_3d_binary: spectrum_x ' // &
                                     'on restart file ignored because' // &
                                     '&total numbers of grid points (nx) ' // &
                                     'do not match'
                         CALL message( 'read_3d_binary', 'PA0293',&
                                                                 0, 1, 0, 6, 0 )
                         READ ( 13 )  rdummy
                      ELSE
                         READ ( 13 )  spectrum_x
                      ENDIF
                   ENDIF

                CASE ( 'spectrum_y' )
                   IF ( k == 1 )  THEN
                      IF ( ny_on_file /= ny )  THEN
                         message_string = 'read_3d_binary: spectrum_y ' //   &
                                     'on restart file ignored because' //    &
                                     '&total numbers of grid points (ny) '// &
                                     'do not match'
                         CALL message( 'read_3d_binary', 'PA0294', &
                                                                 0, 1, 0, 6, 0 )
                      READ ( 13 )  rdummy
                      ELSE
                         READ ( 13 )  spectrum_y
                      ENDIF
                   ENDIF
!#ifdef KPP_CHEM
!
!                CASE ( 'rss' )
!                   IF ( k == 1 )  READ ( 13 )  tmp_2d
!                   rss(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
!                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!                CASE ( 'rssws' )
!                   IF ( k == 1 )  READ ( 13 )  tmp_2d
!                   rssws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
!                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!                CASE ( 'rssws_av' )
!                   IF ( .NOT. ALLOCATED( rssws_av ) )  THEN
!                      ALLOCATE( rssws_av(nysg:nyng,nxlg:nxrg) )
!                   ENDIF  
!                   IF ( k == 1 )  READ ( 13 )  tmp_2d
!                   rssws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
!                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!                CASE ( 'rsswst' )
!                   IF ( k == 1 )  READ ( 13 )  tmp_2d
!                   rsswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
!                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp) 
!
!#endif

                   
                CASE ( 'ss' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ss(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ssws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ssws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ssws_av' )
                   IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                      ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ssws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   sswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp) 
                 
                CASE ( 'ts' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ts(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                     tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ts_av' )
                   IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                      ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ts_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'tswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   tswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u_av' )
                   IF ( .NOT. ALLOCATED( u_av ) )  THEN
                      ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwul(nzb:nzt+1, &
                                          nys_on_file-nbgp:nyn_on_file+nbgp,1:2) )
                      READ ( 13 )  tmp_3dwul
                   ENDIF
                   IF ( outflow_l )  THEN
                      u_m_l(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwul(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'u_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwun(nzb:nzt+1,ny-1:ny, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwun
                   ENDIF
                   IF ( outflow_n )  THEN
                      u_m_n(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwun(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'u_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwur(nzb:nzt+1,&
                                          nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwur
                   ENDIF
                   IF ( outflow_r )  THEN
                      u_m_r(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwur(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'u_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwus(nzb:nzt+1,0:1, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwus
                   ENDIF
                   IF ( outflow_s )  THEN
                      u_m_s(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwus(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'us' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   us(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                     tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'usws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   usws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'uswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   uswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'us_av' )
                   IF ( .NOT. ALLOCATED( us_av ) )  THEN
                      ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   us_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v_av' )
                   IF ( .NOT. ALLOCATED( v_av ) )  THEN
                      ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvl(nzb:nzt+1,&
                                          nys_on_file-nbgp:nyn_on_file+nbgp,0:1) )
                      READ ( 13 )  tmp_3dwvl
                   ENDIF
                   IF ( outflow_l )  THEN
                      v_m_l(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwvl(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'v_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvn(nzb:nzt+1,ny-1:ny, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwvn
                   ENDIF
                   IF ( outflow_n )  THEN
                      v_m_n(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwvn(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'v_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvr(nzb:nzt+1,&
                                          nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwvr
                   ENDIF
                   IF ( outflow_r )  THEN
                      v_m_r(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwvr(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'v_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvs(nzb:nzt+1,1:2, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwvs
                   ENDIF
                   IF ( outflow_s )  THEN
                      v_m_s(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwvs(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'vpt' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   vpt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'vpt_av' )
                   IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                      ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   vpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'vsws' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   vsws(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'vswst' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   vswst(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w_av' )
                   IF ( .NOT. ALLOCATED( w_av ) )  THEN
                      ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwl(nzb:nzt+1,&
                                          nys_on_file-nbgp:nyn_on_file+nbgp,0:1) )
                      READ ( 13 )  tmp_3dwwl
                   ENDIF
                   IF ( outflow_l )  THEN
                      w_m_l(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwwl(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'w_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwn(nzb:nzt+1,ny-1:ny, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwwn
                   ENDIF
                   IF ( outflow_n )  THEN
                      w_m_n(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwwn(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'w_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwr(nzb:nzt+1,&
                                          nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwwr
                   ENDIF
                   IF ( outflow_r )  THEN
                      w_m_r(:,nysc-nbgp:nync+nbgp,:) = tmp_3dwwr(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'w_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwws(nzb:nzt+1,0:1, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwws
                   ENDIF
                   IF ( outflow_s )  THEN
                      w_m_s(:,:,nxlc-nbgp:nxrc+nbgp) = tmp_3dwws(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF
                   DEALLOCATE( tmp_3dwws )

                CASE ( 'z0' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                     tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0_av' )
                   IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                      ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0h' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0h(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                     tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0h_av' )
                   IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                      ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0h_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0q' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0q(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                     tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0q_av' )
                   IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                      ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0q_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE DEFAULT
                   WRITE( message_string, * ) 'unknown field named "', &
                                              TRIM( field_chr ), '" found in', &
                                              '&data from prior run on PE ',myid
                    CALL message( 'read_3d_binary', 'PA0295', 1, 2, 0, 6, 0 )  
                   
             END SELECT

          ENDDO  ! overlap loop

!
!--       Deallocate arrays needed for specific variables only
          IF ( ALLOCATED( tmp_3dwul ) )  DEALLOCATE( tmp_3dwul )
          IF ( ALLOCATED( tmp_3dwun ) )  DEALLOCATE( tmp_3dwun )
          IF ( ALLOCATED( tmp_3dwur ) )  DEALLOCATE( tmp_3dwur )
          IF ( ALLOCATED( tmp_3dwus ) )  DEALLOCATE( tmp_3dwus )
          IF ( ALLOCATED( tmp_3dwvl ) )  DEALLOCATE( tmp_3dwvl )
          IF ( ALLOCATED( tmp_3dwvn ) )  DEALLOCATE( tmp_3dwvn )
          IF ( ALLOCATED( tmp_3dwvr ) )  DEALLOCATE( tmp_3dwvr )
          IF ( ALLOCATED( tmp_3dwvs ) )  DEALLOCATE( tmp_3dwvs )
          IF ( ALLOCATED( tmp_3dwwl ) )  DEALLOCATE( tmp_3dwwl )
          IF ( ALLOCATED( tmp_3dwwn ) )  DEALLOCATE( tmp_3dwwn )
          IF ( ALLOCATED( tmp_3dwwr ) )  DEALLOCATE( tmp_3dwwr )
          IF ( ALLOCATED( tmp_3dwws ) )  DEALLOCATE( tmp_3dwws )
          IF ( ALLOCATED( tmp_4d ) )  DEALLOCATE( tmp_4d )

!
!--       Read next character string
          READ ( 13 )  field_chr

print*,'fm #41@1188 within enddo loop.. i.. ',i,' .. files_to_be_opened .. ',files_to_be_opened
       ENDDO  ! loop over variables

print*,'fm #41@1191 bef calling kchem_restart..'
!
!--    Read land surface restart data
       IF ( land_surface )  THEN
          CALL lsm_read_restart_data( i, nxlfa, nxl_on_file, nxrfa,            &
                                      nxr_on_file, nynfa, nyn_on_file, nysfa,  &
                                      nys_on_file, offset_xa, offset_ya,       &
                                      overlap_count(i), tmp_2d )
       ENDIF

!
!--    Read radiation restart data
       IF ( radiation )  THEN
          CALL radiation_read_restart_data( i, nxlfa, nxl_on_file, nxrfa,      &
                                            nxr_on_file, nynfa, nyn_on_file,   &
                                            nysfa, nys_on_file, offset_xa,     &
                                            offset_ya, overlap_count(i),       &
                                            tmp_2d, tmp_3d )
       ENDIF

#ifdef KPP_CHEM
        print*,'fm 41@1212 bef calling kchem_restart. chemistry is ..', chemistry

       IF ( use_kpp_chemistry )  THEN

          CALL kchem_read_restart_data( i, nxlfa, nxl_on_file, nxrfa,          &
                                            nxr_on_file, nynfa, nyn_on_file,   &
                                            nysfa, nys_on_file, offset_xa,     &
                                            offset_ya, overlap_count(i),tmp_3d )
        print*,'fm 41@1220 aft calling kchem_restart ..'
                                            
       ENDIF
        
!        DO kc_n = 1, NSPEC
!           if (field_chr == chem_species(kc_n)%name) then
!               if (k == 1 ) READ (13)  tmp_3d
!                   chem_species(kc_n)%conc = tmp_3d
!               
!           endif
!        END DO

#endif 



!
!--    Read user-defined restart data
       CALL user_read_restart_data( i, nxlfa, nxl_on_file, nxrfa, nxr_on_file, &
                                    nynfa, nyn_on_file, nysfa, nys_on_file,    &
                                    offset_xa, offset_ya, overlap_count(i),    &
                                    tmp_2d, tmp_3d )

!
!--    Close the restart file
       CALL close_file( 13 )

       DEALLOCATE( tmp_2d, tmp_3d )

    ENDDO  ! loop over restart files


!
!-- Restore the original filename for the restart file to be written
    myid_char = myid_char_save


!
!-- End of time measuring for reading binary data
    CALL cpu_log( log_point_s(14), 'read_3d_binary', 'stop' )

 END SUBROUTINE read_3d_binary
