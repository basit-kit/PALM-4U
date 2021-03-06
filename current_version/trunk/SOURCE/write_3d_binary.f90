!> @file write_3d_binary.f90
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
! $Id: write_3d_binary.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2386 2017-09-04 12:23:03Z basit
! renamed kchem_driver to chemistry_model_mod, use_kpp_chemistry to
! air_chemistry. prefix 'k' is removed from other chem vars
! and subroutine names. KPP_CHEM prep directive replaced with __chem.! 
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho_av to rho_ocean_av
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Radiation actions are now done directly in the respective module.
! Binary version increased to 4.5.
! 
! 1972 2016-07-26 07:52:02Z maronga
! Land surface actions are now done directly in the respective module
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! prr, precipitation_amount moved to arrays_3d
!
! 1833 2016-04-07 14:23:03Z raasch
! statistics module replaced by spectra module
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1788 2016-03-10 11:01:04Z maronga
! Added z0q and z0q_av
! 
! 1709 2015-11-04 14:47:01Z maronga
! Added rad_lw_out_change_0, increased binary_version
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of radiative heating rates for RRTMG. Added output of ol. Removed
! output of rif.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for binary ouput of land surface and radiation model data. 
! 
! 1400 2014-05-09 14:03:54Z knoop
! writing of arrays for random_generator_parallel added
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! Bugfix using cloud_droplets solved. qc, qr*, nr* are no longer written in case
! of cloud_droplets = .TRUE.
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! qr and nr are restricted to precipitation
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
! all actions concerning leapfrog- and ups-scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h, z0h_av
!
! Revision 1.1  1998/03/18 20:20:21  raasch
! Initial revision
!
!
! Description:
! ------------
!> Binary output of variables and arrays for restarts.
!------------------------------------------------------------------------------!
 SUBROUTINE write_3d_binary
 

    USE arrays_3d,                                                             &
        ONLY:  e, kh, km, ol, p, pt, q, ql, qc, nr, nrs, nrsws, nrswst,        &
               prr, precipitation_amount, qr,                                  &
               qrs, qrsws, qrswst, qs, qsws, qswst, rs, rss, rssws, rsswst,    &   !bK added rs, rss, rssws, rsswst
               s, sa, ssws, sswst,                                             &   
               saswsb, saswst, rif_wall, shf, ss, ts, tswst, u, u_m_l, u_m_n,  &
               u_m_r, u_m_s, us, usws, uswst, v, v_m_l, v_m_n, v_m_r, v_m_s,   &
               vpt, vsws, vswst, w, w_m_l, w_m_n, w_m_r, w_m_s, z0, z0h, z0q
        
    USE averaging

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, iran, humidity, passive_scalar, cloud_physics,   &
               cloud_droplets, microphysics_seifert, ocean, topography
                
    USE indices,                                                               &
        ONLY:  nxl, nxr, nys, nyn, nzb, nzt
        
    USE particle_attributes,                                                   &
        ONLY:  iran_part
        
    USE pegrid
    
    USE random_function_mod,                                                   &
        ONLY:  random_iv, random_iy

    USE random_generator_parallel,                                             &
        ONLY:  id_random_array, seq_random_array
        
    USE spectra_mod,                                                           &
        ONLY:  spectrum_x, spectrum_y

#if defined( __chem )         
    USE chemistry_model_mod,                                                   &
        ONLY: chem_last_actions                             
#endif

    IMPLICIT NONE

    CHARACTER (LEN=10) ::  binary_version   !< 


!
!-- Write control parameters and other variables for restart.
    IF ( myid == 0 )  CALL write_var_list
!    if(myid == 0) print*,'fm write_3d_binary, #21 '        !bK debug
!
!-- Write arrays.
    binary_version = '4.5'

    WRITE ( 14 )  binary_version

    WRITE ( 14 )  numprocs, myid, nxl, nxr, nys, nyn, nzb, nzt

!
!-- Attention: After changes to the following output commands the version number
!-- ---------  of the variable binary_version must be changed!
!--            Also, the list of arrays to be read in read_3d_binary must be
!--            adjusted accordingly.
if(myid == 0) print*,'fm 90@167 write bin is: ',air_chemistry

    WRITE ( 14 )  'e                   ';  WRITE ( 14 )  e
    IF ( ALLOCATED( e_av ) )  THEN
       WRITE ( 14 )  'e_av                ';  WRITE ( 14 )  e_av
    ENDIF
    WRITE ( 14 )  'iran                ';  WRITE ( 14 )  iran, iran_part
    WRITE ( 14 )  'kh                  ';  WRITE ( 14 )  kh
    WRITE ( 14 )  'km                  ';  WRITE ( 14 )  km
    IF ( ALLOCATED( lpt_av ) )  THEN
       WRITE ( 14 )  'lpt_av              ';  WRITE ( 14 )  lpt_av
    ENDIF
    IF ( ALLOCATED( lwp_av ) )  THEN
       WRITE ( 14 )  'lwp_av              ';  WRITE ( 14 )  lwp_av
    ENDIF
    WRITE ( 14 )  'ol                  ';  WRITE ( 14 )  ol
    WRITE ( 14 )  'p                   ';  WRITE ( 14 )  p
    IF ( ALLOCATED( p_av ) )  THEN
       WRITE ( 14 )  'p_av                ';  WRITE ( 14 )  p_av
    ENDIF
    IF ( ALLOCATED( pc_av ) )  THEN
       WRITE ( 14 )  'pc_av               ';  WRITE ( 14 )  pc_av
    ENDIF
    IF ( ALLOCATED( pr_av ) )  THEN
       WRITE ( 14 )  'pr_av               ';  WRITE ( 14 )  pr_av
    ENDIF
    IF ( ALLOCATED( prr ) )  THEN
       WRITE ( 14 )  'prr                 ';  WRITE ( 14 )  prr
    ENDIF
    IF ( ALLOCATED( prr_av ) )  THEN
       WRITE ( 14 )  'prr_av              ';  WRITE ( 14 )  prr_av
    ENDIF
    IF ( ALLOCATED( precipitation_amount ) )  THEN
       WRITE ( 14 )  'precipitation_amount';  WRITE ( 14 )  precipitation_amount
    ENDIF
    IF ( ALLOCATED( precipitation_rate_av ) )  THEN
       WRITE ( 14 )  'precipitation_rate_a';  WRITE ( 14 )                     &
                                                           precipitation_rate_av
    ENDIF
    WRITE ( 14 )  'pt                  ';  WRITE ( 14 )  pt
       if(myid==0)  print*,'fm write_3d_binary writing  pt,  #21a '
    IF ( ALLOCATED( pt_av ) )  THEN
       WRITE ( 14 )  'pt_av               ';  WRITE ( 14 )  pt_av
       if(myid==0) print*,'from write_3d_binary writing pt_av #21b '
    ENDIF
    IF ( humidity )  THEN
       WRITE ( 14 )  'q                   ';  WRITE ( 14 )  q 
       IF ( ALLOCATED( q_av ) )  THEN
          WRITE ( 14 )  'q_av                ';  WRITE ( 14 )  q_av
       ENDIF
       IF ( cloud_physics  .OR.  cloud_droplets )  THEN
          WRITE ( 14 )  'ql                  ';  WRITE ( 14 )  ql
          IF ( ALLOCATED( ql_av ) )  THEN
             WRITE ( 14 )  'ql_av               ';  WRITE ( 14 )  ql_av
          ENDIF
          IF ( .NOT. cloud_droplets )  THEN
             WRITE ( 14 )  'qc                  ';  WRITE ( 14 )  qc
             IF ( ALLOCATED( qc_av ) )  THEN
                WRITE ( 14 )  'qc_av               ';  WRITE ( 14 )  qc_av
             ENDIF
             IF ( microphysics_seifert )  THEN
                WRITE ( 14 )  'nr                  ';  WRITE ( 14 )  nr
                IF ( ALLOCATED( nr_av ) )  THEN
                   WRITE ( 14 )  'nr_av               ';  WRITE ( 14 )  nr_av
                ENDIF
                WRITE ( 14 )  'nrs                 ';  WRITE ( 14 )  nrs
                WRITE ( 14 )  'nrsws               ';  WRITE ( 14 )  nrsws
                WRITE ( 14 )  'nrswst              ';  WRITE ( 14 )  nrswst
                WRITE ( 14 )  'qr                  ';  WRITE ( 14 )  qr
                IF ( ALLOCATED( qr_av ) )  THEN
                   WRITE ( 14 )  'qr_av               ';  WRITE ( 14 )  qr_av
                ENDIF
                WRITE ( 14 )  'qrs                 ';  WRITE ( 14 )  qrs
                WRITE ( 14 )  'qrsws               ';  WRITE ( 14 )  qrsws
                WRITE ( 14 )  'qrswst              ';  WRITE ( 14 )  qrswst
             ENDIF
          ENDIF
       ENDIF
       WRITE ( 14 )  'qs                  ';  WRITE ( 14 )  qs
       WRITE ( 14 )  'qsws                ';  WRITE ( 14 )  qsws
       IF ( ALLOCATED( qsws_av ) )  THEN
          WRITE ( 14 )  'qsws_av             ';  WRITE ( 14 )  qsws_av
       ENDIF
       WRITE ( 14 )  'qswst               ';  WRITE ( 14 ) qswst
    ENDIF

    IF ( passive_scalar )  THEN
       WRITE ( 14 )  's                   ';  WRITE ( 14 )  s 
       IF ( ALLOCATED( s_av ) )  THEN
          WRITE ( 14 )  's_av                ';  WRITE ( 14 )  s_av
       ENDIF
       WRITE ( 14 )  'ss                  ';  WRITE ( 14 )  ss
       WRITE ( 14 )  'ssws                ';  WRITE ( 14 )  ssws
       IF ( ALLOCATED( ssws_av ) )  THEN
          WRITE ( 14 )  'ssws_av             ';  WRITE ( 14 )  ssws_av
       ENDIF
       WRITE ( 14 )  'sswst               ';  WRITE ( 14 ) sswst
    ENDIF    

!#if defined( __chem )
!
!    if(myid == 0) print*,'fm 90@267 write bin is: ',air_chemistry
!
!    IF ( air_chemistry )  THEN                                                  !bK added this if block to generate
!        CALL chem_last_actions
!    END IF


!       WRITE ( 14 )  'rs                   ';  WRITE ( 14 )  rs             !   chem_spcs profiles in the
!        if(myid == 0) print*,'rs_av is there or not?    ', allocated(rs_av)
!       IF ( ALLOCATED( rs_av ) )  THEN                                      !   recycle run.
!            if(myid == 0) print*,'i am gonna write rs_av'
!          WRITE ( 14 )  'rs_av                ';  WRITE ( 14 )  rs_av
!       ENDIF
!       WRITE ( 14 )  'rss                  ';  WRITE ( 14 )  rss
!       WRITE ( 14 )  'rssws                ';  WRITE ( 14 )  rssws
!       IF ( ALLOCATED( rssws_av ) )  THEN
!             write(9,*) 'i am gonna write rssws'          
!          WRITE ( 14 )  'rssws_av             ';  WRITE ( 14 )  rssws_av
!       ENDIF
!       WRITE ( 14 )  'rsswst               ';  WRITE ( 14 ) rsswst
!    ENDIF    
!
!#endif

    IF ( ocean )  THEN
       IF ( ALLOCATED( rho_ocean_av ) )  THEN
          WRITE ( 14 )  'rho_ocean_av              ';  WRITE ( 14 )  rho_ocean_av
       ENDIF
       WRITE ( 14 )  'sa                  ';  WRITE ( 14 )  sa
       IF ( ALLOCATED( sa_av ) )  THEN
          WRITE ( 14 )  'sa_av               ';  WRITE ( 14 )  sa_av
       ENDIF
       WRITE ( 14 )  'saswsb              ';  WRITE ( 14 )  saswsb
       WRITE ( 14 )  'saswst              ';  WRITE ( 14 )  saswst
    ENDIF
    IF ( ALLOCATED( ql_c_av ) )  THEN
       WRITE ( 14 )  'ql_c_av             ';  WRITE ( 14 )  ql_c_av
    ENDIF
    IF ( ALLOCATED( ql_v_av ) )  THEN
       WRITE ( 14 )  'ql_v_av             ';  WRITE ( 14 )  ql_v_av
    ENDIF
    IF ( ALLOCATED( ql_vp_av ) )  THEN
       WRITE ( 14 )  'ql_vp_av            ';  WRITE ( 14 )  ql_vp_av
    ENDIF
    IF ( ALLOCATED( qv_av ) )  THEN
       WRITE ( 14 )  'qv_av               ';  WRITE ( 14 )  qv_av
    ENDIF
    WRITE ( 14 )  'random_iv           ';  WRITE ( 14 )  random_iv
                                           WRITE ( 14 )  random_iy
    IF ( ALLOCATED( seq_random_array ) )  THEN
    WRITE ( 14 )  'seq_random_array    ';  WRITE ( 14 )  id_random_array
                                           WRITE ( 14 )  seq_random_array
    ENDIF
    IF ( topography /= 'flat' )  THEN
       WRITE ( 14 )  'rif_wall            ';  WRITE ( 14 )  rif_wall
    ENDIF
    IF ( ALLOCATED( s_av ) )  THEN
       WRITE ( 14 )  's_av                ';  WRITE ( 14 )  s_av
    ENDIF
    WRITE ( 14 )  'shf                 ';  WRITE ( 14 )  shf
    IF ( ALLOCATED( shf_av ) )  THEN
       WRITE ( 14 )  'shf_av              ';  WRITE ( 14 )  shf_av
    ENDIF
    IF ( ALLOCATED( spectrum_x ) )  THEN
       WRITE ( 14 )  'spectrum_x          ';  WRITE ( 14 )  spectrum_x
       WRITE ( 14 )  'spectrum_y          ';  WRITE ( 14 )  spectrum_y
    ENDIF
    WRITE ( 14 )  'ts                  ';  WRITE ( 14 )  ts
    IF ( ALLOCATED( ts_av ) )  THEN
       WRITE ( 14 )  'ts_av               ';  WRITE ( 14 )  ts_av
    ENDIF
    WRITE ( 14 )  'tswst               ';  WRITE ( 14 )  tswst
    WRITE ( 14 )  'u                   ';  WRITE ( 14 )  u
    IF ( ALLOCATED( u_av ) )  THEN
       WRITE ( 14 )  'u_av                ';  WRITE ( 14 )  u_av
    ENDIF
    IF ( ALLOCATED( u_m_l ) )  THEN
       WRITE ( 14 )  'u_m_l               ';  WRITE ( 14 )  u_m_l
    ENDIF
    IF ( ALLOCATED( u_m_n ) )  THEN
       WRITE ( 14 )  'u_m_n               ';  WRITE ( 14 )  u_m_n
    ENDIF
    IF ( ALLOCATED( u_m_r ) )  THEN
       WRITE ( 14 )  'u_m_r               ';  WRITE ( 14 )  u_m_r
    ENDIF
    IF ( ALLOCATED( u_m_s ) )  THEN
       WRITE ( 14 )  'u_m_s               ';  WRITE ( 14 )  u_m_s
    ENDIF
    WRITE ( 14 )  'us                  ';  WRITE ( 14 )  us
    WRITE ( 14 )  'usws                ';  WRITE ( 14 )  usws
    WRITE ( 14 )  'uswst               ';  WRITE ( 14 )  uswst
    IF ( ALLOCATED( us_av ) )  THEN
       WRITE ( 14 )  'us_av               ';  WRITE ( 14 )  us_av
    ENDIF
    WRITE ( 14 )  'v                   ';  WRITE ( 14 )  v
    IF ( ALLOCATED( v_av ) )  THEN
       WRITE ( 14 )  'v_av                ';  WRITE ( 14 )  v_av
    ENDIF
    IF ( ALLOCATED( v_m_l ) )  THEN
       WRITE ( 14 )  'v_m_l               ';  WRITE ( 14 )  v_m_l
    ENDIF
    IF ( ALLOCATED( v_m_n ) )  THEN
       WRITE ( 14 )  'v_m_n               ';  WRITE ( 14 )  v_m_n
    ENDIF
    IF ( ALLOCATED( v_m_r ) )  THEN
       WRITE ( 14 )  'v_m_r               ';  WRITE ( 14 )  v_m_r
    ENDIF
    IF ( ALLOCATED( v_m_s ) )  THEN
       WRITE ( 14 )  'v_m_s               ';  WRITE ( 14 )  v_m_s
    ENDIF
    IF ( humidity )  THEN
       WRITE ( 14 )  'vpt                 ';  WRITE ( 14 )  vpt
       IF ( ALLOCATED( vpt_av ) )  THEN
          WRITE ( 14 )  'vpt_av              ';  WRITE ( 14 )  vpt_av
       ENDIF
    ENDIF
    WRITE ( 14 )  'vsws                ';  WRITE ( 14 )  vsws
    WRITE ( 14 )  'vswst               ';  WRITE ( 14 )  vswst
    WRITE ( 14 )  'w                   ';  WRITE ( 14 )  w
    IF ( ALLOCATED( w_av ) )  THEN
       WRITE ( 14 )  'w_av                ';  WRITE ( 14 )  w_av
    ENDIF
    IF ( ALLOCATED( w_m_l ) )  THEN
       WRITE ( 14 )  'w_m_l               ';  WRITE ( 14 )  w_m_l
    ENDIF
    IF ( ALLOCATED( w_m_n ) )  THEN
       WRITE ( 14 )  'w_m_n               ';  WRITE ( 14 )  w_m_n
    ENDIF
    IF ( ALLOCATED( w_m_r ) )  THEN
       WRITE ( 14 )  'w_m_r               ';  WRITE ( 14 )  w_m_r
    ENDIF
    IF ( ALLOCATED( w_m_s ) )  THEN
       WRITE ( 14 )  'w_m_s               ';  WRITE ( 14 )  w_m_s
    ENDIF
    WRITE ( 14 )  'z0                  ';  WRITE ( 14 )  z0
    IF ( ALLOCATED( z0_av ) )  THEN
       WRITE ( 14 )  'z0_av               ';  WRITE ( 14 )  z0_av
    ENDIF
    WRITE ( 14 )  'z0h                 ';  WRITE ( 14 )  z0h
    IF ( ALLOCATED( z0h_av ) )  THEN
       WRITE ( 14 )  'z0h_av              ';  WRITE ( 14 )  z0h_av
    ENDIF
    WRITE ( 14 )  'z0q                 ';  WRITE ( 14 )  z0q
    IF ( ALLOCATED( z0q_av ) )  THEN
       WRITE ( 14 )  'z0q_av              ';  WRITE ( 14 )  z0q_av
    ENDIF

!
!-- Write end label. Unit 14 is closed in the main program.
    WRITE ( 14 )  '*** end ***         '

 END SUBROUTINE write_3d_binary
 
