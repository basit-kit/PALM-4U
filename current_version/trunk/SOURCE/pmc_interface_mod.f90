MODULE pmc_interface

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
! $Id: pmc_interface_mod.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2019 2016-09-30 13:40:09Z hellstea
! Bugfixes mainly in determining the anterpolation index bounds. These errors 
! were detected when first time tested using 3:1 grid-spacing.
! 
! 2003 2016-08-24 10:22:32Z suehring
! Humidity and passive scalar also separated in nesting mode 
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1938 2016-06-13 15:26:05Z hellstea
! Minor clean-up.
! 
! 1901 2016-05-04 15:39:38Z raasch
! Initial version of purely vertical nesting introduced.
! Code clean up. The words server/client changed to parent/child.
!
! 1900 2016-05-04 15:27:53Z raasch
! unused variables removed
!
! 1894 2016-04-27 09:01:48Z raasch
! bugfix: pt interpolations are omitted in case that the temperature equation is
! switched off
!
! 1892 2016-04-26 13:49:47Z raasch
! bugfix: pt is not set as a data array in case that the temperature equation is
! switched off with neutral = .TRUE.
!
! 1882 2016-04-20 15:24:46Z hellstea
! The factor ijfc for nfc used in anterpolation is redefined as 2-D array 
! and precomputed in pmci_init_anterp_tophat. 
! 
! 1878 2016-04-19 12:30:36Z hellstea
! Synchronization rewritten, logc-array index order changed for cache
! optimization
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1801 2016-04-05 13:12:47Z raasch
! bugfix for r1797: zero setting of temperature within topography does not work
! and has been disabled
!
! 1797 2016-03-21 16:50:28Z raasch
! introduction of different datatransfer modes,
! further formatting cleanup, parameter checks added (including mismatches
! between root and nest model settings),
! +routine pmci_check_setting_mismatches
! comm_world_nesting introduced
!
! 1791 2016-03-11 10:41:25Z raasch
! routine pmci_update_new removed,
! pmc_get_local_model_info renamed pmc_get_model_info, some keywords also
! renamed,
! filling up redundant ghost points introduced,
! some index bound variables renamed,
! further formatting cleanup
!
! 1783 2016-03-06 18:36:17Z raasch
! calculation of nest top area simplified,
! interpolation and anterpolation moved to seperate wrapper subroutines
!
! 1781 2016-03-03 15:12:23Z raasch
! _p arrays are set zero within buildings too, t.._m arrays and respective
! settings within buildings completely removed
!
! 1779 2016-03-03 08:01:28Z raasch
! only the total number of PEs is given for the domains, npe_x and npe_y
! replaced by npe_total, two unused elements removed from array
! define_coarse_grid_real,
! array management changed from linked list to sequential loop
!
! 1766 2016-02-29 08:37:15Z raasch
! modifications to allow for using PALM's pointer version,
! +new routine pmci_set_swaplevel
!
! 1764 2016-02-28 12:45:19Z raasch
! +cpl_parent_id,
! cpp-statements for nesting replaced by __parallel statements,
! errors output with message-subroutine,
! index bugfixes in pmci_interp_tril_all,
! some adjustments to PALM style
!
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by A. Hellsten
!
! Description:
! ------------
! Domain nesting interface routines. The low-level inter-domain communication   
! is conducted by the PMC-library routines.
!-------------------------------------------------------------------------------!

#if defined( __nopointer )
    USE arrays_3d,                                                              &
        ONLY:  dzu, dzw, e, e_p, pt, pt_p, q, q_p, u, u_p, v, v_p, w, w_p, zu,  &
               zw, z0
#else
   USE arrays_3d,                                                               &
        ONLY:  dzu, dzw, e, e_p, e_1, e_2, pt, pt_p, pt_1, pt_2, q, q_p, q_1,   &
               q_2, s, s_2, u, u_p, u_1, u_2, v, v_p, v_1, v_2, w, w_p, w_1,    &
               w_2, zu, zw, z0
#endif

    USE control_parameters,                                                     &
        ONLY:  coupling_char, dt_3d, dz, humidity, message_string,              &
               nest_bound_l, nest_bound_r, nest_bound_s, nest_bound_n,          &
               nest_domain, neutral, passive_scalar, simulated_time,            &
               topography, volume_flow

    USE cpulog,                                                                 &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                         &
        ONLY:  dx, dy

    USE indices,                                                                &
        ONLY:  nbgp, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nyn, nyng, nys, nysg,  &
               nysv, nz, nzb, nzb_s_inner, nzb_u_inner, nzb_u_outer,            &
               nzb_v_inner, nzb_v_outer, nzb_w_inner, nzb_w_outer, nzt

    USE kinds

#if defined( __parallel )
#if defined( __mpifh )
    INCLUDE "mpif.h"
#else
    USE MPI
#endif

    USE pegrid,                                                                 &
        ONLY:  collective_wait, comm1dx, comm1dy, comm2d, myid, myidx, myidy,   &
               numprocs

    USE pmc_child,                                                              &
        ONLY:  pmc_childinit, pmc_c_clear_next_array_list,                      &
               pmc_c_getnextarray, pmc_c_get_2d_index_list, pmc_c_getbuffer,    &
               pmc_c_putbuffer, pmc_c_setind_and_allocmem,                      &
               pmc_c_set_dataarray, pmc_set_dataarray_name

    USE pmc_general,                                                            &
        ONLY:  da_namelen

    USE pmc_handle_communicator,                                                &
        ONLY:  pmc_get_model_info, pmc_init_model, pmc_is_rootmodel,            &
               pmc_no_namelist_found, pmc_parent_for_child

    USE pmc_mpi_wrapper,                                                        &
        ONLY:  pmc_bcast, pmc_recv_from_child, pmc_recv_from_parent,            &
               pmc_send_to_child, pmc_send_to_parent

    USE pmc_parent,                                                             &
        ONLY:  pmc_parentinit, pmc_s_clear_next_array_list, pmc_s_fillbuffer,   &
               pmc_s_getdata_from_buffer, pmc_s_getnextarray,                   &
               pmc_s_setind_and_allocmem, pmc_s_set_active_data_array,          &
               pmc_s_set_dataarray, pmc_s_set_2d_index_list

#endif

    IMPLICIT NONE

    PRIVATE

!
!-- Constants
    INTEGER(iwp), PARAMETER ::  child_to_parent = 2   !:
    INTEGER(iwp), PARAMETER ::  parent_to_child = 1   !:

!
!-- Coupler setup
    INTEGER(iwp), SAVE      ::  comm_world_nesting     !:
    INTEGER(iwp), SAVE      ::  cpl_id  = 1            !:
    CHARACTER(LEN=32), SAVE ::  cpl_name               !:
    INTEGER(iwp), SAVE      ::  cpl_npe_total          !:
    INTEGER(iwp), SAVE      ::  cpl_parent_id          !:

!
!-- Control parameters, will be made input parameters later
    CHARACTER(LEN=7), SAVE ::  nesting_datatransfer_mode = 'mixed'  !: steering
                                                         !: parameter for data-
                                                         !: transfer mode
    CHARACTER(LEN=8), SAVE ::  nesting_mode = 'two-way'  !: steering parameter
                                                         !: for 1- or 2-way nesting

    LOGICAL, SAVE ::  nested_run = .FALSE.  !: general switch

    REAL(wp), SAVE ::  anterp_relax_length_l = -1.0_wp   !:
    REAL(wp), SAVE ::  anterp_relax_length_r = -1.0_wp   !:
    REAL(wp), SAVE ::  anterp_relax_length_s = -1.0_wp   !:
    REAL(wp), SAVE ::  anterp_relax_length_n = -1.0_wp   !:
    REAL(wp), SAVE ::  anterp_relax_length_t = -1.0_wp   !:

!
!-- Geometry
    REAL(wp), SAVE                            ::  area_t               !:
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE ::  coord_x              !:
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE ::  coord_y              !:
    REAL(wp), SAVE                            ::  lower_left_coord_x   !:
    REAL(wp), SAVE                            ::  lower_left_coord_y   !:

!
!-- Child coarse data arrays
    INTEGER(iwp), DIMENSION(5)                  ::  coarse_bound   !:

    REAL(wp), SAVE                              ::  xexl           !:
    REAL(wp), SAVE                              ::  xexr           !:
    REAL(wp), SAVE                              ::  yexs           !:
    REAL(wp), SAVE                              ::  yexn           !:
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE ::  tkefactor_l    !:
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE ::  tkefactor_n    !:
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE ::  tkefactor_r    !:
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE ::  tkefactor_s    !:
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE ::  tkefactor_t    !:

    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ec   !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ptc  !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uc   !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vc   !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wc   !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  qc   !:
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sc   !:

!
!-- Child interpolation coefficients and child-array indices to be 
!-- precomputed and stored.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ico    !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  icu    !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jco    !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jcv    !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kco    !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kcw    !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1xo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2xo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1xu   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2xu   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1yo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2yo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1yv   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2yv   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1zo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2zo   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r1zw   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:)     ::  r2zw   !:

!
!-- Child index arrays and log-ratio arrays for the log-law near-wall
!-- corrections. These are not truly 3-D arrays but multiple 2-D arrays.
    INTEGER(iwp), SAVE :: ncorr  !: 4th dimension of the log_ratio-arrays
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_u_l   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_u_n   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_u_r   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_u_s   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_v_l   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_v_n   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_v_r   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_v_s   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_w_l   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_w_n   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_w_r   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  logc_w_s   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_u_l   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_u_n   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_u_r   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_u_s   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_v_l   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_v_n   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_v_r   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_v_s   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_w_l   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_w_n   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_w_r   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:)   ::  logc_ratio_w_s   !:

!
!-- Upper bounds for k in anterpolation.
    INTEGER(iwp), SAVE ::  kctu   !:
    INTEGER(iwp), SAVE ::  kctw   !:

!
!-- Upper bound for k in log-law correction in interpolation.
    INTEGER(iwp), SAVE ::  nzt_topo_nestbc_l   !:
    INTEGER(iwp), SAVE ::  nzt_topo_nestbc_n   !:
    INTEGER(iwp), SAVE ::  nzt_topo_nestbc_r   !:
    INTEGER(iwp), SAVE ::  nzt_topo_nestbc_s   !:

!
!-- Number of ghost nodes in coarse-grid arrays for i and j in anterpolation.
    INTEGER(iwp), SAVE ::  nhll   !:
    INTEGER(iwp), SAVE ::  nhlr   !:
    INTEGER(iwp), SAVE ::  nhls   !:
    INTEGER(iwp), SAVE ::  nhln   !:

!
!-- Spatial under-relaxation coefficients for anterpolation.
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:) ::  frax   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:) ::  fray   !:
    REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:) ::  fraz   !:

!
!-- Child-array indices to be precomputed and stored for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflu   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuu   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflo   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuo   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflv   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuv   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflo   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuo   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflw   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuw   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflo   !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuo   !:

!
!-- Number of fine-grid nodes inside coarse-grid ij-faces
!-- to be precomputed for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:) ::  ijfc_u        !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:) ::  ijfc_v        !:
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:) ::  ijfc_s        !:
    
    INTEGER(iwp), DIMENSION(3)          ::  define_coarse_grid_int    !:
    REAL(wp), DIMENSION(7)              ::  define_coarse_grid_real   !:

    TYPE coarsegrid_def
       INTEGER(iwp)                        ::  nx
       INTEGER(iwp)                        ::  ny
       INTEGER(iwp)                        ::  nz
       REAL(wp)                            ::  dx
       REAL(wp)                            ::  dy
       REAL(wp)                            ::  dz
       REAL(wp)                            ::  lower_left_coord_x
       REAL(wp)                            ::  lower_left_coord_y
       REAL(wp)                            ::  xend
       REAL(wp)                            ::  yend
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_x
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_y
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzu       
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzw       
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu        
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw        
    END TYPE coarsegrid_def
                                          
    TYPE(coarsegrid_def), SAVE ::  cg   !:


    INTERFACE pmci_check_setting_mismatches
       MODULE PROCEDURE pmci_check_setting_mismatches
    END INTERFACE

    INTERFACE pmci_child_initialize
       MODULE PROCEDURE pmci_child_initialize
    END INTERFACE

    INTERFACE pmci_synchronize
       MODULE PROCEDURE pmci_synchronize
    END INTERFACE

    INTERFACE pmci_datatrans
       MODULE PROCEDURE pmci_datatrans
    END INTERFACE pmci_datatrans

    INTERFACE pmci_ensure_nest_mass_conservation
       MODULE PROCEDURE pmci_ensure_nest_mass_conservation
    END INTERFACE

    INTERFACE pmci_init
       MODULE PROCEDURE pmci_init
    END INTERFACE

    INTERFACE pmci_modelconfiguration
       MODULE PROCEDURE pmci_modelconfiguration
    END INTERFACE

    INTERFACE pmci_parent_initialize
       MODULE PROCEDURE pmci_parent_initialize
    END INTERFACE

    INTERFACE pmci_set_swaplevel
       MODULE PROCEDURE pmci_set_swaplevel
    END INTERFACE pmci_set_swaplevel

    PUBLIC anterp_relax_length_l, anterp_relax_length_r,                        &
           anterp_relax_length_s, anterp_relax_length_n,                        &
           anterp_relax_length_t, child_to_parent, comm_world_nesting,          &
           cpl_id, nested_run, nesting_datatransfer_mode, nesting_mode,         &
           parent_to_child
    PUBLIC pmci_child_initialize
    PUBLIC pmci_datatrans
    PUBLIC pmci_ensure_nest_mass_conservation
    PUBLIC pmci_init
    PUBLIC pmci_modelconfiguration
    PUBLIC pmci_parent_initialize
    PUBLIC pmci_synchronize
    PUBLIC pmci_set_swaplevel


 CONTAINS


 SUBROUTINE pmci_init( world_comm )

    USE control_parameters,                                                     &
        ONLY:  message_string

    IMPLICIT NONE

    INTEGER, INTENT(OUT) ::  world_comm   !:

#if defined( __parallel )

    INTEGER(iwp)         ::  ierr         !:
    INTEGER(iwp)         ::  istat        !:
    INTEGER(iwp)         ::  pmc_status   !:


    CALL pmc_init_model( world_comm, nesting_datatransfer_mode, nesting_mode,   &
                         pmc_status )

    IF ( pmc_status == pmc_no_namelist_found )  THEN
!
!--    This is not a nested run
       world_comm = MPI_COMM_WORLD
       cpl_id     = 1
       cpl_name   = ""

       RETURN

    ENDIF

!
!-- Check steering parameter values
    IF ( TRIM( nesting_mode ) /= 'one-way'  .AND.                               &
         TRIM( nesting_mode ) /= 'two-way'  .AND.                               &
         TRIM( nesting_mode ) /= 'vertical' )                                   &                 
    THEN
       message_string = 'illegal nesting mode: ' // TRIM( nesting_mode )
       CALL message( 'pmci_init', 'PA0417', 3, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( nesting_datatransfer_mode ) /= 'cascade'  .AND.                  &
         TRIM( nesting_datatransfer_mode ) /= 'mixed'    .AND.                  &
         TRIM( nesting_datatransfer_mode ) /= 'overlap' )                       &
    THEN
       message_string = 'illegal nesting datatransfer mode: '                   &
                        // TRIM( nesting_datatransfer_mode )
       CALL message( 'pmci_init', 'PA0418', 3, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the general steering switch which tells PALM that its a nested run
    nested_run = .TRUE.

!
!-- Get some variables required by the pmc-interface (and in some cases in the
!-- PALM code out of the pmci) out of the pmc-core
    CALL pmc_get_model_info( comm_world_nesting = comm_world_nesting,           &
                             cpl_id = cpl_id, cpl_parent_id = cpl_parent_id,    &
                             cpl_name = cpl_name, npe_total = cpl_npe_total,    &
                             lower_left_x = lower_left_coord_x,                 &
                             lower_left_y = lower_left_coord_y )
!
!-- Set the steering switch which tells the models that they are nested (of
!-- course the root domain (cpl_id = 1) is not nested)
    IF ( cpl_id >= 2 )  THEN
       nest_domain = .TRUE.
       WRITE( coupling_char, '(A1,I2.2)') '_', cpl_id
    ENDIF

!
!-- Message that communicators for nesting are initialized.
!-- Attention: myid has been set at the end of pmc_init_model in order to
!-- guarantee that only PE0 of the root domain does the output.
    CALL location_message( 'finished', .TRUE. )
!
!-- Reset myid to its default value
    myid = 0
#else
!
!-- Nesting cannot be used in serial mode. cpl_id is set to root domain (1)
!-- because no location messages would be generated otherwise.
!-- world_comm is given a dummy value to avoid compiler warnings (INTENT(OUT)
!-- should get an explicit value)
    cpl_id     = 1
    nested_run = .FALSE.
    world_comm = 1
#endif

 END SUBROUTINE pmci_init



 SUBROUTINE pmci_modelconfiguration

    IMPLICIT NONE

    CALL location_message( 'setup the nested model configuration', .FALSE. )
!
!-- Compute absolute coordinates for all models
    CALL pmci_setup_coordinates
!
!-- Initialize the child (must be called before pmc_setup_parent)
    CALL pmci_setup_child
!
!-- Initialize PMC parent
    CALL pmci_setup_parent
!
!-- Check for mismatches between settings of master and child variables
!-- (e.g., all children have to follow the end_time settings of the root master)
    CALL pmci_check_setting_mismatches

    CALL location_message( 'finished', .TRUE. )

 END SUBROUTINE pmci_modelconfiguration



 SUBROUTINE pmci_setup_parent

#if defined( __parallel )
    IMPLICIT NONE

    CHARACTER(LEN=32) ::  myname

    INTEGER(iwp) ::  child_id         !:
    INTEGER(iwp) ::  ierr             !:
    INTEGER(iwp) ::  i                !:
    INTEGER(iwp) ::  j                !:
    INTEGER(iwp) ::  k                !:
    INTEGER(iwp) ::  m                !:
    INTEGER(iwp) ::  mm               !:
    INTEGER(iwp) ::  nest_overlap     !:
    INTEGER(iwp) ::  nomatch          !:
    INTEGER(iwp) ::  nx_cl            !:
    INTEGER(iwp) ::  ny_cl            !:
    INTEGER(iwp) ::  nz_cl            !:

    INTEGER(iwp), DIMENSION(5) ::  val    !:


    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ch_xl   !:
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ch_xr   !:    
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ch_ys   !:
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ch_yn   !:
    REAL(wp) ::  dx_cl            !:
    REAL(wp) ::  dy_cl            !:
    REAL(wp) ::  left_limit       !:
    REAL(wp) ::  north_limit      !:
    REAL(wp) ::  right_limit      !:
    REAL(wp) ::  south_limit      !:
    REAL(wp) ::  xez              !:
    REAL(wp) ::  yez              !:

    REAL(wp), DIMENSION(1) ::  fval             !:

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  cl_coord_x   !:
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  cl_coord_y   !:
    

!
!   Initialize the pmc parent
    CALL pmc_parentinit

!
!-- Corners of all children of the present parent
    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 ) .AND. myid == 0 )  THEN 
       ALLOCATE( ch_xl(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( ch_xr(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( ch_ys(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( ch_yn(1:SIZE( pmc_parent_for_child ) - 1) )
    ENDIF

!
!-- Get coordinates from all children
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1

       child_id = pmc_parent_for_child(m)
       IF ( myid == 0 )  THEN        

          CALL pmc_recv_from_child( child_id, val,  size(val),  0, 123, ierr )
          CALL pmc_recv_from_child( child_id, fval, size(fval), 0, 124, ierr )
         
          nx_cl = val(1)
          ny_cl = val(2)
          dx_cl = val(4)
          dy_cl = val(5)

          nz_cl = nz

!
!--       Find the highest nest level in the coarse grid for the reduced z
!--       transfer
          DO  k = 1, nz                 
             IF ( zw(k) > fval(1) )  THEN
                nz_cl = k
                EXIT
             ENDIF
          ENDDO

!    
!--       Get absolute coordinates from the child
          ALLOCATE( cl_coord_x(-nbgp:nx_cl+nbgp) )
          ALLOCATE( cl_coord_y(-nbgp:ny_cl+nbgp) )
         
          CALL pmc_recv_from_child( child_id, cl_coord_x, SIZE( cl_coord_x ),   &
                                     0, 11, ierr )
          CALL pmc_recv_from_child( child_id, cl_coord_y, SIZE( cl_coord_y ),   &
                                     0, 12, ierr )
!          WRITE ( 0, * )  'receive from pmc child ', child_id, nx_cl, ny_cl
         
          define_coarse_grid_real(1) = lower_left_coord_x
          define_coarse_grid_real(2) = lower_left_coord_y
          define_coarse_grid_real(3) = dx
          define_coarse_grid_real(4) = dy
          define_coarse_grid_real(5) = lower_left_coord_x + ( nx + 1 ) * dx
          define_coarse_grid_real(6) = lower_left_coord_y + ( ny + 1 ) * dy
          define_coarse_grid_real(7) = dz

          define_coarse_grid_int(1)  = nx
          define_coarse_grid_int(2)  = ny
          define_coarse_grid_int(3)  = nz_cl

!
!--       Check that the child domain matches parent domain. 
          nomatch = 0
          IF ( nesting_mode == 'vertical' )  THEN
             right_limit = define_coarse_grid_real(5)
             north_limit = define_coarse_grid_real(6)
             IF ( ( cl_coord_x(nx_cl+1) /= right_limit ) .OR.                   &
                  ( cl_coord_y(ny_cl+1) /= north_limit ) )  THEN
                nomatch = 1
             ENDIF
          ELSE
          
!
!--       Check that the children domain is completely inside the parent domain. 
             xez = ( nbgp + 1 ) * dx 
             yez = ( nbgp + 1 ) * dy 
             left_limit  = lower_left_coord_x + xez
             right_limit = define_coarse_grid_real(5) - xez
             south_limit = lower_left_coord_y + yez
             north_limit = define_coarse_grid_real(6) - yez
             IF ( ( cl_coord_x(0) < left_limit )        .OR.                    &
                  ( cl_coord_x(nx_cl+1) > right_limit ) .OR.                    &
                  ( cl_coord_y(0) < south_limit )       .OR.                    &
                  ( cl_coord_y(ny_cl+1) > north_limit ) )  THEN
                nomatch = 1
             ENDIF
          ENDIF

!
!--       Check that parallel nest domains, if any, do not overlap.
          nest_overlap = 0
          IF ( SIZE( pmc_parent_for_child ) - 1 > 0 )  THEN
             ch_xl(m) = cl_coord_x(-nbgp)
             ch_xr(m) = cl_coord_x(nx_cl+nbgp)
             ch_ys(m) = cl_coord_y(-nbgp)
             ch_yn(m) = cl_coord_y(ny_cl+nbgp)

             IF ( m > 1 )  THEN
                DO mm = 1, m-1
                   IF ( ( ch_xl(m) < ch_xr(mm) .OR.                             &
                          ch_xr(m) > ch_xl(mm) )  .AND.                         &
                        ( ch_ys(m) < ch_yn(mm) .OR.                             &
                          ch_yn(m) > ch_ys(mm) ) )  THEN                        
                      nest_overlap = 1
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

          DEALLOCATE( cl_coord_x )
          DEALLOCATE( cl_coord_y )

!
!--       Send coarse grid information to child
          CALL pmc_send_to_child( child_id, define_coarse_grid_real,            &
                                   SIZE( define_coarse_grid_real ), 0, 21,      &
                                   ierr )
          CALL pmc_send_to_child( child_id, define_coarse_grid_int,  3, 0,      &
                                   22, ierr )

!
!--       Send local grid to child
          CALL pmc_send_to_child( child_id, coord_x, nx+1+2*nbgp, 0, 24,        &
                                   ierr )
          CALL pmc_send_to_child( child_id, coord_y, ny+1+2*nbgp, 0, 25,        &
                                   ierr )

!
!--       Also send the dzu-, dzw-, zu- and zw-arrays here
          CALL pmc_send_to_child( child_id, dzu, nz_cl+1, 0, 26, ierr )
          CALL pmc_send_to_child( child_id, dzw, nz_cl+1, 0, 27, ierr )
          CALL pmc_send_to_child( child_id, zu,  nz_cl+2, 0, 28, ierr )
          CALL pmc_send_to_child( child_id, zw,  nz_cl+2, 0, 29, ierr )

       ENDIF

       CALL MPI_BCAST( nomatch, 1, MPI_INTEGER, 0, comm2d, ierr )
       IF ( nomatch /= 0 ) THEN
          WRITE ( message_string, * )  'nested child domain does ',             &
                                       'not fit into its parent domain'
          CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
       ENDIF
  
       CALL MPI_BCAST( nest_overlap, 1, MPI_INTEGER, 0, comm2d, ierr )
       IF ( nest_overlap /= 0 ) THEN
          WRITE ( message_string, * )  'nested parallel child domains overlap'
          CALL message( 'pmci_setup_parent', 'PA0426', 3, 2, 0, 6, 0 )
       ENDIF
      
       CALL MPI_BCAST( nz_cl, 1, MPI_INTEGER, 0, comm2d, ierr )

!
!--    TO_DO: Klaus: please give a comment what is done here
       CALL pmci_create_index_list

!
!--    Include couple arrays into parent content
!--    TO_DO: Klaus: please give a more meaningful comment
       CALL pmc_s_clear_next_array_list
       DO  WHILE ( pmc_s_getnextarray( child_id, myname ) )
          CALL pmci_set_array_pointer( myname, child_id = child_id,             &
                                       nz_cl = nz_cl )
       ENDDO
       CALL pmc_s_setind_and_allocmem( child_id )
    ENDDO

    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 ) .AND. myid == 0 )  THEN 
       DEALLOCATE( ch_xl )
       DEALLOCATE( ch_xr )
       DEALLOCATE( ch_ys )
       DEALLOCATE( ch_yn )
    ENDIF

 CONTAINS


   SUBROUTINE pmci_create_index_list

       IMPLICIT NONE

       INTEGER(iwp) ::  i                  !:
       INTEGER(iwp) ::  ic                 !:
       INTEGER(iwp) ::  ierr               !:
       INTEGER(iwp) ::  j                  !:
       INTEGER(iwp) ::  k                  !:
       INTEGER(iwp) ::  m                  !:
       INTEGER(iwp) ::  n                  !:
       INTEGER(iwp) ::  npx                !:
       INTEGER(iwp) ::  npy                !:
       INTEGER(iwp) ::  nrx                !:
       INTEGER(iwp) ::  nry                !:
       INTEGER(iwp) ::  px                 !:
       INTEGER(iwp) ::  py                 !:
       INTEGER(iwp) ::  parent_pe          !:

       INTEGER(iwp), DIMENSION(2) ::  scoord             !:
       INTEGER(iwp), DIMENSION(2) ::  size_of_array      !:

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  coarse_bound_all   !:
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  index_list         !:

       IF ( myid == 0 )  THEN
!--       TO_DO: Klaus: give more specific comment what size_of_array stands for
          CALL pmc_recv_from_child( child_id, size_of_array, 2, 0, 40, ierr )
          ALLOCATE( coarse_bound_all(size_of_array(1),size_of_array(2)) )
          CALL pmc_recv_from_child( child_id, coarse_bound_all,                 &
                                    SIZE( coarse_bound_all ), 0, 41, ierr )

!
!--       Compute size of index_list.
          ic = 0
          DO  k = 1, size_of_array(2)
             DO  j = coarse_bound_all(3,k), coarse_bound_all(4,k)
                DO  i = coarse_bound_all(1,k), coarse_bound_all(2,k)
                   ic = ic + 1
                ENDDO
             ENDDO
          ENDDO

          ALLOCATE( index_list(6,ic) )

          CALL MPI_COMM_SIZE( comm1dx, npx, ierr )
          CALL MPI_COMM_SIZE( comm1dy, npy, ierr )
!
!--       The +1 in index is because PALM starts with nx=0
          nrx = nxr - nxl + 1
          nry = nyn - nys + 1
          ic  = 0
!
!--       Loop over all children PEs
          DO  k = 1, size_of_array(2)
!
!--          Area along y required by actual child PE
             DO  j = coarse_bound_all(3,k), coarse_bound_all(4,k)
!
!--             Area along x required by actual child PE
                DO  i = coarse_bound_all(1,k), coarse_bound_all(2,k)

                   px = i / nrx
                   py = j / nry
                   scoord(1) = px
                   scoord(2) = py
                   CALL MPI_CART_RANK( comm2d, scoord, parent_pe, ierr )
                  
                   ic = ic + 1
!
!--                First index in parent array
                   index_list(1,ic) = i - ( px * nrx ) + 1 + nbgp
!
!--                Second index in parent array
                   index_list(2,ic) = j - ( py * nry ) + 1 + nbgp
!
!--                x index of child coarse grid
                   index_list(3,ic) = i - coarse_bound_all(1,k) + 1
!
!--                y index of child coarse grid
                   index_list(4,ic) = j - coarse_bound_all(3,k) + 1
!
!--                PE number of child
                   index_list(5,ic) = k - 1
!
!--                PE number of parent
                   index_list(6,ic) = parent_pe

                ENDDO
             ENDDO
          ENDDO
!
!--       TO_DO: Klaus: comment what is done here
          CALL pmc_s_set_2d_index_list( child_id, index_list(:,1:ic) )

       ELSE
!
!--       TO_DO: Klaus: comment why this dummy allocation is required
          ALLOCATE( index_list(6,1) )
          CALL pmc_s_set_2d_index_list( child_id, index_list )
       ENDIF

       DEALLOCATE(index_list)

     END SUBROUTINE pmci_create_index_list

#endif
 END SUBROUTINE pmci_setup_parent



 SUBROUTINE pmci_setup_child

#if defined( __parallel )
    IMPLICIT NONE

    CHARACTER(LEN=da_namelen) ::  myname     !:

    INTEGER(iwp) ::  i          !:
    INTEGER(iwp) ::  ierr       !:
    INTEGER(iwp) ::  icl        !:
    INTEGER(iwp) ::  icr        !:
    INTEGER(iwp) ::  j          !:
    INTEGER(iwp) ::  jcn        !:
    INTEGER(iwp) ::  jcs        !:

    INTEGER(iwp), DIMENSION(5) ::  val        !:
    
    REAL(wp) ::  xcs        !:
    REAL(wp) ::  xce        !:
    REAL(wp) ::  ycs        !:
    REAL(wp) ::  yce        !:

    REAL(wp), DIMENSION(1) ::  fval       !:
                                             
!
!-- TO_DO: describe what is happening in this if-clause
!-- Root model does not have a parent and is not a child
    IF ( .NOT. pmc_is_rootmodel() )  THEN

       CALL pmc_childinit
!
!--    Here AND ONLY HERE the arrays are defined, which actualy will be
!--    exchanged between child and parent.
!--    If a variable is removed, it only has to be removed from here.
!--    Please check, if the arrays are in the list of POSSIBLE exchange arrays
!--    in subroutines:
!--    pmci_set_array_pointer (for parent arrays)
!--    pmci_create_child_arrays (for child arrays)
       CALL pmc_set_dataarray_name( 'coarse', 'u'  ,'fine', 'u',  ierr )
       CALL pmc_set_dataarray_name( 'coarse', 'v'  ,'fine', 'v',  ierr )
       CALL pmc_set_dataarray_name( 'coarse', 'w'  ,'fine', 'w',  ierr )
       CALL pmc_set_dataarray_name( 'coarse', 'e'  ,'fine', 'e',  ierr )
       IF ( .NOT. neutral )  THEN
          CALL pmc_set_dataarray_name( 'coarse', 'pt' ,'fine', 'pt', ierr )
       ENDIF
       IF ( humidity )  THEN
          CALL pmc_set_dataarray_name( 'coarse', 'q'  ,'fine', 'q',  ierr )
       ENDIF
       IF ( passive_scalar )  THEN
          CALL pmc_set_dataarray_name( 'coarse', 's'  ,'fine', 's',  ierr )
       ENDIF

       CALL pmc_set_dataarray_name( lastentry = .TRUE. )

!
!--    Send grid to parent
       val(1)  = nx
       val(2)  = ny
       val(3)  = nz
       val(4)  = dx
       val(5)  = dy
       fval(1) = zw(nzt+1)

       IF ( myid == 0 )  THEN

          CALL pmc_send_to_parent( val, SIZE( val ), 0, 123, ierr )
          CALL pmc_send_to_parent( fval, SIZE( fval ), 0, 124, ierr )
          CALL pmc_send_to_parent( coord_x, nx + 1 + 2 * nbgp, 0, 11, ierr )
          CALL pmc_send_to_parent( coord_y, ny + 1 + 2 * nbgp, 0, 12, ierr )

!
!--       Receive Coarse grid information.
!--       TO_DO: find shorter and more meaningful name for  define_coarse_grid_real
          CALL pmc_recv_from_parent( define_coarse_grid_real,                  &
                                     SIZE(define_coarse_grid_real), 0, 21, ierr )
          CALL pmc_recv_from_parent( define_coarse_grid_int,  3, 0, 22, ierr )
!
!--        Debug-printouts - keep them
!          WRITE(0,*) 'Coarse grid from parent '
!          WRITE(0,*) 'startx_tot    = ',define_coarse_grid_real(1)
!          WRITE(0,*) 'starty_tot    = ',define_coarse_grid_real(2)
!          WRITE(0,*) 'endx_tot      = ',define_coarse_grid_real(5)
!          WRITE(0,*) 'endy_tot      = ',define_coarse_grid_real(6)
!          WRITE(0,*) 'dx            = ',define_coarse_grid_real(3)
!          WRITE(0,*) 'dy            = ',define_coarse_grid_real(4)
!          WRITE(0,*) 'dz            = ',define_coarse_grid_real(7)
!          WRITE(0,*) 'nx_coarse     = ',define_coarse_grid_int(1)
!          WRITE(0,*) 'ny_coarse     = ',define_coarse_grid_int(2)
!          WRITE(0,*) 'nz_coarse     = ',define_coarse_grid_int(3)
       ENDIF

       CALL MPI_BCAST( define_coarse_grid_real, SIZE(define_coarse_grid_real),  &
                       MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( define_coarse_grid_int, 3, MPI_INTEGER, 0, comm2d, ierr )

       cg%dx = define_coarse_grid_real(3)
       cg%dy = define_coarse_grid_real(4)
       cg%dz = define_coarse_grid_real(7)
       cg%nx = define_coarse_grid_int(1)
       cg%ny = define_coarse_grid_int(2)
       cg%nz = define_coarse_grid_int(3)

!
!--    Get parent coordinates on coarse grid
       ALLOCATE( cg%coord_x(-nbgp:cg%nx+nbgp) )
       ALLOCATE( cg%coord_y(-nbgp:cg%ny+nbgp) )
      
       ALLOCATE( cg%dzu(1:cg%nz+1) )
       ALLOCATE( cg%dzw(1:cg%nz+1) )
       ALLOCATE( cg%zu(0:cg%nz+1) )
       ALLOCATE( cg%zw(0:cg%nz+1) )

!
!--    Get coarse grid coordinates and values of the z-direction from the parent
       IF ( myid == 0)  THEN

          CALL pmc_recv_from_parent( cg%coord_x, cg%nx+1+2*nbgp, 0, 24, ierr )
          CALL pmc_recv_from_parent( cg%coord_y, cg%ny+1+2*nbgp, 0, 25, ierr )
          CALL pmc_recv_from_parent( cg%dzu, cg%nz + 1, 0, 26, ierr )
          CALL pmc_recv_from_parent( cg%dzw, cg%nz + 1, 0, 27, ierr )
          CALL pmc_recv_from_parent( cg%zu, cg%nz + 2, 0, 28, ierr )
          CALL pmc_recv_from_parent( cg%zw, cg%nz + 2, 0, 29, ierr )

       ENDIF

!
!--    Broadcast this information
       CALL MPI_BCAST( cg%coord_x, cg%nx+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( cg%coord_y, cg%ny+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( cg%dzu, cg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( cg%dzw, cg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( cg%zu, cg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( cg%zw, cg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       
!
!--    Find the index bounds for the nest domain in the coarse-grid index space
       CALL pmci_map_fine_to_coarse_grid
!
!--    TO_DO: Klaus give a comment what is happening here
       CALL pmc_c_get_2d_index_list

!
!--    Include couple arrays into child content
!--    TO_DO: Klaus: better explain the above comment (what is child content?)
       CALL  pmc_c_clear_next_array_list
       DO  WHILE ( pmc_c_getnextarray( myname ) )
!--       TO_DO: Klaus, why the child-arrays are still up to cg%nz??
          CALL pmci_create_child_arrays ( myname, icl, icr, jcs, jcn, cg%nz )
       ENDDO
       CALL pmc_c_setind_and_allocmem

!
!--    Precompute interpolation coefficients and child-array indices
       CALL pmci_init_interp_tril

!
!--    Precompute the log-law correction index- and ratio-arrays 
       CALL pmci_init_loglaw_correction 

!
!--    Define the SGS-TKE scaling factor based on the grid-spacing ratio
       CALL pmci_init_tkefactor

!
!--    Two-way coupling for general and vertical nesting.
!--    Precompute the index arrays and relaxation functions for the
!--    anterpolation
       IF ( TRIM( nesting_mode ) == 'two-way' .OR.                              &
                  nesting_mode == 'vertical' )  THEN
          CALL pmci_init_anterp_tophat
       ENDIF

!
!--    Finally, compute the total area of the top-boundary face of the domain.
!--    This is needed in the pmc_ensure_nest_mass_conservation      
       area_t = ( nx + 1 ) * (ny + 1 ) * dx * dy

    ENDIF

 CONTAINS

    SUBROUTINE pmci_map_fine_to_coarse_grid
!
!--    Determine index bounds of interpolation/anterpolation area in the coarse
!--    grid index space
       IMPLICIT NONE

       INTEGER(iwp), DIMENSION(5,numprocs) ::  coarse_bound_all   !:
       INTEGER(iwp), DIMENSION(2)          ::  size_of_array      !:
                                              
       REAL(wp) ::  loffset     !:
       REAL(wp) ::  noffset     !:
       REAL(wp) ::  roffset     !:
       REAL(wp) ::  soffset     !:

!
!--    If the fine- and coarse grid nodes do not match:
       loffset = MOD( coord_x(nxl), cg%dx )
       xexl    = cg%dx + loffset
!
!--    This is needed in the anterpolation phase
       nhll = CEILING( xexl / cg%dx )
       xcs  = coord_x(nxl) - xexl
       DO  i = 0, cg%nx
          IF ( cg%coord_x(i) > xcs )  THEN
             icl = MAX( -1, i-1 )
             EXIT
          ENDIF
       ENDDO
!
!--    If the fine- and coarse grid nodes do not match
       roffset = MOD( coord_x(nxr+1), cg%dx )
       xexr    = cg%dx + roffset
!
!--    This is needed in the anterpolation phase
       nhlr = CEILING( xexr / cg%dx )
       xce  = coord_x(nxr+1) + xexr
!--    One "extra" layer is taken behind the right boundary 
!--    because it may be needed in cases of non-integer grid-spacing ratio
       DO  i = cg%nx, 0 , -1
          IF ( cg%coord_x(i) < xce )  THEN
             icr = MIN( cg%nx+1, i+1 )
             EXIT
          ENDIF
       ENDDO
!
!--    If the fine- and coarse grid nodes do not match
       soffset = MOD( coord_y(nys), cg%dy )
       yexs    = cg%dy + soffset
!
!--    This is needed in the anterpolation phase
       nhls = CEILING( yexs / cg%dy )
       ycs  = coord_y(nys) - yexs
       DO  j = 0, cg%ny
          IF ( cg%coord_y(j) > ycs )  THEN
             jcs = MAX( -nbgp, j-1 )
             EXIT
          ENDIF
       ENDDO
!
!--    If the fine- and coarse grid nodes do not match
       noffset = MOD( coord_y(nyn+1), cg%dy )
       yexn    = cg%dy + noffset
!
!--    This is needed in the anterpolation phase
       nhln = CEILING( yexn / cg%dy )
       yce  = coord_y(nyn+1) + yexn
!--    One "extra" layer is taken behind the north boundary 
!--    because it may be needed in cases of non-integer grid-spacing ratio
       DO  j = cg%ny, 0, -1
          IF ( cg%coord_y(j) < yce )  THEN
             jcn = MIN( cg%ny + nbgp, j+1 )
             EXIT
          ENDIF
       ENDDO

       coarse_bound(1) = icl
       coarse_bound(2) = icr
       coarse_bound(3) = jcs
       coarse_bound(4) = jcn
       coarse_bound(5) = myid
!
!--    Note that MPI_Gather receives data from all processes in the rank order
!--    TO_DO: refer to the line where this fact becomes important
       CALL MPI_GATHER( coarse_bound, 5, MPI_INTEGER, coarse_bound_all, 5,      &
                        MPI_INTEGER, 0, comm2d, ierr )

       IF ( myid == 0 )  THEN
          size_of_array(1) = SIZE( coarse_bound_all, 1 )
          size_of_array(2) = SIZE( coarse_bound_all, 2 )
          CALL pmc_send_to_parent( size_of_array, 2, 0, 40, ierr )
          CALL pmc_send_to_parent( coarse_bound_all, SIZE( coarse_bound_all ),  &
                                   0, 41, ierr )
       ENDIF

    END SUBROUTINE pmci_map_fine_to_coarse_grid



    SUBROUTINE pmci_init_interp_tril
!
!--    Precomputation of the interpolation coefficients and child-array indices
!--    to be used by the interpolation routines interp_tril_lr, interp_tril_ns
!--    and interp_tril_t.

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !:
       INTEGER(iwp) ::  i1      !:
       INTEGER(iwp) ::  j       !:
       INTEGER(iwp) ::  j1      !:
       INTEGER(iwp) ::  k       !:
       INTEGER(iwp) ::  kc      !:

       REAL(wp) ::  xb          !:
       REAL(wp) ::  xcsu        !:
       REAL(wp) ::  xfso        !:
       REAL(wp) ::  xcso        !:
       REAL(wp) ::  xfsu        !:
       REAL(wp) ::  yb          !:
       REAL(wp) ::  ycso        !:
       REAL(wp) ::  ycsv        !:
       REAL(wp) ::  yfso        !:
       REAL(wp) ::  yfsv        !:
       REAL(wp) ::  zcso        !:
       REAL(wp) ::  zcsw        !:
       REAL(wp) ::  zfso        !:
       REAL(wp) ::  zfsw        !:
      

       xb = nxl * dx
       yb = nys * dy
      
       ALLOCATE( icu(nxlg:nxrg) )
       ALLOCATE( ico(nxlg:nxrg) )
       ALLOCATE( jcv(nysg:nyng) )
       ALLOCATE( jco(nysg:nyng) )
       ALLOCATE( kcw(nzb:nzt+1) )
       ALLOCATE( kco(nzb:nzt+1) )
       ALLOCATE( r1xu(nxlg:nxrg) )
       ALLOCATE( r2xu(nxlg:nxrg) )
       ALLOCATE( r1xo(nxlg:nxrg) )
       ALLOCATE( r2xo(nxlg:nxrg) )
       ALLOCATE( r1yv(nysg:nyng) )
       ALLOCATE( r2yv(nysg:nyng) )
       ALLOCATE( r1yo(nysg:nyng) )
       ALLOCATE( r2yo(nysg:nyng) )
       ALLOCATE( r1zw(nzb:nzt+1) )
       ALLOCATE( r2zw(nzb:nzt+1) )
       ALLOCATE( r1zo(nzb:nzt+1) )
       ALLOCATE( r2zo(nzb:nzt+1) )

!
!--    Note that the node coordinates xfs... and xcs... are relative to the
!--    lower-left-bottom corner of the fc-array, not the actual child domain
!--    corner
       DO  i = nxlg, nxrg
          xfsu    = coord_x(i) - ( lower_left_coord_x + xb - xexl )
          xfso    = coord_x(i) + 0.5_wp * dx - ( lower_left_coord_x + xb - xexl )
          icu(i)  = icl + FLOOR( xfsu / cg%dx )
          ico(i)  = icl + FLOOR( ( xfso - 0.5_wp * cg%dx ) / cg%dx )
          xcsu    = ( icu(i) - icl ) * cg%dx
          xcso    = ( ico(i) - icl ) * cg%dx + 0.5_wp * cg%dx
          r2xu(i) = ( xfsu - xcsu ) / cg%dx
          r2xo(i) = ( xfso - xcso ) / cg%dx
          r1xu(i) = 1.0_wp - r2xu(i)
          r1xo(i) = 1.0_wp - r2xo(i)
       ENDDO

       DO  j = nysg, nyng
          yfsv    = coord_y(j) - ( lower_left_coord_y + yb - yexs )
          yfso    = coord_y(j) + 0.5_wp * dy - ( lower_left_coord_y + yb - yexs )
          jcv(j)  = jcs + FLOOR( yfsv / cg%dy )
          jco(j)  = jcs + FLOOR( ( yfso -0.5_wp * cg%dy ) / cg%dy )
          ycsv    = ( jcv(j) - jcs ) * cg%dy
          ycso    = ( jco(j) - jcs ) * cg%dy + 0.5_wp * cg%dy
          r2yv(j) = ( yfsv - ycsv ) / cg%dy
          r2yo(j) = ( yfso - ycso ) / cg%dy
          r1yv(j) = 1.0_wp - r2yv(j)
          r1yo(j) = 1.0_wp - r2yo(j)
       ENDDO

       DO  k = nzb, nzt + 1
          zfsw = zw(k)
          zfso = zu(k)

          kc = 0
          DO  WHILE ( cg%zw(kc) <= zfsw )
             kc = kc + 1
          ENDDO
          kcw(k) = kc - 1
         
          kc = 0
          DO  WHILE ( cg%zu(kc) <= zfso )
             kc = kc + 1
          ENDDO
          kco(k) = kc - 1

          zcsw    = cg%zw(kcw(k))
          zcso    = cg%zu(kco(k))
          r2zw(k) = ( zfsw - zcsw ) / cg%dzw(kcw(k)+1)
          r2zo(k) = ( zfso - zcso ) / cg%dzu(kco(k)+1)
          r1zw(k) = 1.0_wp - r2zw(k)
          r1zo(k) = 1.0_wp - r2zo(k)
       ENDDO
      
    END SUBROUTINE pmci_init_interp_tril



    SUBROUTINE pmci_init_loglaw_correction
!
!--    Precomputation of the index and log-ratio arrays for the log-law
!--    corrections for near-wall nodes after the nest-BC interpolation.
!--    These are used by the interpolation routines interp_tril_lr and
!--    interp_tril_ns.

       IMPLICIT NONE

       INTEGER(iwp) ::  direction    !:  Wall normal index: 1=k, 2=j, 3=i.
       INTEGER(iwp) ::  i            !:
       INTEGER(iwp) ::  icorr        !:
       INTEGER(iwp) ::  inc          !:  Wall outward-normal index increment -1
                                     !: or 1, for direction=1, inc=1 always
       INTEGER(iwp) ::  iw           !:
       INTEGER(iwp) ::  j            !:
       INTEGER(iwp) ::  jcorr        !:
       INTEGER(iwp) ::  jw           !:
       INTEGER(iwp) ::  k            !:
       INTEGER(iwp) ::  kb           !:
       INTEGER(iwp) ::  kcorr        !:
       INTEGER(iwp) ::  lc           !:
       INTEGER(iwp) ::  ni           !:
       INTEGER(iwp) ::  nj           !:
       INTEGER(iwp) ::  nk           !:
       INTEGER(iwp) ::  nzt_topo_max !:
       INTEGER(iwp) ::  wall_index   !:  Index of the wall-node coordinate

       REAL(wp), ALLOCATABLE, DIMENSION(:) ::  lcr   !:

!
!--    First determine the maximum k-index needed for the near-wall corrections.
!--    This maximum is individual for each boundary to minimize the storage
!--    requirements and to minimize the corresponding loop k-range in the
!--    interpolation routines.
       nzt_topo_nestbc_l = nzb
       IF ( nest_bound_l )  THEN
          DO  i = nxl-1, nxl
             DO  j = nys, nyn
                nzt_topo_nestbc_l = MAX( nzt_topo_nestbc_l, nzb_u_inner(j,i),   &
                                         nzb_v_inner(j,i), nzb_w_inner(j,i) )
             ENDDO
          ENDDO
          nzt_topo_nestbc_l = nzt_topo_nestbc_l + 1
       ENDIF
      
       nzt_topo_nestbc_r = nzb
       IF ( nest_bound_r )  THEN
          i = nxr + 1
          DO  j = nys, nyn
             nzt_topo_nestbc_r = MAX( nzt_topo_nestbc_r, nzb_u_inner(j,i),      &
                                      nzb_v_inner(j,i), nzb_w_inner(j,i) )
          ENDDO
          nzt_topo_nestbc_r = nzt_topo_nestbc_r + 1
       ENDIF

       nzt_topo_nestbc_s = nzb
       IF ( nest_bound_s )  THEN
          DO  j = nys-1, nys
             DO  i = nxl, nxr
                nzt_topo_nestbc_s = MAX( nzt_topo_nestbc_s, nzb_u_inner(j,i),   &
                                         nzb_v_inner(j,i), nzb_w_inner(j,i) )
             ENDDO
          ENDDO
          nzt_topo_nestbc_s = nzt_topo_nestbc_s + 1
       ENDIF

       nzt_topo_nestbc_n = nzb
       IF ( nest_bound_n )  THEN
          j = nyn + 1
          DO  i = nxl, nxr
             nzt_topo_nestbc_n = MAX( nzt_topo_nestbc_n, nzb_u_inner(j,i),      &
                                      nzb_v_inner(j,i), nzb_w_inner(j,i) )
          ENDDO
          nzt_topo_nestbc_n = nzt_topo_nestbc_n + 1
       ENDIF

!
!--    Then determine the maximum number of near-wall nodes per wall point based
!--    on the grid-spacing ratios.
       nzt_topo_max = MAX( nzt_topo_nestbc_l, nzt_topo_nestbc_r,                &
                           nzt_topo_nestbc_s, nzt_topo_nestbc_n )

!
!--    Note that the outer division must be integer division.
       ni = CEILING( cg%dx / dx ) / 2
       nj = CEILING( cg%dy / dy ) / 2
       nk = 1
       DO  k = 1, nzt_topo_max
          nk = MAX( nk, CEILING( cg%dzu(k) / dzu(k) ) )
       ENDDO
       nk = nk / 2   !  Note that this must be integer division.
       ncorr =  MAX( ni, nj, nk )

       ALLOCATE( lcr(0:ncorr-1) )
       lcr = 1.0_wp

!
!--    First horizontal walls. Note that also logc_w_? and logc_ratio_w_? need to 
!--    be allocated and initialized here.
!--    Left boundary
       IF ( nest_bound_l )  THEN

          ALLOCATE( logc_u_l(1:2,nzb:nzt_topo_nestbc_l,nys:nyn) )
          ALLOCATE( logc_v_l(1:2,nzb:nzt_topo_nestbc_l,nys:nyn) )
          ALLOCATE( logc_w_l(1:2,nzb:nzt_topo_nestbc_l,nys:nyn) )
          ALLOCATE( logc_ratio_u_l(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_l,nys:nyn) )
          ALLOCATE( logc_ratio_v_l(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_l,nys:nyn) )
          ALLOCATE( logc_ratio_w_l(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_l,nys:nyn) )
          logc_u_l       = 0
          logc_v_l       = 0
          logc_w_l       = 0
          logc_ratio_u_l = 1.0_wp
          logc_ratio_v_l = 1.0_wp
          logc_ratio_w_l = 1.0_wp
          direction      = 1
          inc            = 1

          DO  j = nys, nyn
!
!--          Left boundary for u
             i   = 0
             kb  = nzb_u_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_u_l(1,k,j) = lc
             logc_ratio_u_l(1,0:ncorr-1,k,j) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp
!
!--          Left boundary for v
             i   = -1
             kb  =  nzb_v_inner(j,i)
             k   =  kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_v_l(1,k,j) = lc
             logc_ratio_v_l(1,0:ncorr-1,k,j) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp

          ENDDO

       ENDIF

!
!--    Right boundary
       IF ( nest_bound_r )  THEN
           
          ALLOCATE( logc_u_r(1:2,nzb:nzt_topo_nestbc_r,nys:nyn) )
          ALLOCATE( logc_v_r(1:2,nzb:nzt_topo_nestbc_r,nys:nyn) )
          ALLOCATE( logc_w_r(1:2,nzb:nzt_topo_nestbc_r,nys:nyn) )
          ALLOCATE( logc_ratio_u_r(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_r,nys:nyn) )
          ALLOCATE( logc_ratio_v_r(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_r,nys:nyn) )
          ALLOCATE( logc_ratio_w_r(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_r,nys:nyn) )
          logc_u_r       = 0
          logc_v_r       = 0
          logc_w_r       = 0
          logc_ratio_u_r = 1.0_wp
          logc_ratio_v_r = 1.0_wp
          logc_ratio_w_r = 1.0_wp
          direction      = 1
          inc            = 1

          DO  j = nys, nyn
!
!--          Right boundary for u
             i   = nxr + 1
             kb  = nzb_u_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_u_r(1,k,j) = lc
             logc_ratio_u_r(1,0:ncorr-1,k,j) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp
!
!--          Right boundary for v
             i   = nxr + 1
             kb  = nzb_v_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_v_r(1,k,j) = lc
             logc_ratio_v_r(1,0:ncorr-1,k,j) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp

          ENDDO

       ENDIF

!
!--    South boundary
       IF ( nest_bound_s )  THEN

          ALLOCATE( logc_u_s(1:2,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          ALLOCATE( logc_v_s(1:2,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          ALLOCATE( logc_w_s(1:2,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          ALLOCATE( logc_ratio_u_s(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          ALLOCATE( logc_ratio_v_s(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          ALLOCATE( logc_ratio_w_s(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_s,nxl:nxr) )
          logc_u_s       = 0
          logc_v_s       = 0
          logc_w_s       = 0
          logc_ratio_u_s = 1.0_wp
          logc_ratio_v_s = 1.0_wp
          logc_ratio_w_s = 1.0_wp
          direction      = 1
          inc            = 1

          DO  i = nxl, nxr
!
!--          South boundary for u
             j   = -1
             kb  =  nzb_u_inner(j,i)
             k   =  kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_u_s(1,k,i) = lc
             logc_ratio_u_s(1,0:ncorr-1,k,i) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp
!
!--          South boundary for v
             j   = 0
             kb  = nzb_v_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_v_s(1,k,i) = lc
             logc_ratio_v_s(1,0:ncorr-1,k,i) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp

          ENDDO

       ENDIF

!
!--    North boundary
       IF ( nest_bound_n )  THEN

          ALLOCATE( logc_u_n(1:2,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          ALLOCATE( logc_v_n(1:2,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          ALLOCATE( logc_w_n(1:2,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          ALLOCATE( logc_ratio_u_n(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          ALLOCATE( logc_ratio_v_n(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          ALLOCATE( logc_ratio_w_n(1:2,0:ncorr-1,nzb:nzt_topo_nestbc_n,nxl:nxr) )
          logc_u_n       = 0
          logc_v_n       = 0
          logc_w_n       = 0
          logc_ratio_u_n = 1.0_wp
          logc_ratio_v_n = 1.0_wp
          logc_ratio_w_n = 1.0_wp
          direction      = 1
          inc            = 1

          DO  i = nxl, nxr
!
!--          North boundary for u
             j   = nyn + 1
             kb  = nzb_u_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_u_n(1,k,i) = lc
             logc_ratio_u_n(1,0:ncorr-1,k,i) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp
!
!--          North boundary for v
             j   = nyn + 1
             kb  = nzb_v_inner(j,i)
             k   = kb + 1
             wall_index = kb
             CALL pmci_define_loglaw_correction_parameters( lc, lcr, k, j,      &
                                inc, wall_index, z0(j,i), kb, direction, ncorr )
             logc_v_n(1,k,i) = lc
             logc_ratio_v_n(1,0:ncorr-1,k,i) = lcr(0:ncorr-1)
             lcr(0:ncorr-1) = 1.0_wp

          ENDDO

       ENDIF

!       
!--    Then vertical walls and corners if necessary
       IF ( topography /= 'flat' )  THEN

          kb = 0       ! kb is not used when direction > 1
!        
!--       Left boundary
          IF ( nest_bound_l )  THEN

             direction  = 2

             DO  j = nys, nyn
                DO  k = nzb, nzt_topo_nestbc_l
!
!--                Wall for u on the south side, but not on the north side
                   i  = 0
                   IF ( ( nzb_u_outer(j,i) > nzb_u_outer(j+1,i) ) .AND.         &
                        ( nzb_u_outer(j,i) == nzb_u_outer(j-1,i) ) )            &
                   THEN
                      inc        =  1
                      wall_index =  j
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_u_l(2,k,j) = inc * lc
                      logc_ratio_u_l(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for u on the north side, but not on the south side
                   i  = 0
                   IF ( ( nzb_u_outer(j,i) > nzb_u_outer(j-1,i) ) .AND.         &
                        ( nzb_u_outer(j,i) == nzb_u_outer(j+1,i) ) )  THEN
                      inc        = -1
                      wall_index =  j + 1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_u_l(2,k,j) = inc * lc
                      logc_ratio_u_l(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the south side, but not on the north side.
                   i  = -1
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j+1,i) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j-1,i) ) )  THEN
                      inc        =  1
                      wall_index =  j
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_l(2,k,j) = inc * lc
                      logc_ratio_w_l(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the north side, but not on the south side.
                   i  = -1
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j-1,i) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j+1,i) ) )  THEN
                      inc        = -1
                      wall_index =  j+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_l(2,k,j) = inc * lc
                      logc_ratio_w_l(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

                ENDDO
             ENDDO

          ENDIF   !  IF ( nest_bound_l )

!        
!--       Right boundary
          IF ( nest_bound_r )  THEN

             direction  = 2
             i  = nxr + 1

             DO  j = nys, nyn
                DO  k = nzb, nzt_topo_nestbc_r
!
!--                Wall for u on the south side, but not on the north side
                   IF ( ( nzb_u_outer(j,i) > nzb_u_outer(j+1,i) )  .AND.        &
                        ( nzb_u_outer(j,i) == nzb_u_outer(j-1,i) ) )  THEN
                      inc        = 1
                      wall_index = j
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_u_r(2,k,j) = inc * lc
                      logc_ratio_u_r(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for u on the north side, but not on the south side
                   IF ( ( nzb_u_outer(j,i) > nzb_u_outer(j-1,i) )  .AND.        &
                        ( nzb_u_outer(j,i) == nzb_u_outer(j+1,i) ) )  THEN
                      inc        = -1
                      wall_index =  j+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_u_r(2,k,j) = inc * lc
                      logc_ratio_u_r(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the south side, but not on the north side
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j+1,i) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j-1,i) ) )  THEN
                      inc        =  1
                      wall_index =  j
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_r(2,k,j) = inc * lc
                      logc_ratio_w_r(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the north side, but not on the south side
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j-1,i) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j+1,i) ) )  THEN
                      inc        = -1
                      wall_index =  j+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, j, inc, wall_index, z0(j,i), kb, direction, ncorr )

!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_r(2,k,j) = inc * lc
                      logc_ratio_w_r(2,0:ncorr-1,k,j) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

                ENDDO
             ENDDO

          ENDIF   !  IF ( nest_bound_r )

!        
!--       South boundary
          IF ( nest_bound_s )  THEN

             direction  = 3

             DO  i = nxl, nxr
                DO  k = nzb, nzt_topo_nestbc_s
!
!--                Wall for v on the left side, but not on the right side
                   j  = 0
                   IF ( ( nzb_v_outer(j,i) > nzb_v_outer(j,i+1) )  .AND.        &
                        ( nzb_v_outer(j,i) == nzb_v_outer(j,i-1) ) )  THEN
                      inc        =  1
                      wall_index =  i
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_v_s(2,k,i) = inc * lc
                      logc_ratio_v_s(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for v on the right side, but not on the left side
                   j  = 0
                   IF ( ( nzb_v_outer(j,i) > nzb_v_outer(j,i-1) )  .AND.        &
                        ( nzb_v_outer(j,i) == nzb_v_outer(j,i+1) ) )  THEN
                      inc        = -1
                      wall_index =  i+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_v_s(2,k,i) = inc * lc
                      logc_ratio_v_s(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the left side, but not on the right side
                   j  = -1
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j,i+1) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j,i-1) ) )  THEN
                      inc        =  1
                      wall_index =  i
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_s(2,k,i) = inc * lc
                      logc_ratio_w_s(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the right side, but not on the left side
                   j  = -1
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j,i-1) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j,i+1) ) )  THEN
                      inc        = -1
                      wall_index =  i+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_s(2,k,i) = inc * lc
                      logc_ratio_w_s(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

                ENDDO
             ENDDO

          ENDIF   !  IF (nest_bound_s )

!        
!--       North boundary
          IF ( nest_bound_n )  THEN

             direction  = 3
             j  = nyn + 1

             DO  i = nxl, nxr
                DO  k = nzb, nzt_topo_nestbc_n
!
!--                Wall for v on the left side, but not on the right side
                   IF ( ( nzb_v_outer(j,i) > nzb_v_outer(j,i+1) )  .AND.        &
                        ( nzb_v_outer(j,i) == nzb_v_outer(j,i-1) ) )  THEN
                      inc        = 1
                      wall_index = i
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_v_n(2,k,i) = inc * lc
                      logc_ratio_v_n(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for v on the right side, but not on the left side
                   IF ( ( nzb_v_outer(j,i) > nzb_v_outer(j,i-1) )  .AND.        &
                        ( nzb_v_outer(j,i) == nzb_v_outer(j,i+1) ) )  THEN
                      inc        = -1
                      wall_index =  i + 1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_v_n(2,k,i) = inc * lc
                      logc_ratio_v_n(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the left side, but not on the right side
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j,i+1) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j,i-1) ) )  THEN
                      inc        = 1
                      wall_index = i
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_n(2,k,i) = inc * lc
                      logc_ratio_w_n(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

!
!--                Wall for w on the right side, but not on the left side
                   IF ( ( nzb_w_outer(j,i) > nzb_w_outer(j,i-1) )  .AND.        &
                        ( nzb_w_outer(j,i) == nzb_w_outer(j,i+1) ) )  THEN
                      inc        = -1
                      wall_index =  i+1
                      CALL pmci_define_loglaw_correction_parameters( lc, lcr,   &
                          k, i, inc, wall_index, z0(j,i), kb, direction, ncorr )
!
!--                   The direction of the wall-normal index is stored as the
!--                   sign of the logc-element.
                      logc_w_n(2,k,i) = inc * lc
                      logc_ratio_w_n(2,0:ncorr-1,k,i) = lcr(0:ncorr-1)
                      lcr(0:ncorr-1) = 1.0_wp
                   ENDIF

                ENDDO
             ENDDO

          ENDIF   !  IF ( nest_bound_n )

       ENDIF   !  IF ( topography /= 'flat' )

    END SUBROUTINE pmci_init_loglaw_correction



    SUBROUTINE pmci_define_loglaw_correction_parameters( lc, lcr, k, ij, inc,   &
                                        wall_index, z0_l, kb, direction, ncorr )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)  ::  direction                 !:
       INTEGER(iwp), INTENT(IN)  ::  ij                        !:
       INTEGER(iwp), INTENT(IN)  ::  inc                       !:
       INTEGER(iwp), INTENT(IN)  ::  k                         !:
       INTEGER(iwp), INTENT(IN)  ::  kb                        !:
       INTEGER(iwp), INTENT(OUT) ::  lc                        !:
       INTEGER(iwp), INTENT(IN)  ::  ncorr                     !:
       INTEGER(iwp), INTENT(IN)  ::  wall_index                !:

       INTEGER(iwp) ::  alcorr       !:
       INTEGER(iwp) ::  corr_index   !:
       INTEGER(iwp) ::  lcorr        !:

       LOGICAL      ::  more         !:

       REAL(wp), DIMENSION(0:ncorr-1), INTENT(OUT) ::  lcr     !:
       REAL(wp), INTENT(IN)      ::  z0_l                      !:
      
       REAL(wp)     ::  logvelc1     !:
      

       SELECT CASE ( direction )

          CASE (1)   !  k
             more = .TRUE.
             lcorr = 0
             DO  WHILE ( more .AND. lcorr <= ncorr-1 )
                corr_index = k + lcorr
                IF ( lcorr == 0 )  THEN
                   CALL pmci_find_logc_pivot_k( lc, logvelc1, z0_l, kb )
                ENDIF
               
                IF ( corr_index < lc )  THEN
                   lcr(lcorr) = LOG( ( zu(k) - zw(kb) ) / z0_l ) / logvelc1
                   more = .TRUE.
                ELSE
                   lcr(lcorr) = 1.0
                   more = .FALSE.
                ENDIF
                lcorr = lcorr + 1
             ENDDO

          CASE (2)   !  j
             more = .TRUE.
             lcorr  = 0
             alcorr = 0
             DO  WHILE ( more  .AND.  alcorr <= ncorr-1 )
                corr_index = ij + lcorr   ! In this case (direction = 2) ij is j
                IF ( lcorr == 0 )  THEN
                   CALL pmci_find_logc_pivot_j( lc, logvelc1, ij, wall_index,   &
                                                z0_l, inc )
                ENDIF

!
!--             The role of inc here is to make the comparison operation "<"
!--             valid in both directions
                IF ( inc * corr_index < inc * lc )  THEN
                   lcr(alcorr) = LOG( ABS( coord_y(corr_index) + 0.5_wp * dy    &
                                         - coord_y(wall_index) ) / z0_l )       &
                                 / logvelc1
                   more = .TRUE.
                ELSE
                   lcr(alcorr) = 1.0_wp
                   more = .FALSE.
                ENDIF
                lcorr  = lcorr + inc
                alcorr = ABS( lcorr )
             ENDDO

          CASE (3)   !  i
             more = .TRUE.
             lcorr  = 0
             alcorr = 0
             DO  WHILE ( more  .AND.  alcorr <= ncorr-1 )
                corr_index = ij + lcorr   ! In this case (direction = 3) ij is i
                IF ( lcorr == 0 )  THEN
                   CALL pmci_find_logc_pivot_i( lc, logvelc1, ij, wall_index,   &
                                                z0_l, inc )
                ENDIF
!
!--             The role of inc here is to make the comparison operation "<"
!--             valid in both directions
                IF ( inc * corr_index < inc * lc )  THEN
                   lcr(alcorr) = LOG( ABS( coord_x(corr_index) + 0.5_wp * dx    &
                                         - coord_x(wall_index) ) / z0_l )       &
                                 / logvelc1
                   more = .TRUE.
                ELSE
                   lcr(alcorr) = 1.0_wp
                   more = .FALSE.
                ENDIF
                lcorr  = lcorr + inc
                alcorr = ABS( lcorr )
             ENDDO

       END SELECT

    END SUBROUTINE pmci_define_loglaw_correction_parameters



    SUBROUTINE pmci_find_logc_pivot_k( lc, logzc1, z0_l, kb )
!
!--    Finds the pivot node and te log-law factor for near-wall nodes for
!--    which the wall-parallel velocity components will be log-law corrected
!--    after interpolation. This subroutine is only for horizontal walls.

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)  ::  kb   !:
       INTEGER(iwp), INTENT(OUT) ::  lc   !:

       INTEGER(iwp) ::  kbc    !:
       INTEGER(iwp) ::  k1     !:

       REAL(wp), INTENT(OUT) ::  logzc1     !:
       REAL(wp), INTENT(IN) ::  z0_l       !:

       REAL(wp) ::  zuc1   !:


       kbc = nzb + 1
!
!--    kbc is the first coarse-grid point above the surface
       DO  WHILE ( cg%zu(kbc) < zu(kb) )
          kbc = kbc + 1
       ENDDO
       zuc1  = cg%zu(kbc)
       k1    = kb + 1
       DO  WHILE ( zu(k1) < zuc1 )  !  Important: must be <, not <=
          k1 = k1 + 1
       ENDDO
       logzc1 = LOG( (zu(k1) - zw(kb) ) / z0_l )
       lc = k1

    END SUBROUTINE pmci_find_logc_pivot_k



    SUBROUTINE pmci_find_logc_pivot_j( lc, logyc1, j, jw, z0_l, inc )
!
!--    Finds the pivot node and te log-law factor for near-wall nodes for
!--    which the wall-parallel velocity components will be log-law corrected
!--    after interpolation. This subroutine is only for vertical walls on
!--    south/north sides of the node.

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)  ::  inc    !:  increment must be 1 or -1.
       INTEGER(iwp), INTENT(IN)  ::  j      !:
       INTEGER(iwp), INTENT(IN)  ::  jw     !:
       INTEGER(iwp), INTENT(OUT) ::  lc     !:

       INTEGER(iwp) ::  j1       !:

       REAL(wp), INTENT(IN) ::  z0_l   !:

       REAL(wp) ::  logyc1   !:
       REAL(wp) ::  yc1      !:

!
!--    yc1 is the y-coordinate of the first coarse-grid u- and w-nodes out from
!--    the wall
       yc1  = coord_y(jw) + 0.5_wp * inc * cg%dy
!
!--    j1 is the first fine-grid index further away from the wall than yc1
       j1 = j
!
!--    Important: must be <, not <=
       DO  WHILE ( inc * ( coord_y(j1) + 0.5_wp * dy ) < inc * yc1 )
          j1 = j1 + inc
       ENDDO

       logyc1 = LOG( ABS( coord_y(j1) + 0.5_wp * dy - coord_y(jw) ) / z0_l )
       lc = j1

    END SUBROUTINE pmci_find_logc_pivot_j



    SUBROUTINE pmci_find_logc_pivot_i( lc, logxc1, i, iw, z0_l, inc )
!
!--    Finds the pivot node and the log-law factor for near-wall nodes for
!--    which the wall-parallel velocity components will be log-law corrected
!--    after interpolation. This subroutine is only for vertical walls on
!--    south/north sides of the node.

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)  ::  i      !:
       INTEGER(iwp), INTENT(IN)  ::  inc    !: increment must be 1 or -1.
       INTEGER(iwp), INTENT(IN)  ::  iw     !:
       INTEGER(iwp), INTENT(OUT) ::  lc     !:

       INTEGER(iwp) ::  i1       !:

       REAL(wp), INTENT(IN) ::  z0_l   !:

       REAL(wp) ::  logxc1   !:
       REAL(wp) ::  xc1      !:

!
!--    xc1 is the x-coordinate of the first coarse-grid v- and w-nodes out from
!--    the wall
       xc1  = coord_x(iw) + 0.5_wp * inc * cg%dx
!
!--    i1 is the first fine-grid index futher away from the wall than xc1.
       i1 = i
!
!--    Important: must be <, not <=
       DO  WHILE ( inc * ( coord_x(i1) + 0.5_wp *dx ) < inc * xc1 )
          i1 = i1 + inc
       ENDDO
      
       logxc1 = LOG( ABS( coord_x(i1) + 0.5_wp*dx - coord_x(iw) ) / z0_l )
       lc = i1

    END SUBROUTINE pmci_find_logc_pivot_i




    SUBROUTINE pmci_init_anterp_tophat
!
!--    Precomputation of the child-array indices for
!--    corresponding coarse-grid array index and the
!--    Under-relaxation coefficients to be used by anterp_tophat.

       IMPLICIT NONE

       INTEGER(iwp) ::  i        !: Fine-grid index
       INTEGER(iwp) ::  ifc_o    !:
       INTEGER(iwp) ::  ifc_u    !:
       INTEGER(iwp) ::  ii       !: Coarse-grid index
       INTEGER(iwp) ::  istart   !:
       INTEGER(iwp) ::  j        !: Fine-grid index
       INTEGER(iwp) ::  jj       !: Coarse-grid index
       INTEGER(iwp) ::  jstart   !:
       INTEGER(iwp) ::  k        !: Fine-grid index
       INTEGER(iwp) ::  kk       !: Coarse-grid index
       INTEGER(iwp) ::  kstart   !:
       REAL(wp)     ::  xi       !:
       REAL(wp)     ::  eta      !:
       REAL(wp)     ::  zeta     !:
     
!
!--    Default values:
       IF ( anterp_relax_length_l < 0.0_wp )  THEN
          anterp_relax_length_l = 0.1_wp * ( nx + 1 ) * dx
       ENDIF
       IF ( anterp_relax_length_r < 0.0_wp )  THEN
          anterp_relax_length_r = 0.1_wp * ( nx + 1 ) * dx
       ENDIF
       IF ( anterp_relax_length_s < 0.0_wp )  THEN
          anterp_relax_length_s = 0.1_wp * ( ny + 1 ) * dy
       ENDIF
       IF ( anterp_relax_length_n < 0.0_wp )  THEN
          anterp_relax_length_n = 0.1_wp * ( ny + 1 ) * dy
       ENDIF
       IF ( anterp_relax_length_t < 0.0_wp )  THEN
          anterp_relax_length_t = 0.1_wp * zu(nzt)
       ENDIF

!
!--    First determine kctu and kctw that are the coarse-grid upper bounds for
!--    index k
       kk = 0
       DO  WHILE ( cg%zu(kk) <= zu(nzt) )
          kk = kk + 1
       ENDDO
       kctu = kk - 1

       kk = 0
       DO  WHILE ( cg%zw(kk) <= zw(nzt-1) )
          kk = kk + 1
       ENDDO
       kctw = kk - 1

       ALLOCATE( iflu(icl:icr) )
       ALLOCATE( iflo(icl:icr) )
       ALLOCATE( ifuu(icl:icr) )
       ALLOCATE( ifuo(icl:icr) )
       ALLOCATE( jflv(jcs:jcn) )
       ALLOCATE( jflo(jcs:jcn) )
       ALLOCATE( jfuv(jcs:jcn) )
       ALLOCATE( jfuo(jcs:jcn) )
       ALLOCATE( kflw(0:kctw) )
       ALLOCATE( kflo(0:kctu) )
       ALLOCATE( kfuw(0:kctw) )
       ALLOCATE( kfuo(0:kctu) )

       ALLOCATE( ijfc_u(jcs:jcn,icl:icr) )
       ALLOCATE( ijfc_v(jcs:jcn,icl:icr) )
       ALLOCATE( ijfc_s(jcs:jcn,icl:icr) )

!
!--    i-indices of u for each ii-index value
!--    ii=icr is redundant for anterpolation
       istart = nxlg
       DO  ii = icl, icr-1
          i = istart
          DO  WHILE ( ( coord_x(i) < cg%coord_x(ii) - 0.5_wp * cg%dx )  .AND.   &
                      ( i < nxrg ) )
             i = i + 1
          ENDDO
          iflu(ii) = MIN( MAX( i, nxlg ), nxrg )
          DO  WHILE ( ( coord_x(i) <= cg%coord_x(ii) + 0.5_wp * cg%dx )  .AND.  &
                      ( i < nxrg+1 ) )
             i = i + 1
          ENDDO
          ifuu(ii) = MIN( MAX( i-1, iflu(ii) ), nxrg )
          istart = iflu(ii)
       ENDDO
       iflu(icr) = nxrg
       ifuu(icr) = nxrg

!
!--    i-indices of others for each ii-index value
!--    ii=icr is redundant for anterpolation
       istart = nxlg
       DO  ii = icl, icr-1
          i = istart
          DO  WHILE ( ( coord_x(i) + 0.5_wp * dx < cg%coord_x(ii) )  .AND.      &
                      ( i < nxrg ) )
             i = i + 1
          ENDDO
          iflo(ii) = MIN( MAX( i, nxlg ), nxrg )
          DO  WHILE ( ( coord_x(i) + 0.5_wp * dx <= cg%coord_x(ii) + cg%dx )    &
                      .AND.  ( i < nxrg+1 ) )
             i = i + 1
          ENDDO
          ifuo(ii) = MIN( MAX( i-1, iflo(ii) ), nxrg )
          istart = iflo(ii)
       ENDDO
       iflo(icr) = nxrg
       ifuo(icr) = nxrg

!
!--    j-indices of v for each jj-index value
!--    jj=jcn is redundant for anterpolation
       jstart = nysg
       DO  jj = jcs, jcn-1
          j = jstart
          DO  WHILE ( ( coord_y(j) < cg%coord_y(jj) - 0.5_wp * cg%dy )  .AND.   &
                      ( j < nyng ) )
             j = j + 1
          ENDDO
          jflv(jj) = MIN( MAX( j, nysg ), nyng )
          DO  WHILE ( ( coord_y(j) <= cg%coord_y(jj) + 0.5_wp * cg%dy )  .AND.  &
                      ( j < nyng+1 ) )
             j = j + 1
          ENDDO
          jfuv(jj) = MIN( MAX( j-1, jflv(jj) ), nyng )
          jstart = jflv(jj)
       ENDDO
       jflv(jcn) = nyng
       jfuv(jcn) = nyng
!
!--    j-indices of others for each jj-index value
!--    jj=jcn is redundant for anterpolation
       jstart = nysg
       DO  jj = jcs, jcn-1
          j = jstart
          DO  WHILE ( ( coord_y(j) + 0.5_wp * dy < cg%coord_y(jj) )  .AND.      &
                      ( j < nyng ) )
             j = j + 1
          ENDDO
          jflo(jj) = MIN( MAX( j, nysg ), nyng )
          DO  WHILE ( ( coord_y(j) + 0.5_wp * dy <= cg%coord_y(jj) + cg%dy )    &
                      .AND.  ( j < nyng+1 ) )
             j = j + 1
          ENDDO
          jfuo(jj) = MIN( MAX( j-1, jflo(jj) ), nyng )
          jstart = jflo(jj)
       ENDDO
       jflo(jcn) = nyng
       jfuo(jcn) = nyng

!
!--    k-indices of w for each kk-index value
       kstart  = 0
       kflw(0) = 0
       kfuw(0) = 0
       DO  kk = 1, kctw
          k = kstart
          DO  WHILE ( ( zw(k) < cg%zu(kk) )  .AND.  ( k < nzt ) )
             k = k + 1
          ENDDO
          kflw(kk) = MIN( MAX( k, 1 ), nzt + 1 )
          DO  WHILE ( ( zw(k) <= cg%zu(kk+1) )  .AND.  ( k < nzt+1 ) )
             k = k + 1
          ENDDO
          kfuw(kk) = MIN( MAX( k-1, kflw(kk) ), nzt + 1 )
          kstart = kflw(kk)
       ENDDO

!
!--    k-indices of others for each kk-index value
       kstart  = 0
       kflo(0) = 0
       kfuo(0) = 0
       DO  kk = 1, kctu
          k = kstart
          DO  WHILE ( ( zu(k) < cg%zw(kk-1) )  .AND.  ( k < nzt ) )
             k = k + 1
          ENDDO
          kflo(kk) = MIN( MAX( k, 1 ), nzt + 1 )
          DO  WHILE ( ( zu(k) <= cg%zw(kk) )  .AND.  ( k < nzt+1 ) )
             k = k + 1
          ENDDO
          kfuo(kk) = MIN( MAX( k-1, kflo(kk) ), nzt + 1 )
          kstart = kflo(kk)
       ENDDO
 
!
!--    Precomputation of number of fine-grid nodes inside coarse-grid ij-faces.
!--    Note that ii, jj, and kk are coarse-grid indices. 
!--    This information is needed in anterpolation.
       DO  ii = icl, icr
          ifc_u = ifuu(ii) - iflu(ii) + 1
          ifc_o = ifuo(ii) - iflo(ii) + 1
          DO  jj = jcs, jcn
             ijfc_u(jj,ii) = ifc_u * ( jfuo(jj) - jflo(jj) + 1 )
             ijfc_v(jj,ii) = ifc_o * ( jfuv(jj) - jflv(jj) + 1 )
             ijfc_s(jj,ii) = ifc_o * ( jfuo(jj) - jflo(jj) + 1 )             
          ENDDO
       ENDDO
    
!
!--    Spatial under-relaxation coefficients
       ALLOCATE( frax(icl:icr) )
       ALLOCATE( fray(jcs:jcn) )
       
       frax(icl:icr) = 1.0_wp
       fray(jcs:jcn) = 1.0_wp

       IF ( nesting_mode /= 'vertical' )  THEN
          DO  ii = icl, icr
             IF ( nest_bound_l )  THEN
                xi = ( MAX( 0.0_wp, ( cg%coord_x(ii) -                          &
                       lower_left_coord_x ) ) / anterp_relax_length_l )**4
             ELSEIF ( nest_bound_r )  THEN
                xi = ( MAX( 0.0_wp, ( lower_left_coord_x + ( nx + 1 ) * dx -    &
                                      cg%coord_x(ii) ) ) /                      &
                       anterp_relax_length_r )**4
             ELSE
                xi = 999999.9_wp
             ENDIF
             frax(ii) = xi / ( 1.0_wp + xi )
          ENDDO


          DO  jj = jcs, jcn
             IF ( nest_bound_s )  THEN
                eta = ( MAX( 0.0_wp, ( cg%coord_y(jj) -                         &
                        lower_left_coord_y ) ) / anterp_relax_length_s )**4
             ELSEIF ( nest_bound_n )  THEN
                eta = ( MAX( 0.0_wp, ( lower_left_coord_y + ( ny + 1 ) * dy -   &
                                       cg%coord_y(jj)) ) /                      &
                        anterp_relax_length_n )**4
             ELSE
                eta = 999999.9_wp
             ENDIF
             fray(jj) = eta / ( 1.0_wp + eta )
          ENDDO
       ENDIF
      
       ALLOCATE( fraz(0:kctu) )
       DO  kk = 0, kctu
          zeta = ( ( zu(nzt) - cg%zu(kk) ) / anterp_relax_length_t )**4
          fraz(kk) = zeta / ( 1.0_wp + zeta )
       ENDDO

    END SUBROUTINE pmci_init_anterp_tophat



    SUBROUTINE pmci_init_tkefactor

!
!--    Computes the scaling factor for the SGS TKE from coarse grid to be used
!--    as BC for the fine grid. Based on the Kolmogorov energy spectrum
!--    for the inertial subrange and assumption of sharp cut-off of the resolved
!--    energy spectrum. Near the surface, the reduction of TKE is made
!--    smaller than further away from the surface.

       IMPLICIT NONE
       REAL(wp), PARAMETER ::  cfw = 0.2_wp          !:
       REAL(wp), PARAMETER ::  c_tkef = 0.6_wp       !:
       REAL(wp)            ::  fw                    !:
       REAL(wp), PARAMETER ::  fw0 = 0.9_wp          !:
       REAL(wp)            ::  glsf                  !:
       REAL(wp)            ::  glsc                  !:
       REAL(wp)            ::  height                !:
       REAL(wp), PARAMETER ::  p13 = 1.0_wp/3.0_wp   !:
       REAL(wp), PARAMETER ::  p23 = 2.0_wp/3.0_wp   !:
       INTEGER(iwp)        ::  k                     !:
       INTEGER(iwp)        ::  kc                    !:
        

       IF ( nest_bound_l )  THEN
          ALLOCATE( tkefactor_l(nzb:nzt+1,nysg:nyng) )
          tkefactor_l = 0.0_wp
          i = nxl - 1
          DO  j = nysg, nyng
             DO  k = nzb_s_inner(j,i) + 1, nzt
                kc     = kco(k+1)
                glsf   = ( dx * dy * dzu(k) )**p13
                glsc   = ( cg%dx * cg%dy *cg%dzu(kc) )**p13
                height = zu(k) - zu(nzb_s_inner(j,i))
                fw     = EXP( -cfw * height / glsf )
                tkefactor_l(k,j) = c_tkef * ( fw0 * fw + ( 1.0_wp - fw ) *      &
                                              ( glsf / glsc )**p23 )
             ENDDO
             tkefactor_l(nzb_s_inner(j,i),j) = c_tkef * fw0
          ENDDO
       ENDIF

       IF ( nest_bound_r )  THEN
          ALLOCATE( tkefactor_r(nzb:nzt+1,nysg:nyng) )
          tkefactor_r = 0.0_wp
          i = nxr + 1
          DO  j = nysg, nyng
             DO  k = nzb_s_inner(j,i) + 1, nzt
                kc     = kco(k+1)
                glsf   = ( dx * dy * dzu(k) )**p13
                glsc   = ( cg%dx * cg%dy * cg%dzu(kc) )**p13
                height = zu(k) - zu(nzb_s_inner(j,i))
                fw     = EXP( -cfw * height / glsf )
                tkefactor_r(k,j) = c_tkef * ( fw0 * fw + ( 1.0_wp - fw ) *      &
                                              ( glsf / glsc )**p23 )
             ENDDO
             tkefactor_r(nzb_s_inner(j,i),j) = c_tkef * fw0
          ENDDO
       ENDIF

      IF ( nest_bound_s )  THEN
          ALLOCATE( tkefactor_s(nzb:nzt+1,nxlg:nxrg) )
          tkefactor_s = 0.0_wp
          j = nys - 1
          DO  i = nxlg, nxrg
             DO  k = nzb_s_inner(j,i) + 1, nzt
                kc     = kco(k+1)
                glsf   = ( dx * dy * dzu(k) )**p13
                glsc   = ( cg%dx * cg%dy * cg%dzu(kc) ) ** p13
                height = zu(k) - zu(nzb_s_inner(j,i))
                fw     = EXP( -cfw*height / glsf )
                tkefactor_s(k,i) = c_tkef * ( fw0 * fw + ( 1.0_wp - fw ) *      &
                                              ( glsf / glsc )**p23 )
             ENDDO
             tkefactor_s(nzb_s_inner(j,i),i) = c_tkef * fw0
          ENDDO
       ENDIF

       IF ( nest_bound_n )  THEN
          ALLOCATE( tkefactor_n(nzb:nzt+1,nxlg:nxrg) )
          tkefactor_n = 0.0_wp
          j = nyn + 1
          DO  i = nxlg, nxrg
             DO  k = nzb_s_inner(j,i)+1, nzt
                kc     = kco(k+1)
                glsf   = ( dx * dy * dzu(k) )**p13
                glsc   = ( cg%dx * cg%dy * cg%dzu(kc) )**p13
                height = zu(k) - zu(nzb_s_inner(j,i))
                fw     = EXP( -cfw * height / glsf )
                tkefactor_n(k,i) = c_tkef * ( fw0 * fw + ( 1.0_wp - fw ) *      &
                                              ( glsf / glsc )**p23 )
             ENDDO
             tkefactor_n(nzb_s_inner(j,i),i) = c_tkef * fw0
          ENDDO
       ENDIF

       ALLOCATE( tkefactor_t(nysg:nyng,nxlg:nxrg) )
       k = nzt
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             kc     = kco(k+1)
             glsf   = ( dx * dy * dzu(k) )**p13
             glsc   = ( cg%dx * cg%dy * cg%dzu(kc) )**p13
             height = zu(k) - zu(nzb_s_inner(j,i))
             fw     = EXP( -cfw * height / glsf )
             tkefactor_t(j,i) = c_tkef * ( fw0 * fw + ( 1.0_wp - fw ) *         &
                                           ( glsf / glsc )**p23 )
          ENDDO
       ENDDO
      
    END SUBROUTINE pmci_init_tkefactor

#endif
 END SUBROUTINE pmci_setup_child



 SUBROUTINE pmci_setup_coordinates

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  i   !:
    INTEGER(iwp) ::  j   !:

!
!-- Create coordinate arrays.
    ALLOCATE( coord_x(-nbgp:nx+nbgp) )
    ALLOCATE( coord_y(-nbgp:ny+nbgp) )
     
    DO  i = -nbgp, nx + nbgp
       coord_x(i) = lower_left_coord_x + i * dx
    ENDDO
     
    DO  j = -nbgp, ny + nbgp
       coord_y(j) = lower_left_coord_y + j * dy
    ENDDO

#endif
 END SUBROUTINE pmci_setup_coordinates




 SUBROUTINE pmci_set_array_pointer( name, child_id, nz_cl )

    IMPLICIT NONE

    INTEGER, INTENT(IN)          ::  child_id    !:
    INTEGER, INTENT(IN)          ::  nz_cl       !:
    CHARACTER(LEN=*), INTENT(IN) ::  name        !:

#if defined( __parallel )
    INTEGER(iwp) ::  ierr        !:
    INTEGER(iwp) ::  istat       !:

    REAL(wp), POINTER, DIMENSION(:,:)   ::  p_2d        !:
    REAL(wp), POINTER, DIMENSION(:,:)   ::  p_2d_sec    !:
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d        !:
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d_sec    !:


    NULLIFY( p_3d )
    NULLIFY( p_2d )

!
!-- List of array names, which can be coupled.
!-- In case of 3D please change also the second array for the pointer version
    IF ( TRIM(name) == "u"  )  p_3d => u
    IF ( TRIM(name) == "v"  )  p_3d => v
    IF ( TRIM(name) == "w"  )  p_3d => w
    IF ( TRIM(name) == "e"  )  p_3d => e
    IF ( TRIM(name) == "pt" )  p_3d => pt
    IF ( TRIM(name) == "q"  )  p_3d => q
    IF ( TRIM(name) == "s"  )  p_3d => s
!
!-- Next line is just an example for a 2D array (not active for coupling!)
!-- Please note, that z0 has to be declared as TARGET array in modules.f90
!    IF ( TRIM(name) == "z0" )    p_2d => z0

#if defined( __nopointer )
    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_3d, nz_cl, nz )
    ELSEIF ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_2d )
    ELSE
!
!--    Give only one message for the root domain
       IF ( myid == 0  .AND.  cpl_id == 1 )  THEN

          message_string = 'pointer for array "' // TRIM( name ) //             &
                           '" can''t be associated'
          CALL message( 'pmci_set_array_pointer', 'PA0117', 3, 2, 0, 6, 0 )
       ELSE
!
!--       Avoid others to continue
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF
    ENDIF
#else
    IF ( TRIM(name) == "u"  )  p_3d_sec => u_2
    IF ( TRIM(name) == "v"  )  p_3d_sec => v_2
    IF ( TRIM(name) == "w"  )  p_3d_sec => w_2
    IF ( TRIM(name) == "e"  )  p_3d_sec => e_2
    IF ( TRIM(name) == "pt" )  p_3d_sec => pt_2
    IF ( TRIM(name) == "q"  )  p_3d_sec => q_2
    IF ( TRIM(name) == "s"  )  p_3d_sec => s_2

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_3d, nz_cl, nz,                     &
                                 array_2 = p_3d_sec )
    ELSEIF ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_2d )
    ELSE
!
!--    Give only one message for the root domain
       IF ( myid == 0  .AND.  cpl_id == 1 )  THEN

          message_string = 'pointer for array "' // TRIM( name ) //             &
                           '" can''t be associated'
          CALL message( 'pmci_set_array_pointer', 'PA0117', 3, 2, 0, 6, 0 )
       ELSE
!
!--       Avoid others to continue
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF

   ENDIF
#endif

#endif
 END SUBROUTINE pmci_set_array_pointer



 SUBROUTINE pmci_create_child_arrays( name, is, ie, js, je, nzc  )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name    !:

    INTEGER(iwp), INTENT(IN) ::  ie      !:
    INTEGER(iwp), INTENT(IN) ::  is      !:
    INTEGER(iwp), INTENT(IN) ::  je      !:
    INTEGER(iwp), INTENT(IN) ::  js      !:
    INTEGER(iwp), INTENT(IN) ::  nzc     !:  Note that nzc is cg%nz

#if defined( __parallel )
    INTEGER(iwp) ::  ierr    !:
    INTEGER(iwp) ::  istat   !:

    REAL(wp), POINTER,DIMENSION(:,:)   ::  p_2d    !:
    REAL(wp), POINTER,DIMENSION(:,:,:) ::  p_3d    !:


    NULLIFY( p_3d )
    NULLIFY( p_2d )

!
!-- List of array names, which can be coupled
    IF ( TRIM( name ) == "u" )  THEN
       IF ( .NOT. ALLOCATED( uc ) )  ALLOCATE( uc(0:nzc+1, js:je, is:ie) )
       p_3d => uc
    ELSEIF ( TRIM( name ) == "v" )  THEN
       IF ( .NOT. ALLOCATED( vc ) )  ALLOCATE( vc(0:nzc+1, js:je, is:ie) )
       p_3d => vc
    ELSEIF ( TRIM( name ) == "w" )  THEN
       IF ( .NOT. ALLOCATED( wc ) )  ALLOCATE( wc(0:nzc+1, js:je, is:ie) )
       p_3d => wc
    ELSEIF ( TRIM( name ) == "e" )  THEN
       IF ( .NOT. ALLOCATED( ec ) )  ALLOCATE( ec(0:nzc+1, js:je, is:ie) )
       p_3d => ec
    ELSEIF ( TRIM( name ) == "pt")  THEN
       IF ( .NOT. ALLOCATED( ptc ) ) ALLOCATE( ptc(0:nzc+1, js:je, is:ie) )
       p_3d => ptc
    ELSEIF ( TRIM( name ) == "q")  THEN
       IF ( .NOT. ALLOCATED( qc ) ) ALLOCATE( qc(0:nzc+1, js:je, is:ie) )
       p_3d => qc
    ELSEIF ( TRIM( name ) == "s")  THEN
       IF ( .NOT. ALLOCATED( sc ) ) ALLOCATE( sc(0:nzc+1, js:je, is:ie) )
       p_3d => sc
    !ELSEIF (trim(name) == "z0") then
       !IF (.not.allocated(z0c))  allocate(z0c(js:je, is:ie))
       !p_2d => z0c
    ENDIF

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_c_set_dataarray( p_3d )
    ELSEIF ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_c_set_dataarray( p_2d )
    ELSE
!
!--    Give only one message for the first child domain
       IF ( myid == 0  .AND.  cpl_id == 2 )  THEN

          message_string = 'pointer for array "' // TRIM( name ) //             &
                           '" can''t be associated'
          CALL message( 'pmci_create_child_arrays', 'PA0170', 3, 2, 0, 6, 0 )
       ELSE
!
!--       Prevent others from continuing
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF
    ENDIF

#endif
 END SUBROUTINE pmci_create_child_arrays



 SUBROUTINE pmci_parent_initialize

!
!-- Send data for the children in order to let them create initial 
!-- conditions by interpolating the parent-domain fields.
#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  child_id    !:
    INTEGER(iwp) ::  m           !:

    REAL(wp) ::  waittime        !:


    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_fillbuffer( child_id, waittime=waittime )
    ENDDO

#endif
 END SUBROUTINE pmci_parent_initialize



 SUBROUTINE pmci_child_initialize

!
!-- Create initial conditions for the current child domain by interpolating 
!-- the parent-domain fields.
#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  i          !:
    INTEGER(iwp) ::  icl        !:
    INTEGER(iwp) ::  icr        !:
    INTEGER(iwp) ::  j          !:
    INTEGER(iwp) ::  jcn        !:
    INTEGER(iwp) ::  jcs        !:

    REAL(wp) ::  waittime       !:

!
!-- Root id is never a child
    IF ( cpl_id > 1 )  THEN

!
!--    Child domain boundaries in the parent index space
       icl = coarse_bound(1)
       icr = coarse_bound(2)
       jcs = coarse_bound(3)
       jcn = coarse_bound(4)

!
!--    Get data from the parent
       CALL pmc_c_getbuffer( waittime = waittime )

!
!--    The interpolation.
       CALL pmci_interp_tril_all ( u,  uc,  icu, jco, kco, r1xu, r2xu, r1yo,    &
                                   r2yo, r1zo, r2zo, nzb_u_inner, 'u' )
       CALL pmci_interp_tril_all ( v,  vc,  ico, jcv, kco, r1xo, r2xo, r1yv,    &
                                   r2yv, r1zo, r2zo, nzb_v_inner, 'v' )
       CALL pmci_interp_tril_all ( w,  wc,  ico, jco, kcw, r1xo, r2xo, r1yo,    &
                                   r2yo, r1zw, r2zw, nzb_w_inner, 'w' )
       CALL pmci_interp_tril_all ( e,  ec,  ico, jco, kco, r1xo, r2xo, r1yo,    &
                                   r2yo, r1zo, r2zo, nzb_s_inner, 'e' )
       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_tril_all ( pt, ptc, ico, jco, kco, r1xo, r2xo,       &
                                      r1yo, r2yo, r1zo, r2zo, nzb_s_inner, 's' )
       ENDIF
       IF ( humidity )  THEN
          CALL pmci_interp_tril_all ( q, qc, ico, jco, kco, r1xo, r2xo, r1yo,   &
                                      r2yo, r1zo, r2zo, nzb_s_inner, 's' )
       ENDIF
       IF ( passive_scalar )  THEN
          CALL pmci_interp_tril_all ( s, sc, ico, jco, kco, r1xo, r2xo, r1yo,   &
                                      r2yo, r1zo, r2zo, nzb_s_inner, 's' )
       ENDIF

       IF ( topography /= 'flat' )  THEN
!
!--       Inside buildings set velocities and TKE back to zero.
!--       Other scalars (pt, q, s, km, kh, p, sa, ...) are ignored at present,
!--       maybe revise later.
          DO   i = nxlg, nxrg
             DO   j = nysg, nyng
                u(nzb:nzb_u_inner(j,i),j,i)   = 0.0_wp
                v(nzb:nzb_v_inner(j,i),j,i)   = 0.0_wp
                w(nzb:nzb_w_inner(j,i),j,i)   = 0.0_wp
                e(nzb:nzb_s_inner(j,i),j,i)   = 0.0_wp
                u_p(nzb:nzb_u_inner(j,i),j,i) = 0.0_wp
                v_p(nzb:nzb_v_inner(j,i),j,i) = 0.0_wp
                w_p(nzb:nzb_w_inner(j,i),j,i) = 0.0_wp
                e_p(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDIF
    ENDIF


 CONTAINS


    SUBROUTINE pmci_interp_tril_all( f, fc, ic, jc, kc, r1x, r2x, r1y, r2y,     &
                                     r1z, r2z, kb, var )
!
!--    Interpolation of the internal values for the child-domain initialization
!--    This subroutine is based on trilinear interpolation.
!--    Coding based on interp_tril_lr/sn/t
       IMPLICIT NONE

       CHARACTER(LEN=1), INTENT(IN) :: var  !:

       INTEGER(iwp), DIMENSION(nxlg:nxrg), INTENT(IN)           ::  ic    !:
       INTEGER(iwp), DIMENSION(nysg:nyng), INTENT(IN)           ::  jc    !:
       INTEGER(iwp), DIMENSION(nzb:nzt+1), INTENT(IN)           ::  kc    !:
       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  kb    !:

       INTEGER(iwp) ::  i      !:
       INTEGER(iwp) ::  ib     !:
       INTEGER(iwp) ::  ie     !:
       INTEGER(iwp) ::  j      !:
       INTEGER(iwp) ::  jb     !:
       INTEGER(iwp) ::  je     !:
       INTEGER(iwp) ::  k      !:
       INTEGER(iwp) ::  k1     !:
       INTEGER(iwp) ::  kbc    !:
       INTEGER(iwp) ::  l      !:
       INTEGER(iwp) ::  m      !:
       INTEGER(iwp) ::  n      !:

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) :: f !:
       REAL(wp), DIMENSION(0:cg%nz+1,jcs:jcn,icl:icr), INTENT(IN) :: fc       !:
       REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN) :: r1x   !:
       REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN) :: r2x   !:
       REAL(wp), DIMENSION(nysg:nyng), INTENT(IN) :: r1y   !:
       REAL(wp), DIMENSION(nysg:nyng), INTENT(IN) :: r2y   !:
       REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN) :: r1z   !:
       REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN) :: r2z   !:

       REAL(wp) ::  fk         !:
       REAL(wp) ::  fkj        !:
       REAL(wp) ::  fkjp       !:
       REAL(wp) ::  fkp        !:
       REAL(wp) ::  fkpj       !:
       REAL(wp) ::  fkpjp      !:
       REAL(wp) ::  logratio   !:
       REAL(wp) ::  logzuc1    !:
       REAL(wp) ::  zuc1       !:


       ib = nxl
       ie = nxr
       jb = nys
       je = nyn
       IF ( nesting_mode /= 'vertical' )  THEN
          IF ( nest_bound_l )  THEN
             ib = nxl - 1
!
!--          For u, nxl is a ghost node, but not for the other variables
             IF ( var == 'u' )  THEN
                ib = nxl
             ENDIF
          ENDIF
          IF ( nest_bound_s )  THEN
             jb = nys - 1
!
!--          For v, nys is a ghost node, but not for the other variables
             IF ( var == 'v' )  THEN
                jb = nys
             ENDIF
          ENDIF
          IF ( nest_bound_r )  THEN
             ie = nxr + 1
          ENDIF
          IF ( nest_bound_n )  THEN
             je = nyn + 1
          ENDIF
       ENDIF
!
!--    Trilinear interpolation.
       DO  i = ib, ie
          DO  j = jb, je
             DO  k = kb(j,i), nzt + 1
                l = ic(i)
                m = jc(j)
                n = kc(k)
                fkj      = r1x(i) * fc(n,m,l)     + r2x(i) * fc(n,m,l+1)
                fkjp     = r1x(i) * fc(n,m+1,l)   + r2x(i) * fc(n,m+1,l+1)
                fkpj     = r1x(i) * fc(n+1,m,l)   + r2x(i) * fc(n+1,m,l+1)
                fkpjp    = r1x(i) * fc(n+1,m+1,l) + r2x(i) * fc(n+1,m+1,l+1)
                fk       = r1y(j) * fkj  + r2y(j) * fkjp
                fkp      = r1y(j) * fkpj + r2y(j) * fkpjp
                f(k,j,i) = r1z(k) * fk   + r2z(k) * fkp
             ENDDO
          ENDDO
       ENDDO

!
!--    Correct the interpolated values of u and v in near-wall nodes, i.e. in
!--    the nodes below the coarse-grid nodes with k=1. The corrction is only
!--    made over horizontal wall surfaces in this phase. For the nest boundary
!--    conditions, a corresponding correction is made for all vertical walls,
!--    too.
       IF ( var == 'u' .OR. var == 'v' )  THEN
          DO  i = ib, nxr
             DO  j = jb, nyn
                kbc = 1
!
!--             kbc is the first coarse-grid point above the surface
                DO  WHILE ( cg%zu(kbc) < zu(kb(j,i)) )
                   kbc = kbc + 1
                ENDDO
                zuc1 = cg%zu(kbc)
                k1   = kb(j,i) + 1
                DO  WHILE ( zu(k1) < zuc1 )
                   k1 = k1 + 1
                ENDDO
                logzuc1 = LOG( ( zu(k1) - zu(kb(j,i)) ) / z0(j,i) )

                k = kb(j,i) + 1
                DO  WHILE ( zu(k) < zuc1 )
                   logratio = ( LOG( ( zu(k) - zu(kb(j,i)) ) / z0(j,i)) ) /     &
                                logzuc1
                   f(k,j,i) = logratio * f(k1,j,i)
                   k  = k + 1
                ENDDO
                f(kb(j,i),j,i) = 0.0_wp
             ENDDO
          ENDDO

       ELSEIF ( var == 'w' )  THEN

          DO  i = ib, nxr
              DO  j = jb, nyn
                f(kb(j,i),j,i) = 0.0_wp
             ENDDO
          ENDDO

       ENDIF

    END SUBROUTINE pmci_interp_tril_all

#endif
 END SUBROUTINE pmci_child_initialize



 SUBROUTINE pmci_check_setting_mismatches
!
!-- Check for mismatches between settings of master and child variables
!-- (e.g., all children have to follow the end_time settings of the root model).
!-- The root model overwrites variables in the other models, so these variables
!-- only need to be set once in file PARIN.

#if defined( __parallel )

    USE control_parameters,                                                     &
        ONLY:  dt_restart, end_time, message_string, restart_time, time_restart

    IMPLICIT NONE

    INTEGER ::  ierr

    REAL(wp) ::  dt_restart_root
    REAL(wp) ::  end_time_root
    REAL(wp) ::  restart_time_root
    REAL(wp) ::  time_restart_root

!
!-- Check the time to be simulated.
!-- Here, and in the following, the root process communicates the respective
!-- variable to all others, and its value will then be compared with the local
!-- values.
    IF ( pmc_is_rootmodel() )  end_time_root = end_time
    CALL MPI_BCAST( end_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( end_time /= end_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',       &
               'child settings &   end_time(root) = ', end_time_root,           &
               ' &   end_time(child) = ', end_time, ' & child value is set',    &
               ' to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6,  &
                        0 )
          end_time = end_time_root
       ENDIF
    ENDIF

!
!-- Same for restart time
    IF ( pmc_is_rootmodel() )  restart_time_root = restart_time
    CALL MPI_BCAST( restart_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( restart_time /= restart_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',       &
               'child settings &   restart_time(root) = ', restart_time_root,   &
               ' &   restart_time(child) = ', restart_time, ' & child ',        &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6,  &
                        0 )
          restart_time = restart_time_root
       ENDIF
    ENDIF

!
!-- Same for dt_restart
    IF ( pmc_is_rootmodel() )  dt_restart_root = dt_restart
    CALL MPI_BCAST( dt_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( dt_restart /= dt_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',       &
               'child settings &   dt_restart(root) = ', dt_restart_root,       &
               ' &   dt_restart(child) = ', dt_restart, ' & child ',            &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6,  &
                        0 )
          dt_restart = dt_restart_root
       ENDIF
    ENDIF

!
!-- Same for time_restart
    IF ( pmc_is_rootmodel() )  time_restart_root = time_restart
    CALL MPI_BCAST( time_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( time_restart /= time_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',       &
               'child settings &   time_restart(root) = ', time_restart_root,   &
               ' &   time_restart(child) = ', time_restart, ' & child ',        &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6,  &
                        0 )
          time_restart = time_restart_root
       ENDIF
    ENDIF

#endif

 END SUBROUTINE pmci_check_setting_mismatches



 SUBROUTINE pmci_ensure_nest_mass_conservation

!
!-- Adjust the volume-flow rate through the top boundary so that the net volume
!-- flow through all boundaries of the current nest domain becomes zero.
    IMPLICIT NONE

    INTEGER(iwp) ::  i                           !:
    INTEGER(iwp) ::  ierr                        !:
    INTEGER(iwp) ::  j                           !:
    INTEGER(iwp) ::  k                           !:

    REAL(wp) ::  dxdy                            !:
    REAL(wp) ::  innor                           !:
    REAL(wp) ::  w_lt                            !:
    REAL(wp), DIMENSION(1:3) ::  volume_flow_l   !:

!
!-- Sum up the volume flow through the left/right boundaries
    volume_flow(1)   = 0.0_wp
    volume_flow_l(1) = 0.0_wp

    IF ( nest_bound_l )  THEN
       i = 0
       innor = dy
       DO   j = nys, nyn
          DO   k = nzb_u_inner(j,i)+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + innor * u(k,j,i) * dzw(k)
          ENDDO
       ENDDO
    ENDIF

    IF ( nest_bound_r )  THEN
       i = nx + 1
       innor = -dy
       DO   j = nys, nyn
          DO   k = nzb_u_inner(j,i)+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + innor * u(k,j,i) * dzw(k)
          ENDDO
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flow_l(1), volume_flow(1), 1, MPI_REAL,          &
                        MPI_SUM, comm2d, ierr )
#else
    volume_flow(1) = volume_flow_l(1)
#endif

!
!-- Sum up the volume flow through the south/north boundaries
    volume_flow(2)   = 0.0_wp
    volume_flow_l(2) = 0.0_wp

    IF ( nest_bound_s )  THEN
       j = 0
       innor = dx
       DO   i = nxl, nxr
          DO   k = nzb_v_inner(j,i)+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + innor * v(k,j,i) * dzw(k)
          ENDDO
       ENDDO
    ENDIF

    IF ( nest_bound_n )  THEN
       j = ny + 1
       innor = -dx
       DO   i = nxl, nxr
          DO   k = nzb_v_inner(j,i)+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + innor * v(k,j,i) * dzw(k)
          ENDDO
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flow_l(2), volume_flow(2), 1, MPI_REAL,          &
                        MPI_SUM, comm2d, ierr )
#else
    volume_flow(2) = volume_flow_l(2)
#endif

!
!-- Sum up the volume flow through the top boundary
    volume_flow(3)   = 0.0_wp
    volume_flow_l(3) = 0.0_wp
    dxdy = dx * dy
    k = nzt
    DO   i = nxl, nxr
       DO   j = nys, nyn
          volume_flow_l(3) = volume_flow_l(3) - w(k,j,i) * dxdy
       ENDDO
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flow_l(3), volume_flow(3), 1, MPI_REAL,          &
                        MPI_SUM, comm2d, ierr )
#else
    volume_flow(3) = volume_flow_l(3)
#endif

!
!-- Correct the top-boundary value of w
    w_lt = (volume_flow(1) + volume_flow(2) + volume_flow(3)) / area_t
    DO   i = nxl, nxr
       DO   j = nys, nyn
          DO  k = nzt, nzt + 1
             w(k,j,i) = w(k,j,i) + w_lt
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE pmci_ensure_nest_mass_conservation



 SUBROUTINE pmci_synchronize

#if defined( __parallel )
!
!-- Unify the time steps for each model and synchronize using 
!-- MPI_ALLREDUCE with the MPI_MIN operator over all processes using 
!-- the global communicator MPI_COMM_WORLD.
   IMPLICIT NONE

   INTEGER(iwp)           :: ierr  !:
   REAL(wp), DIMENSION(1) :: dtl   !:
   REAL(wp), DIMENSION(1) :: dtg   !:

   
   dtl(1) = dt_3d 
   CALL MPI_ALLREDUCE( dtl, dtg, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr )
   dt_3d  = dtg(1)

#endif
 END SUBROUTINE pmci_synchronize
                


 SUBROUTINE pmci_set_swaplevel( swaplevel )
!
!-- After each Runge-Kutta sub-timestep, alternately set buffer one or buffer
!-- two active

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  swaplevel  !: swaplevel (1 or 2) of PALM's
                                            !: timestep

    INTEGER(iwp)            ::  child_id    !:
    INTEGER(iwp)            ::  m           !:

    DO  m = 1, SIZE( pmc_parent_for_child )-1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_set_active_data_array( child_id, swaplevel )
    ENDDO

 END SUBROUTINE pmci_set_swaplevel



 SUBROUTINE pmci_datatrans( local_nesting_mode )
!
!-- This subroutine controls the nesting according to the nestpar 
!-- parameter nesting_mode (two-way (default) or one-way) and the 
!-- order of anterpolations according to the nestpar parameter 
!-- nesting_datatransfer_mode (cascade, overlap or mixed (default)).
!-- Although nesting_mode is a variable of this model, pass it as 
!-- an argument to allow for example to force one-way initialization 
!-- phase.

    IMPLICIT NONE

    INTEGER(iwp)           ::  ierr   !:
    INTEGER(iwp)           ::  istat  !:

    CHARACTER(LEN=*), INTENT(IN) ::  local_nesting_mode

    IF ( TRIM( local_nesting_mode ) == 'one-way' )  THEN

       CALL pmci_child_datatrans( parent_to_child )
       CALL pmci_parent_datatrans( parent_to_child )

    ELSE

       IF( nesting_datatransfer_mode == 'cascade' )  THEN

          CALL pmci_child_datatrans( parent_to_child )
          CALL pmci_parent_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ELSEIF( nesting_datatransfer_mode == 'overlap')  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_child_datatrans( child_to_parent )
          CALL pmci_parent_datatrans( child_to_parent )

       ELSEIF( TRIM( nesting_datatransfer_mode ) == 'mixed' )  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ENDIF

    ENDIF

 END SUBROUTINE pmci_datatrans




 SUBROUTINE pmci_parent_datatrans( direction )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction   !:

#if defined( __parallel )
    INTEGER(iwp) ::  child_id    !:
    INTEGER(iwp) ::  i           !:
    INTEGER(iwp) ::  j           !:
    INTEGER(iwp) ::  ierr        !:
    INTEGER(iwp) ::  m           !:

    REAL(wp)               ::  waittime    !:
    REAL(wp), DIMENSION(1) ::  dtc         !:
    REAL(wp), DIMENSION(1) ::  dtl         !:


    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       
       IF ( direction == parent_to_child )  THEN
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'start' )
          CALL pmc_s_fillbuffer( child_id )
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'stop' )
       ELSE
!
!--       Communication from child to parent
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'start' )
          child_id = pmc_parent_for_child(m)
          CALL pmc_s_getdata_from_buffer( child_id )
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'stop' )

!
!--       The anterpolated data is now available in u etc
          IF ( topography /= 'flat' )  THEN

!
!--          Inside buildings/topography reset velocities and TKE back to zero.
!--          Other scalars (pt, q, s, km, kh, p, sa, ...) are ignored at
!--          present, maybe revise later.
             DO   i = nxlg, nxrg
                DO   j = nysg, nyng
                   u(nzb:nzb_u_inner(j,i),j,i)  = 0.0_wp
                   v(nzb:nzb_v_inner(j,i),j,i)  = 0.0_wp
                   w(nzb:nzb_w_inner(j,i),j,i)  = 0.0_wp
                   e(nzb:nzb_s_inner(j,i),j,i)  = 0.0_wp
!
!--                TO_DO: zero setting of temperature within topography creates
!--                       wrong results
!                   pt(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   IF ( humidity  .OR.  passive_scalar )  THEN
!                      q(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO

#endif
 END SUBROUTINE pmci_parent_datatrans



 SUBROUTINE pmci_child_datatrans( direction )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction   !:

#if defined( __parallel )
    INTEGER(iwp) ::  ierr        !:
    INTEGER(iwp) ::  icl         !:
    INTEGER(iwp) ::  icr         !:
    INTEGER(iwp) ::  jcs         !:
    INTEGER(iwp) ::  jcn         !:
    
    REAL(wp), DIMENSION(1) ::  dtl         !:
    REAL(wp), DIMENSION(1) ::  dts         !:


    dtl = dt_3d
    IF ( cpl_id > 1 )  THEN
!
!--    Child domain boundaries in the parent indice space.
       icl = coarse_bound(1)
       icr = coarse_bound(2)
       jcs = coarse_bound(3)
       jcn = coarse_bound(4)

       IF ( direction == parent_to_child )  THEN

          CALL cpu_log( log_point_s(73), 'pmc child recv', 'start' )
          CALL pmc_c_getbuffer( )
          CALL cpu_log( log_point_s(73), 'pmc child recv', 'stop' )

          CALL cpu_log( log_point_s(75), 'pmc interpolation', 'start' )
          CALL pmci_interpolation
          CALL cpu_log( log_point_s(75), 'pmc interpolation', 'stop' )

       ELSE
!
!--       direction == child_to_parent
          CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'start' )
          CALL pmci_anterpolation
          CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'stop' )

          CALL cpu_log( log_point_s(74), 'pmc child send', 'start' )
          CALL pmc_c_putbuffer( )
          CALL cpu_log( log_point_s(74), 'pmc child send', 'stop' )

       ENDIF
    ENDIF

 CONTAINS

    SUBROUTINE pmci_interpolation

!
!--    A wrapper routine for all interpolation and extrapolation actions
       IMPLICIT NONE

!
!--    In case of vertical nesting no interpolation is needed for the
!--    horizontal boundaries
       IF ( nesting_mode /= 'vertical' )  THEN
       
!
!--    Left border pe:
          IF ( nest_bound_l )  THEN
             CALL pmci_interp_tril_lr( u,  uc,  icu, jco, kco, r1xu, r2xu,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_u_inner,     &
                                       logc_u_l, logc_ratio_u_l,                &
                                       nzt_topo_nestbc_l, 'l', 'u' )
             CALL pmci_interp_tril_lr( v,  vc,  ico, jcv, kco, r1xo, r2xo,      &
                                       r1yv, r2yv, r1zo, r2zo, nzb_v_inner,     &
                                       logc_v_l, logc_ratio_v_l,                &
                                       nzt_topo_nestbc_l, 'l', 'v' )
             CALL pmci_interp_tril_lr( w,  wc,  ico, jco, kcw, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zw, r2zw, nzb_w_inner,     &
                                       logc_w_l, logc_ratio_w_l,                &
                                       nzt_topo_nestbc_l, 'l', 'w' )
             CALL pmci_interp_tril_lr( e,  ec,  ico, jco, kco, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_s_inner,     &
                                       logc_u_l, logc_ratio_u_l,                &
                                       nzt_topo_nestbc_l, 'l', 'e' )
             IF ( .NOT. neutral )  THEN
                CALL pmci_interp_tril_lr( pt, ptc, ico, jco, kco, r1xo, r2xo,   &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_l, logc_ratio_u_l,             &
                                          nzt_topo_nestbc_l, 'l', 's' )
             ENDIF
             IF ( humidity )  THEN
                CALL pmci_interp_tril_lr( q, qc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_l, logc_ratio_u_l,             &
                                          nzt_topo_nestbc_l, 'l', 's' )
             ENDIF
             IF ( passive_scalar )  THEN
                CALL pmci_interp_tril_lr( s, sc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_l, logc_ratio_u_l,             &
                                          nzt_topo_nestbc_l, 'l', 's' )
             ENDIF

             IF ( TRIM( nesting_mode ) == 'one-way' )  THEN
                CALL pmci_extrap_ifoutflow_lr( u, nzb_u_inner, 'l', 'u' )
                CALL pmci_extrap_ifoutflow_lr( v, nzb_v_inner, 'l', 'v' )
                CALL pmci_extrap_ifoutflow_lr( w, nzb_w_inner, 'l', 'w' )
                CALL pmci_extrap_ifoutflow_lr( e, nzb_s_inner, 'l', 'e' )
                IF ( .NOT. neutral )  THEN
                   CALL pmci_extrap_ifoutflow_lr( pt,nzb_s_inner, 'l', 's' )
                ENDIF
                IF ( humidity )  THEN
                   CALL pmci_extrap_ifoutflow_lr( q, nzb_s_inner, 'l', 's' )
                ENDIF
                IF ( passive_scalar )  THEN
                   CALL pmci_extrap_ifoutflow_lr( s, nzb_s_inner, 'l', 's' )
                ENDIF
             ENDIF

          ENDIF

   !
   !--    Right border pe
          IF ( nest_bound_r )  THEN
             CALL pmci_interp_tril_lr( u,  uc,  icu, jco, kco, r1xu, r2xu,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_u_inner,     &
                                       logc_u_r, logc_ratio_u_r,                &
                                       nzt_topo_nestbc_r, 'r', 'u' )
             CALL pmci_interp_tril_lr( v,  vc,  ico, jcv, kco, r1xo, r2xo,      &
                                       r1yv, r2yv, r1zo, r2zo, nzb_v_inner,     &
                                       logc_v_r, logc_ratio_v_r,                &
                                       nzt_topo_nestbc_r, 'r', 'v' )
             CALL pmci_interp_tril_lr( w,  wc,  ico, jco, kcw, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zw, r2zw, nzb_w_inner,     &
                                       logc_w_r, logc_ratio_w_r,                &
                                       nzt_topo_nestbc_r, 'r', 'w' )
             CALL pmci_interp_tril_lr( e,  ec,  ico, jco, kco, r1xo, r2xo,      &
                                       r1yo,r2yo, r1zo, r2zo, nzb_s_inner,      &
                                       logc_u_r, logc_ratio_u_r,                &
                                       nzt_topo_nestbc_r, 'r', 'e' )
             IF ( .NOT. neutral )  THEN
                CALL pmci_interp_tril_lr( pt, ptc, ico, jco, kco, r1xo, r2xo,   &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_r, logc_ratio_u_r,             &
                                          nzt_topo_nestbc_r, 'r', 's' )
             ENDIF
             IF ( humidity )  THEN
                CALL pmci_interp_tril_lr( q, qc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_r, logc_ratio_u_r,             &
                                          nzt_topo_nestbc_r, 'r', 's' )
             ENDIF
             IF ( passive_scalar )  THEN
                CALL pmci_interp_tril_lr( s, sc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_r, logc_ratio_u_r,             &
                                          nzt_topo_nestbc_r, 'r', 's' )
             ENDIF

             IF ( TRIM( nesting_mode ) == 'one-way' )  THEN
                CALL pmci_extrap_ifoutflow_lr( u, nzb_u_inner, 'r', 'u' )
                CALL pmci_extrap_ifoutflow_lr( v, nzb_v_inner, 'r', 'v' )
                CALL pmci_extrap_ifoutflow_lr( w, nzb_w_inner, 'r', 'w' )
                CALL pmci_extrap_ifoutflow_lr( e, nzb_s_inner, 'r', 'e' )
                IF ( .NOT. neutral )  THEN
                   CALL pmci_extrap_ifoutflow_lr( pt,nzb_s_inner, 'r', 's' )
                ENDIF
                IF ( humidity )  THEN
                   CALL pmci_extrap_ifoutflow_lr( q, nzb_s_inner, 'r', 's' )
                ENDIF
                IF ( passive_scalar )  THEN
                   CALL pmci_extrap_ifoutflow_lr( s, nzb_s_inner, 'r', 's' )
                ENDIF
             ENDIF

          ENDIF

   !
   !--    South border pe
          IF ( nest_bound_s )  THEN
             CALL pmci_interp_tril_sn( u,  uc,  icu, jco, kco, r1xu, r2xu,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_u_inner,     &
                                       logc_u_s, logc_ratio_u_s,                &
                                       nzt_topo_nestbc_s, 's', 'u' )
             CALL pmci_interp_tril_sn( v,  vc,  ico, jcv, kco, r1xo, r2xo,      &
                                       r1yv, r2yv, r1zo, r2zo, nzb_v_inner,     &
                                       logc_v_s, logc_ratio_v_s,                &
                                       nzt_topo_nestbc_s, 's', 'v' )
             CALL pmci_interp_tril_sn( w,  wc,  ico, jco, kcw, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zw, r2zw, nzb_w_inner,     &
                                       logc_w_s, logc_ratio_w_s,                &
                                       nzt_topo_nestbc_s, 's','w' )
             CALL pmci_interp_tril_sn( e,  ec,  ico, jco, kco, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_s_inner,     &
                                       logc_u_s, logc_ratio_u_s,                &
                                       nzt_topo_nestbc_s, 's', 'e' )
             IF ( .NOT. neutral )  THEN
                CALL pmci_interp_tril_sn( pt, ptc, ico, jco, kco, r1xo, r2xo,   &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_s, logc_ratio_u_s,             &
                                          nzt_topo_nestbc_s, 's', 's' )
             ENDIF
             IF ( humidity )  THEN
                CALL pmci_interp_tril_sn( q, qc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo,r2yo, r1zo, r2zo, nzb_s_inner,   &
                                          logc_u_s, logc_ratio_u_s,             &
                                          nzt_topo_nestbc_s, 's', 's' )
             ENDIF
             IF ( passive_scalar )  THEN
                CALL pmci_interp_tril_sn( s, sc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo,r2yo, r1zo, r2zo, nzb_s_inner,   &
                                          logc_u_s, logc_ratio_u_s,             &
                                          nzt_topo_nestbc_s, 's', 's' )
             ENDIF

             IF ( TRIM( nesting_mode ) == 'one-way' )  THEN
                CALL pmci_extrap_ifoutflow_sn( u, nzb_u_inner, 's', 'u' )
                CALL pmci_extrap_ifoutflow_sn( v, nzb_v_inner, 's', 'v' )
                CALL pmci_extrap_ifoutflow_sn( w, nzb_w_inner, 's', 'w' )
                CALL pmci_extrap_ifoutflow_sn( e, nzb_s_inner, 's', 'e' )
                IF ( .NOT. neutral )  THEN
                   CALL pmci_extrap_ifoutflow_sn( pt,nzb_s_inner, 's', 's' )
                ENDIF
                IF ( humidity )  THEN
                   CALL pmci_extrap_ifoutflow_sn( q, nzb_s_inner, 's', 's' )
                ENDIF
                IF ( passive_scalar )  THEN
                   CALL pmci_extrap_ifoutflow_sn( s, nzb_s_inner, 's', 's' )
                ENDIF
             ENDIF

          ENDIF

   !
   !--    North border pe
          IF ( nest_bound_n )  THEN
             CALL pmci_interp_tril_sn( u,  uc,  icu, jco, kco, r1xu, r2xu,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_u_inner,     &
                                       logc_u_n, logc_ratio_u_n,                &
                                       nzt_topo_nestbc_n, 'n', 'u' )
             CALL pmci_interp_tril_sn( v,  vc,  ico, jcv, kco, r1xo, r2xo,      &
                                       r1yv, r2yv, r1zo, r2zo, nzb_v_inner,     &
                                       logc_v_n, logc_ratio_v_n,                & 
                                       nzt_topo_nestbc_n, 'n', 'v' )
             CALL pmci_interp_tril_sn( w,  wc,  ico, jco, kcw, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zw, r2zw, nzb_w_inner,     &
                                       logc_w_n, logc_ratio_w_n,                &
                                       nzt_topo_nestbc_n, 'n', 'w' )
             CALL pmci_interp_tril_sn( e,  ec,  ico, jco, kco, r1xo, r2xo,      &
                                       r1yo, r2yo, r1zo, r2zo, nzb_s_inner,     &
                                       logc_u_n, logc_ratio_u_n,                &
                                       nzt_topo_nestbc_n, 'n', 'e' )
             IF ( .NOT. neutral )  THEN
                CALL pmci_interp_tril_sn( pt, ptc, ico, jco, kco, r1xo, r2xo,   &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_n, logc_ratio_u_n,             &
                                          nzt_topo_nestbc_n, 'n', 's' )
             ENDIF
             IF ( humidity )  THEN
                CALL pmci_interp_tril_sn( q, qc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_n, logc_ratio_u_n,             &
                                          nzt_topo_nestbc_n, 'n', 's' )
             ENDIF
             IF ( passive_scalar )  THEN
                CALL pmci_interp_tril_sn( s, sc, ico, jco, kco, r1xo, r2xo,     &
                                          r1yo, r2yo, r1zo, r2zo, nzb_s_inner,  &
                                          logc_u_n, logc_ratio_u_n,             &
                                          nzt_topo_nestbc_n, 'n', 's' )
             ENDIF

             IF ( TRIM( nesting_mode ) == 'one-way' )  THEN
                CALL pmci_extrap_ifoutflow_sn( u, nzb_u_inner, 'n', 'u' )
                CALL pmci_extrap_ifoutflow_sn( v, nzb_v_inner, 'n', 'v' )
                CALL pmci_extrap_ifoutflow_sn( w, nzb_w_inner, 'n', 'w' )
                CALL pmci_extrap_ifoutflow_sn( e, nzb_s_inner, 'n', 'e' )
                IF ( .NOT. neutral )  THEN
                   CALL pmci_extrap_ifoutflow_sn( pt,nzb_s_inner, 'n', 's' )
                ENDIF
                IF ( humidity )  THEN
                   CALL pmci_extrap_ifoutflow_sn( q, nzb_s_inner, 'n', 's' )
                ENDIF
                IF ( passive_scalar )  THEN
                   CALL pmci_extrap_ifoutflow_sn( s, nzb_s_inner, 'n', 's' )
                ENDIF

             ENDIF

          ENDIF

       ENDIF       !: IF ( nesting_mode /= 'vertical' )

!
!--    All PEs are top-border PEs
       CALL pmci_interp_tril_t( u,  uc,  icu, jco, kco, r1xu, r2xu, r1yo,       &
                                r2yo, r1zo, r2zo, 'u' )
       CALL pmci_interp_tril_t( v,  vc,  ico, jcv, kco, r1xo, r2xo, r1yv,       &
                                r2yv, r1zo, r2zo, 'v' )
       CALL pmci_interp_tril_t( w,  wc,  ico, jco, kcw, r1xo, r2xo, r1yo,       &
                                r2yo, r1zw, r2zw, 'w' )
       CALL pmci_interp_tril_t( e,  ec,  ico, jco, kco, r1xo, r2xo, r1yo,       &
                                r2yo, r1zo, r2zo, 'e' )
       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_tril_t( pt, ptc, ico, jco, kco, r1xo, r2xo, r1yo,    &
                                   r2yo, r1zo, r2zo, 's' )
       ENDIF
       IF ( humidity )  THEN
          CALL pmci_interp_tril_t( q, qc, ico, jco, kco, r1xo, r2xo, r1yo,      &
                                   r2yo, r1zo, r2zo, 's' )
       ENDIF
       IF ( passive_scalar )  THEN
          CALL pmci_interp_tril_t( s, sc, ico, jco, kco, r1xo, r2xo, r1yo,      &
                                   r2yo, r1zo, r2zo, 's' )
       ENDIF

       IF ( TRIM( nesting_mode ) == 'one-way' )  THEN
          CALL pmci_extrap_ifoutflow_t( u,  'u' )
          CALL pmci_extrap_ifoutflow_t( v,  'v' )
          CALL pmci_extrap_ifoutflow_t( w,  'w' )
          CALL pmci_extrap_ifoutflow_t( e,  'e' )
          IF ( .NOT. neutral )  THEN
             CALL pmci_extrap_ifoutflow_t( pt, 's' )
          ENDIF
          IF ( humidity )  THEN
             CALL pmci_extrap_ifoutflow_t( q, 's' )
          ENDIF
          IF ( passive_scalar )  THEN
             CALL pmci_extrap_ifoutflow_t( s, 's' )
          ENDIF
       ENDIF

   END SUBROUTINE pmci_interpolation



   SUBROUTINE pmci_anterpolation

!
!--   A wrapper routine for all anterpolation actions.
!--   Note that TKE is not anterpolated.
      IMPLICIT NONE

      CALL pmci_anterp_tophat( u,  uc,  kctu, iflu, ifuu, jflo, jfuo, kflo,     &
                               kfuo, ijfc_u, 'u' )
      CALL pmci_anterp_tophat( v,  vc,  kctu, iflo, ifuo, jflv, jfuv, kflo,     &
                               kfuo, ijfc_v, 'v' )
      CALL pmci_anterp_tophat( w,  wc,  kctw, iflo, ifuo, jflo, jfuo, kflw,     &
                               kfuw, ijfc_s, 'w' )
      IF ( .NOT. neutral )  THEN
         CALL pmci_anterp_tophat( pt, ptc, kctu, iflo, ifuo, jflo, jfuo, kflo,  &
                                  kfuo, ijfc_s, 's' )
      ENDIF
      IF ( humidity )  THEN
         CALL pmci_anterp_tophat( q, qc, kctu, iflo, ifuo, jflo, jfuo, kflo,    &
                                  kfuo, ijfc_s, 's' )
      ENDIF
      IF ( passive_scalar )  THEN
         CALL pmci_anterp_tophat( s, sc, kctu, iflo, ifuo, jflo, jfuo, kflo,    &
                                  kfuo, ijfc_s, 's' )
      ENDIF

   END SUBROUTINE pmci_anterpolation



   SUBROUTINE pmci_interp_tril_lr( f, fc, ic, jc, kc, r1x, r2x, r1y, r2y, r1z,  &
                                   r2z, kb, logc, logc_ratio, nzt_topo_nestbc,  &
                                   edge, var )
!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the left and right boundaries. It is
!--   based on trilinear interpolation.

      IMPLICIT NONE

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                                      INTENT(INOUT) ::  f       !:
      REAL(wp), DIMENSION(0:cg%nz+1,jcs:jcn,icl:icr),                           &
                                      INTENT(IN)    ::  fc      !:
      REAL(wp), DIMENSION(1:2,0:ncorr-1,nzb:nzt_topo_nestbc,nys:nyn),           &
                                      INTENT(IN)    ::  logc_ratio   !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r1x     !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r2x     !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r1y     !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r2y     !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r1z     !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r2z     !:
      
      INTEGER(iwp), DIMENSION(nxlg:nxrg), INTENT(IN)           ::  ic     !:
      INTEGER(iwp), DIMENSION(nysg:nyng), INTENT(IN)           ::  jc     !:
      INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  kb     !:
      INTEGER(iwp), DIMENSION(nzb:nzt+1), INTENT(IN)           ::  kc     !:
      INTEGER(iwp), DIMENSION(1:2,nzb:nzt_topo_nestbc,nys:nyn),                 &
                                          INTENT(IN)           ::  logc   !:
      INTEGER(iwp) ::  nzt_topo_nestbc   !:

      CHARACTER(LEN=1), INTENT(IN) ::  edge   !:
      CHARACTER(LEN=1), INTENT(IN) ::  var    !:

      INTEGER(iwp) ::  i       !:
      INTEGER(iwp) ::  ib      !:
      INTEGER(iwp) ::  ibgp    !:
      INTEGER(iwp) ::  iw      !:
      INTEGER(iwp) ::  j       !:
      INTEGER(iwp) ::  jco     !:
      INTEGER(iwp) ::  jcorr   !:
      INTEGER(iwp) ::  jinc    !:
      INTEGER(iwp) ::  jw      !:
      INTEGER(iwp) ::  j1      !:
      INTEGER(iwp) ::  k       !:
      INTEGER(iwp) ::  kco     !:
      INTEGER(iwp) ::  kcorr   !:
      INTEGER(iwp) ::  k1      !:
      INTEGER(iwp) ::  l       !:
      INTEGER(iwp) ::  m       !:
      INTEGER(iwp) ::  n       !:
      INTEGER(iwp) ::  kbc     !:
      
      REAL(wp) ::  coarse_dx   !:
      REAL(wp) ::  coarse_dy   !:
      REAL(wp) ::  coarse_dz   !:
      REAL(wp) ::  fkj         !:
      REAL(wp) ::  fkjp        !:
      REAL(wp) ::  fkpj        !:
      REAL(wp) ::  fkpjp       !:
      REAL(wp) ::  fk          !:
      REAL(wp) ::  fkp         !:
      
! 
!--   Check which edge is to be handled
      IF ( edge == 'l' )  THEN
!
!--      For u, nxl is a ghost node, but not for the other variables
         IF ( var == 'u' )  THEN
            i  = nxl
            ib = nxl - 1 
         ELSE
            i  = nxl - 1
            ib = nxl - 2
         ENDIF
      ELSEIF ( edge == 'r' )  THEN
         i  = nxr + 1
         ib = nxr + 2
      ENDIF
      
      DO  j = nys, nyn+1
         DO  k = kb(j,i), nzt+1
            l = ic(i)
            m = jc(j)
            n = kc(k)
            fkj      = r1x(i) * fc(n,m,l)     + r2x(i) * fc(n,m,l+1)
            fkjp     = r1x(i) * fc(n,m+1,l)   + r2x(i) * fc(n,m+1,l+1)
            fkpj     = r1x(i) * fc(n+1,m,l)   + r2x(i) * fc(n+1,m,l+1)
            fkpjp    = r1x(i) * fc(n+1,m+1,l) + r2x(i) * fc(n+1,m+1,l+1)
            fk       = r1y(j) * fkj  + r2y(j) * fkjp
            fkp      = r1y(j) * fkpj + r2y(j) * fkpjp
            f(k,j,i) = r1z(k) * fk   + r2z(k) * fkp
         ENDDO
      ENDDO

!
!--   Generalized log-law-correction algorithm.
!--   Doubly two-dimensional index arrays logc(1:2,:,:) and log-ratio arrays 
!--   logc_ratio(1:2,0:ncorr-1,:,:) have been precomputed in subroutine
!--   pmci_init_loglaw_correction.
!
!--   Solid surface below the node 
      IF ( var == 'u' .OR. var == 'v' )  THEN           
         DO  j = nys, nyn
            k = kb(j,i)+1
            IF ( ( logc(1,k,j) /= 0 )  .AND.  ( logc(2,k,j) == 0 ) )  THEN
               k1 = logc(1,k,j)
               DO  kcorr = 0, ncorr - 1
                  kco = k + kcorr
                  f(kco,j,i) = logc_ratio(1,kcorr,k,j) * f(k1,j,i)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!
!--   In case of non-flat topography, also vertical walls and corners need to be
!--   treated. Only single and double wall nodes are corrected. Triple and
!--   higher-multiple wall nodes are not corrected as the log law would not be
!--   valid anyway in such locations.
      IF ( topography /= 'flat' )  THEN
         IF ( var == 'u' .OR. var == 'w' )  THEN                 

!
!--         Solid surface only on south/north side of the node                   
            DO  j = nys, nyn
               DO  k = kb(j,i)+1, nzt_topo_nestbc
                  IF ( ( logc(2,k,j) /= 0 )  .AND.  ( logc(1,k,j) == 0 ) )  THEN

!
!--                  Direction of the wall-normal index is carried in as the
!--                  sign of logc
                     jinc = SIGN( 1, logc(2,k,j) )
                     j1   = ABS( logc(2,k,j) )
                     DO  jcorr = 0, ncorr-1
                        jco = j + jinc * jcorr
                        f(k,jco,i) = logc_ratio(2,jcorr,k,j) * f(k,j1,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

!
!--      Solid surface on both below and on south/north side of the node           
         IF ( var == 'u' )  THEN
            DO  j = nys, nyn
               k = kb(j,i) + 1
               IF ( ( logc(2,k,j) /= 0 )  .AND.  ( logc(1,k,j) /= 0 ) )  THEN
                  k1   = logc(1,k,j)                 
                  jinc = SIGN( 1, logc(2,k,j) )
                  j1   = ABS( logc(2,k,j) )                 
                  DO  jcorr = 0, ncorr-1
                     jco = j + jinc * jcorr
                     DO  kcorr = 0, ncorr-1
                        kco = k + kcorr
                        f(kco,jco,i) = 0.5_wp * ( logc_ratio(1,kcorr,k,j) *     &
                                                  f(k1,j,i)                     &
                                                + logc_ratio(2,jcorr,k,j) *     &
                                                  f(k,j1,i) )
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

      ENDIF  ! ( topography /= 'flat' )

!
!--   Rescale if f is the TKE.
      IF ( var == 'e')  THEN
         IF ( edge == 'l' )  THEN
            DO  j = nys, nyn + 1
               DO  k = kb(j,i), nzt + 1
                  f(k,j,i) = tkefactor_l(k,j) * f(k,j,i)
               ENDDO
            ENDDO
         ELSEIF ( edge == 'r' )  THEN           
            DO  j = nys, nyn+1
               DO  k = kb(j,i), nzt+1
                  f(k,j,i) = tkefactor_r(k,j) * f(k,j,i)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!
!--   Store the boundary values also into the other redundant ghost node layers
      IF ( edge == 'l' )  THEN
         DO  ibgp = -nbgp, ib
            f(0:nzt+1,nysg:nyng,ibgp) = f(0:nzt+1,nysg:nyng,i)
         ENDDO
      ELSEIF ( edge == 'r' )  THEN
         DO  ibgp = ib, nx+nbgp
            f(0:nzt+1,nysg:nyng,ibgp) = f(0:nzt+1,nysg:nyng,i)
         ENDDO
      ENDIF

   END SUBROUTINE pmci_interp_tril_lr



   SUBROUTINE pmci_interp_tril_sn( f, fc, ic, jc, kc, r1x, r2x, r1y, r2y, r1z,  &
                                   r2z, kb, logc, logc_ratio,                   &
                                   nzt_topo_nestbc, edge, var )

!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the south and north boundaries. 
!--   This subroutine is based on trilinear interpolation.

      IMPLICIT NONE

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                                      INTENT(INOUT) ::  f             !:
      REAL(wp), DIMENSION(0:cg%nz+1,jcs:jcn,icl:icr),                           &
                                      INTENT(IN)    ::  fc            !:
      REAL(wp), DIMENSION(1:2,0:ncorr-1,nzb:nzt_topo_nestbc,nxl:nxr),           &
                                      INTENT(IN)    ::  logc_ratio    !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r1x           !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r2x           !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r1y           !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r2y           !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r1z           !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r2z           !:
      
      INTEGER(iwp), DIMENSION(nxlg:nxrg), INTENT(IN)           ::  ic    !:
      INTEGER(iwp), DIMENSION(nysg:nyng), INTENT(IN)           ::  jc    !:
      INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  kb    !:
      INTEGER(iwp), DIMENSION(nzb:nzt+1), INTENT(IN)           ::  kc    !:
      INTEGER(iwp), DIMENSION(1:2,nzb:nzt_topo_nestbc,nxl:nxr),                 &
                                          INTENT(IN)           ::  logc  !:
      INTEGER(iwp) ::  nzt_topo_nestbc   !:

      CHARACTER(LEN=1), INTENT(IN) ::  edge   !:
      CHARACTER(LEN=1), INTENT(IN) ::  var    !:
      
      INTEGER(iwp) ::  i       !:
      INTEGER(iwp) ::  iinc    !:
      INTEGER(iwp) ::  icorr   !:
      INTEGER(iwp) ::  ico     !:
      INTEGER(iwp) ::  i1      !:
      INTEGER(iwp) ::  j       !:
      INTEGER(iwp) ::  jb      !:
      INTEGER(iwp) ::  jbgp    !:
      INTEGER(iwp) ::  k       !:
      INTEGER(iwp) ::  kcorr   !:
      INTEGER(iwp) ::  kco     !:
      INTEGER(iwp) ::  k1      !:
      INTEGER(iwp) ::  l       !:
      INTEGER(iwp) ::  m       !:
      INTEGER(iwp) ::  n       !:
                            
      REAL(wp) ::  coarse_dx   !:
      REAL(wp) ::  coarse_dy   !:
      REAL(wp) ::  coarse_dz   !:
      REAL(wp) ::  fk          !:
      REAL(wp) ::  fkj         !:
      REAL(wp) ::  fkjp        !:
      REAL(wp) ::  fkpj        !:
      REAL(wp) ::  fkpjp       !:
      REAL(wp) ::  fkp         !:
      
!
!--   Check which edge is to be handled: south or north
      IF ( edge == 's' )  THEN
!
!--      For v, nys is a ghost node, but not for the other variables
         IF ( var == 'v' )  THEN
            j  = nys
            jb = nys - 1 
         ELSE
            j  = nys - 1
            jb = nys - 2
         ENDIF
      ELSEIF ( edge == 'n' )  THEN
         j  = nyn + 1
         jb = nyn + 2
      ENDIF

      DO  i = nxl, nxr+1
         DO  k = kb(j,i), nzt+1
            l = ic(i)
            m = jc(j)
            n = kc(k)              
            fkj      = r1x(i) * fc(n,m,l)     + r2x(i) * fc(n,m,l+1)
            fkjp     = r1x(i) * fc(n,m+1,l)   + r2x(i) * fc(n,m+1,l+1)
            fkpj     = r1x(i) * fc(n+1,m,l)   + r2x(i) * fc(n+1,m,l+1)
            fkpjp    = r1x(i) * fc(n+1,m+1,l) + r2x(i) * fc(n+1,m+1,l+1)
            fk       = r1y(j) * fkj  + r2y(j) * fkjp
            fkp      = r1y(j) * fkpj + r2y(j) * fkpjp
            f(k,j,i) = r1z(k) * fk   + r2z(k) * fkp
         ENDDO
      ENDDO

!
!--   Generalized log-law-correction algorithm.
!--   Multiply two-dimensional index arrays logc(1:2,:,:) and log-ratio arrays 
!--   logc_ratio(1:2,0:ncorr-1,:,:) have been precomputed in subroutine
!--   pmci_init_loglaw_correction.
!
!--   Solid surface below the node 
      IF ( var == 'u'  .OR.  var == 'v' )  THEN           
         DO  i = nxl, nxr
            k = kb(j,i) + 1
            IF ( ( logc(1,k,i) /= 0 )  .AND.  ( logc(2,k,i) == 0 ) )  THEN
               k1 = logc(1,k,i)
               DO  kcorr = 0, ncorr-1
                  kco = k + kcorr
                  f(kco,j,i) = logc_ratio(1,kcorr,k,i) * f(k1,j,i)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!
!--   In case of non-flat topography, also vertical walls and corners need to be
!--   treated. Only single and double wall nodes are corrected.
!--   Triple and higher-multiple wall nodes are not corrected as it would be
!--   extremely complicated and the log law would not be valid anyway in such
!--   locations.
      IF ( topography /= 'flat' )  THEN
         IF ( var == 'v' .OR. var == 'w' )  THEN
            DO  i = nxl, nxr
               DO  k = kb(j,i), nzt_topo_nestbc

!
!--               Solid surface only on left/right side of the node           
                  IF ( ( logc(2,k,i) /= 0 )  .AND.  ( logc(1,k,i) == 0 ) )  THEN

!
!--                  Direction of the wall-normal index is carried in as the
!--                  sign of logc
                     iinc = SIGN( 1, logc(2,k,i) )
                     i1  = ABS( logc(2,k,i) )
                     DO  icorr = 0, ncorr-1
                        ico = i + iinc * icorr
                        f(k,j,ico) = logc_ratio(2,icorr,k,i) * f(k,j,i1)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

!
!--      Solid surface on both below and on left/right side of the node           
         IF ( var == 'v' )  THEN
            DO  i = nxl, nxr
               k = kb(j,i) + 1
               IF ( ( logc(2,k,i) /= 0 )  .AND.  ( logc(1,k,i) /= 0 ) )  THEN
                  k1   = logc(1,k,i)         
                  iinc = SIGN( 1, logc(2,k,i) )
                  i1   = ABS( logc(2,k,i) )
                  DO  icorr = 0, ncorr-1
                     ico = i + iinc * icorr
                     DO  kcorr = 0, ncorr-1
                        kco = k + kcorr
                        f(kco,i,ico) = 0.5_wp * ( logc_ratio(1,kcorr,k,i) *     &
                                                  f(k1,j,i)  &
                                                + logc_ratio(2,icorr,k,i) *     &
                                                  f(k,j,i1) )
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         
      ENDIF  ! ( topography /= 'flat' )

!
!--   Rescale if f is the TKE.
      IF ( var == 'e')  THEN
         IF ( edge == 's' )  THEN
            DO  i = nxl, nxr + 1
               DO  k = kb(j,i), nzt+1
                  f(k,j,i) = tkefactor_s(k,i) * f(k,j,i)
               ENDDO
            ENDDO
         ELSEIF ( edge == 'n' )  THEN
            DO  i = nxl, nxr + 1
               DO  k = kb(j,i), nzt+1
                  f(k,j,i) = tkefactor_n(k,i) * f(k,j,i)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!
!--   Store the boundary values also into the other redundant ghost node layers
      IF ( edge == 's' )  THEN
         DO  jbgp = -nbgp, jb
            f(0:nzt+1,jbgp,nxlg:nxrg) = f(0:nzt+1,j,nxlg:nxrg)
         ENDDO
      ELSEIF ( edge == 'n' )  THEN
         DO  jbgp = jb, ny+nbgp
            f(0:nzt+1,jbgp,nxlg:nxrg) = f(0:nzt+1,j,nxlg:nxrg)
         ENDDO
      ENDIF

   END SUBROUTINE pmci_interp_tril_sn

 

   SUBROUTINE pmci_interp_tril_t( f, fc, ic, jc, kc, r1x, r2x, r1y, r2y, r1z,   &
                                  r2z, var )

!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the top boundary. 
!--   This subroutine is based on trilinear interpolation.

      IMPLICIT NONE

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                                      INTENT(INOUT) ::  f     !:
      REAL(wp), DIMENSION(0:cg%nz+1,jcs:jcn,icl:icr),                           &
                                      INTENT(IN)    ::  fc    !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r1x   !:
      REAL(wp), DIMENSION(nxlg:nxrg), INTENT(IN)    ::  r2x   !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r1y   !:
      REAL(wp), DIMENSION(nysg:nyng), INTENT(IN)    ::  r2y   !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r1z   !:
      REAL(wp), DIMENSION(nzb:nzt+1), INTENT(IN)    ::  r2z   !:
      
      INTEGER(iwp), DIMENSION(nxlg:nxrg), INTENT(IN) ::  ic    !:
      INTEGER(iwp), DIMENSION(nysg:nyng), INTENT(IN) ::  jc    !:
      INTEGER(iwp), DIMENSION(nzb:nzt+1), INTENT(IN) ::  kc    !:
      
      CHARACTER(LEN=1), INTENT(IN) :: var   !:

      INTEGER(iwp) ::  i   !:
      INTEGER(iwp) ::  j   !:
      INTEGER(iwp) ::  k   !:
      INTEGER(iwp) ::  l   !:
      INTEGER(iwp) ::  m   !:
      INTEGER(iwp) ::  n   !:
      
      REAL(wp) ::  coarse_dx   !:
      REAL(wp) ::  coarse_dy   !:
      REAL(wp) ::  coarse_dz   !:
      REAL(wp) ::  fk          !:
      REAL(wp) ::  fkj         !:
      REAL(wp) ::  fkjp        !:
      REAL(wp) ::  fkpj        !:
      REAL(wp) ::  fkpjp       !:
      REAL(wp) ::  fkp         !:

      
      IF ( var == 'w' )  THEN
         k  = nzt
      ELSE
         k  = nzt + 1
      ENDIF
     
      DO  i = nxl-1, nxr+1
         DO  j = nys-1, nyn+1
            l = ic(i)
            m = jc(j)
            n = kc(k)              
            fkj      = r1x(i) * fc(n,m,l)     + r2x(i) * fc(n,m,l+1)
            fkjp     = r1x(i) * fc(n,m+1,l)   + r2x(i) * fc(n,m+1,l+1)
            fkpj     = r1x(i) * fc(n+1,m,l)   + r2x(i) * fc(n+1,m,l+1)
            fkpjp    = r1x(i) * fc(n+1,m+1,l) + r2x(i) * fc(n+1,m+1,l+1)
            fk       = r1y(j) * fkj  + r2y(j) * fkjp
            fkp      = r1y(j) * fkpj + r2y(j) * fkpjp
            f(k,j,i) = r1z(k) * fk   + r2z(k) * fkp
         ENDDO
      ENDDO

!
!--   Just fill up the second ghost-node layer for w.
      IF ( var == 'w' )  THEN
         f(nzt+1,:,:) = f(nzt,:,:)
      ENDIF

!
!--   Rescale if f is the TKE.
!--   It is assumed that the bottom surface never reaches the top boundary of a
!--   nest domain.
      IF ( var == 'e' )  THEN
         DO  i = nxl, nxr
            DO  j = nys, nyn
               f(k,j,i) = tkefactor_t(j,i) * f(k,j,i)
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE pmci_interp_tril_t



    SUBROUTINE pmci_extrap_ifoutflow_lr( f, kb, edge, var )
!
!--    After the interpolation of ghost-node values for the child-domain
!--    boundary conditions, this subroutine checks if there is a local outflow
!--    through the boundary. In that case this subroutine overwrites the
!--    interpolated values by values extrapolated from the domain. This
!--    subroutine handles the left and right boundaries. However, this operation
!--    is only needed in case of one-way coupling.

       IMPLICIT NONE

       CHARACTER(LEN=1), INTENT(IN) ::  edge   !:
       CHARACTER(LEN=1), INTENT(IN) ::  var    !:

       INTEGER(iwp) ::  i     !:
       INTEGER(iwp) ::  ib    !:
       INTEGER(iwp) ::  ibgp  !:
       INTEGER(iwp) ::  ied   !:
       INTEGER(iwp) ::  j     !:
       INTEGER(iwp) ::  k     !:
      
       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  kb   !:

       REAL(wp) ::  outnor    !:
       REAL(wp) ::  vdotnor   !:

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) :: f !:

!
!--    Check which edge is to be handled: left or right
       IF ( edge == 'l' )  THEN
          IF ( var == 'u' )  THEN
             i   = nxl
             ib  = nxl - 1
             ied = nxl + 1
          ELSE
             i   = nxl - 1
             ib  = nxl - 2
             ied = nxl
          ENDIF
          outnor = -1.0_wp
       ELSEIF ( edge == 'r' )  THEN
          i      = nxr + 1
          ib     = nxr + 2
          ied    = nxr
          outnor = 1.0_wp
       ENDIF

       DO  j = nys, nyn+1
          DO  k = kb(j,i), nzt+1
             vdotnor = outnor * u(k,j,ied)
!
!--          Local outflow
             IF ( vdotnor > 0.0_wp )  THEN
                f(k,j,i) = f(k,j,ied)
             ENDIF
          ENDDO
          IF ( (var == 'u' )  .OR.  (var == 'v' )  .OR.  (var == 'w') )  THEN
             f(kb(j,i),j,i) = 0.0_wp
          ENDIF
       ENDDO

!
!--    Store the boundary values also into the redundant ghost node layers.
       IF ( edge == 'l' )  THEN
          DO ibgp = -nbgp, ib
             f(0:nzt+1,nysg:nyng,ibgp) = f(0:nzt+1,nysg:nyng,i)
          ENDDO
       ELSEIF ( edge == 'r' )  THEN
          DO ibgp = ib, nx+nbgp
             f(0:nzt+1,nysg:nyng,ibgp) = f(0:nzt+1,nysg:nyng,i)
          ENDDO
       ENDIF

    END SUBROUTINE pmci_extrap_ifoutflow_lr



    SUBROUTINE pmci_extrap_ifoutflow_sn( f, kb, edge, var )
!
!--    After  the interpolation of ghost-node values for the child-domain
!--    boundary conditions, this subroutine checks if there is a local outflow
!--    through the boundary. In that case this subroutine overwrites the
!--    interpolated values by values extrapolated from the domain. This
!--    subroutine handles the south and north boundaries.

       IMPLICIT NONE

       CHARACTER(LEN=1), INTENT(IN) ::  edge   !:
       CHARACTER(LEN=1), INTENT(IN) ::  var    !:
      
       INTEGER(iwp) ::  i         !:
       INTEGER(iwp) ::  j         !:
       INTEGER(iwp) ::  jb        !:
       INTEGER(iwp) ::  jbgp      !:
       INTEGER(iwp) ::  jed       !:
       INTEGER(iwp) ::  k         !:

       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  kb   !:

       REAL(wp)     ::  outnor    !:
       REAL(wp)     ::  vdotnor   !:

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) :: f !:

!
!--    Check which edge is to be handled: left or right
       IF ( edge == 's' )  THEN
          IF ( var == 'v' )  THEN
             j   = nys
             jb  = nys - 1
             jed = nys + 1
          ELSE
             j   = nys - 1
             jb  = nys - 2
             jed = nys
          ENDIF
          outnor = -1.0_wp
       ELSEIF ( edge == 'n' )  THEN
          j      = nyn + 1
          jb     = nyn + 2
          jed    = nyn
          outnor = 1.0_wp
       ENDIF

       DO  i = nxl, nxr+1
          DO  k = kb(j,i), nzt+1
             vdotnor = outnor * v(k,jed,i)
!
!--          Local outflow
             IF ( vdotnor > 0.0_wp )  THEN
                f(k,j,i) = f(k,jed,i)
             ENDIF
          ENDDO
          IF ( (var == 'u' )  .OR.  (var == 'v' )  .OR.  (var == 'w') )  THEN
             f(kb(j,i),j,i) = 0.0_wp
          ENDIF
       ENDDO

!
!--    Store the boundary values also into the redundant ghost node layers.
       IF ( edge == 's' )  THEN
          DO  jbgp = -nbgp, jb
             f(0:nzt+1,jbgp,nxlg:nxrg) = f(0:nzt+1,j,nxlg:nxrg)
          ENDDO
       ELSEIF ( edge == 'n' )  THEN
          DO  jbgp = jb, ny+nbgp
             f(0:nzt+1,jbgp,nxlg:nxrg) = f(0:nzt+1,j,nxlg:nxrg)
          ENDDO
       ENDIF

    END SUBROUTINE pmci_extrap_ifoutflow_sn

 

    SUBROUTINE pmci_extrap_ifoutflow_t( f, var )
!
!--    Interpolation of ghost-node values used as the child-domain boundary
!--    conditions. This subroutine handles the top boundary. It is based on
!--    trilinear interpolation.

       IMPLICIT NONE

       CHARACTER(LEN=1), INTENT(IN) ::  var   !:
      
       INTEGER(iwp) ::  i     !:
       INTEGER(iwp) ::  j     !:
       INTEGER(iwp) ::  k     !:
       INTEGER(iwp) ::  ked   !:

       REAL(wp) ::  vdotnor   !:

       REAL(wp), DIMENSION(nzb:nzt+1,nys-nbgp:nyn+nbgp,nxl-nbgp:nxr+nbgp),      &
                 INTENT(INOUT) ::  f   !:
     

       IF ( var == 'w' )  THEN
          k    = nzt
          ked  = nzt - 1
       ELSE
          k    = nzt + 1
          ked  = nzt
       ENDIF

       DO  i = nxl, nxr
          DO  j = nys, nyn
             vdotnor = w(ked,j,i)
!
!--          Local outflow
             IF ( vdotnor > 0.0_wp )  THEN
                f(k,j,i) = f(ked,j,i)
             ENDIF
          ENDDO
       ENDDO

!
!--    Just fill up the second ghost-node layer for w
       IF ( var == 'w' )  THEN
          f(nzt+1,:,:) = f(nzt,:,:)
       ENDIF

    END SUBROUTINE pmci_extrap_ifoutflow_t



    SUBROUTINE pmci_anterp_tophat( f, fc, kct, ifl, ifu, jfl, jfu, kfl, kfu,    &
                                   ijfc, var )
!
!--    Anterpolation of internal-node values to be used as the parent-domain
!--    values. This subroutine is based on the first-order numerical
!--    integration of the fine-grid values contained within the coarse-grid
!--    cell.

       IMPLICIT NONE

       CHARACTER(LEN=1), INTENT(IN) ::  var   !:

       INTEGER(iwp) ::  i         !: Fine-grid index
       INTEGER(iwp) ::  ii        !: Coarse-grid index
       INTEGER(iwp) ::  iclp      !:
       INTEGER(iwp) ::  icrm      !:
       INTEGER(iwp) ::  j         !: Fine-grid index
       INTEGER(iwp) ::  jj        !: Coarse-grid index
       INTEGER(iwp) ::  jcnm      !:
       INTEGER(iwp) ::  jcsp      !:
       INTEGER(iwp) ::  k         !: Fine-grid index
       INTEGER(iwp) ::  kk        !: Coarse-grid index
       INTEGER(iwp) ::  kcb       !:
       INTEGER(iwp) ::  nfc       !:

       INTEGER(iwp), INTENT(IN) ::  kct   !:

       INTEGER(iwp), DIMENSION(icl:icr), INTENT(IN) ::  ifl   !:
       INTEGER(iwp), DIMENSION(icl:icr), INTENT(IN) ::  ifu   !:
       INTEGER(iwp), DIMENSION(jcs:jcn), INTENT(IN) ::  jfl   !:
       INTEGER(iwp), DIMENSION(jcs:jcn), INTENT(IN) ::  jfu   !:
       INTEGER(iwp), DIMENSION(0:kct), INTENT(IN)   ::  kfl   !:
       INTEGER(iwp), DIMENSION(0:kct), INTENT(IN)   ::  kfu   !:

       INTEGER(iwp), DIMENSION(jcs:jcn,icl:icr), INTENT(IN) ::  ijfc !:

       REAL(wp) ::  cellsum   !:
       REAL(wp) ::  f1f       !:
       REAL(wp) ::  fra       !:

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  f   !:
       REAL(wp), DIMENSION(0:cg%nz+1,jcs:jcn,icl:icr), INTENT(INOUT)  ::  fc  !:
  

!
!--    Initialize the index bounds for anterpolation 
       iclp = icl 
       icrm = icr 
       jcsp = jcs 
       jcnm = jcn 

!
!--    Define the index bounds iclp, icrm, jcsp and jcnm.
!--    Note that kcb is simply zero and kct enters here as a parameter and it is
!--    determined in pmci_init_anterp_tophat

       IF ( nesting_mode == 'vertical' )  THEN
          IF ( nest_bound_l )  THEN
             iclp = icl + nhll
          ENDIF
          IF ( nest_bound_r ) THEN
             icrm = icr - nhlr
          ENDIF
          IF ( nest_bound_s )  THEN
             jcsp = jcs + nhls
          ENDIF
          IF ( nest_bound_n )  THEN
             jcnm = jcn - nhln
          ENDIF
       ELSE
          IF ( nest_bound_l )  THEN
             IF ( var == 'u' )  THEN
                iclp = icl + nhll + 1
             ELSE
                iclp = icl + nhll
             ENDIF
          ENDIF
          IF ( nest_bound_r )  THEN
             icrm = icr - nhlr
          ENDIF

          IF ( nest_bound_s )  THEN
             IF ( var == 'v' )  THEN
                jcsp = jcs + nhls + 1
             ELSE
                jcsp = jcs + nhls
             ENDIF
          ENDIF
          IF ( nest_bound_n )  THEN
             jcnm = jcn - nhln
          ENDIF
          kcb = 0
       ENDIF
       
!
!--    Note that ii, jj, and kk are coarse-grid indices and i,j, and k 
!--    are fine-grid indices.
       DO  ii = iclp, icrm
          DO  jj = jcsp, jcnm
!
!--          For simplicity anterpolate within buildings too
             DO  kk = kcb, kct
!
!--             ijfc is precomputed in pmci_init_anterp_tophat
                nfc =  ijfc(jj,ii) * ( kfu(kk) - kfl(kk) + 1 )
                cellsum = 0.0_wp
                DO  i = ifl(ii), ifu(ii)
                   DO  j = jfl(jj), jfu(jj)
                      DO  k = kfl(kk), kfu(kk)
                         cellsum = cellsum + f(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
!
!--             Spatial under-relaxation.
                fra  = frax(ii) * fray(jj) * fraz(kk)
!
!--             Block out the fine-grid corner patches from the anterpolation
                IF ( ( ifl(ii) < nxl ) .OR. ( ifu(ii) > nxr ) )  THEN
                   IF ( ( jfl(jj) < nys ) .OR. ( jfu(jj) > nyn ) )  THEN
                      fra = 0.0_wp
                   ENDIF
                ENDIF

                fc(kk,jj,ii) = ( 1.0_wp - fra ) * fc(kk,jj,ii) +                &
                               fra * cellsum / REAL( nfc, KIND = wp )

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE pmci_anterp_tophat

#endif
 END SUBROUTINE pmci_child_datatrans

END MODULE pmc_interface
