!> @file modules.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: modules.f90 2055 2016-11-09 09:30:24Z ketelsen $
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2011 2016-09-19 17:29:57Z kanani
! +urban_surface, +lsf_exception, +varnamelength
! 
! 2007 2016-08-24 15:47:17Z kanani
! Increased DIMENSION of data_output, data_output_user, do2d, do3d
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! +constant_top_scalarflux, top_scalarflux
! default of bc_s_t adjusted
! 
! 1968 2016-07-18 12:01:49Z suehring
! Changed dimension for MPI-datatypes type_x_int and type_y_int
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! +bc_s_t_val, diss_s_s, diss_l_s, flux_s_s, flux_l_s, s, sp, s1, s2, s3, ssws_av, 
!  s_init, s_surf, sums_wsss_ws_l, ss, ssws, sswst, ts_m, wall_sflux
! +constant_scalarflux, ibc_s_b, ibc_s_t, s_vertical_gradient_level_ind
!
! Unused variables removed
! -gamma_x, gamma_y, gamma_z, var_x, var_y, var_z
!
! Change initial values (in order to unify gradient calculation): 
! pt_vertical_gradient_level, sa_vertical_gradient_level
!
! 1957 2016-07-07 10:43:48Z suehring
! +fl_max, num_leg, num_var_fl, num_var_fl_user, var_fl_max, virtual_flight
! 
! 1918 2016-05-27 14:35:57Z raasch
! default timestep switched from -1.0 to +1.0 in order to avoid wrong sign of
! initially calculated divergence
!
! 1906 2016-05-24 14:38:08Z suehring
! default value of mg_switch_to_pe0_level changed to -1
!
! 1849 2016-04-08 11:33:18Z hoffmann
! bfactor, mass_of_solute, molecular_weight_of_solute, molecular_weight_of_water,
! vanthoff moved to mod_particle_attributes.
! dt_micro and several cloud_parameters moved to microphysics_mod.
! 1d-microphysics profiles moved to microphysics_mod.
!
! 1845 2016-04-08 08:29:13Z raasch
! -nzb_2d
!
! 1833 2016-04-07 14:23:03Z raasch
! spectra parameter moved to spectra module
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects removed
! turbulence renamed collision_turbulence, drizzle renamed 
! cloud_water_sedimentation
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphysics_sat_adjust, microphysics_kessler,
! microphysics_seifert added.
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-directive for decalpha removed
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
! 
! 1788 2016-03-10 11:01:04Z maronga
! Added roughness length for moisture (z0q)
!
! 1786 2016-03-08 05:49:27Z raasch
! module spectrum moved to new separate module
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables moved to the netcdf-interface module
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char extended to LEN=3
!
! 1764 2016-02-28 12:45:19Z raasch
! some reformatting
!
! 1762 2016-02-25 12:31:13Z hellstea
! +nest_* variables, size of volume_flow arrays increased by one element
!
! 1738 2015-12-18 13:56:05Z raasch
! +mean_surface_level_height
!
! 1695 2015-10-27 10:03:11Z maronga
! Removed rif (forgotten in last revision)
! 
! 1693 2015-10-27 08:35:45Z maronga
! Renamed zp -> z_mo
! 
! 1691 2015-10-26 16:17:44Z maronga
! Renamed Obukhov length. Added ol, removed rif. Increased number of profiles 
! (pr_palm). Changed default values for dissipation_1d to 'detering' and 
! (mixing_length_1d to 'blackadar'. Added most_method. rif_min and rif_max
! renamed to zeta_min and zeta_max and new values assigned.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1677 2015-10-02 13:25:23Z boeske
! +ngp_yz_int, type_xz_int, type_yz_int
! 
! 1666 2015-09-23 07:31:10Z raasch
! +user_interface_current_revision, user_interface_required_revision in
! control_parameters
!
! 1639 2015-08-31 14:46:48Z knoop
! Bugfix: string 'unknown' extended to match LEN=13 
!
! 1575 2015-03-27 09:56:27Z raasch
! +ngp_xz
!
! 1560 2015-03-06 10:48:54Z keck
! +recycling_yshift
! 
! 1557 2015-03-05 16:43:04Z suehring
! +monotonic_adjustment
! 
! 1551 2015-03-03 14:18:16Z maronga
! Increased pr_palm to 120. Increased length of dots_unit and dots_label to 13
! digits. Increased length of domask, do2d, and do3d to 20 digits.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Renamed "radiation" -> "cloud_top_radiation"
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   canopy-model related parameters/variables moved to module 
!   plant_canopy_model_mod
! 
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores.
! Increased identifier string length for user-defined quantities to 20.
! 
! 1450 2014-08-21 07:31:51Z heinze
! ntnudge from 100 to 1000 increased to allow longer simulations
! 
! 1431 2014-07-15 14:47:17Z suehring
! +var_d
! 
! 1429 2014-07-15 12:53:45Z knoop
! +ensemble_member_nr to prepare the random_generator for ensemble runs
!
! 1382 2014-04-30 12:15:41Z boeske
! Renamed variables which store large scale forcing tendencies
! pt_lsa -> td_lsa_lpt, pt_subs -> td_sub_lpt, 
! q_lsa  -> td_lsa_q,   q_subs  -> td_sub_q
! 
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing enabled:
! increase pr_palm from 90 to 100 to allow for more standard profiles
! + ngp_sums_ls, pt_lsa, pt_subs, q_lsa, q_subs, tmp_tnudge, sums_ls_l, 
! use_subsidence_tendencies
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! tend_* removed
! call_microphysics_at_all_substeps added
! default of drizzle set to true
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! particle_attributes moved to mod_particle_attributes.f90
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1327 2014-03-21 11:00:16Z raasch
! REAL constants defined as wp-kind
! -avs_output, data_output_format, do3d_compress, iso2d_output, netcdf_output
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
! module cpulog moved to new separate module-file
! interface for cpu_log removed
!
! 1314 2014-03-14 18:25:17Z suehring
! + log_z_z0, number_of_sublayers, z0_av_global 
! 1308 2014-03-13 14:58:42Z fricke
! +ntdim_2d_xy, ntdim_2d_xz, ntdim_2d_yz, ntdim_3d
!
! 1257 2013-11-08 15:18:40Z raasch
! set default values for grid indices of maximum velocity components
! u|v|w_max_ijk
!
! 1241 2013-10-30 11:36:58Z heinze
! Usage of nudging enabled 
! +nudging, ntnudge, ptnudge, qnudge, tnudge, unudge, vnudge, wnudge
! increase pr_palm from 80 to 90 to allow for more standard profiles
!
! Enable prescribed time depenend surface fluxes and geostrophic wind read in
! from external file LSF_DATA
! +large_scale_forcing, lsf_surf, lsf_vert, nlsf, time_surf, shf_surf, qsws_surf,
!  pt_surf, q_surf, p_surf, time_vert, ug_vert, vg_vert, wsubs_vert
!
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_0 changed to 32bit int, +wall_flags_00,
! +rflags_s_inner, rflags_invers
!
! 1216 2013-08-26 09:31:42Z raasch
! +transpose_compute_overlap,
! several variables are now defined in the serial (non-parallel) case also
!
! 1212 2013-08-15 08:46:27Z raasch
! +tri
!
! 1179 2013-06-14 05:57:58Z raasch
! +reference_state, ref_state, use_initial_profile_as_reference, vpt_reference,
! use_reference renamed use_single_reference_value
!
! 1159 2013-05-21 11:58:22Z fricke
! -bc_lr_dirneu, bc_lr_neudir, bc_ns_dirneu, bc_ns_neudir
! +use_cmax
!
! 1128 2013-04-12 06:19:32Z raasch
! +background_communication, i_left, i_right, j_north, j_south, req, req_count,
! send_receive, sendrecv_in_background, wait_stat
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1113 2013-03-10 02:48:14Z raasch
! +on_device
!
! 1111 2013-03-08 23:54:10Z raasch
! +tric, nr_timesteps_this_run
!
! 1106 2013-03-04 05:31:38Z raasch
! array_kind renamed precision_kind, pdims defined in serial code
! bugfix: default value assigned to coupling_start_time
!
! 1095 2013-02-03 02:21:01Z raasch
! FORTRAN error in r1092 removed
!
! 1092 2013-02-02 11:24:22Z raasch
! character length in some derived type changed for better alignment
!
! 1065 2012-11-22 17:42:36Z hoffmann
! + c_sedimentation, limiter_sedimentation, turbulence, a_1, a_2, a_3, b_1, b_2, 
! + b_3, c_1, c_2, c_3, beta_cc
!
! bottom boundary condition of qr, nr changed from Dirichlet to Neumann
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr) 
! of the two-moment cloud physics scheme:
! +*_init, flux_l_*, diss_l_*, flux_s_*, diss_s_*, *sws, *swst, tend_*, *, *_p
! +t*_m, *_1, *_2, *_3, *_av, bc_*_b, bc_*_t, ibc_*_b, ibc_*_t, bc_*_t_val, 
! +*_vertical_gradient, *_surface_initial_change, *_vertical_gradient_level, 
! +*_vertical_gradient_level_ind, *_surface, constant_waterflux_*,  
! +cloud_scheme, icloud_scheme, surface_waterflux_*, sums_ws*s_ws_l, wall_*flux
!
! constants for the two-moment scheme:
! +a_vent, a_term, b_vent, b_term, c_evap, c_term, cof, eps_sb, k_cc, k_cr, k_rr, 
! +k_br, kappa_rr, kin_vis_air, mu_constant_value, nc, pirho_l, dpirho_l, rho_1,
! +schmidt, schmidt_p_1d3, stp, x0, xmin, xmax, dt_precipitation, w_precipitation
!
! steering parameters for the two_moment scheme:
! +mu_constant, ventilation_effect
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! +output_format_netcdf
!
! 1015 2012-09-27 09:23:24Z raasch
! +acc_rank, num_acc_per_node,
! -adjust_mixing_length
!
! 1010 2012-09-20 07:59:54Z raasch
! pointer free version can be generated with cpp switch __nopointer
!
! 1003 2012-09-14 14:35:53Z raasch
! -grid_matching, nxa, nya, etc., nnx_pe, nny_pe, spl_*
!
! 1001 2012-09-13 14:08:46Z raasch
! -asselin_filter_factor, cut_spline_overshoot, dt_changed, last_dt_change,
! last_dt_change_1d, long_filter_factor, overshoot_limit_*, ups_limit_*
! several pointer/target arrays converted to normal ones
!
! 996 2012-09-07 10:41:47Z raasch
! -use_prior_plot1d_parameters
!
! 978 2012-08-09 08:28:32Z fricke
! +c_u_m, c_u_m_l, c_v_m, c_v_m_l, c_w_m, c_w_m_l,
! +bc_lr_dirneu, bc_lr_neudir, bc_ns_dirneu, bc_ns_neudir
! -km_damp_x, km_damp_y, km_damp_max, outflow_damping_width
! +z0h, z0h_av, z0h_factor, z0h1d
! +ptdf_x, ptdf_y, pt_damping_width, pt_damping_factor 
!
! 964 2012-07-26 09:14:24Z raasch
! -cross_linecolors, cross_linestyles, cross_normalized_x, cross_normx_factor,
! cross_normalized_y, cross_normy_factor, cross_pnc_local,
! cross_profile_numbers, cross_profile_number_counter, cross_uxmax,
! cross_uxmax_computed, cross_uxmax_normalized,
! cross_uxmax_normalized_computed, cross_uxmin, cross_uxmin_computed,
! cross_uxmin_normalized, cross_uxmin_normalized_computed, cross_uymax,
! cross_uymin, cross_xtext, dopr_crossindex, dopr_label, linecolors, linestyles,
! nz_do1d, profil_output, z_max_do1d, z_max_do1d_normalized
!
! 951 2012-07-19 14:22:52Z hoffmann
! changing profile_columns and profile_rows
!
! 940 2012-07-09 14:31:00Z raasch
! +neutral
!
! 927 2012-06-06 19:15:04Z raasch
! +masking_method
!
! 880 2012-04-13 06:28:59Z raasch
! gathered_size, subdomain_size moved to control_parameters
!
! 866 2012-03-28 06:44:41Z raasch
! interface for global_min_max changed
!
! 861 2012-03-26 14:18:34Z suehring
! +wall_flags_0
! -boundary_flags
! +nzb_max
! +adv_sca_1, +adv_mom_1
!
! 849 2012-03-15 10:35:09Z raasch
! +deleted_particles, deleted_tails, tr.._count_sum, tr.._count_recv_sum in
! particle_attributes,
! +de_dx, de_dy, de_dz in arrays_3d,
! first_call_advec_particles renamed first_call_lpm
!
! 828 2012-02-21 12:00:36Z raasch
! +dissipation_classes, radius_classes, use_kernel_tables,
! particle feature color renamed class
!
! 825 2012-02-19 03:03:44Z raasch
! +bfactor, curvature_solution_effects, eps_ros, molecular_weight_of_water,
! vanthoff, -b_cond in cloud_parameters,
! dopts_num increased to 29, particle attributes speed_x|y|z_sgs renamed
! rvar1|2|3
! wang_collision_kernel and turbulence_effects_on_collision replaced by
! collision_kernel, hall_kernel, palm_kernel, wang_kernel
!
! 809 2012-01-30 13:32:58Z marongas
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files.
! New parameter check_restart has been defined which is needed by 
! check_namelist_files only.
!
! 805 2012-01-17 15:53:28Z franke
! Bugfix collective_wait must be out of parallel branch for runs in serial mode
!
! 801 2012-01-10 17:30:36Z suehring
! Dimesion of sums_wsus_ws_l, ! sums_wsvs_ws_l, sums_us2_ws_l, sums_vs2_ws_l,
! sums_ws2_ws_l, sums_wspts_ws_l, sums_wsqs_ws_l, sums_wssas_ws_l increased.
! for thread-safe summation in advec_ws.
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.95  2007/02/11 13:18:30  raasch
! version 3.1b (last under RCS control)
!
! Revision 1.1  1997/07/24 11:21:26  raasch
! Initial revision
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of all variables
!>
!> @todo Add description for each variable
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables for special advection schemes.
!------------------------------------------------------------------------------!
 MODULE advection
 
    USE kinds

    REAL(wp), DIMENSION(:), ALLOCATABLE   ::  aex, bex, dex, eex

    SAVE

 END MODULE advection


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of all arrays defined on the computational grid.
!------------------------------------------------------------------------------!
 MODULE arrays_3d

    USE kinds

    REAL(wp), DIMENSION(:), ALLOCATABLE ::                                     &
          c_u_m, c_u_m_l, c_v_m, c_v_m_l, c_w_m, c_w_m_l, ddzu, ddzu_pres,     &
          dd2zu, dzu, ddzw, dzw, hyp, inflow_damping_factor, l_grid,           &
          ptdf_x, ptdf_y, p_surf, pt_surf, pt_init, qsws_surf, q_init, q_surf, &
          rdf, rdf_sc, ref_state, rs_init, rs_surf, s_init, s_surf, sa_init, shf_surf,           &      !bK added rs_init, rs_surf
          timenudge, time_surf, time_vert, tmp_tnudge, ug, u_init,             &
          u_nzb_p1_for_vfc, vg, v_init, v_nzb_p1_for_vfc, w_subs, zu, zw

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::                                   &
          c_u, c_v, c_w, diss_s_e, diss_s_nr, diss_s_pt, diss_s_q,             &
          diss_s_qr, diss_s_s, diss_s_sa, diss_s_u, diss_s_v, diss_s_w, dzu_mg,&
          dzw_mg, flux_s_e, flux_s_nr, flux_s_pt, flux_s_q, flux_s_qr,         &
          flux_s_s, flux_s_sa, flux_s_u, flux_s_v, flux_s_w, f1_mg, f2_mg,     &
          f3_mg, mean_inflow_profiles, nrs, nrsws, nrswst,                     &
          ol,                                                                  & !< Obukhov length
          precipitation_amount, precipitation_rate, ptnudge, pt_slope_ref,     &
          qnudge, qs, qsws, qswst, qswst_remote, qrs, qrsws, qrswst, rss, rssws, rsswst,   &  !bK added rss, rssws, rsswst
          saswsb, saswst, shf, ss, ssws, sswst, tnudge, td_lsa_lpt, td_lsa_q,  &
          td_sub_lpt,                                                          &
          td_sub_q, total_2d_a, total_2d_o, ts, tswst, ug_vert, unudge, us,    &
          usws, uswst, vnudge, vg_vert, vsws, vswst, wnudge, wsubs_vert,       &
          z0,                                                                  & !< roughness length for momentum
          z0h,                                                                 & !< roughness length for heat
          z0q                                                                    !< roughness length for moisture

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::                                 &
          d, de_dx, de_dy, de_dz, diss, diss_l_e,                              &
          diss_l_nr, diss_l_pt, diss_l_q, diss_l_qr, diss_l_s, diss_l_sa,      &
          diss_l_u, diss_l_v, diss_l_w, flux_l_e, flux_l_nr, flux_l_pt,        &
          flux_l_q, flux_l_qr, flux_l_s, flux_l_sa, flux_l_u, flux_l_v,        &
          flux_l_w, kh, km, l_wall, prr, p_loc, tend, tric,                    &
          u_m_l, u_m_n, u_m_r, u_m_s, v_m_l, v_m_n, v_m_r, v_m_s, w_m_l,       &
          w_m_n, w_m_r, w_m_s

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                         &
          e, e_p, nr, nr_p, p, prho, pt, pt_p, q, q_p, qc, ql, ql_c, ql_v,     &  
          ql_vp, qr, qr_p, rho_ocean, s, s_p, sa, sa_p, te_m, tnr_m, tpt_m, tq_m,   &
          tqr_m, ts_m, tsa_m, tu_m, tv_m, tw_m, u, u_p, v, v_p, vpt, w, w_p 
#else
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                         &
          e_1, e_2, e_3, p, prho_1, nr_1, nr_2, nr_3, pt_1, pt_2, pt_3, q_1,   &        
          q_2, q_3, qc_1, ql_v, ql_vp, ql_1, ql_2, qr_1, qr_2, qr_3, rho_1, rs_1, rs_2, rs_3,    & !bK added rs_1, rs_2, rs_3
          s_1, s_2, s_3, sa_1, sa_2, sa_3, u_1, u_2, u_3, v_1, v_2, v_3, vpt_1,&
          w_1, w_2, w_3

    REAL(wp), DIMENSION(:,:,:), POINTER ::                                     &
          e, e_p,  nr, nr_p, prho, pt, pt_p, q, q_p, qc, ql, ql_c, qr, qr_p,   &
          rho_ocean, s, s_p, sa, sa_p, te_m, tnr_m, tpt_m, tq_m,               &
          tqr_m, ts_m, tsa_m, tu_m, tv_m, tw_m, u, u_p, v, v_p,               &
          vpt, w, w_p
#endif

#ifdef KPP_CHEM
    REAL(wp), DIMENSION(:,:,:), POINTER ::      rs_p, rs, trs_m                       !bK added pointer arrays  copies addresses
                                                                                        ! from pointer array s_p, s and ts_m
 
#endif

    
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  rif_wall, tri

    REAL(wp), DIMENSION(:), ALLOCATABLE :: rho_air     !< air density profile on the uv grid
    REAL(wp), DIMENSION(:), ALLOCATABLE :: rho_air_zw  !< air density profile on the w grid
    REAL(wp), DIMENSION(:), ALLOCATABLE :: drho_air    !< inverse air density profile on the uv grid
    REAL(wp), DIMENSION(:), ALLOCATABLE :: drho_air_zw !< inverse air density profile on the w grid

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rho_air_mg    !< air density profiles on the uv grid for multigrid
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rho_air_zw_mg !< air density profiles on the w grid for multigrid

    REAL(wp), DIMENSION(:), ALLOCATABLE :: heatflux_input_conversion      !< conversion factor array for heatflux input
    REAL(wp), DIMENSION(:), ALLOCATABLE :: waterflux_input_conversion     !< conversion factor array for waterflux input
    REAL(wp), DIMENSION(:), ALLOCATABLE :: momentumflux_input_conversion  !< conversion factor array for momentumflux input
    REAL(wp), DIMENSION(:), ALLOCATABLE :: heatflux_output_conversion     !< conversion factor array for heatflux output
    REAL(wp), DIMENSION(:), ALLOCATABLE :: waterflux_output_conversion    !< conversion factor array for waterflux output
    REAL(wp), DIMENSION(:), ALLOCATABLE :: momentumflux_output_conversion !< conversion factor array for momentumflux output


    SAVE

 END MODULE arrays_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables needed for time-averaging of 2d/3d data.
!------------------------------------------------------------------------------!
 MODULE averaging
 
    USE kinds

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lwp_av,                & !< Avg. liquid water path
                                              precipitation_rate_av, & !< Avg. of precipitation rate
                                              ol_av,                 & !< Avg. of Obukhov length
                                              qsws_av,               & !< Avg. of surface moisture flux
                                              rssws_av,              & !< Avg. of surface rscalar flux      !bK added
                                              ssws_av,               & !< Avg. of surface scalar flux
                                              shf_av,                & !< Avg. of surface heat flux
                                              ts_av,                 & !< Avg. of characteristic temperature scale
                                              us_av,                 & !< Avg. of friction velocity
                                              z0_av,                 & !< Avg. of roughness length for momentum
                                              z0h_av,                & !< Avg. of roughness length for heat
                                              z0q_av                   !< Avg. of roughness length for moisture

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                         &
          e_av, lpt_av, nr_av, p_av, pc_av, pr_av, prr_av, pt_av, q_av, qc_av, &
          ql_av, ql_c_av, ql_v_av, ql_vp_av, qr_av, qv_av, rho_ocean_av, rs_av, s_av, sa_av,&                !bK added rs_av
          u_av, v_av, vpt_av, w_av
 

!    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: rs_av
 
 END MODULE averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables and constants for cloud physics.
!------------------------------------------------------------------------------!
 MODULE cloud_parameters
 
    USE kinds

    REAL(wp) :: cp = 1005.0_wp,       & !< heat capacity of dry air (J kg-1 K-1)
                l_v = 2.5E+06_wp,     & !< latent heat of vaporization (J kg-1)
                l_d_cp, l_d_r, l_d_rv, & !< l_v / cp, l_v / r_d, l_v / r_v
                rho_l = 1.0E3_wp,     & !< density of water (kg m-3)
                r_d = 287.0_wp,       & !< sp. gas const. dry air (J kg-1 K-1)
                r_v = 461.51_wp         !< sp. gas const. water vapor (J kg-1 K-1)


    REAL(wp), DIMENSION(:), ALLOCATABLE     ::  hyrho, pt_d_t, t_d_pt  

    SAVE

 END MODULE cloud_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of general constants.
!------------------------------------------------------------------------------!
 MODULE constants
 
    USE kinds

    REAL(wp)  ::  pi = 3.141592654_wp
    REAL(wp)  ::  adv_mom_1, adv_mom_3, adv_mom_5, adv_sca_1, adv_sca_3, adv_sca_5
    

    SAVE

 END MODULE constants

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of parameters for program control
!------------------------------------------------------------------------------!
 MODULE control_parameters

    USE kinds

    TYPE plot_precision
       CHARACTER (LEN=8) ::  variable
       INTEGER(iwp)      ::  precision
    END TYPE plot_precision

    TYPE(plot_precision), DIMENSION(100) ::  plot_3d_precision =               &
                        (/ plot_precision( 'u', 2 ), plot_precision( 'v', 2 ), &
                           plot_precision( 'w', 2 ), plot_precision( 'p', 5 ), &
                           plot_precision( 'pt', 2 ),                          &
                           ( plot_precision( ' ', 1 ), i9 = 1,95 ) /)

    TYPE file_status
       LOGICAL ::  opened, opened_before
    END TYPE file_status
    
    INTEGER, PARAMETER :: mask_xyz_dimension = 100, max_masks = 50
    INTEGER(iwp), PARAMETER ::  varnamelength = 30  !< length of output variable names

    TYPE(file_status), DIMENSION(200+2*max_masks) ::                &
                             openfile = file_status(.FALSE.,.FALSE.)

    CHARACTER (LEN=1)    ::  cycle_mg = 'w', timestep_reason = ' '
    CHARACTER (LEN=3)    ::  coupling_char = ''
    CHARACTER (LEN=5)    ::  write_binary = 'false'
    CHARACTER (LEN=8)    ::  most_method = 'lookup', & !< NAMELIST parameter defining method to be used to calculate Okukhov length, 
                             run_date,               & !< 
                             run_time                  !< 
    CHARACTER (LEN=9)    ::  simulated_time_chr
    CHARACTER (LEN=11)   ::  topography_grid_convention = ' '
    CHARACTER (LEN=12)   ::  version = ' ', revision = ' ', &
                             user_interface_current_revision = ' ', &
                             user_interface_required_revision = ' '
    CHARACTER (LEN=16)   ::  conserve_volume_flow_mode = 'default', &
                             loop_optimization = 'default', &
                             momentum_advec = 'ws-scheme', &
                             psolver = 'poisfft', &
                             scalar_advec = 'ws-scheme'
    CHARACTER (LEN=20)   ::  approximation = 'boussinesq'
    CHARACTER (LEN=40)   ::  flux_input_mode = 'approximation-specific'
    CHARACTER (LEN=40)   ::  flux_output_mode = 'approximation-specific'
    CHARACTER (LEN=20)   ::  bc_e_b = 'neumann', bc_lr = 'cyclic', &
                             bc_ns = 'cyclic', bc_p_b = 'neumann', &
                             bc_p_t = 'dirichlet', bc_pt_b = 'dirichlet', &
                             bc_pt_t = 'initial_gradient', &
                             bc_q_b = 'dirichlet', bc_q_t = 'neumann', &
                             bc_s_b = 'dirichlet', bc_s_t = 'initial_gradient', &
                             bc_rs_b = 'dirichlet', bc_rs_t = 'initial_gradient', &   ! bK
                             bc_sa_t = 'neumann', &
                             bc_uv_b = 'dirichlet', bc_uv_t = 'dirichlet', &
                             cloud_scheme = 'saturation_adjust', &
                             coupling_mode = 'uncoupled', &
                             coupling_mode_remote = 'uncoupled', &
                             dissipation_1d = 'detering', &
                             fft_method = 'system-specific', &
                             mixing_length_1d = 'blackadar', &
                             random_generator = 'numerical-recipes', &
                             reference_state = 'initial_profile', &
                             return_addres, return_username, &
                             timestep_scheme = 'runge-kutta-3'                             
    CHARACTER (LEN=40)   ::  avs_data_file, topography = 'flat'
    CHARACTER (LEN=64)   ::  host = ' '
    CHARACTER (LEN=80)   ::  log_message, run_identifier
    CHARACTER (LEN=100)  ::  initializing_actions = ' ' 
    CHARACTER (LEN=110)  ::  run_description_header
    CHARACTER (LEN=1000) ::  message_string = ' '

    CHARACTER (LEN=varnamelength), DIMENSION(500) ::  data_output = ' ',    &
                                           data_output_user = ' ', doav = ' '
    CHARACTER (LEN=varnamelength), DIMENSION(max_masks,100) ::  &
         data_output_masks = ' ', data_output_masks_user = ' '

    CHARACTER (LEN=varnamelength), DIMENSION(300) ::  data_output_pr = ' '
    CHARACTER (LEN=varnamelength), DIMENSION(200) ::  data_output_pr_user = ' '
    CHARACTER (LEN=varnamelength), DIMENSION(max_masks,0:1,100) ::  domask = ' '
    CHARACTER (LEN=varnamelength), DIMENSION(0:1,500) ::  do2d = ' ', do3d = ' '

    INTEGER(iwp), PARAMETER :: fl_max = 100, var_fl_max = 20
    
    INTEGER(iwp) ::  abort_mode = 1, average_count_pr = 0, &
                     average_count_3d = 0, current_timestep_number = 0, &
                     coupling_topology = 0, &
                     dist_range = 0, disturbance_level_ind_b, &
                     disturbance_level_ind_t, doav_n = 0, dopr_n = 0, &
                     dopr_time_count = 0, dopts_time_count = 0, &
                     dots_time_count = 0, &
                     do2d_xy_n = 0, do2d_xz_n = 0, do2d_yz_n = 0, do3d_avs_n = 0, &
                     dp_level_ind_b = 0, dvrp_filecount = 0, &
                     dz_stretch_level_index, ensemble_member_nr = 0, gamma_mg, gathered_size, &
                     grid_level, ibc_e_b, ibc_p_b, ibc_p_t, &
                     ibc_pt_b, ibc_pt_t, ibc_q_b, ibc_q_t, ibc_s_b, ibc_s_t, ibc_rs_b, ibc_rs_t, &   ! bK added  ibc_rs_s and t
                     ibc_sa_t, ibc_uv_b, ibc_uv_t, &
                     inflow_disturbance_begin = -1, inflow_disturbance_end = -1, &
                     intermediate_timestep_count, intermediate_timestep_count_max, &
                     io_group = 0, io_blocks = 1, iran = -1234567, &
                     masks = 0, maximum_grid_level, &
                     maximum_parallel_io_streams = -1, max_pr_user = 0, &
                     mgcycles = 0, mg_cycles = -1, mg_switch_to_pe0_level = -1, mid, &
                     nlsf = 1000, ntnudge = 1000, ngsrb = 2, &
                     nr_timesteps_this_run = 0, &
                     nsor = 20, nsor_ini = 100, n_sor, normalizing_region = 0, &
                     num_leg=0, num_var_fl, num_var_fl_user=0, &
                     nz_do3d = -9999, prt_time_count = 0, &
                     recycling_plane, runnr = 0, &
                     skip_do_avs = 0, subdomain_size, terminate_coupled = 0, &
                     terminate_coupled_remote = 0, timestep_count = 0

    INTEGER(iwp) ::  dist_nxl(0:1), dist_nxr(0:1), dist_nyn(0:1), dist_nys(0:1), &
                     do2d_no(0:1) = 0, do2d_xy_time_count(0:1), &
                     do2d_xz_time_count(0:1), do2d_yz_time_count(0:1), &
                     do3d_no(0:1) = 0, do3d_time_count(0:1), &
                     domask_no(max_masks,0:1) = 0, domask_time_count(max_masks,0:1),&
                     mask_size(max_masks,3) = -1, mask_size_l(max_masks,3) = -1, &
                     mask_start_l(max_masks,3) = -1, &
                     pt_vertical_gradient_level_ind(10) = -9999, &
                     q_vertical_gradient_level_ind(10) = -9999, &
                     rs_vertical_gradient_level_ind(10) = -9999, &      !bK added
                     s_vertical_gradient_level_ind(10) = -9999, &                     
                     sa_vertical_gradient_level_ind(10) = -9999, &
                     section(100,3), section_xy(100) = -9999, &
                     section_xz(100) = -9999, section_yz(100) = -9999, &
                     ug_vertical_gradient_level_ind(10) = -9999, &
                     vg_vertical_gradient_level_ind(10) = -9999, &
                     subs_vertical_gradient_level_i(10) = -9999


    INTEGER(iwp), DIMENSION(0:1) :: ntdim_2d_xy, ntdim_2d_xz, ntdim_2d_yz, ntdim_3d

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  grid_level_count

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  mask_i, mask_j, mask_k
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  &
                    mask_i_global, mask_j_global, mask_k_global

    LOGICAL ::  bc_lr_cyc =.TRUE., bc_lr_dirrad = .FALSE., &
                bc_lr_raddir = .FALSE., bc_ns_cyc = .TRUE., &
                bc_ns_dirrad = .FALSE., bc_ns_raddir = .FALSE.,&
                call_microphysics_at_all_substeps = .FALSE., &
                call_psolver_at_all_substeps = .TRUE., &
                cloud_droplets = .FALSE., cloud_physics = .FALSE., &
                cloud_top_radiation = .FALSE., &
                conserve_volume_flow = .FALSE., constant_diffusion = .FALSE., &
                constant_chemistryflux = .TRUE., &                                       !bK added constant_chemistryflux
                constant_flux_layer = .TRUE., &
                constant_heatflux = .TRUE., constant_top_heatflux = .TRUE., &
                constant_top_momentumflux = .FALSE., &
                constant_top_salinityflux = .TRUE., &
                constant_top_scalarflux = .TRUE., &
                constant_top_chemistryflux= .TRUE., &                                   !bK added constant_top_chemistryflux
                constant_scalarflux = .TRUE., &              
                constant_waterflux = .TRUE., create_disturbances = .TRUE., &
                data_output_2d_on_each_pe = .TRUE., & 
                dissipation_control = .FALSE., disturbance_created = .FALSE., &
                do2d_at_begin = .FALSE., do3d_at_begin = .FALSE., &
                do_sum = .FALSE., &
                dp_external = .FALSE., dp_smooth = .FALSE., &
                dt_fixed = .FALSE., &
                dt_3d_reached, dt_3d_reached_l, exchange_mg = .FALSE., &
                first_call_lpm = .TRUE., &
                force_print_header = .FALSE., galilei_transformation = .FALSE.,&
                humidity = .FALSE., humidity_remote = .FALSE., &
                inflow_l = .FALSE., inflow_n = .FALSE., &
                inflow_r = .FALSE., inflow_s = .FALSE., &
                large_scale_forcing = .FALSE., &
                large_scale_subsidence = .FALSE.
    LOGICAL ::  lsf_exception = .FALSE.  !< temporary flag for use of lsf with buildings on flat terrain
    LOGICAL ::  lsf_surf = .TRUE., &
                lsf_vert = .TRUE., lptnudge = .FALSE., lqnudge = .FALSE., &
                lunudge = .FALSE., lvnudge = .FALSE., lwnudge = .FALSE., &
                masking_method = .FALSE., &
                microphysics_sat_adjust = .FALSE., &
                microphysics_kessler = .FALSE., &
                microphysics_seifert = .FALSE., &
                mg_switch_to_pe0 = .FALSE., &
                monotonic_adjustment = .FALSE.
    LOGICAL ::  nest_bound_l = .FALSE. !< nested boundary on left side
    LOGICAL ::  nest_bound_n = .FALSE. !< nested boundary on north side
    LOGICAL ::  nest_bound_r = .FALSE. !< nested boundary on right side
    LOGICAL ::  nest_bound_s = .FALSE. !< nested boundary on south side
    LOGICAL ::  nest_domain  = .FALSE. !< domain is nested into a parent domain
    LOGICAL ::  neutral = .FALSE., nudging = .FALSE., &
                ocean = .FALSE., on_device = .FALSE., &
                outflow_l = .FALSE., outflow_n = .FALSE., outflow_r = .FALSE., &
                outflow_s = .FALSE., passive_scalar = .FALSE., chemistry = .FALSE., &       ! bK added reactive_scalar
                precipitation = .FALSE., &
                random_heatflux = .FALSE., recycling_yshift = .FALSE.,&
                run_control_header = .FALSE., run_coupled = .TRUE., &
                scalar_rayleigh_damping = .TRUE., sloping_surface = .FALSE., &
                stop_dt = .FALSE., synchronous_exchange = .FALSE., &
                terminate_run = .FALSE., transpose_compute_overlap = .FALSE., &
                turbulent_inflow = .FALSE.
    LOGICAL ::  turbulent_outflow = .FALSE. !< flag for turbulent outflow condition
    LOGICAL ::  urban_surface = .FALSE.  !< flag for urban surface model
    LOGICAL ::  use_cmax = .TRUE., use_initial_profile_as_reference = .FALSE., &
                use_prescribed_profile_data = .FALSE., &
                use_single_reference_value = .FALSE., &
                use_subsidence_tendencies = .FALSE., &
                use_surface_fluxes = .FALSE., use_top_fluxes = .FALSE., &
                use_ug_for_galilei_tr = .TRUE., use_upstream_for_tke = .FALSE.
    LOGICAL ::  virtual_flight = .FALSE.  !< flag for virtual flight model
    LOGICAL ::  wall_adjustment = .TRUE., ws_scheme_sca = .FALSE., &
                ws_scheme_mom = .FALSE.

    LOGICAL ::  data_output_xy(0:1) = .FALSE., data_output_xz(0:1) = .FALSE., &
                data_output_yz(0:1) = .FALSE.

    REAL(wp) ::  advected_distance_x = 0.0_wp, advected_distance_y = 0.0_wp, &
                 alpha_surface = 0.0_wp, atmos_ocean_sign = 1.0_wp, &
                 averaging_interval = 0.0_wp, averaging_interval_pr = 9999999.9_wp, &
                 bc_rs_t_val, bc_pt_t_val, bc_q_t_val, bc_s_t_val, bottom_salinityflux = 0.0_wp, &    ! bK added bc_rs_t_val
                 building_height = 50.0_wp, building_length_x = 50.0_wp, &
                 building_length_y = 50.0_wp, building_wall_left = 9999999.9_wp, &
                 building_wall_south = 9999999.9_wp, canyon_height = 50.0_wp, &
                 canyon_width_x = 9999999.9_wp, canyon_width_y = 9999999.9_wp, &
                 canyon_wall_left = 9999999.9_wp, canyon_wall_south = 9999999.9_wp, &
                 cfl_factor = -1.0_wp, cos_alpha_surface, &
                 coupling_start_time = 0.0_wp, disturbance_amplitude = 0.25_wp, &
                 disturbance_energy_limit = 0.01_wp, &
                 disturbance_level_b = -9999999.9_wp, &
                 disturbance_level_t = -9999999.9_wp, &
                 dp_level_b = 0.0_wp, &
                 dt = -1.0_wp, dt_averaging_input = 0.0_wp, &
                 dt_averaging_input_pr = 9999999.9_wp, dt_coupling = 9999999.9_wp, &
                 dt_data_output = 9999999.9_wp, &
                 dt_data_output_av = 9999999.9_wp, dt_disturb = 9999999.9_wp, &
                 dt_dopr = 9999999.9_wp, dt_dopr_listing = 9999999.9_wp, &
                 dt_dopts = 9999999.9_wp, dt_dots = 9999999.9_wp, &
                 dt_do2d_xy = 9999999.9_wp, dt_do2d_xz = 9999999.9_wp, &
                 dt_do2d_yz = 9999999.9_wp, dt_do3d = 9999999.9_wp, dt_dvrp = 9999999.9_wp, &
                 dt_max = 20.0_wp, &
                 dt_restart = 9999999.9_wp, &
                 dt_run_control = 60.0_wp, dt_3d = 1.0_wp, dz = -1.0_wp, &
                 dz_max = 9999999.9_wp, dz_stretch_factor = 1.08_wp, &
                 dz_stretch_level = 100000.0_wp, e_init = 0.0_wp, e_min = 0.0_wp, &
                 end_time = 0.0_wp, &
                 f = 0.0_wp, fs = 0.0_wp, g = 9.81_wp, inflow_damping_height = 9999999.9_wp, &
                 inflow_damping_width = 9999999.9_wp, kappa = 0.4_wp, km_constant = -1.0_wp,&
                 mask_scale_x = 1.0_wp, mask_scale_y = 1.0_wp, mask_scale_z = 1.0_wp, &
                 maximum_cpu_time_allowed = 0.0_wp,  &
                 molecular_viscosity = 1.461E-5_wp, &
                 old_dt = 1.0E-10_wp, omega = 7.29212E-5_wp, omega_sor = 1.8_wp, &
                 outflow_source_plane = -9999999.9_wp, & !< x-position of outflow-source plane (turbulent outflow method)
                 particle_maximum_age = 9999999.9_wp, &
                 phi = 55.0_wp, prandtl_number = 1.0_wp, &
                 precipitation_amount_interval = 9999999.9_wp, prho_reference, &
                 pt_damping_factor = 0.0_wp, pt_damping_width = 0.0_wp, &
                 pt_reference = 9999999.9_wp, pt_slope_offset = 0.0_wp, &
                 pt_surface = 300.0_wp, pt_surface_initial_change = 0.0_wp, &
                 q_surface = 0.0_wp, q_surface_initial_change = 0.0_wp, &
                 rayleigh_damping_factor = -1.0_wp, rayleigh_damping_height = -1.0_wp, &
                 recycling_width = 9999999.9_wp, residual_limit = 1.0E-4_wp, &
                 restart_time = 9999999.9_wp, rho_reference, rho_surface, &
                 roughness_length = 0.1_wp, &
                 sa_surface = 35.0_wp, &
                 simulated_time = 0.0_wp, simulated_time_at_begin, sin_alpha_surface, &
                 skip_time_data_output = 0.0_wp, skip_time_data_output_av = 9999999.9_wp,&
                 skip_time_dopr = 9999999.9_wp, &
                 skip_time_do2d_xy = 9999999.9_wp, skip_time_do2d_xz = 9999999.9_wp, &
                 skip_time_do2d_yz = 9999999.9_wp, skip_time_do3d = 9999999.9_wp, &
                 surface_chemistryflux= 9999999.9_wp, surface_heatflux = 9999999.9_wp, &        ! bK added surface_chemistryflux
                 surface_pressure = 1013.25_wp, surface_scalarflux = 9999999.9_wp, &
                 surface_waterflux = 9999999.9_wp,  s_surface = 0.0_wp,  rs_surface = 0.0_wp, &  ! bK added rs_surface
                 rs_surface_initial_change = 0.0_wp, s_surface_initial_change = 0.0_wp, termination_time_needed = -1.0_wp, &     !bk added rs_surface_initial_change
                 time_coupling = 0.0_wp, time_disturb = 0.0_wp, time_dopr = 0.0_wp, &
                 time_dopr_av = 0.0_wp, time_dopr_listing = 0.0_wp, time_dopts = 0.0_wp, &
                 time_dosp = 0.0_wp, time_dosp_av = 0.0_wp, time_dots = 0.0_wp, & 
                 time_do2d_xy = 0.0_wp, time_do2d_xz = 0.0_wp, time_do2d_yz = 0.0_wp, &
                 time_do3d = 0.0_wp, &
                 time_do_av = 0.0_wp, time_do_sla = 0.0_wp, time_dvrp = 0.0_wp, &
                 time_restart = 9999999.9_wp, time_run_control = 0.0_wp,&
                 time_since_reference_point, top_heatflux = 9999999.9_wp, &
                 top_momentumflux_u = 9999999.9_wp, &
                 top_momentumflux_v = 9999999.9_wp, top_salinityflux = 9999999.9_wp, &
                 top_scalarflux = 9999999.9_wp, top_chemistryflux,  &                              !bK added top_chemistry_flux
                 ug_surface = 0.0_wp, u_bulk = 0.0_wp, u_gtrans = 0.0_wp, &
                 vg_surface = 0.0_wp, vpt_reference = 9999999.9_wp, &
                 v_bulk = 0.0_wp, v_gtrans = 0.0_wp, wall_adjustment_factor = 1.8_wp, &
                 zeta_max = 20.0_wp,    & !< Maximum value of zeta = z/L
                 zeta_min = -9990.0_wp, & !< Minimum value of zeta = z/L
                 z_max_do2d = -1.0_wp, z0h_factor = 1.0_wp

    REAL(wp) ::  do2d_xy_last_time(0:1) = -1.0_wp, do2d_xz_last_time(0:1) = -1.0_wp, &
                 do2d_yz_last_time(0:1) = -1.0_wp, dpdxy(1:2) = 0.0_wp, &
                 dt_domask(max_masks) = 9999999.9_wp, &
                 mask_scale(3), &
                 pt_vertical_gradient(10) = 0.0_wp, &
                 pt_vertical_gradient_level(10) = -1.0_wp, &
                 q_vertical_gradient(10) = 0.0_wp, &
                 q_vertical_gradient_level(10) = -1.0_wp, &
                 rs_vertical_gradient(10) = 0.0_wp, &               !bK added
                 rs_vertical_gradient_level(10) = -1.0_wp, &        !bK added
                 s_vertical_gradient(10) = 0.0_wp, &
                 s_vertical_gradient_level(10) = -1.0_wp, &
                 sa_vertical_gradient(10) = 0.0_wp, &
                 sa_vertical_gradient_level(10) = -1.0_wp, &
                 skip_time_domask(max_masks) = 9999999.9_wp, threshold(20) = 0.0_wp, &
                 time_domask(max_masks) = 0.0_wp, &
                 tsc(10) = (/ 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /), &
                 u_profile(100) = 9999999.9_wp, uv_heights(100) = 9999999.9_wp, &
                 v_profile(100) = 9999999.9_wp, &
                 ug_vertical_gradient(10) = 0.0_wp, &
                 ug_vertical_gradient_level(10) = -9999999.9_wp, &
                 vg_vertical_gradient(10) = 0.0_wp, &
                 vg_vertical_gradient_level(10) = -9999999.9_wp, &
                 volume_flow(1:3) = 0.0_wp, volume_flow_area(1:3) = 0.0_wp, &
                 volume_flow_initial(1:3) = 0.0_wp, wall_heatflux(0:4) = 0.0_wp, &
                 wall_humidityflux(0:4) = 0.0_wp, wall_nrflux(0:4) = 0.0_wp, &
                 wall_qflux(0:4) = 0.0_wp, wall_qrflux(0:4) = 0.0_wp, &
                 wall_salinityflux(0:4) = 0.0_wp, wall_sflux(0:4) = 0.0_wp, wall_scalarflux(0:4) = 0.0_wp, &
                 subs_vertical_gradient(10) = 0.0_wp, &
                 subs_vertical_gradient_level(10) = -9999999.9_wp

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dp_smooth_factor

    REAL(wp), DIMENSION(max_masks,mask_xyz_dimension) :: &
            mask_x = -1.0_wp, mask_y = -1.0_wp, mask_z = -1.0_wp
    REAL(wp), DIMENSION(max_masks,3) ::                  &
            mask_x_loop = -1.0_wp, mask_y_loop = -1.0_wp, mask_z_loop = -1.0_wp
    
!
!--    internal mask arrays ("mask,dimension,selection")
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  mask, mask_loop

    SAVE

 END MODULE control_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables used with dvrp-software.
!------------------------------------------------------------------------------!
 MODULE dvrp_variables

     USE kinds

    CHARACTER (LEN=10) ::  dvrp_output = 'rtsp', particle_color = 'none', &
                           particle_dvrpsize = 'none'

    CHARACTER (LEN=20), DIMENSION(10) ::  mode_dvrp = &
                                     (/ ( '                    ', i9 = 1,10 ) /)

    CHARACTER (LEN=80) ::  dvrp_directory = 'default',                    &
                           dvrp_file      = 'default',                    &
                           dvrp_host      = 'origin.rvs.uni-hannover.de', &
                           dvrp_password  = '********',                   &
                           dvrp_username  = ' '

    INTEGER(iwp) ::  cluster_size = 1, dvrp_colortable_entries = 4,                 &
                     dvrp_colortable_entries_prt = 22, islice_dvrp,                 &
                     nx_dvrp, nxl_dvrp, nxr_dvrp, ny_dvrp, nyn_dvrp, nys_dvrp,      &
                     nz_dvrp, pathlines_fadeintime = 5, pathlines_fadeouttime = 5,  &
                     pathlines_linecount = 1000, pathlines_maxhistory = 40,         &
                     pathlines_wavecount = 10, pathlines_wavetime = 50,             &
                     vc_gradient_normals = 0, vc_mode = 0, vc_size_x = 2,           &
                     vc_size_y = 2, vc_size_z = 2

    INTEGER(iwp), DIMENSION(10) ::  slicer_position_dvrp

    LOGICAL ::  cyclic_dvrp = .FALSE., dvrp_overlap, dvrp_total_overlap, &
                local_dvrserver_running, lock_steering_update = .FALSE., &
                use_seperate_pe_for_dvrp_output = .FALSE.

    REAL(wp)    ::  clip_dvrp_l = 9999999.9_wp, clip_dvrp_n = 9999999.9_wp, &
                    clip_dvrp_r = 9999999.9_wp, clip_dvrp_s = 9999999.9_wp, &
                    superelevation = 1.0_wp, superelevation_x = 1.0_wp,     &
                    superelevation_y = 1.0_wp, vc_alpha = 38.0_wp

    REAL(wp), DIMENSION(2) ::  color_interval = (/ 0.0_wp, 1.0_wp /), &
                               dvrpsize_interval = (/ 0.0_wp, 1.0_wp /)

    REAL(wp), DIMENSION(3) ::  groundplate_color = (/ 0.0_wp, 0.6_wp, 0.0_wp /), &
                               topography_color = (/ 0.8_wp, 0.7_wp, 0.6_wp /)

    REAL(wp), DIMENSION(2,10)     ::  slicer_range_limits_dvrp

    REAL(wp), DIMENSION(3,10)     ::  isosurface_color

    REAL(sp), DIMENSION(2,100) ::  interval_values_dvrp,                       &
                                   interval_values_dvrp_prt, interval_h_dvrp,  &
                                   interval_h_dvrp_prt, interval_l_dvrp = 0.5_sp, &
                                   interval_l_dvrp_prt = 0.5_sp, interval_s_dvrp = 1.0_sp, &
                                   interval_s_dvrp_prt = 1.0_sp, interval_a_dvrp = 0.0_sp, &
                                   interval_a_dvrp_prt = 0.0_sp

    DATA  slicer_range_limits_dvrp / -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, &
                                     -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, &
                                     -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, &
                                     -1.0_wp, 1.0_wp /

    DATA  isosurface_color / 0.9_wp, 0.9_wp, 0.9_wp,  0.8_wp, 0.1_wp, 0.1_wp,  0.1_wp, 0.1_wp, 0.8_wp, &
                             0.1_wp, 0.8_wp, 0.1_wp,  0.6_wp, 0.1_wp, 0.1_wp,  0.1_wp, 0.1_wp, 0.6_wp, &
                             0.1_wp, 0.6_wp, 0.1_wp,  0.4_wp, 0.1_wp, 0.1_wp,  0.1_wp, 0.1_wp, 0.4_wp, &
                             0.1_wp, 0.4_wp, 0.1_wp /

    DATA  interval_h_dvrp / 270.0_wp, 225.0_wp, 225.0_wp, 180.0_wp, 70.0_wp, 25.0_wp, &
                            25.0_wp, -25.0_wp, 192 * 0.0_wp /

    DATA  interval_h_dvrp_prt / 270.0_wp, 225.0_wp, 225.0_wp, 180.0_wp, 70.0_wp, 25.0_wp, &
                                25.0_wp, -25.0_wp, 192 * 0.0_wp /

    REAL(sp), DIMENSION(:), ALLOCATABLE ::  xcoor_dvrp, ycoor_dvrp, zcoor_dvrp

    TYPE steering
       CHARACTER (LEN=24) ::  name
       REAL(sp)           ::  min, max
       INTEGER(iwp)       ::  imin, imax
    END TYPE steering

    TYPE(steering), DIMENSION(:), ALLOCATABLE ::  steering_dvrp

    SAVE

 END MODULE dvrp_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of grid spacings.
!------------------------------------------------------------------------------!
 MODULE grid_variables

    USE kinds

    REAL(wp) ::  ddx, ddx2, dx = 1.0_wp, dx2, ddy, ddy2, dy = 1.0_wp, dy2

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddx2_mg, ddy2_mg

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fwxm, fwxp, fwym, fwyp, fxm, fxp,   &
                                              fym, fyp, wall_e_x, wall_e_y,       &
                                              wall_u, wall_v, wall_w_x, wall_w_y, &
                                              zu_s_inner, zw_w_inner

    SAVE

 END MODULE grid_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of array bounds, number of gridpoints, and wall flag arrays.
!------------------------------------------------------------------------------!
 MODULE indices

    USE kinds

    INTEGER(iwp) ::  i_left, i_right, j_north, j_south, nbgp = 3, ngp_sums,         &
                     ngp_sums_ls, nnx, nx = 0, nx_a, nx_o, nxl, nxlg, nxlu, nxr,    &
                     nxrg, nx_on_file, nny, ny = 0, ny_a, ny_o, nyn, nyng, nys,     &
                     nysg, nysv, ny_on_file, nnz, nz = 0, nzb, nzb_diff, nzb_max,   &
                     nzt, nzt_diff


    INTEGER(idp), DIMENSION(:), ALLOCATABLE ::                                    &
                ngp_3d, ngp_3d_inner   ! need to have 64 bit for grids > 2E9

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::                                    &
                   ngp_2dh, nxl_mg, nxr_mg, nyn_mg, nys_mg, nzt_mg


    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: ngp_2dh_outer, ngp_2dh_s_inner,  &
                   mg_loc_ind, nzb_diff_s_inner, nzb_diff_s_outer, nzb_diff_u,    &
                   nzb_diff_v, nzb_inner, nzb_outer, nzb_s_inner, nzb_s_outer,    &
                   nzb_u_inner, nzb_u_outer, nzb_v_inner, nzb_v_outer,            &
                   nzb_w_inner, nzb_w_outer

    INTEGER(iwp), DIMENSION(:,:,:), POINTER ::  flags

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  wall_flags_0, wall_flags_00

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE,  TARGET ::                       &
                   wall_flags_1, wall_flags_2, wall_flags_3, wall_flags_4,        &
                   wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,        &
                   wall_flags_9, wall_flags_10

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rflags_s_inner, rflags_invers

    SAVE

 END MODULE indices


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interfaces for special subroutines which use optional parameters.
!------------------------------------------------------------------------------!
 MODULE interfaces

    INTERFACE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
       SUBROUTINE global_min_max ( i1, i2, j1, j2, k1, k2, feld, mode, offset, &
                                   wert, wert_ijk, wert1, wert1_ijk )

          USE kinds

          CHARACTER (LEN=*), INTENT(IN)      ::  mode
          INTEGER(iwp), INTENT(IN)           ::  i1, i2, j1, j2, k1, k2
          INTEGER(iwp)                       ::  wert_ijk(3)
          INTEGER(iwp), OPTIONAL             ::  wert1_ijk(3)
          REAL(wp)                           ::  offset, wert
          REAL(wp), OPTIONAL                 ::  wert1
          REAL(wp), INTENT(IN)               ::  feld(i1:i2,j1:j2,k1:k2)

       END SUBROUTINE global_min_max

    END INTERFACE

    SAVE

 END MODULE interfaces


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interfaces for subroutines with pointer arguments called in
!> prognostic_equations.
!------------------------------------------------------------------------------!
 MODULE pointer_interfaces

    INTERFACE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
       SUBROUTINE advec_s_bc( sk, sk_char )

          USE kinds

          CHARACTER (LEN=*), INTENT(IN)   ::  sk_char
#if defined( __nopointer )
          REAL(wp), DIMENSION(:,:,:) ::  sk
#else
          REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif
       END SUBROUTINE advec_s_bc

    END INTERFACE


    SAVE

 END MODULE pointer_interfaces


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables for the 1D-model.
!------------------------------------------------------------------------------!
 MODULE model_1d

    USE kinds

    INTEGER(iwp) ::  current_timestep_number_1d = 0, damp_level_ind_1d

    LOGICAL ::  run_control_header_1d = .FALSE., stop_dt_1d = .FALSE. 

    REAL(wp) ::     damp_level_1d = -1.0_wp, dt_1d = 60.0_wp, dt_max_1d = 300.0_wp, &
                    dt_pr_1d = 9999999.9_wp, dt_run_control_1d = 60.0_wp, &
                    end_time_1d = 864000.0_wp, old_dt_1d = 1.0E-10_wp, &
                    qs1d, simulated_time_1d = 0.0_wp, time_pr_1d = 0.0_wp, &
                    time_run_control_1d = 0.0_wp, ts1d, us1d, usws1d, &
                    vsws1d, z01d, z0h1d


    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e1d, e1d_p, kh1d, km1d, l_black, l1d,  &
                                            rif1d, te_e, te_em, te_u, te_um, te_v, &
                                            te_vm, u1d, u1d_p, v1d, v1d_p

    SAVE

 END MODULE model_1d




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables which define processor topology and the exchange of
!> ghost point layers. This modules must be placed in all routines which contain
!> MPI-calls.
!------------------------------------------------------------------------------!
 MODULE pegrid

    USE kinds

#if defined( __parallel )
#if defined( __mpifh )
    INCLUDE "mpif.h"
#else
    USE MPI
#endif
#endif
    CHARACTER(LEN=2) ::  send_receive = 'al'
    CHARACTER(LEN=7) ::  myid_char = ''
    INTEGER(iwp)     ::  acc_rank, comm1dx, comm1dy, comm2d, comm_inter,       &
                         comm_palm, id_inflow = 0,                             &
                         id_outflow = 0,        & !< myidx of procs at outflow (turbulent outflow method)
                         id_outflow_source = 0, & !< myidx of procs including ouflow source plane (turbulent outflow method)
                         id_recycling = 0, ierr,                               &
                         myid = 0, myidx = 0, myidy = 0, ndim = 2, ngp_a,      &
                         ngp_o, ngp_xy, ngp_y, npex = -1, npey = -1,           &
                         numprocs = 1, numprocs_previous_run = -1,             &
                         num_acc_per_node = 0, pleft, pnorth, pright, psouth,  &
                         req_count = 0, sendrecvcount_xy, sendrecvcount_yz,    &
                         sendrecvcount_zx, sendrecvcount_zyd,                  &
                         sendrecvcount_yxd, target_id, tasks_per_node = -9999, &
                         threads_per_task = 1, type_x, type_xy,    &
                         type_y

    INTEGER(iwp)          ::  pdims(2) = 1, req(100)

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  hor_index_bounds, &
                                                  hor_index_bounds_previous_run

    LOGICAL ::  background_communication =.FALSE., collective_wait = .FALSE., &
                sendrecv_in_background = .FALSE.

#if defined( __parallel )
#if defined( __mpi2 )
    CHARACTER (LEN=MPI_MAX_PORT_NAME) ::  port_name
#endif

    INTEGER(iwp) ::  ibuf(12), pcoord(2)
    INTEGER(iwp) ::  status(MPI_STATUS_SIZE)
    INTEGER(iwp), DIMENSION(MPI_STATUS_SIZE,100) ::  wait_stat

    INTEGER(iwp) :: ngp_yz_int, type_xz_int, type_yz_int
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ngp_xz, ngp_yz, type_x_int,    &
                                                type_xz, type_y_int, type_yz

    LOGICAL ::  left_border_pe  = .FALSE., north_border_pe = .FALSE., &
                reorder = .TRUE., right_border_pe = .FALSE.,          &
                south_border_pe = .FALSE.

    LOGICAL, DIMENSION(2) ::  cyclic = (/ .TRUE. , .TRUE. /), &
                              remain_dims
#endif

    SAVE

 END MODULE pegrid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables which control PROFIL-output.
!------------------------------------------------------------------------------!
 MODULE profil_parameter

    USE kinds

    INTEGER(iwp), PARAMETER ::  crmax = 100

    CHARACTER (LEN=20), DIMENSION(20) ::  cross_ts_profiles = &
                           (/ ' E E*               ', ' dt                 ', &
                              ' u* w*              ', ' th*                ', &
                              ' umax vmax wmax     ', ' div_old div_new    ', &
                              ' z_i_wpt z_i_pt     ', ' w"pt"0 w"pt" wpt   ', &
                              ' pt(0) pt(zp)       ', ' splux spluy spluz  ', &
                              ' L                  ',                         &
                            ( '                    ', i9 = 1, 9 ) /)

    CHARACTER (LEN=100), DIMENSION(crmax) ::  cross_profiles = &
                           (/ ' u v                           ', &
                              ' pt                            ', &
                              ' w"pt" w*pt* w*pt*BC wpt wptBC ', &
                              ' w"u" w*u* wu w"v" w*v* wv     ', &
                              ' km kh                         ', &
                              ' l                             ', &
                         ( '                               ', i9 = 1, 94 ) /)

    INTEGER(iwp) ::  profile_columns = 2, profile_rows = 3, profile_number = 0

    INTEGER(iwp) ::  cross_ts_numbers(crmax,crmax) = 0, &
                     cross_ts_number_count(crmax) = 0, &
                     dopr_index(300) = 0, dopr_initial_index(300) = 0, &
                     dots_crossindex(100) = 0, dots_index(100) = 0
                

    REAL(wp) ::  cross_ts_uymax(20) = &
                             (/ 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp,   &
                                999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp,   &
                                999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp,   &
                                999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp /),&
                 cross_ts_uymax_computed(20) = 999.999_wp, &
                 cross_ts_uymin(20) = &
                             (/ 999.999_wp, 999.999_wp, 999.999_wp,  -5.000_wp, 999.999_wp,   &
                                999.999_wp,   0.000_wp, 999.999_wp, 999.999_wp, 999.999_wp,   &
                                999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp,   &
                                999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp, 999.999_wp /),&
                 cross_ts_uymin_computed(20) = 999.999_wp

    SAVE

 END MODULE profil_parameter

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of statistical quantities, e.g. global sums.
!------------------------------------------------------------------------------!
 MODULE statistics

    USE kinds

    CHARACTER (LEN=40) ::  region(0:9)
    INTEGER(iwp) ::  pr_palm = 130, statistic_regions = 0
    INTEGER(iwp) ::  u_max_ijk(3) = -1, v_max_ijk(3) = -1, w_max_ijk(3) = -1
    LOGICAL ::  flow_statistics_called = .FALSE.
    REAL(wp) ::     u_max, v_max, w_max
    REAL(wp), DIMENSION(:), ALLOCATABLE       ::  mean_surface_level_height,      &
                                                  sums_divnew_l, sums_divold_l,   &
                                                  weight_substep, weight_pres
    REAL(wp), DIMENSION(:,:), ALLOCATABLE     ::  sums, sums_wsts_bc_l, ts_value, &
                                                  sums_wsus_ws_l, sums_wsvs_ws_l, &
                                                  sums_us2_ws_l, sums_vs2_ws_l,   &
                                                  sums_ws2_ws_l,                  &
                                                  sums_wsnrs_ws_l,                &
                                                  sums_wspts_ws_l,                &
                                                  sums_wsqs_ws_l,                 &   
                                                  sums_wsqrs_ws_l,                &
                                                  sums_wssas_ws_l,                &
                                                  sums_wsss_ws_l,                 &
                                                  sums_ls_l
                                              
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  hom_sum, rmask, sums_l,      &
                                                  sums_l_l, sums_up_fraction_l

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  hom, hom_kchem

    SAVE

 END MODULE statistics



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of indices for transposed arrays.
!------------------------------------------------------------------------------!
 MODULE transpose_indices

    USE kinds

    INTEGER(iwp) ::  nxl_y, nxl_yd, nxl_z, nxr_y, nxr_yd, nxr_z, nyn_x, nyn_z, &
                     nys_x, nys_z, nzb_x, nzb_y, nzb_yd, nzt_x, nzt_y, nzt_yd
                
    SAVE

 END MODULE transpose_indices