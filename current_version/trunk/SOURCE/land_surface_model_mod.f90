!> @file land_surface_model_mod.f90
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
! $Id: land_surface_model_mod.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1978 2016-07-29 12:08:31Z maronga
! Bugfix: initial values of pave_surface and water_surface were not set.
! 
! 1976 2016-07-27 13:28:04Z maronga
! Parts of the code have been reformatted. Use of radiation model output is
! generalized and simplified. Added more output quantities due to modularization
! 
! 1972 2016-07-26 07:52:02Z maronga
! Further modularization: output of cross sections and 3D data is now done in this
! module. Moreover, restart data is written and read directly within this module.
!
! 
! 1966 2016-07-18 11:54:18Z maronga
! Bugfix: calculation of m_total in soil model was not set to zero at model start
! 
! 1949 2016-06-17 07:19:16Z maronga
! Bugfix: calculation of qsws_soil_eb with precipitation = .TRUE. gave
! qsws_soil_eb = 0 due to a typo
! 
! 1856 2016-04-13 12:56:17Z maronga
! Bugfix: for water surfaces, the initial water surface temperature is set equal
! to the intital skin temperature. Moreover, the minimum value of r_a is now
! 1.0 to avoid too large fluxes at the first model time step
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! prr moved to arrays_3d
!
! 1826 2016-04-07 12:01:39Z maronga
! Cleanup after modularization
! 
! 1817 2016-04-06 15:44:20Z maronga
! Added interface for lsm_init_arrays. Added subroutines for check_parameters,
! header, and parin. Renamed some subroutines.
! 
! 1788 2016-03-10 11:01:04Z maronga
! Bugfix: calculate lambda_surface based on temperature gradient between skin
! layer and soil layer instead of Obukhov length
! Changed: moved calculation of surface specific humidity to energy balance solver
! New: water surfaces are available by using a fixed sea surface temperature.
! The roughness lengths are calculated dynamically using the Charnock 
! parameterization. This involves the new roughness length for moisture z0q.
! New: modified solution of the energy balance solver and soil model for
! paved surfaces (i.e. asphalt concrete).
! Syntax layout improved.
! Changed: parameter dewfall removed.
! 
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables moved to netcdf module
!
! 1757 2016-02-22 15:49:32Z maronga
! Bugfix: set tm_soil_m to zero after allocation. Added parameter 
! unscheduled_radiation_calls to control calls of the radiation model based on
! the skin temperature change during one time step (preliminary version). Set 
! qsws_soil_eb to zero at model start (previously set to qsws_eb). Removed MAX
! function as it cannot be vectorized.
! 
! 1709 2015-11-04 14:47:01Z maronga
! Renamed pt_1 and qv_1 to pt1 and qv1.
! Bugfix: set initial values for t_surface_p in case of restart runs
! Bugfix: zero resistance caused crash when using radiation_scheme = 'clear-sky'
! Bugfix: calculation of rad_net when using radiation_scheme = 'clear-sky'
! Added todo action
!
! 1697 2015-10-28 17:14:10Z raasch
! bugfix: misplaced cpp-directive
!
! 1695 2015-10-27 10:03:11Z maronga
! Bugfix: REAL constants provided with KIND-attribute in call of 
! Replaced rif with ol
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added skip_time_do_lsm to allow for spin-ups without LSM. Various bugfixes:
! Soil temperatures are now defined at the edges of the layers, calculation of
! shb_eb corrected, prognostic equation for skin temperature corrected. Surface
! fluxes are now directly transfered to atmosphere
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: definition of character strings requires same length for all elements
! 
! 1585 2015-04-30 07:05:52Z maronga
! Modifications for RRTMG. Changed tables to PARAMETER type. 
! 
! 1571 2015-03-12 16:12:49Z maronga
! Removed upper-case variable names. Corrected distribution of precipitation to 
! the liquid water reservoir and the bare soil fractions.
! 
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
! 
! 1553 2015-03-03 17:33:54Z maronga
! Improved better treatment of roughness lengths. Added default soil temperature
! profile
! 
! 1551 2015-03-03 14:18:16Z maronga
! Flux calculation is now done in prandtl_fluxes. Added support for data output.
! Vertical indices have been replaced. Restart runs are now possible. Some
! variables have beem renamed. Bugfix in the prognostic equation for the surface
! temperature. Introduced z0_eb and z0h_eb, which overwrite the setting of
! roughness_length and z0_factor. Added Clapp & Hornberger parametrization for
! the hydraulic conductivity. Bugfix for root fraction and extraction
! calculation
! 
! intrinsic function MAX and MIN
!
! 1500 2014-12-03 17:42:41Z maronga
! Corrected calculation of aerodynamic resistance (r_a).
! Precipitation is now added to liquid water reservoir using LE_liq.
! Added support for dry runs.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
!
! Description:
! ------------
!> Land surface model, consisting of a solver for the energy balance at the
!> surface and a four layer soil scheme. The scheme is similar to the TESSEL
!> scheme implemented in the ECMWF IFS model, with modifications according to
!> H-TESSEL. The implementation is based on the formulation implemented in the 
!> DALES and UCLA-LES models.
!>
!> @todo Consider partial absorption of the net shortwave radiation by the 
!>       skin layer.
!> @todo Improve surface water parameterization
!> @todo Invert indices (running from -3 to 0. Currently: nzb_soil=0, 
!>       nzt_soil=3)).
!> @todo Implement surface runoff model (required when performing long-term LES 
!>       with considerable precipitation.
!> @todo Fix crashes with radiation_scheme == 'constant'
!>
!> @note No time step criterion is required as long as the soil layers do not
!>       become too thin.
!------------------------------------------------------------------------------!
 MODULE land_surface_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  hyp, ol, pt, pt_p, prr, q, q_p, ql, qsws, shf, ts, us, vpt, z0, &
               z0h, z0q

    USE cloud_parameters,                                                      &
        ONLY:  cp, hyrho, l_d_cp, l_d_r, l_v, pt_d_t, rho_l, r_d, r_v

    USE control_parameters,                                                    &
        ONLY:  cloud_physics, dt_3d, humidity, intermediate_timestep_count,    &
               initializing_actions, intermediate_timestep_count_max,          &
               max_masks, precipitation, pt_surface,                           &
               rho_surface, roughness_length, surface_pressure,                &
               timestep_scheme, tsc, z0h_factor, time_since_reference_point

    USE indices,                                                               &
        ONLY:  nbgp, nxlg, nxrg, nyng, nysg, nzb, nzb_s_inner 

    USE kinds

    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  force_radiation_call, rad_net, rad_sw_in, rad_lw_out,           &
               rad_lw_out_change_0, unscheduled_radiation_calls
        
    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    IMPLICIT NONE

!
!-- LSM model constants
    INTEGER(iwp), PARAMETER :: nzb_soil = 0, & !< bottom of the soil model (to be switched)
                               nzt_soil = 3, & !< top of the soil model (to be switched)
                               nzs = 4         !< number of soil layers (fixed for now)

    REAL(wp), PARAMETER ::                     &
              b_ch               = 6.04_wp,    & ! Clapp & Hornberger exponent
              lambda_h_dry       = 0.19_wp,    & ! heat conductivity for dry soil   
              lambda_h_sm        = 3.44_wp,    & ! heat conductivity of the soil matrix
              lambda_h_water     = 0.57_wp,    & ! heat conductivity of water
              psi_sat            = -0.388_wp,  & ! soil matrix potential at saturation
              rho_c_soil         = 2.19E6_wp,  & ! volumetric heat capacity of soil
              rho_c_water        = 4.20E6_wp,  & ! volumetric heat capacity of water
              m_max_depth        = 0.0002_wp     ! Maximum capacity of the water reservoir (m)


!
!-- LSM variables
    INTEGER(iwp) :: veg_type  = 2, & !< NAMELIST veg_type_2d
                    soil_type = 3    !< NAMELIST soil_type_2d

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  soil_type_2d, &  !< soil type, 0: user-defined, 1-7: generic (see list)
                                                  veg_type_2d      !< vegetation type, 0: user-defined, 1-19: generic (see list)

    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: water_surface,     & !< flag parameter for water surfaces (classes 14+15)
                                            pave_surface,      & !< flag parameter for pavements (asphalt etc.) (class 20)
                                            building_surface     !< flag parameter indicating that the surface element is covered by buildings (no LSM actions, not implemented yet)

    LOGICAL :: conserve_water_content = .TRUE.,  & !< open or closed bottom surface for the soil model
               force_radiation_call_l = .FALSE., & !< flag parameter for unscheduled radiation model calls
               land_surface = .FALSE.              !< flag parameter indicating wheather the lsm is used

!   value 9999999.9_wp -> generic available or user-defined value must be set
!   otherwise -> no generic variable and user setting is optional
    REAL(wp) :: alpha_vangenuchten = 9999999.9_wp,      & !< NAMELIST alpha_vg
                canopy_resistance_coefficient = 9999999.9_wp, & !< NAMELIST g_d
                c_surface   = 20000.0_wp,               & !< Surface (skin) heat capacity
                drho_l_lv,                              & !< (rho_l * l_v)**-1
                exn,                                    & !< value of the Exner function
                e_s = 0.0_wp,                           & !< saturation water vapour pressure
                field_capacity = 9999999.9_wp,          & !< NAMELIST m_fc
                f_shortwave_incoming = 9999999.9_wp,    & !< NAMELIST f_sw_in
                hydraulic_conductivity = 9999999.9_wp,  & !< NAMELIST gamma_w_sat
                ke = 0.0_wp,                            & !< Kersten number
                lambda_h_sat = 0.0_wp,                  & !< heat conductivity for saturated soil
                lambda_surface_stable = 9999999.9_wp,   & !< NAMELIST lambda_surface_s
                lambda_surface_unstable = 9999999.9_wp, & !< NAMELIST lambda_surface_u
                leaf_area_index = 9999999.9_wp,         & !< NAMELIST lai
                l_vangenuchten = 9999999.9_wp,          & !< NAMELIST l_vg
                min_canopy_resistance = 9999999.9_wp,   & !< NAMELIST r_canopy_min
                min_soil_resistance = 50.0_wp,          & !< NAMELIST r_soil_min
                m_total = 0.0_wp,                       & !< weighted total water content of the soil (m3/m3)
                n_vangenuchten = 9999999.9_wp,          & !< NAMELIST n_vg
                pave_depth = 9999999.9_wp,              & !< depth of the pavement
                pave_heat_capacity = 1.94E6_wp,         & !< volumetric heat capacity of pavement (e.g. roads)
                pave_heat_conductivity = 1.00_wp,       & !< heat conductivity for pavements (e.g. roads)
                q_s = 0.0_wp,                           & !< saturation specific humidity
                residual_moisture = 9999999.9_wp,       & !< NAMELIST m_res
                rho_cp,                                 & !< rho_surface * cp
                rho_lv,                                 & !< rho_ocean * l_v
                rd_d_rv,                                & !< r_d / r_v
                saturation_moisture = 9999999.9_wp,     & !< NAMELIST m_sat
                skip_time_do_lsm = 0.0_wp,              & !< LSM is not called before this time
                vegetation_coverage = 9999999.9_wp,     & !< NAMELIST c_veg
                wilting_point = 9999999.9_wp,           & !< NAMELIST m_wilt
                z0_eb  = 9999999.9_wp,                  & !< NAMELIST z0 (lsm_par)
                z0h_eb = 9999999.9_wp,                  & !< NAMELIST z0h (lsm_par)
                z0q_eb = 9999999.9_wp                     !< NAMELIST z0q (lsm_par)

    REAL(wp), DIMENSION(nzb_soil:nzt_soil) :: &
              ddz_soil_stag,                  & !< 1/dz_soil_stag
              dz_soil_stag,                   & !< soil grid spacing (center-center)
              root_extr = 0.0_wp,             & !< root extraction
              root_fraction = (/9999999.9_wp, 9999999.9_wp,    &
                                9999999.9_wp, 9999999.9_wp /), & !< distribution of root surface area to the individual soil layers
              zs = (/0.07_wp, 0.28_wp, 1.00_wp,  2.89_wp/),    & !< soil layer depths (m)
              soil_moisture = 0.0_wp          !< soil moisture content (m3/m3)

    REAL(wp), DIMENSION(nzb_soil:nzt_soil+1) ::   &
              soil_temperature = (/290.0_wp, 287.0_wp, 285.0_wp,  283.0_wp,    & !< soil temperature (K)
                                   283.0_wp /),                                &                                   
              ddz_soil,                                                        & !< 1/dz_soil
              dz_soil                                                            !< soil grid spacing (edge-edge)

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: t_surface,   & !< surface temperature (K)
                                                     t_surface_p, & !< progn. surface temperature (K)
                                                     m_liq_eb,    & !< liquid water reservoir (m)
                                                     m_liq_eb_av, & !< liquid water reservoir (m)
                                                     m_liq_eb_p     !< progn. liquid water reservoir (m)
#else
    REAL(wp), DIMENSION(:,:), POINTER :: t_surface,      &
                                         t_surface_p,    & 
                                         m_liq_eb,       & 
                                         m_liq_eb_p

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: t_surface_1, t_surface_2, &
                                                     m_liq_eb_av,              &
                                                     m_liq_eb_1, m_liq_eb_2
#endif

!
!-- Temporal tendencies for time stepping            
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: tt_surface_m,  & !< surface temperature tendency (K)
                                             tm_liq_eb_m      !< liquid water reservoir tendency (m)

!
!-- Energy balance variables                
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
              alpha_vg,         & !< coef. of Van Genuchten
              c_liq,            & !< liquid water coverage (of vegetated area)
              c_liq_av,         & !< average of c_liq
              c_soil_av,        & !< average of c_soil
              c_veg,            & !< vegetation coverage
              c_veg_av,         & !< average of c_veg
              f_sw_in,          & !< fraction of absorbed shortwave radiation by the surface layer (not implemented yet)
              ghf_eb,           & !< ground heat flux
              ghf_eb_av,        & !< average of ghf_eb
              gamma_w_sat,      & !< hydraulic conductivity at saturation
              g_d,              & !< coefficient for dependence of r_canopy on water vapour pressure deficit
              lai,              & !< leaf area index
              lai_av,           & !< average of lai
              lambda_surface_s, & !< coupling between surface and soil (depends on vegetation type)
              lambda_surface_u, & !< coupling between surface and soil (depends on vegetation type)
              l_vg,             & !< coef. of Van Genuchten
              m_fc,             & !< soil moisture at field capacity (m3/m3)
              m_res,            & !< residual soil moisture
              m_sat,            & !< saturation soil moisture (m3/m3)
              m_wilt,           & !< soil moisture at permanent wilting point (m3/m3)
              n_vg,             & !< coef. Van Genuchten  
              qsws_eb,          & !< surface flux of latent heat (total)
              qsws_eb_av,       & !< average of qsws_eb
              qsws_liq_eb,      & !< surface flux of latent heat (liquid water portion)
              qsws_liq_eb_av,   & !< average of qsws_liq_eb
              qsws_soil_eb,     & !< surface flux of latent heat (soil portion)
              qsws_soil_eb_av,  & !< average of qsws_soil_eb
              qsws_veg_eb,      & !< surface flux of latent heat (vegetation portion)
              qsws_veg_eb_av,   & !< average of qsws_veg_eb
              rad_net_l,        & !< local copy of rad_net (net radiation at surface)
              r_a,              & !< aerodynamic resistance 
              r_a_av,           & !< average of r_a
              r_canopy,         & !< canopy resistance
              r_soil,           & !< soil resistance
              r_soil_min,       & !< minimum soil resistance
              r_s,              & !< total surface resistance (combination of r_soil and r_canopy)
              r_s_av,           & !< average of r_s
              r_canopy_min,     & !< minimum canopy (stomatal) resistance
              shf_eb,           & !< surface flux of sensible heat
              shf_eb_av           !< average of shf_eb


    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::                                 &
              lambda_h, &   !< heat conductivity of soil (W/m/K)                            
              lambda_w, &   !< hydraulic diffusivity of soil (?)
              gamma_w,  &   !< hydraulic conductivity of soil (W/m/K)
              rho_c_total   !< volumetric heat capacity of the actual soil matrix (?) 

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                         &
              t_soil,    & !< Soil temperature (K)
              t_soil_av, & !< Average of t_soil
              t_soil_p,  & !< Prog. soil temperature (K)
              m_soil,    & !< Soil moisture (m3/m3)
              m_soil_av, & !< Average of m_soil
              m_soil_p     !< Prog. soil moisture (m3/m3)
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::                                     &
              t_soil, t_soil_p, &
              m_soil, m_soil_p    

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                         &
              t_soil_av, t_soil_1, t_soil_2,                                   &
              m_soil_av, m_soil_1, m_soil_2
#endif


    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::                                 &
              tt_soil_m, & !< t_soil storage array
              tm_soil_m, & !< m_soil storage array
              root_fr      !< root fraction (sum=1) 


!
!-- Predefined Land surface classes (veg_type)
    CHARACTER(26), DIMENSION(0:20), PARAMETER :: veg_type_name = (/ &
                                   'user defined              ',    & ! 0 
                                   'crops, mixed farming      ',    & !  1
                                   'short grass               ',    & !  2
                                   'evergreen needleleaf trees',    & !  3
                                   'deciduous needleleaf trees',    & !  4
                                   'evergreen broadleaf trees ',    & !  5
                                   'deciduous broadleaf trees ',    & !  6
                                   'tall grass                ',    & !  7
                                   'desert                    ',    & !  8
                                   'tundra                    ',    & !  9
                                   'irrigated crops           ',    & ! 10
                                   'semidesert                ',    & ! 11
                                   'ice caps and glaciers     ',    & ! 12
                                   'bogs and marshes          ',    & ! 13
                                   'inland water              ',    & ! 14
                                   'ocean                     ',    & ! 15
                                   'evergreen shrubs          ',    & ! 16
                                   'deciduous shrubs          ',    & ! 17
                                   'mixed forest/woodland     ',    & ! 18
                                   'interrupted forest        ',    & ! 19
                                   'pavements/roads           '     & ! 20
                                                                 /)

!
!-- Soil model classes (soil_type)
    CHARACTER(12), DIMENSION(0:7), PARAMETER :: soil_type_name = (/ &
                                   'user defined',                  & ! 0 
                                   'coarse      ',                  & ! 1
                                   'medium      ',                  & ! 2
                                   'medium-fine ',                  & ! 3
                                   'fine        ',                  & ! 4
                                   'very fine   ',                  & ! 5
                                   'organic     ',                  & ! 6
                                   'loamy (CH)  '                   & ! 7
                                                                 /)
!
!-- Land surface parameters according to the respective classes (veg_type)

!
!-- Land surface parameters I
!--                          r_canopy_min,     lai,   c_veg,     g_d
    REAL(wp), DIMENSION(0:3,1:20), PARAMETER :: veg_pars = RESHAPE( (/ &
                                 180.0_wp, 3.00_wp, 0.90_wp, 0.00_wp,  & !  1
                                 110.0_wp, 2.00_wp, 0.85_wp, 0.00_wp,  & !  2
                                 500.0_wp, 5.00_wp, 0.90_wp, 0.03_wp,  & !  3
                                 500.0_wp, 5.00_wp, 0.90_wp, 0.03_wp,  & !  4
                                 175.0_wp, 5.00_wp, 0.90_wp, 0.03_wp,  & !  5
                                 240.0_wp, 6.00_wp, 0.99_wp, 0.13_wp,  & !  6
                                 100.0_wp, 2.00_wp, 0.70_wp, 0.00_wp,  & !  7
                                 250.0_wp, 0.05_wp, 0.00_wp, 0.00_wp,  & !  8
                                  80.0_wp, 1.00_wp, 0.50_wp, 0.00_wp,  & !  9
                                 180.0_wp, 3.00_wp, 0.90_wp, 0.00_wp,  & ! 10
                                 150.0_wp, 0.50_wp, 0.10_wp, 0.00_wp,  & ! 11
                                   0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp,  & ! 12
                                 240.0_wp, 4.00_wp, 0.60_wp, 0.00_wp,  & ! 13
                                   0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp,  & ! 14
                                   0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp,  & ! 15
                                 225.0_wp, 3.00_wp, 0.50_wp, 0.00_wp,  & ! 16
                                 225.0_wp, 1.50_wp, 0.50_wp, 0.00_wp,  & ! 17
                                 250.0_wp, 5.00_wp, 0.90_wp, 0.03_wp,  & ! 18
                                 175.0_wp, 2.50_wp, 0.90_wp, 0.03_wp,  & ! 19
                                   0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp   & ! 20
                                 /), (/ 4, 20 /) )

!
!-- Land surface parameters II          z0,         z0h,         z0q
    REAL(wp), DIMENSION(0:2,1:20), PARAMETER :: roughness_par = RESHAPE( (/ & 
                                   0.25_wp,  0.25E-2_wp,  0.25E-2_wp,       & !  1
                                   0.20_wp,  0.20E-2_wp,  0.20E-2_wp,       & !  2
                                   2.00_wp,     2.00_wp,     2.00_wp,       & !  3
                                   2.00_wp,     2.00_wp,     2.00_wp,       & !  4
                                   2.00_wp,     2.00_wp,     2.00_wp,       & !  5
                                   2.00_wp,     2.00_wp,     2.00_wp,       & !  6
                                   0.47_wp,  0.47E-2_wp,  0.47E-2_wp,       & !  7
                                  0.013_wp, 0.013E-2_wp, 0.013E-2_wp,       & !  8
                                  0.034_wp, 0.034E-2_wp, 0.034E-2_wp,       & !  9
                                    0.5_wp,  0.50E-2_wp,  0.50E-2_wp,       & ! 10
                                   0.17_wp,  0.17E-2_wp,  0.17E-2_wp,       & ! 11
                                 1.3E-3_wp,   1.3E-4_wp,   1.3E-4_wp,       & ! 12
                                   0.83_wp,  0.83E-2_wp,  0.83E-2_wp,       & ! 13
                                   0.00_wp,     0.00_wp,     0.00_wp,       & ! 14
                                   0.00_wp,     0.00_wp,     0.00_wp,       & ! 15
                                   0.10_wp,  0.10E-2_wp,  0.10E-2_wp,       & ! 16
                                   0.25_wp,  0.25E-2_wp,  0.25E-2_wp,       & ! 17
                                   2.00_wp,  2.00E-2_wp,  2.00E-2_wp,       & ! 18
                                   1.10_wp,  1.10E-2_wp,  1.10E-2_wp,       & ! 19
                                 1.0E-4_wp,   1.0E-5_wp,   1.0E-5_wp        & ! 20
                                 /), (/ 3, 20 /) )

!
!-- Land surface parameters III lambda_surface_s, lambda_surface_u, f_sw_in
    REAL(wp), DIMENSION(0:2,1:20), PARAMETER :: surface_pars = RESHAPE( (/ &
                                      10.0_wp,       10.0_wp, 0.05_wp,     & !  1
                                      10.0_wp,       10.0_wp, 0.05_wp,     & !  2
                                      20.0_wp,       15.0_wp, 0.03_wp,     & !  3
                                      20.0_wp,       15.0_wp, 0.03_wp,     & !  4
                                      20.0_wp,       15.0_wp, 0.03_wp,     & !  5
                                      20.0_wp,       15.0_wp, 0.03_wp,     & !  6
                                      10.0_wp,       10.0_wp, 0.05_wp,     & !  7
                                      15.0_wp,       15.0_wp, 0.00_wp,     & !  8
                                      10.0_wp,       10.0_wp, 0.05_wp,     & !  9
                                      10.0_wp,       10.0_wp, 0.05_wp,     & ! 10 
                                      10.0_wp,       10.0_wp, 0.05_wp,     & ! 11
                                      58.0_wp,       58.0_wp, 0.00_wp,     & ! 12
                                      10.0_wp,       10.0_wp, 0.05_wp,     & ! 13 
                                    1.0E10_wp,     1.0E10_wp, 0.00_wp,     & ! 14
                                    1.0E10_wp,     1.0E10_wp, 0.00_wp,     & ! 15
                                      10.0_wp,       10.0_wp, 0.05_wp,     & ! 16
                                      10.0_wp,       10.0_wp, 0.05_wp,     & ! 17
                                      20.0_wp,       15.0_wp, 0.03_wp,     & ! 18
                                      20.0_wp,       15.0_wp, 0.03_wp,     & ! 19
                                       0.0_wp,        0.0_wp, 0.00_wp      & ! 20
                                      /), (/ 3, 20 /) )

!
!-- Root distribution (sum = 1)  level 1, level 2, level 3, level 4,
    REAL(wp), DIMENSION(0:3,1:20), PARAMETER :: root_distribution = RESHAPE( (/ &
                                 0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,            & !  1
                                 0.35_wp, 0.38_wp, 0.23_wp, 0.04_wp,            & !  2
                                 0.26_wp, 0.39_wp, 0.29_wp, 0.06_wp,            & !  3
                                 0.26_wp, 0.38_wp, 0.29_wp, 0.07_wp,            & !  4
                                 0.24_wp, 0.38_wp, 0.31_wp, 0.07_wp,            & !  5
                                 0.25_wp, 0.34_wp, 0.27_wp, 0.14_wp,            & !  6
                                 0.27_wp, 0.27_wp, 0.27_wp, 0.09_wp,            & !  7
                                 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & !  8
                                 0.47_wp, 0.45_wp, 0.08_wp, 0.00_wp,            & !  9
                                 0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,            & ! 10
                                 0.17_wp, 0.31_wp, 0.33_wp, 0.19_wp,            & ! 11
                                 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & ! 12
                                 0.25_wp, 0.34_wp, 0.27_wp, 0.11_wp,            & ! 13
                                 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & ! 14
                                 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & ! 15
                                 0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,            & ! 16 
                                 0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,            & ! 17 
                                 0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp,            & ! 18
                                 0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp,            & ! 19
                                 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp             & ! 20
                                 /), (/ 4, 20 /) )

!
!-- Soil parameters according to the following porosity classes (soil_type)

!
!-- Soil parameters I           alpha_vg,      l_vg,    n_vg, gamma_w_sat
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER :: soil_pars = RESHAPE( (/     &
                                 3.83_wp,  1.250_wp, 1.38_wp,  6.94E-6_wp, & ! 1
                                 3.14_wp, -2.342_wp, 1.28_wp,  1.16E-6_wp, & ! 2
                                 0.83_wp, -0.588_wp, 1.25_wp,  0.26E-6_wp, & ! 3
                                 3.67_wp, -1.977_wp, 1.10_wp,  2.87E-6_wp, & ! 4
                                 2.65_wp,  2.500_wp, 1.10_wp,  1.74E-6_wp, & ! 5
                                 1.30_wp,  0.400_wp, 1.20_wp,  0.93E-6_wp, & ! 6
                                 0.00_wp,  0.00_wp,  0.00_wp,  0.57E-6_wp  & ! 7
                                 /), (/ 4, 7 /) )

!
!-- Soil parameters II              m_sat,     m_fc,   m_wilt,    m_res  
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER :: m_soil_pars = RESHAPE( (/ &
                                 0.403_wp, 0.244_wp, 0.059_wp, 0.025_wp, & ! 1
                                 0.439_wp, 0.347_wp, 0.151_wp, 0.010_wp, & ! 2
                                 0.430_wp, 0.383_wp, 0.133_wp, 0.010_wp, & ! 3
                                 0.520_wp, 0.448_wp, 0.279_wp, 0.010_wp, & ! 4
                                 0.614_wp, 0.541_wp, 0.335_wp, 0.010_wp, & ! 5
                                 0.766_wp, 0.663_wp, 0.267_wp, 0.010_wp, & ! 6
                                 0.472_wp, 0.323_wp, 0.171_wp, 0.000_wp  & ! 7
                                 /), (/ 4, 7 /) )


    SAVE


    PRIVATE

    
!
!-- Public functions
    PUBLIC lsm_check_data_output, lsm_check_data_output_pr,                    &
           lsm_check_parameters, lsm_define_netcdf_grid, lsm_3d_data_averaging,& 
           lsm_data_output_2d, lsm_data_output_3d, lsm_energy_balance,         &
           lsm_header, lsm_init, lsm_init_arrays, lsm_parin, lsm_soil_model,   &
           lsm_swap_timelevel, lsm_read_restart_data, lsm_last_actions
!
!-- Public parameters, constants and initial values
    PUBLIC land_surface, skip_time_do_lsm

!
!-- Public grid variables
    PUBLIC nzb_soil, nzs, nzt_soil, zs

!
!-- Public 2D output variables
    PUBLIC ghf_eb, qsws_eb, qsws_liq_eb, qsws_soil_eb,qsws_veg_eb, r_a, r_s,   &
           shf_eb

!
!-- Public prognostic variables
    PUBLIC m_soil, t_soil


    INTERFACE lsm_check_data_output
       MODULE PROCEDURE lsm_check_data_output
    END INTERFACE lsm_check_data_output
    
    INTERFACE lsm_check_data_output_pr
       MODULE PROCEDURE lsm_check_data_output_pr
    END INTERFACE lsm_check_data_output_pr
    
    INTERFACE lsm_check_parameters
       MODULE PROCEDURE lsm_check_parameters
    END INTERFACE lsm_check_parameters
    
    INTERFACE lsm_3d_data_averaging
       MODULE PROCEDURE lsm_3d_data_averaging
    END INTERFACE lsm_3d_data_averaging

    INTERFACE lsm_data_output_2d
       MODULE PROCEDURE lsm_data_output_2d
    END INTERFACE lsm_data_output_2d

    INTERFACE lsm_data_output_3d
       MODULE PROCEDURE lsm_data_output_3d
    END INTERFACE lsm_data_output_3d

    INTERFACE lsm_define_netcdf_grid
       MODULE PROCEDURE lsm_define_netcdf_grid
    END INTERFACE lsm_define_netcdf_grid

    INTERFACE lsm_energy_balance
       MODULE PROCEDURE lsm_energy_balance
    END INTERFACE lsm_energy_balance

    INTERFACE lsm_header
       MODULE PROCEDURE lsm_header
    END INTERFACE lsm_header
    
    INTERFACE lsm_init
       MODULE PROCEDURE lsm_init
    END INTERFACE lsm_init

    INTERFACE lsm_init_arrays
       MODULE PROCEDURE lsm_init_arrays
    END INTERFACE lsm_init_arrays
    
    INTERFACE lsm_parin
       MODULE PROCEDURE lsm_parin
    END INTERFACE lsm_parin
    
    INTERFACE lsm_soil_model
       MODULE PROCEDURE lsm_soil_model
    END INTERFACE lsm_soil_model

    INTERFACE lsm_swap_timelevel
       MODULE PROCEDURE lsm_swap_timelevel
    END INTERFACE lsm_swap_timelevel

    INTERFACE lsm_read_restart_data
       MODULE PROCEDURE lsm_read_restart_data
    END INTERFACE lsm_read_restart_data

    INTERFACE lsm_last_actions
       MODULE PROCEDURE lsm_last_actions
    END INTERFACE lsm_last_actions

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_data_output( var, unit, i, ilen, k )
 
 
    USE control_parameters,                                                 &
        ONLY:  data_output, message_string

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< 
    CHARACTER (LEN=*) ::  var !<

    INTEGER(iwp) :: i
    INTEGER(iwp) :: ilen   
    INTEGER(iwp) :: k

    SELECT CASE ( TRIM( var ) )

       CASE ( 'm_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                      'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'm3/m3'
           
       CASE ( 't_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                      'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'K'   
             
       CASE ( 'lai*', 'c_liq*', 'c_soil*', 'c_veg*', 'ghf_eb*', 'm_liq_eb*',&
              'qsws_eb*', 'qsws_liq_eb*', 'qsws_soil_eb*', 'qsws_veg_eb*',  &
              'r_a*', 'r_s*', 'shf_eb*' )
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //         &
                              TRIM( var ) // '" & only 2d-horizontal ' //   &
                              'cross sections are allowed for this value'
             CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'lai*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_liq*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_soil*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_veg*'  .AND.  .NOT. land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'ghf_eb*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'm_liq_eb*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_eb*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_liq_eb*'  .AND.  .NOT. land_surface )   &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_soil_eb*'  .AND.  .NOT.  land_surface ) &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_veg_eb*'  .AND.  .NOT. land_surface )   &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'r_a*'  .AND.  .NOT.  land_surface ) &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'r_s*'  .AND.  .NOT.  land_surface ) &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res land_surface = .TRUE.'
             CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( TRIM( var ) == 'lai*'   )  unit = 'none' 
          IF ( TRIM( var ) == 'c_liq*' )  unit = 'none'
          IF ( TRIM( var ) == 'c_soil*')  unit = 'none'
          IF ( TRIM( var ) == 'c_veg*' )  unit = 'none'
          IF ( TRIM( var ) == 'ghf_eb*')  unit = 'W/m2'
          IF ( TRIM( var ) == 'm_liq_eb*'     )  unit = 'm'
          IF ( TRIM( var ) == 'qsws_eb*'      ) unit = 'W/m2'
          IF ( TRIM( var ) == 'qsws_liq_eb*'  ) unit = 'W/m2'
          IF ( TRIM( var ) == 'qsws_soil_eb*' ) unit = 'W/m2'
          IF ( TRIM( var ) == 'qsws_veg_eb*'  ) unit = 'W/m2'
          IF ( TRIM( var ) == 'r_a*')     unit = 's/m'     
          IF ( TRIM( var ) == 'r_s*')     unit = 's/m' 
          IF ( TRIM( var ) == 'shf_eb*')  unit = 'W/m2'
             
       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE lsm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_data_output_pr( variable, var_count, unit, dopr_unit )
 
    USE control_parameters,                                                 &
        ONLY:  data_output_pr, message_string

    USE indices

    USE profil_parameter

    USE statistics

    IMPLICIT NONE
   
    CHARACTER (LEN=*) ::  unit      !< 
    CHARACTER (LEN=*) ::  variable  !< 
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
 
    INTEGER(iwp) ::  user_pr_index !< 
    INTEGER(iwp) ::  var_count     !< 

    SELECT CASE ( TRIM( variable ) )
       
       CASE ( 't_soil', '#t_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'data_output_pr = ' //                        &
                              TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not implemented for land_surface = .FALSE.'
             CALL message( 'check_parameters', 'PA0402', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 89
             dopr_unit     = 'K'
             hom(0:nzs-1,2,89,:)  = SPREAD( - zs, 2, statistic_regions+1 )
             IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
                dopr_initial_index(var_count) = 90
                hom(0:nzs-1,2,90,:)   = SPREAD( - zs, 2, statistic_regions+1 )
                data_output_pr(var_count)     = data_output_pr(var_count)(2:)
             ENDIF
             unit = dopr_unit
          ENDIF

       CASE ( 'm_soil', '#m_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'data_output_pr = ' //                        &
                              TRIM( data_output_pr(var_count) ) // ' is' // &
                              ' not implemented for land_surface = .FALSE.'
             CALL message( 'check_parameters', 'PA0402', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 91
             dopr_unit     = 'm3/m3'
             hom(0:nzs-1,2,91,:)  = SPREAD( - zs, 2, statistic_regions+1 )
             IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
                dopr_initial_index(var_count) = 92
                hom(0:nzs-1,2,92,:)   = SPREAD( - zs, 2, statistic_regions+1 )
                data_output_pr(var_count)     = data_output_pr(var_count)(2:)
             ENDIF
             unit = dopr_unit
          ENDIF


       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE lsm_check_data_output_pr
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_parameters

    USE control_parameters,                                                    &
        ONLY:  bc_pt_b, bc_q_b, constant_flux_layer, message_string,           &
               most_method, topography
                 
    USE radiation_model_mod,                                                   &
        ONLY:  radiation
    
    
    IMPLICIT NONE

 
!
!-- Dirichlet boundary conditions are required as the surface fluxes are
!-- calculated from the temperature/humidity gradients in the land surface
!-- model
    IF ( bc_pt_b == 'neumann'  .OR.  bc_q_b == 'neumann' )  THEN
       message_string = 'lsm requires setting of'//                         &
                        'bc_pt_b = "dirichlet" and '//                      &
                        'bc_q_b  = "dirichlet"'
       CALL message( 'check_parameters', 'PA0399', 1, 2, 0, 6, 0 )
    ENDIF

    IF (  .NOT.  constant_flux_layer )  THEN
       message_string = 'lsm requires '//                                   &
                        'constant_flux_layer = .T.'
       CALL message( 'check_parameters', 'PA0400', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( topography /= 'flat' )  THEN
       message_string = 'lsm cannot be used ' //                            & 
                        'in combination with  topography /= "flat"'
       CALL message( 'check_parameters', 'PA0415', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ( veg_type == 14  .OR.  veg_type == 15 ) .AND.                     &
           most_method == 'lookup' )  THEN
        WRITE( message_string, * ) 'veg_type = ', veg_type, ' is not ',     &
                                   'allowed in combination with ',          &
                                   'most_method = ', most_method
       CALL message( 'check_parameters', 'PA0417', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( veg_type == 0 )  THEN
       IF ( SUM( root_fraction ) /= 1.0_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of root_fraction(0:3)'//       &
                           '/= 9999999.9 and SUM(root_fraction) = 1'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF
 
       IF ( min_canopy_resistance == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user defined)'//                  &
                           'requires setting of min_canopy_resistance'//    &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( leaf_area_index == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of leaf_area_index'//          &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( vegetation_coverage == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of vegetation_coverage'//      &
                           '/= 9999999.9'
             CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( canopy_resistance_coefficient == 9999999.9_wp)  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of'//                          &
                           'canopy_resistance_coefficient /= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( lambda_surface_stable == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of lambda_surface_stable'//    &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( lambda_surface_unstable == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of lambda_surface_unstable'//  &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( f_shortwave_incoming == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of f_shortwave_incoming'//     &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( z0_eb == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of z0_eb'//                    &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( z0h_eb == 9999999.9_wp )  THEN
          message_string = 'veg_type = 0 (user_defined)'//                  &
                           'requires setting of z0h_eb'//                   &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF


    ENDIF

    IF ( soil_type == 0 )  THEN

       IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of alpha_vangenuchten'//       &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( l_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of l_vangenuchten'//           &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( n_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of n_vangenuchten'//           &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of hydraulic_conductivity'//   &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( saturation_moisture == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of saturation_moisture'//      &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( field_capacity == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of field_capacity'//           &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( wilting_point == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of wilting_point'//            &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( residual_moisture == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                 &
                           'requires setting of residual_moisture'//        &
                           '/= 9999999.9'
          CALL message( 'check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    IF (  .NOT.  radiation )  THEN
       message_string = 'lsm requires '//                                   &
                        'radiation = .T.'
       CALL message( 'check_parameters', 'PA0400', 1, 2, 0, 6, 0 )
    ENDIF
       
       
 END SUBROUTINE lsm_check_parameters
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the surface.
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_energy_balance


    IMPLICIT NONE

    INTEGER(iwp) ::  i         !< running index
    INTEGER(iwp) ::  j         !< running index
    INTEGER(iwp) ::  k, ks     !< running index

    REAL(wp) :: c_surface_tmp,& !< temporary variable for storing the volumetric heat capacity of the surface
                f1,          & !< resistance correction term 1
                f2,          & !< resistance correction term 2
                f3,          & !< resistance correction term 3
                m_min,       & !< minimum soil moisture
                e,           & !< water vapour pressure
                e_s,         & !< water vapour saturation pressure
                e_s_dt,      & !< derivate of e_s with respect to T
                tend,        & !< tendency
                dq_s_dt,     & !< derivate of q_s with respect to T
                coef_1,      & !< coef. for prognostic equation
                coef_2,      & !< coef. for prognostic equation
                f_qsws,      & !< factor for qsws_eb
                f_qsws_veg,  & !< factor for qsws_veg_eb
                f_qsws_soil, & !< factor for qsws_soil_eb
                f_qsws_liq,  & !< factor for qsws_liq_eb
                f_shf,       & !< factor for shf_eb
                lambda_surface, & !< Current value of lambda_surface
                m_liq_eb_max,   & !< maxmimum value of the liq. water reservoir
                pt1,         & !< potential temperature at first grid level
                qv1            !< specific humidity at first grid level

!
!-- Calculate the exner function for the current time step
    exn = ( surface_pressure / 1000.0_wp )**0.286_wp

    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          k = nzb_s_inner(j,i)

!
!--       Set lambda_surface according to stratification between skin layer and soil
          IF (  .NOT.  pave_surface(j,i) )  THEN

             c_surface_tmp = c_surface

             IF ( t_surface(j,i) >= t_soil(nzb_soil,j,i))  THEN
                lambda_surface = lambda_surface_s(j,i)
             ELSE
                lambda_surface = lambda_surface_u(j,i)
             ENDIF
          ELSE

             c_surface_tmp = pave_heat_capacity * dz_soil(nzb_soil) * 0.5_wp
             lambda_surface = pave_heat_conductivity * ddz_soil(nzb_soil)

          ENDIF

!
!--       First step: calculate aerodyamic resistance. As pt, us, ts
!--       are not available for the prognostic time step, data from the last
!--       time step is used here. Note that this formulation is the
!--       equivalent to the ECMWF formulation using drag coefficients
          IF ( cloud_physics )  THEN
             pt1 = pt(k+1,j,i) + l_d_cp * pt_d_t(k+1) * ql(k+1,j,i)
             qv1 = q(k+1,j,i) - ql(k+1,j,i)
          ELSE
             pt1 = pt(k+1,j,i)
             qv1 = q(k+1,j,i)
          ENDIF

          r_a(j,i) = (pt1 - pt(k,j,i)) / (ts(j,i) * us(j,i) + 1.0E-20_wp)

!
!--       Make sure that the resistance does not drop to zero for neutral 
!--       stratification
          IF ( ABS(r_a(j,i)) < 1.0_wp )  r_a(j,i) = 1.0_wp

!
!--       Second step: calculate canopy resistance r_canopy
!--       f1-f3 here are defined as 1/f1-f3 as in ECMWF documentation
 
!--       f1: correction for incoming shortwave radiation (stomata close at 
!--       night)
          f1 = MIN( 1.0_wp, ( 0.004_wp * rad_sw_in(k,j,i) + 0.05_wp ) /        &
                           (0.81_wp * (0.004_wp * rad_sw_in(k,j,i)             &
                            + 1.0_wp)) )



!
!--       f2: correction for soil moisture availability to plants (the 
!--       integrated soil moisture must thus be considered here)
!--       f2 = 0 for very dry soils
          m_total = 0.0_wp
          DO  ks = nzb_soil, nzt_soil
              m_total = m_total + root_fr(ks,j,i)                              &
                        * MAX(m_soil(ks,j,i),m_wilt(j,i))
          ENDDO 

          IF ( m_total > m_wilt(j,i)  .AND.  m_total < m_fc(j,i) )  THEN
             f2 = ( m_total - m_wilt(j,i) ) / (m_fc(j,i) - m_wilt(j,i) )
          ELSEIF ( m_total >= m_fc(j,i) )  THEN
             f2 = 1.0_wp
          ELSE
             f2 = 1.0E-20_wp
          ENDIF

!
!--       Calculate water vapour pressure at saturation 
          e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp * ( t_surface(j,i)        &
                        - 273.16_wp ) / ( t_surface(j,i) - 35.86_wp ) )

!
!--       f3: correction for vapour pressure deficit
          IF ( g_d(j,i) /= 0.0_wp )  THEN
!
!--          Calculate vapour pressure
             e  = qv1 * surface_pressure / 0.622_wp
             f3 = EXP ( -g_d(j,i) * (e_s - e) )
          ELSE
             f3 = 1.0_wp
          ENDIF

!
!--       Calculate canopy resistance. In case that c_veg is 0 (bare soils),
!--       this calculation is obsolete, as r_canopy is not used below.
!--       To do: check for very dry soil -> r_canopy goes to infinity
          r_canopy(j,i) = r_canopy_min(j,i) / (lai(j,i) * f1 * f2 * f3         &
                                          + 1.0E-20_wp)

!
!--       Third step: calculate bare soil resistance r_soil. The Clapp & 
!--       Hornberger parametrization does not consider c_veg.
          IF ( soil_type_2d(j,i) /= 7 )  THEN
             m_min = c_veg(j,i) * m_wilt(j,i) + (1.0_wp - c_veg(j,i)) *        &
                     m_res(j,i)
          ELSE
             m_min = m_wilt(j,i)
          ENDIF

          f2 = ( m_soil(nzb_soil,j,i) - m_min ) / ( m_fc(j,i) - m_min )
          f2 = MAX(f2,1.0E-20_wp)
          f2 = MIN(f2,1.0_wp)

          r_soil(j,i) = r_soil_min(j,i) / f2

!
!--       Calculate the maximum possible liquid water amount on plants and
!--       bare surface. For vegetated surfaces, a maximum depth of 0.2 mm is
!--       assumed, while paved surfaces might hold up 1 mm of water. The 
!--       liquid water fraction for paved surfaces is calculated after 
!--       Noilhan & Planton (1989), while the ECMWF formulation is used for
!--       vegetated surfaces and bare soils.
          IF ( pave_surface(j,i) )  THEN
             m_liq_eb_max = m_max_depth * 5.0_wp
             c_liq(j,i) = MIN( 1.0_wp, (m_liq_eb(j,i) / m_liq_eb_max)**0.67 )
          ELSE
             m_liq_eb_max = m_max_depth * ( c_veg(j,i) * lai(j,i)              &
                            + (1.0_wp - c_veg(j,i)) )
             c_liq(j,i) = MIN( 1.0_wp, m_liq_eb(j,i) / m_liq_eb_max )
          ENDIF

!
!--       Calculate saturation specific humidity
          q_s = 0.622_wp * e_s / surface_pressure

!
!--       In case of dewfall, set evapotranspiration to zero
!--       All super-saturated water is then removed from the air
          IF ( humidity  .AND.  q_s <= qv1 )  THEN
             r_canopy(j,i) = 0.0_wp
             r_soil(j,i)   = 0.0_wp
          ENDIF

!
!--       Calculate coefficients for the total evapotranspiration 
!--       In case of water surface, set vegetation and soil fluxes to zero.
!--       For pavements, only evaporation of liquid water is possible.
          IF ( water_surface(j,i) )  THEN
             f_qsws_veg  = 0.0_wp
             f_qsws_soil = 0.0_wp
             f_qsws_liq  = rho_lv / r_a(j,i)
          ELSEIF ( pave_surface (j,i) )  THEN
             f_qsws_veg  = 0.0_wp
             f_qsws_soil = 0.0_wp
             f_qsws_liq  = rho_lv * c_liq(j,i) / r_a(j,i)
          ELSE
             f_qsws_veg  = rho_lv * c_veg(j,i) * (1.0_wp - c_liq(j,i))/        &
                           (r_a(j,i) + r_canopy(j,i))
             f_qsws_soil = rho_lv * (1.0_wp - c_veg(j,i)) / (r_a(j,i) +        &
                                                             r_soil(j,i))
             f_qsws_liq  = rho_lv * c_veg(j,i) * c_liq(j,i) / r_a(j,i)
          ENDIF
!
!--       If soil moisture is below wilting point, plants do no longer
!--       transpirate.
!           IF ( m_soil(k,j,i) < m_wilt(j,i) )  THEN
!              f_qsws_veg = 0.0_wp
!           ENDIF

          f_shf  = rho_cp / r_a(j,i)
          f_qsws = f_qsws_veg + f_qsws_soil + f_qsws_liq

!
!--       Calculate derivative of q_s for Taylor series expansion
          e_s_dt = e_s * ( 17.269_wp / (t_surface(j,i) - 35.86_wp) -           &
                           17.269_wp*(t_surface(j,i) - 273.16_wp)              &
                           / (t_surface(j,i) - 35.86_wp)**2 )

          dq_s_dt = 0.622_wp * e_s_dt / surface_pressure

!
!--       Add LW up so that it can be removed in prognostic equation
          rad_net_l(j,i) = rad_net(j,i) + rad_lw_out(nzb,j,i)

!
!--       Calculate new skin temperature
          IF ( humidity )  THEN

!
!--          Numerator of the prognostic equation
             coef_1 = rad_net_l(j,i) + rad_lw_out_change_0(j,i)                &
                      * t_surface(j,i) - rad_lw_out(nzb,j,i)                   &
                      + f_shf * pt1 + f_qsws * ( qv1 - q_s                     &
                      + dq_s_dt * t_surface(j,i) ) + lambda_surface            &
                      * t_soil(nzb_soil,j,i)

!
!--          Denominator of the prognostic equation
             coef_2 = rad_lw_out_change_0(j,i) + f_qsws * dq_s_dt              &
                      + lambda_surface + f_shf / exn
          ELSE

!
!--          Numerator of the prognostic equation
             coef_1 = rad_net_l(j,i) + rad_lw_out_change_0(j,i)                &
                      * t_surface(j,i) - rad_lw_out(nzb,j,i)                   &
                      + f_shf * pt1  + lambda_surface                          &
                      * t_soil(nzb_soil,j,i)

!
!--          Denominator of the prognostic equation
             coef_2 = rad_lw_out_change_0(j,i) + lambda_surface + f_shf / exn

          ENDIF

          tend = 0.0_wp

!
!--       Implicit solution when the surface layer has no heat capacity,
!--       otherwise use RK3 scheme.
          t_surface_p(j,i) = ( coef_1 * dt_3d * tsc(2) + c_surface_tmp *       &
                             t_surface(j,i) ) / ( c_surface_tmp + coef_2       &
                                * dt_3d * tsc(2) ) 

!
!--       Add RK3 term
          IF ( c_surface_tmp /= 0.0_wp )  THEN

             t_surface_p(j,i) = t_surface_p(j,i) + dt_3d * tsc(3)              &
                                * tt_surface_m(j,i)

!
!--          Calculate true tendency
             tend = (t_surface_p(j,i) - t_surface(j,i) - dt_3d * tsc(3)        &
                    * tt_surface_m(j,i)) / (dt_3d  * tsc(2))
!
!--          Calculate t_surface tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   tt_surface_m(j,i) = tend
                ELSEIF ( intermediate_timestep_count <                         &
                         intermediate_timestep_count_max )  THEN
                   tt_surface_m(j,i) = -9.5625_wp * tend + 5.3125_wp           &
                                       * tt_surface_m(j,i)
                ENDIF
             ENDIF
          ENDIF

!
!--       In case of fast changes in the skin temperature, it is possible to
!--       update the radiative fluxes independently from the prescribed
!--       radiation call frequency. This effectively prevents oscillations,
!--       especially when setting skip_time_do_radiation /= 0. The threshold
!--       value of 0.2 used here is just a first guess. This method should be
!--       revised in the future as tests have shown that the threshold is 
!--       often reached, when no oscillations would occur (causes immense
!--       computing time for the radiation code).
          IF ( ABS( t_surface_p(j,i) - t_surface(j,i) ) > 0.2_wp  .AND.        &
               unscheduled_radiation_calls )  THEN
             force_radiation_call_l = .TRUE.
          ENDIF

          pt(k,j,i) = t_surface_p(j,i) / exn

!
!--       Calculate fluxes
          rad_net_l(j,i)   = rad_net_l(j,i) + rad_lw_out_change_0(j,i)         &
                             * t_surface(j,i) - rad_lw_out(nzb,j,i)            &
                             - rad_lw_out_change_0(j,i) * t_surface_p(j,i)

          rad_net(j,i) = rad_net_l(j,i)
          rad_lw_out(nzb,j,i) = rad_lw_out(nzb,j,i) + rad_lw_out_change_0(j,i) &
                                * ( t_surface_p(j,i) - t_surface(j,i) )

          ghf_eb(j,i)    = lambda_surface * (t_surface_p(j,i)                  &
                           - t_soil(nzb_soil,j,i))

          shf_eb(j,i)    = - f_shf * ( pt1 - pt(k,j,i) )

          shf(j,i) = shf_eb(j,i) / rho_cp

          IF ( humidity )  THEN
             qsws_eb(j,i)  = - f_qsws    * ( qv1 - q_s + dq_s_dt               &
                             * t_surface(j,i) - dq_s_dt * t_surface_p(j,i) )

             qsws(j,i) = qsws_eb(j,i) / rho_lv

             qsws_veg_eb(j,i)  = - f_qsws_veg  * ( qv1 - q_s                   &
                                 + dq_s_dt * t_surface(j,i) - dq_s_dt          &
                                 * t_surface_p(j,i) )

             qsws_soil_eb(j,i) = - f_qsws_soil * ( qv1 - q_s                   &
                                 + dq_s_dt * t_surface(j,i) - dq_s_dt          &
                                 * t_surface_p(j,i) )

             qsws_liq_eb(j,i)  = - f_qsws_liq  * ( qv1 - q_s                   &
                                 + dq_s_dt * t_surface(j,i) - dq_s_dt          &
                                 * t_surface_p(j,i) )
          ENDIF

!
!--       Calculate the true surface resistance
          IF ( qsws_eb(j,i) == 0.0_wp )  THEN
             r_s(j,i) = 1.0E10_wp
          ELSE
             r_s(j,i) = - rho_lv * ( qv1 - q_s + dq_s_dt                       &
                        * t_surface(j,i) - dq_s_dt * t_surface_p(j,i) )        &
                        / qsws_eb(j,i) - r_a(j,i)
          ENDIF

!
!--       Calculate change in liquid water reservoir due to dew fall or 
!--       evaporation of liquid water
          IF ( humidity )  THEN
!
!--          If precipitation is activated, add rain water to qsws_liq_eb
!--          and qsws_soil_eb according the the vegetation coverage.
!--          precipitation_rate is given in mm.
             IF ( precipitation )  THEN

!
!--             Add precipitation to liquid water reservoir, if possible.
!--             Otherwise, add the water to soil. In case of
!--             pavements, the exceeding water amount is implicitely removed 
!--             as runoff as qsws_soil_eb is then not used in the soil model
                IF ( m_liq_eb(j,i) /= m_liq_eb_max )  THEN
                   qsws_liq_eb(j,i) = qsws_liq_eb(j,i)                         &
                                    + c_veg(j,i) * prr(k,j,i) * hyrho(k)       &
                                    * 0.001_wp * rho_l * l_v
                ELSE
                   qsws_soil_eb(j,i) = qsws_soil_eb(j,i)                       &
                                     + c_veg(j,i) * prr(k,j,i) * hyrho(k)      &
                                     * 0.001_wp * rho_l * l_v
                ENDIF

!--             Add precipitation to bare soil according to the bare soil
!--             coverage.
                qsws_soil_eb(j,i) = qsws_soil_eb(j,i) + (1.0_wp                &
                                    - c_veg(j,i)) * prr(k,j,i) * hyrho(k)      &
                                    * 0.001_wp * rho_l * l_v
             ENDIF

!
!--          If the air is saturated, check the reservoir water level
             IF ( qsws_eb(j,i) < 0.0_wp )  THEN

!
!--             Check if reservoir is full (avoid values > m_liq_eb_max)
!--             In that case, qsws_liq_eb goes to qsws_soil_eb. In this 
!--             case qsws_veg_eb is zero anyway (because c_liq = 1),       
!--             so that tend is zero and no further check is needed
                IF ( m_liq_eb(j,i) == m_liq_eb_max )  THEN
                   qsws_soil_eb(j,i) = qsws_soil_eb(j,i)                       &
                                        + qsws_liq_eb(j,i)

                   qsws_liq_eb(j,i)  = 0.0_wp
                ENDIF

!
!--             In case qsws_veg_eb becomes negative (unphysical behavior), 
!--             let the water enter the liquid water reservoir as dew on the
!--             plant
                IF ( qsws_veg_eb(j,i) < 0.0_wp )  THEN
                   qsws_liq_eb(j,i) = qsws_liq_eb(j,i) + qsws_veg_eb(j,i)
                   qsws_veg_eb(j,i) = 0.0_wp
                ENDIF
             ENDIF                    
  
             tend = - qsws_liq_eb(j,i) * drho_l_lv
             m_liq_eb_p(j,i) = m_liq_eb(j,i) + dt_3d * ( tsc(2) * tend         &
                                                + tsc(3) * tm_liq_eb_m(j,i) )

!
!--          Check if reservoir is overfull -> reduce to maximum
!--          (conservation of water is violated here)
             m_liq_eb_p(j,i) = MIN(m_liq_eb_p(j,i),m_liq_eb_max)

!
!--          Check if reservoir is empty (avoid values < 0.0)
!--          (conservation of water is violated here)
             m_liq_eb_p(j,i) = MAX(m_liq_eb_p(j,i),0.0_wp)


!
!--          Calculate m_liq_eb tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   tm_liq_eb_m(j,i) = tend
                ELSEIF ( intermediate_timestep_count <                         &
                         intermediate_timestep_count_max )  THEN
                   tm_liq_eb_m(j,i) = -9.5625_wp * tend + 5.3125_wp            &
                                    * tm_liq_eb_m(j,i)
                ENDIF
             ENDIF

          ENDIF

       ENDDO
    ENDDO

!
!-- Make a logical OR for all processes. Force radiation call if at
!-- least one processor reached the threshold change in skin temperature
    IF ( unscheduled_radiation_calls  .AND.  intermediate_timestep_count       &
         == intermediate_timestep_count_max-1 )  THEN
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( force_radiation_call_l, force_radiation_call,    &
                           1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#else
       force_radiation_call = force_radiation_call_l
#endif
       force_radiation_call_l = .FALSE.
    ENDIF

!
!-- Calculate surface specific humidity
    IF ( humidity )  THEN
       CALL calc_q_surface
    ENDIF

!
!-- Calculate new roughness lengths (for water surfaces only)
    CALL calc_z0_water_surface


 END SUBROUTINE lsm_energy_balance


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_header ( io )


       IMPLICIT NONE

       CHARACTER (LEN=86) ::  t_soil_chr          !< String for soil temperature profile
       CHARACTER (LEN=86) ::  roots_chr           !< String for root profile
       CHARACTER (LEN=86) ::  vertical_index_chr  !< String for the vertical index
       CHARACTER (LEN=86) ::  m_soil_chr          !< String for soil moisture
       CHARACTER (LEN=86) ::  soil_depth_chr      !< String for soil depth
       CHARACTER (LEN=10) ::  coor_chr            !< Temporary string
    
       INTEGER(iwp) ::  i                         !< Loop index over soil layers
 
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file
 
       t_soil_chr = ''
       m_soil_chr    = ''
       soil_depth_chr  = '' 
       roots_chr        = '' 
       vertical_index_chr   = ''

       i = 1
       DO i = nzb_soil, nzt_soil
          WRITE (coor_chr,'(F10.2,7X)') soil_temperature(i)
          t_soil_chr = TRIM( t_soil_chr ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)') soil_moisture(i)
          m_soil_chr = TRIM( m_soil_chr ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)')  - zs(i)
          soil_depth_chr = TRIM( soil_depth_chr ) // ' '  // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)')  root_fraction(i)
          roots_chr = TRIM( roots_chr ) // ' '  // TRIM( coor_chr )

          WRITE (coor_chr,'(I10,7X)')  i
          vertical_index_chr = TRIM( vertical_index_chr ) // ' '  //           &
                               TRIM( coor_chr )
       ENDDO

!
!--    Write land surface model header
       WRITE( io,  1 )
       IF ( conserve_water_content )  THEN
          WRITE( io, 2 )
       ELSE
          WRITE( io, 3 )
       ENDIF

       WRITE( io, 4 ) TRIM( veg_type_name(veg_type) ),                         &
                        TRIM (soil_type_name(soil_type) )
       WRITE( io, 5 ) TRIM( soil_depth_chr ), TRIM( t_soil_chr ),              &
                        TRIM( m_soil_chr ), TRIM( roots_chr ),                 &
                        TRIM( vertical_index_chr )

1   FORMAT (//' Land surface model information:'/                              &
              ' ------------------------------'/)
2   FORMAT ('    --> Soil bottom is closed (water content is conserved',       &
            ', default)')
3   FORMAT ('    --> Soil bottom is open (water content is not conserved)')         
4   FORMAT ('    --> Land surface type  : ',A,/                                &
            '    --> Soil porosity type : ',A)
5   FORMAT (/'    Initial soil temperature and moisture profile:'//            &
            '       Height:        ',A,'  m'/                                  &
            '       Temperature:   ',A,'  K'/                                  &
            '       Moisture:      ',A,'  m**3/m**3'/                          &
            '       Root fraction: ',A,'  '/                                   &
            '       Grid point:    ',A)

    END SUBROUTINE lsm_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_init
    

       IMPLICIT NONE

       INTEGER(iwp) ::  i !< running index
       INTEGER(iwp) ::  j !< running index
       INTEGER(iwp) ::  k !< running index

       REAL(wp) :: pt1   !< potential temperature at first grid level


!
!--    Calculate Exner function
       exn = ( surface_pressure / 1000.0_wp )**0.286_wp


!
!--    If no cloud physics is used, rho_surface has not been calculated before
       IF (  .NOT.  cloud_physics )  THEN
          rho_surface = surface_pressure * 100.0_wp / ( r_d * pt_surface * exn )
       ENDIF

!
!--    Calculate frequently used parameters
       rho_cp    = cp * rho_surface
       rd_d_rv   = r_d / r_v
       rho_lv    = rho_surface * l_v
       drho_l_lv = 1.0_wp / (rho_l * l_v)

!
!--    Set inital values for prognostic quantities
       tt_surface_m = 0.0_wp
       tt_soil_m    = 0.0_wp
       tm_soil_m    = 0.0_wp
       tm_liq_eb_m  = 0.0_wp
       c_liq        = 0.0_wp

       ghf_eb = 0.0_wp
       shf_eb = rho_cp * shf

       IF ( humidity )  THEN
          qsws_eb = rho_lv * qsws
       ELSE
          qsws_eb = 0.0_wp
       ENDIF

       qsws_liq_eb  = 0.0_wp
       qsws_soil_eb = 0.0_wp
       qsws_veg_eb  = 0.0_wp

       r_a        = 50.0_wp
       r_s        = 50.0_wp
       r_canopy   = 0.0_wp
       r_soil     = 0.0_wp

!
!--    Allocate 3D soil model arrays
       ALLOCATE ( root_fr(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( lambda_h(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( rho_c_total(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )

       lambda_h = 0.0_wp
!
!--    If required, allocate humidity-related variables for the soil model
       IF ( humidity )  THEN
          ALLOCATE ( lambda_w(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( gamma_w(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )   

          lambda_w = 0.0_wp 
       ENDIF

!
!--    Calculate grid spacings. Temperature and moisture are defined at
!--    the edges of the soil layers (_stag), whereas gradients/fluxes are defined
!--    at the centers
       dz_soil(nzb_soil) = zs(nzb_soil)

       DO  k = nzb_soil+1, nzt_soil
          dz_soil(k) = zs(k) - zs(k-1)
       ENDDO
       dz_soil(nzt_soil+1) = dz_soil(nzt_soil)

       DO  k = nzb_soil, nzt_soil-1
          dz_soil_stag(k) = 0.5_wp * (dz_soil(k+1) + dz_soil(k))
       ENDDO
       dz_soil_stag(nzt_soil) = dz_soil(nzt_soil)

       ddz_soil      = 1.0_wp / dz_soil
       ddz_soil_stag = 1.0_wp / dz_soil_stag

!
!--    Initialize standard soil types. It is possible to overwrite each
!--    parameter by setting the respecticy NAMELIST variable to a 
!--    value /= 9999999.9. 
       IF ( soil_type /= 0 )  THEN  
 
          IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
             alpha_vangenuchten = soil_pars(0,soil_type)
          ENDIF

          IF ( l_vangenuchten == 9999999.9_wp )  THEN
             l_vangenuchten = soil_pars(1,soil_type)
          ENDIF

          IF ( n_vangenuchten == 9999999.9_wp )  THEN
             n_vangenuchten = soil_pars(2,soil_type)            
          ENDIF

          IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
             hydraulic_conductivity = soil_pars(3,soil_type)            
          ENDIF

          IF ( saturation_moisture == 9999999.9_wp )  THEN
             saturation_moisture = m_soil_pars(0,soil_type)           
          ENDIF

          IF ( field_capacity == 9999999.9_wp )  THEN
             field_capacity = m_soil_pars(1,soil_type)           
          ENDIF

          IF ( wilting_point == 9999999.9_wp )  THEN
             wilting_point = m_soil_pars(2,soil_type)            
          ENDIF

          IF ( residual_moisture == 9999999.9_wp )  THEN
             residual_moisture = m_soil_pars(3,soil_type)       
          ENDIF

       ENDIF    

!
!--    Map values to the respective 2D arrays
       alpha_vg      = alpha_vangenuchten
       l_vg          = l_vangenuchten
       n_vg          = n_vangenuchten 
       gamma_w_sat   = hydraulic_conductivity
       m_sat         = saturation_moisture
       m_fc          = field_capacity
       m_wilt        = wilting_point
       m_res         = residual_moisture
       r_soil_min    = min_soil_resistance

!
!--    Initial run actions
       IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

          t_soil    = 0.0_wp
          m_liq_eb  = 0.0_wp
          m_soil    = 0.0_wp

!
!--       Map user settings of T and q for each soil layer
!--       (make sure that the soil moisture does not drop below the permanent
!--       wilting point) -> problems with devision by zero)
          DO  k = nzb_soil, nzt_soil
             t_soil(k,:,:)    = soil_temperature(k)
             m_soil(k,:,:)    = MAX(soil_moisture(k),m_wilt(:,:))
             soil_moisture(k) = MAX(soil_moisture(k),wilting_point)
          ENDDO
          t_soil(nzt_soil+1,:,:) = soil_temperature(nzt_soil+1)

!
!--       Calculate surface temperature
          t_surface   = pt_surface * exn

!
!--       Set artifical values for ts and us so that r_a has its initial value 
!--       for the first time step
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                k = nzb_s_inner(j,i)

                IF ( cloud_physics )  THEN
                   pt1 = pt(k+1,j,i) + l_d_cp * pt_d_t(k+1) * ql(k+1,j,i)
                ELSE
                   pt1 = pt(k+1,j,i)
                ENDIF

!
!--             Assure that r_a cannot be zero at model start
                IF ( pt1 == pt(k,j,i) )  pt1 = pt1 + 1.0E-10_wp

                us(j,i)  = 0.1_wp
                ts(j,i)  = (pt1 - pt(k,j,i)) / r_a(j,i)
                shf(j,i) = - us(j,i) * ts(j,i)
             ENDDO
          ENDDO

!
!--    Actions for restart runs
       ELSE

          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                k = nzb_s_inner(j,i)                
                t_surface(j,i) = pt(k,j,i) * exn
             ENDDO
          ENDDO

       ENDIF

       DO  k = nzb_soil, nzt_soil
          root_fr(k,:,:) = root_fraction(k)
       ENDDO

       IF ( veg_type /= 0 )  THEN
          IF ( min_canopy_resistance == 9999999.9_wp )  THEN
             min_canopy_resistance = veg_pars(0,veg_type)
          ENDIF
          IF ( leaf_area_index == 9999999.9_wp )  THEN
             leaf_area_index = veg_pars(1,veg_type)         
          ENDIF
          IF ( vegetation_coverage == 9999999.9_wp )  THEN
             vegetation_coverage = veg_pars(2,veg_type)     
          ENDIF
          IF ( canopy_resistance_coefficient == 9999999.9_wp )  THEN
              canopy_resistance_coefficient= veg_pars(3,veg_type)     
          ENDIF
          IF ( lambda_surface_stable == 9999999.9_wp )  THEN
             lambda_surface_stable = surface_pars(0,veg_type)         
          ENDIF
          IF ( lambda_surface_unstable == 9999999.9_wp )  THEN
             lambda_surface_unstable = surface_pars(1,veg_type)        
          ENDIF
          IF ( f_shortwave_incoming == 9999999.9_wp )  THEN
             f_shortwave_incoming = surface_pars(2,veg_type)        
          ENDIF
          IF ( z0_eb == 9999999.9_wp )  THEN
             roughness_length = roughness_par(0,veg_type) 
             z0_eb            = roughness_par(0,veg_type) 
          ENDIF
          IF ( z0h_eb == 9999999.9_wp )  THEN
             z0h_eb = roughness_par(1,veg_type)
          ENDIF
          IF ( z0q_eb == 9999999.9_wp )  THEN
             z0q_eb = roughness_par(2,veg_type)
          ENDIF
          z0h_factor = z0h_eb / ( z0_eb + 1.0E-20_wp )

          IF ( ANY( root_fraction == 9999999.9_wp ) )  THEN
             DO  k = nzb_soil, nzt_soil
                root_fr(k,:,:) = root_distribution(k,veg_type)
                root_fraction(k) = root_distribution(k,veg_type)
             ENDDO
          ENDIF

       ELSE

          IF ( z0_eb == 9999999.9_wp )  THEN
             z0_eb = roughness_length
          ENDIF
          IF ( z0h_eb == 9999999.9_wp )  THEN
             z0h_eb = z0_eb * z0h_factor
          ENDIF
          IF ( z0q_eb == 9999999.9_wp )  THEN
             z0q_eb = z0_eb * z0h_factor
          ENDIF

       ENDIF

!
!--    For surfaces covered with pavement, set depth of the pavement (with dry
!--    soil below). The depth must be greater than the first soil layer depth
       IF ( veg_type == 20 )  THEN
          IF ( pave_depth == 9999999.9_wp )  THEN
             pave_depth = zs(nzb_soil)  
          ELSE
             pave_depth = MAX( zs(nzb_soil), pave_depth )
          ENDIF
       ENDIF

!
!--    Map vegetation and soil types to 2D array to allow for heterogeneous 
!--    surfaces via user interface see below
       veg_type_2d = veg_type
       soil_type_2d = soil_type

!
!--    Map vegetation parameters to the respective 2D arrays
       r_canopy_min         = min_canopy_resistance
       lai                  = leaf_area_index
       c_veg                = vegetation_coverage
       g_d                  = canopy_resistance_coefficient
       lambda_surface_s     = lambda_surface_stable
       lambda_surface_u     = lambda_surface_unstable
       f_sw_in              = f_shortwave_incoming
       z0                   = z0_eb
       z0h                  = z0h_eb
       z0q                  = z0q_eb

!
!--    Possibly do user-defined actions (e.g. define heterogeneous land surface)
       CALL user_init_land_surface

!
!--    Set flag parameter if vegetation type was set to a water surface. Also
!--    set temperature to a constant value in all "soil" layers.
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             IF ( veg_type_2d(j,i) == 14  .OR.  veg_type_2d(j,i) == 15 )  THEN
                water_surface(j,i) = .TRUE.
                t_soil(:,j,i) = t_surface(j,i)
             ELSEIF ( veg_type_2d(j,i) == 20 )  THEN
                pave_surface(j,i) = .TRUE.
                m_soil(:,j,i) = 0.0_wp
             ENDIF

          ENDDO
       ENDDO

!
!--    Calculate new roughness lengths (for water surfaces only)
       CALL calc_z0_water_surface

       t_soil_p    = t_soil
       m_soil_p    = m_soil
       m_liq_eb_p  = m_liq_eb
       t_surface_p = t_surface



!--    Store initial profiles of t_soil and m_soil (assuming they are 
!--    horizontally homogeneous on this PE)
       hom(nzb_soil:nzt_soil,1,90,:)  = SPREAD( t_soil(nzb_soil:nzt_soil,      &
                                                nysg,nxlg), 2,                 &
                                                statistic_regions+1 )
       hom(nzb_soil:nzt_soil,1,92,:)  = SPREAD( m_soil(nzb_soil:nzt_soil,      &
                                                nysg,nxlg), 2,                 &
                                                statistic_regions+1 )

    END SUBROUTINE lsm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate land surface model arrays and define pointers
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_init_arrays
    

       IMPLICIT NONE

!
!--    Allocate surface and soil temperature / humidity
#if defined( __nopointer )
       ALLOCATE ( m_liq_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_liq_eb_p(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_soil(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_soil_p(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_surface(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_surface_p(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_soil(nzb_soil:nzt_soil+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_soil_p(nzb_soil:nzt_soil+1,nysg:nyng,nxlg:nxrg) )
#else
       ALLOCATE ( m_liq_eb_1(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_liq_eb_2(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_soil_1(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_soil_2(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_surface_1(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_surface_2(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_soil_1(nzb_soil:nzt_soil+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( t_soil_2(nzb_soil:nzt_soil+1,nysg:nyng,nxlg:nxrg) )
#endif

!
!--    Allocate intermediate timestep arrays
       ALLOCATE ( tm_liq_eb_m(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( tm_soil_m(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( tt_surface_m(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( tt_soil_m(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )

!
!--    Allocate 2D vegetation model arrays
       ALLOCATE ( alpha_vg(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( building_surface(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( c_liq(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( c_veg(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( f_sw_in(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( ghf_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( gamma_w_sat(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( g_d(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( lai(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( l_vg(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( lambda_surface_u(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( lambda_surface_s(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_fc(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_res(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_sat(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( m_wilt(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( n_vg(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( pave_surface(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( qsws_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( qsws_soil_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( qsws_liq_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( qsws_veg_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( rad_net_l(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_a(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_canopy(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_soil(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_soil_min(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_s(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( r_canopy_min(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( shf_eb(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( soil_type_2d(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( veg_type_2d(nysg:nyng,nxlg:nxrg) )
       ALLOCATE ( water_surface(nysg:nyng,nxlg:nxrg) )

       water_surface = .FALSE.
       pave_surface  = .FALSE.

#if ! defined( __nopointer )
!
!--    Initial assignment of the pointers
       t_soil    => t_soil_1;    t_soil_p    => t_soil_2
       t_surface => t_surface_1; t_surface_p => t_surface_2
       m_soil    => m_soil_1;    m_soil_p    => m_soil_2
       m_liq_eb  => m_liq_eb_1;  m_liq_eb_p  => m_liq_eb_2
#endif


    END SUBROUTINE lsm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &lsmpar for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_parin


       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /lsm_par/         alpha_vangenuchten, c_surface,               &
                                  canopy_resistance_coefficient,               &
                                  conserve_water_content,                      &
                                  f_shortwave_incoming, field_capacity,        & 
                                  hydraulic_conductivity,                      &
                                  lambda_surface_stable,                       &
                                  lambda_surface_unstable, leaf_area_index,    &
                                  l_vangenuchten, min_canopy_resistance,       &
                                  min_soil_resistance, n_vangenuchten,         &
                                  pave_depth, pave_heat_capacity,              &
                                  pave_heat_conductivity,                      &
                                  residual_moisture, root_fraction,            &
                                  saturation_moisture, skip_time_do_lsm,       &
                                  soil_moisture, soil_temperature, soil_type,  &
                                  vegetation_coverage, veg_type, wilting_point,& 
                                  zs, z0_eb, z0h_eb, z0q_eb
       
       line = ' '
       
!
!--    Try to find land surface model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&lsm_par' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, lsm_par )

!
!--    Set flag that indicates that the land surface model is switched on
       land_surface = .TRUE.

 10    CONTINUE
       

    END SUBROUTINE lsm_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Soil model as part of the land surface model. The model predicts soil
!> temperature and water content.
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_soil_model


       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index
       INTEGER(iwp) ::  j   !< running index
       INTEGER(iwp) ::  k   !< running index

       REAL(wp)     :: h_vg !< Van Genuchten coef. h

       REAL(wp), DIMENSION(nzb_soil:nzt_soil) :: gamma_temp,  & !< temp. gamma
                                                 lambda_temp, & !< temp. lambda
                                                 tend           !< tendency

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng

             IF ( pave_surface(j,i) )  THEN
                rho_c_total(nzb_soil,j,i) = pave_heat_capacity
                lambda_temp(nzb_soil)     = pave_heat_conductivity
             ENDIF

             IF (  .NOT.  water_surface(j,i) )  THEN
                DO  k = nzb_soil, nzt_soil


                   IF ( pave_surface(j,i)  .AND.  zs(k) <= pave_depth )  THEN
                    
                      rho_c_total(k,j,i) = pave_heat_capacity
                      lambda_temp(k)     = pave_heat_conductivity   

                   ELSE            
!
!--                   Calculate volumetric heat capacity of the soil, taking 
!--                   into account water content 
                      rho_c_total(k,j,i) = (rho_c_soil * (1.0_wp - m_sat(j,i)) &
                                           + rho_c_water * m_soil(k,j,i))

!
!--                   Calculate soil heat conductivity at the center of the soil
!--                   layers
                      lambda_h_sat = lambda_h_sm ** (1.0_wp - m_sat(j,i)) *    &
                                     lambda_h_water ** m_soil(k,j,i)

                      ke = 1.0_wp + LOG10(MAX(0.1_wp,m_soil(k,j,i)             &
                                                     / m_sat(j,i)))

                      lambda_temp(k) = ke * (lambda_h_sat - lambda_h_dry) +    &
                                       lambda_h_dry
                   ENDIF

                ENDDO

!
!--             Calculate soil heat conductivity (lambda_h) at the _stag level 
!--             using linear interpolation. For pavement surface, the
!--             true pavement depth is considered
                DO  k = nzb_soil, nzt_soil-1
                   IF ( pave_surface(j,i)  .AND.  zs(k)   < pave_depth         &
                                           .AND.  zs(k+1) > pave_depth )  THEN
                      lambda_h(k,j,i) = ( pave_depth - zs(k) ) / dz_soil(k+1)  &
                                        * lambda_temp(k)                       &
                                        + ( 1.0_wp - ( pave_depth - zs(k) )    &
                                        / dz_soil(k+1) ) * lambda_temp(k+1)
                   ELSE
                      lambda_h(k,j,i) = ( lambda_temp(k+1) + lambda_temp(k) )  &
                                        * 0.5_wp
                   ENDIF
                ENDDO
                lambda_h(nzt_soil,j,i) = lambda_temp(nzt_soil)




!
!--             Prognostic equation for soil temperature t_soil
                tend(:) = 0.0_wp

                tend(nzb_soil) = (1.0_wp/rho_c_total(nzb_soil,j,i)) *          &
                          ( lambda_h(nzb_soil,j,i) * ( t_soil(nzb_soil+1,j,i)  &
                            - t_soil(nzb_soil,j,i) ) * ddz_soil(nzb_soil+1)    &
                            + ghf_eb(j,i) ) * ddz_soil_stag(nzb_soil)

                DO  k = nzb_soil+1, nzt_soil
                   tend(k) = (1.0_wp/rho_c_total(k,j,i))                       &
                             * (   lambda_h(k,j,i)                             &
                                 * ( t_soil(k+1,j,i) - t_soil(k,j,i) )         &
                                 * ddz_soil(k+1)                               &
                                 - lambda_h(k-1,j,i)                           &
                                 * ( t_soil(k,j,i) - t_soil(k-1,j,i) )         &
                                 * ddz_soil(k)                                 &
                               ) * ddz_soil_stag(k)

                ENDDO

                t_soil_p(nzb_soil:nzt_soil,j,i) = t_soil(nzb_soil:nzt_soil,j,i)&
                                                  + dt_3d * ( tsc(2)           &
                                                  * tend(nzb_soil:nzt_soil)    & 
                                                  + tsc(3)                     &
                                                  * tt_soil_m(:,j,i) )   

!
!--             Calculate t_soil tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb_soil, nzt_soil
                         tt_soil_m(k,j,i) = tend(k)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count <                      &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb_soil, nzt_soil
                         tt_soil_m(k,j,i) = -9.5625_wp * tend(k) + 5.3125_wp   &
                                         * tt_soil_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF


                DO  k = nzb_soil, nzt_soil

!
!--                Calculate soil diffusivity at the center of the soil layers
                   lambda_temp(k) = (- b_ch * gamma_w_sat(j,i) * psi_sat       &
                                     / m_sat(j,i) ) * ( MAX( m_soil(k,j,i),    &
                                     m_wilt(j,i) ) / m_sat(j,i) )**(           &
                                     b_ch + 2.0_wp )

!
!--                Parametrization of Van Genuchten
                   IF ( soil_type /= 7 )  THEN
!
!--                   Calculate the hydraulic conductivity after Van Genuchten 
!--                   (1980)
                      h_vg = ( ( (m_res(j,i) - m_sat(j,i)) / ( m_res(j,i) -    &
                                 MAX( m_soil(k,j,i), m_wilt(j,i) ) ) )**(      &
                                 n_vg(j,i) / (n_vg(j,i) - 1.0_wp ) ) - 1.0_wp  &
                             )**( 1.0_wp / n_vg(j,i) ) / alpha_vg(j,i)


                      gamma_temp(k) = gamma_w_sat(j,i) * ( ( (1.0_wp +         &
                                      ( alpha_vg(j,i) * h_vg )**n_vg(j,i))**(  &
                                      1.0_wp - 1.0_wp / n_vg(j,i) ) - (        &
                                      alpha_vg(j,i) * h_vg )**( n_vg(j,i)      &
                                      - 1.0_wp) )**2 )                         &
                                      / ( ( 1.0_wp + ( alpha_vg(j,i) * h_vg    &
                                      )**n_vg(j,i) )**( ( 1.0_wp  - 1.0_wp     &
                                      / n_vg(j,i) ) *( l_vg(j,i) + 2.0_wp) ) )

!
!--                Parametrization of Clapp & Hornberger
                   ELSE
                      gamma_temp(k) = gamma_w_sat(j,i) * ( m_soil(k,j,i)       &
                                      / m_sat(j,i) )**(2.0_wp * b_ch + 3.0_wp)
                   ENDIF

                ENDDO

!
!--             Prognostic equation for soil moisture content. Only performed,
!--             when humidity is enabled in the atmosphere and the surface type
!--             is not pavement (implies dry soil below).
                IF ( humidity  .AND.  .NOT.  pave_surface(j,i) )  THEN
!
!--                Calculate soil diffusivity (lambda_w) at the _stag level 
!--                using linear interpolation. To do: replace this with
!--                ECMWF-IFS Eq. 8.81
                   DO  k = nzb_soil, nzt_soil-1
                     
                      lambda_w(k,j,i) = ( lambda_temp(k+1) + lambda_temp(k) )  &
                                        * 0.5_wp
                      gamma_w(k,j,i)  = ( gamma_temp(k+1) + gamma_temp(k) )    &
                                        * 0.5_wp

                   ENDDO

!
!
!--                In case of a closed bottom (= water content is conserved), 
!--                set hydraulic conductivity to zero to that no water will be 
!--                lost in the bottom layer.
                   IF ( conserve_water_content )  THEN
                      gamma_w(nzt_soil,j,i) = 0.0_wp
                   ELSE
                      gamma_w(nzt_soil,j,i) = gamma_temp(nzt_soil)
                   ENDIF     

!--                The root extraction (= root_extr * qsws_veg_eb / (rho_l     
!--                * l_v)) ensures the mass conservation for water. The         
!--                transpiration of plants equals the cumulative withdrawals by 
!--                the roots in the soil. The scheme takes into account the 
!--                availability of water in the soil layers as well as the root 
!--                fraction in the respective layer. Layer with moisture below 
!--                wilting point will not contribute, which reflects the 
!--                preference of plants to take water from moister layers.

!
!--                Calculate the root extraction (ECMWF 7.69, the sum of 
!--                root_extr = 1). The energy balance solver guarantees a 
!--                positive transpiration, so that there is no need for an 
!--                additional check.
                   m_total = 0.0_wp
                   DO  k = nzb_soil, nzt_soil
                       IF ( m_soil(k,j,i) > m_wilt(j,i) )  THEN
                          m_total = m_total + root_fr(k,j,i) * m_soil(k,j,i)
                       ENDIF
                   ENDDO  

                   IF ( m_total > 0.0_wp )  THEN
                      DO  k = nzb_soil, nzt_soil
                         IF ( m_soil(k,j,i) > m_wilt(j,i) )  THEN
                            root_extr(k) = root_fr(k,j,i) * m_soil(k,j,i)      &
                                                            / m_total
                         ELSE
                            root_extr(k) = 0.0_wp
                         ENDIF
                      ENDDO
                   ENDIF

!
!--                Prognostic equation for soil water content m_soil.
                   tend(:) = 0.0_wp

                   tend(nzb_soil) = ( lambda_w(nzb_soil,j,i) * (               &
                            m_soil(nzb_soil+1,j,i) - m_soil(nzb_soil,j,i) )    &
                            * ddz_soil(nzb_soil+1) - gamma_w(nzb_soil,j,i) - ( &
                               root_extr(nzb_soil) * qsws_veg_eb(j,i)          &
                               + qsws_soil_eb(j,i) ) * drho_l_lv )             &
                               * ddz_soil_stag(nzb_soil)

                   DO  k = nzb_soil+1, nzt_soil-1
                      tend(k) = ( lambda_w(k,j,i) * ( m_soil(k+1,j,i)          &
                                - m_soil(k,j,i) ) * ddz_soil(k+1)              &
                                - gamma_w(k,j,i)                               &
                                - lambda_w(k-1,j,i) * (m_soil(k,j,i) -         &
                                m_soil(k-1,j,i)) * ddz_soil(k)                 &
                                + gamma_w(k-1,j,i) - (root_extr(k)             &
                                * qsws_veg_eb(j,i) * drho_l_lv)                &
                                ) * ddz_soil_stag(k)

                   ENDDO
                   tend(nzt_soil) = ( - gamma_w(nzt_soil,j,i)                  &
                                           - lambda_w(nzt_soil-1,j,i)          &
                                           * (m_soil(nzt_soil,j,i)             &
                                           - m_soil(nzt_soil-1,j,i))           &
                                           * ddz_soil(nzt_soil)                &
                                           + gamma_w(nzt_soil-1,j,i) - (       &
                                             root_extr(nzt_soil)               &
                                           * qsws_veg_eb(j,i) * drho_l_lv  )   &
                                     ) * ddz_soil_stag(nzt_soil)             

                   m_soil_p(nzb_soil:nzt_soil,j,i) = m_soil(nzb_soil:nzt_soil,j,i)&
                                                   + dt_3d * ( tsc(2) * tend(:)   &
                                                   + tsc(3) * tm_soil_m(:,j,i) )   
   
!
!--                Account for dry soils (find a better solution here!)
                   DO  k = nzb_soil, nzt_soil
                      IF ( m_soil_p(k,j,i) < 0.0_wp )  m_soil_p(k,j,i) = 0.0_wp
                   ENDDO

!
!--                Calculate m_soil tendencies for the next Runge-Kutta step
                   IF ( timestep_scheme(1:5) == 'runge' )  THEN
                      IF ( intermediate_timestep_count == 1 )  THEN
                         DO  k = nzb_soil, nzt_soil
                            tm_soil_m(k,j,i) = tend(k)
                         ENDDO
                      ELSEIF ( intermediate_timestep_count <                   &
                               intermediate_timestep_count_max )  THEN
                         DO  k = nzb_soil, nzt_soil
                            tm_soil_m(k,j,i) = -9.5625_wp * tend(k) + 5.3125_wp&
                                     * tm_soil_m(k,j,i)
                         ENDDO
                      ENDIF
                   ENDIF

                ENDIF

             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE lsm_soil_model

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: mod_count

#if defined( __nopointer )

       t_surface    = t_surface_p
       t_soil       = t_soil_p
       IF ( humidity )  THEN
          m_soil    = m_soil_p
          m_liq_eb  = m_liq_eb_p
       ENDIF

#else
    
       SELECT CASE ( mod_count )

          CASE ( 0 )

             t_surface  => t_surface_1; t_surface_p  => t_surface_2
             t_soil     => t_soil_1;    t_soil_p     => t_soil_2
             IF ( humidity )  THEN
                m_soil    => m_soil_1;   m_soil_p    => m_soil_2
                m_liq_eb  => m_liq_eb_1; m_liq_eb_p  => m_liq_eb_2
             ENDIF


          CASE ( 1 )

             t_surface  => t_surface_2; t_surface_p  => t_surface_1
             t_soil     => t_soil_2;    t_soil_p     => t_soil_1
             IF ( humidity )  THEN
                m_soil    => m_soil_2;   m_soil_p    => m_soil_1
                m_liq_eb  => m_liq_eb_2; m_liq_eb_p  => m_liq_eb_1
             ENDIF

       END SELECT
#endif

    END SUBROUTINE lsm_swap_timelevel




!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
SUBROUTINE lsm_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

             CASE ( 'c_liq*' )
                IF ( .NOT. ALLOCATED( c_liq_av ) )  THEN
                   ALLOCATE( c_liq_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_liq_av = 0.0_wp

             CASE ( 'c_soil*' )
                IF ( .NOT. ALLOCATED( c_soil_av ) )  THEN
                   ALLOCATE( c_soil_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_soil_av = 0.0_wp

             CASE ( 'c_veg*' )
                IF ( .NOT. ALLOCATED( c_veg_av ) )  THEN
                   ALLOCATE( c_veg_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_veg_av = 0.0_wp

             CASE ( 'ghf_eb*' )
                IF ( .NOT. ALLOCATED( ghf_eb_av ) )  THEN
                   ALLOCATE( ghf_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ghf_eb_av = 0.0_wp

             CASE ( 'lai*' )
                IF ( .NOT. ALLOCATED( lai_av ) )  THEN
                   ALLOCATE( lai_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                lai_av = 0.0_wp

             CASE ( 'm_liq_eb*' )
                IF ( .NOT. ALLOCATED( m_liq_eb_av ) )  THEN
                   ALLOCATE( m_liq_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                m_liq_eb_av = 0.0_wp

             CASE ( 'm_soil' )
                IF ( .NOT. ALLOCATED( m_soil_av ) )  THEN
                   ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                ENDIF
                m_soil_av = 0.0_wp

             CASE ( 'qsws_eb*' )
                IF ( .NOT. ALLOCATED( qsws_eb_av ) )  THEN
                   ALLOCATE( qsws_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_eb_av = 0.0_wp

             CASE ( 'qsws_liq_eb*' )
                IF ( .NOT. ALLOCATED( qsws_liq_eb_av ) )  THEN
                   ALLOCATE( qsws_liq_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_liq_eb_av = 0.0_wp

             CASE ( 'qsws_soil_eb*' )
                IF ( .NOT. ALLOCATED( qsws_soil_eb_av ) )  THEN
                   ALLOCATE( qsws_soil_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_soil_eb_av = 0.0_wp

             CASE ( 'qsws_veg_eb*' )
                IF ( .NOT. ALLOCATED( qsws_veg_eb_av ) )  THEN
                   ALLOCATE( qsws_veg_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_veg_eb_av = 0.0_wp

             CASE ( 'r_a*' )
                IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                   ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_a_av = 0.0_wp

             CASE ( 'r_s*' )
                IF ( .NOT. ALLOCATED( r_s_av ) )  THEN
                   ALLOCATE( r_s_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_s_av = 0.0_wp

             CASE ( 'shf_eb*' )
                IF ( .NOT. ALLOCATED( shf_eb_av ) )  THEN
                   ALLOCATE( shf_eb_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                shf_eb_av = 0.0_wp

             CASE ( 't_soil' )
                IF ( .NOT. ALLOCATED( t_soil_av ) )  THEN
                   ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                ENDIF
                t_soil_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'c_liq*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_liq_av(j,i) = c_liq_av(j,i) + c_liq(j,i)
                ENDDO
             ENDDO

          CASE ( 'c_soil*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_soil_av(j,i) = c_soil_av(j,i) + (1.0 - c_veg(j,i))
                ENDDO
             ENDDO

          CASE ( 'c_veg*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_veg_av(j,i) = c_veg_av(j,i) + c_veg(j,i)
                ENDDO
             ENDDO

          CASE ( 'ghf_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ghf_eb_av(j,i) = ghf_eb_av(j,i) + ghf_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'lai*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lai_av(j,i) = lai_av(j,i) + lai(j,i)
                ENDDO
             ENDDO

          CASE ( 'm_liq_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   m_liq_eb_av(j,i) = m_liq_eb_av(j,i) + m_liq_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'm_soil' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_soil, nzt_soil
                      m_soil_av(k,j,i) = m_soil_av(k,j,i) + m_soil(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qsws_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_eb_av(j,i) = qsws_eb_av(j,i) + qsws_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'qsws_liq_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_liq_eb_av(j,i) = qsws_liq_eb_av(j,i) + qsws_liq_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'qsws_soil_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_soil_eb_av(j,i) = qsws_soil_eb_av(j,i) + qsws_soil_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'qsws_veg_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_veg_eb_av(j,i) = qsws_veg_eb_av(j,i) + qsws_veg_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 'r_a*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   r_a_av(j,i) = r_a_av(j,i) + r_a(j,i)
                ENDDO
             ENDDO

          CASE ( 'r_s*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   r_s_av(j,i) = r_s_av(j,i) + r_s(j,i)
                ENDDO
             ENDDO

          CASE ( 'shf_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   shf_eb_av(j,i) = shf_eb_av(j,i) + shf_eb(j,i)
                ENDDO
             ENDDO

          CASE ( 't_soil' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_soil, nzt_soil
                      t_soil_av(k,j,i) = t_soil_av(k,j,i) + t_soil(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'c_liq*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_liq_av(j,i) = c_liq_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'c_soil*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_soil_av(j,i) = c_soil_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'c_veg*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   c_veg_av(j,i) = c_veg_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'ghf_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   ghf_eb_av(j,i) = ghf_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

         CASE ( 'lai*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lai_av(j,i) = lai_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'm_liq_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   m_liq_eb_av(j,i) = m_liq_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'm_soil' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_soil, nzt_soil
                      m_soil_av(k,j,i) = m_soil_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'qsws_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_eb_av(j,i) = qsws_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'qsws_liq_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_liq_eb_av(j,i) = qsws_liq_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'qsws_soil_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_soil_eb_av(j,i) = qsws_soil_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'qsws_veg_eb*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   qsws_veg_eb_av(j,i) = qsws_veg_eb_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'r_a*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   r_a_av(j,i) = r_a_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'r_s*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   r_s_av(j,i) = r_s_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 't_soil' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_soil, nzt_soil
                      t_soil_av(k,j,i) = t_soil_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

    ENDIF

END SUBROUTINE lsm_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
     IMPLICIT NONE

     CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
     LOGICAL, INTENT(OUT)           ::  found       !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

     found  = .TRUE.

!
!--  Check for the grid
     SELECT CASE ( TRIM( var ) )

        CASE ( 'm_soil', 't_soil', 'm_soil_xy', 't_soil_xy', 'm_soil_xz',      &
               't_soil_xz', 'm_soil_yz', 't_soil_yz' )
           grid_x = 'x'
           grid_y = 'y'
           grid_z = 'zs'

        CASE DEFAULT
           found  = .FALSE.
           grid_x = 'none'
           grid_y = 'none'
           grid_z = 'none'
     END SELECT

 END SUBROUTINE lsm_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_data_output_2d( av, variable, found, grid, mode, local_pf,     &
                                two_d, nzb_do, nzt_do )
 
    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av !< 
    INTEGER(iwp) ::  i  !< 
    INTEGER(iwp) ::  j  !< 
    INTEGER(iwp) ::  k  !< 
    INTEGER(iwp) ::  nzb_do  !< 
    INTEGER(iwp) ::  nzt_do  !< 

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !< 

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'c_liq*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = c_liq(j,i) * c_veg(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = c_liq_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'c_soil*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = 1.0_wp - c_veg(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = c_soil_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'c_veg*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = c_veg(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = c_veg_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'ghf_eb*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = ghf_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = ghf_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'lai*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = lai(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = lai_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'm_liq_eb*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = m_liq_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = m_liq_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'm_soil_xy', 'm_soil_xz', 'm_soil_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = m_soil(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = m_soil_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          nzb_do = nzb_soil
          nzt_do = nzt_soil

          IF ( mode == 'xy' ) grid = 'zs'

       CASE ( 'qsws_eb*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  qsws_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   local_pf(i,j,nzb+1) =  qsws_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'qsws_liq_eb*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  qsws_liq_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   local_pf(i,j,nzb+1) =  qsws_liq_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'qsws_soil_eb*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  qsws_soil_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   local_pf(i,j,nzb+1) =  qsws_soil_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'qsws_veg_eb*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  qsws_veg_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   local_pf(i,j,nzb+1) =  qsws_veg_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'


       CASE ( 'r_a*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = r_a(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = r_a_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'r_s*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = r_s(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = r_s_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'shf_eb*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  shf_eb(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) =  shf_eb_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 't_soil_xy', 't_soil_xz', 't_soil_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = t_soil(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = t_soil_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          nzb_do = nzb_soil
          nzt_do = nzt_soil

          IF ( mode == 'xy' )  grid = 'zs'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT
 
 END SUBROUTINE lsm_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_data_output_3d( av, variable, found, local_pf )
 

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av    !< 
    INTEGER(iwp) ::  i     !< 
    INTEGER(iwp) ::  j     !< 
    INTEGER(iwp) ::  k     !< 

    LOGICAL      ::  found !< 

    REAL(sp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb_soil:nzt_soil) ::  local_pf !< 


    found = .TRUE.


    SELECT CASE ( TRIM( variable ) )


      CASE ( 'm_soil' )

         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = m_soil(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = m_soil_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 't_soil' )

         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = t_soil(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = t_soil_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF


       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE lsm_data_output_3d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Write restart data for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_last_actions
 

    USE control_parameters
        
    USE kinds

    IMPLICIT NONE

    IF ( write_binary(1:4) == 'true' )  THEN
       IF ( ALLOCATED( c_liq_av ) )  THEN
          WRITE ( 14 )  'c_liq_av            ';  WRITE ( 14 ) c_liq_av
       ENDIF
       IF ( ALLOCATED( c_soil_av ) )  THEN
          WRITE ( 14 )  'c_soil_av           ';  WRITE ( 14 ) c_soil_av
       ENDIF
       IF ( ALLOCATED( c_veg_av ) )  THEN
          WRITE ( 14 )  'c_veg_av            ';  WRITE ( 14 ) c_veg_av
       ENDIF
       IF ( ALLOCATED( ghf_eb_av ) )  THEN
          WRITE ( 14 )  'ghf_eb_av           ';  WRITE ( 14 )  ghf_eb_av
       ENDIF
       IF ( ALLOCATED( lai_av ) )  THEN
          WRITE ( 14 )  'lai_av              ';  WRITE ( 14 )  lai_av
       ENDIF
       WRITE ( 14 )  'm_liq_eb            ';  WRITE ( 14 )  m_liq_eb
       IF ( ALLOCATED( m_liq_eb_av ) )  THEN
          WRITE ( 14 )  'm_liq_eb_av         ';  WRITE ( 14 )  m_liq_eb_av
       ENDIF
       WRITE ( 14 )  'm_soil              ';  WRITE ( 14 )  m_soil
       IF ( ALLOCATED( m_soil_av ) )  THEN
          WRITE ( 14 )  'm_soil_av           ';  WRITE ( 14 )  m_soil_av
       ENDIF
       IF ( ALLOCATED( qsws_eb_av ) )  THEN
          WRITE ( 14 )  'qsws_eb_av          ';  WRITE ( 14 )  qsws_eb_av
       ENDIF   
       IF ( ALLOCATED( qsws_liq_eb_av ) )  THEN
          WRITE ( 14 )  'qsws_liq_eb_av      ';  WRITE ( 14 )  qsws_liq_eb_av
       ENDIF  
       IF ( ALLOCATED( qsws_soil_eb_av ) )  THEN
          WRITE ( 14 )  'qsws_soil_eb_av     ';  WRITE ( 14 )  qsws_soil_eb_av
       ENDIF 
       IF ( ALLOCATED( qsws_veg_eb_av ) )  THEN
          WRITE ( 14 )  'qsws_veg_eb_av      ';  WRITE ( 14 )  qsws_veg_eb_av
       ENDIF 
       IF ( ALLOCATED( shf_eb_av ) )  THEN
          WRITE ( 14 )  'shf_eb_av           ';  WRITE ( 14 )  shf_eb_av
       ENDIF
       WRITE ( 14 )  't_soil              ';  WRITE ( 14 )  t_soil
       IF ( ALLOCATED( t_soil_av ) )  THEN
          WRITE ( 14 )  't_soil_av           ';  WRITE ( 14 )  t_soil_av
       ENDIF

       WRITE ( 14 )  '*** end lsm ***     '

    ENDIF

 END SUBROUTINE lsm_last_actions


SUBROUTINE lsm_read_restart_data( i, nxlfa, nxl_on_file, nxrfa, nxr_on_file,   &
                                     nynfa, nyn_on_file, nysfa, nys_on_file,   &
                                     offset_xa, offset_ya, overlap_count,      &
                                     tmp_2d )
 

    USE control_parameters
        
    USE indices
    
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=20) :: field_char   !< 

    INTEGER(iwp) ::  i               !< 
    INTEGER(iwp) ::  k               !< 
    INTEGER(iwp) ::  nxlc            !< 
    INTEGER(iwp) ::  nxlf            !< 
    INTEGER(iwp) ::  nxl_on_file     !< 
    INTEGER(iwp) ::  nxrc            !< 
    INTEGER(iwp) ::  nxrf            !< 
    INTEGER(iwp) ::  nxr_on_file     !< 
    INTEGER(iwp) ::  nync            !< 
    INTEGER(iwp) ::  nynf            !< 
    INTEGER(iwp) ::  nyn_on_file     !< 
    INTEGER(iwp) ::  nysc            !< 
    INTEGER(iwp) ::  nysf            !< 
    INTEGER(iwp) ::  nys_on_file     !< 
    INTEGER(iwp) ::  overlap_count   !< 

    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxlfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxrfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nynfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nysfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_xa   !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_ya   !< 

    REAL(wp),                                                                  &
       DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_2d   !< 

    REAL(wp),                                                                  &
       DIMENSION(nzb_soil:nzt_soil+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_3d   !< 

    REAL(wp),                                                                  &
       DIMENSION(nzb_soil:nzt_soil,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_3d2   !< 


   IF ( initializing_actions == 'read_restart_data' )  THEN
      READ ( 13 )  field_char

      DO  WHILE ( TRIM( field_char ) /= '*** end lsm ***' )

         DO  k = 1, overlap_count

            nxlf = nxlfa(i,k)
            nxlc = nxlfa(i,k) + offset_xa(i,k)
            nxrf = nxrfa(i,k)
            nxrc = nxrfa(i,k) + offset_xa(i,k)
            nysf = nysfa(i,k)
            nysc = nysfa(i,k) + offset_ya(i,k)
            nynf = nynfa(i,k)
            nync = nynfa(i,k) + offset_ya(i,k)


            SELECT CASE ( TRIM( field_char ) )

                CASE ( 'c_liq_av' )
                   IF ( .NOT. ALLOCATED( c_liq_av ) )  THEN
                      ALLOCATE( c_liq_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   c_liq_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'c_soil_av' )
                   IF ( .NOT. ALLOCATED( c_soil_av ) )  THEN
                      ALLOCATE( c_soil_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   c_soil_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'c_veg_av' )
                   IF ( .NOT. ALLOCATED( c_veg_av ) )  THEN
                      ALLOCATE( c_veg_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   c_veg_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ghf_eb_av' )
                   IF ( .NOT. ALLOCATED( ghf_eb_av ) )  THEN
                      ALLOCATE( ghf_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ghf_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'm_liq_eb' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   m_liq_eb(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =        &
                                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'lai_av' )
                   IF ( .NOT. ALLOCATED( lai_av ) )  THEN
                      ALLOCATE( lai_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   lai_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'm_liq_eb_av' )
                   IF ( .NOT. ALLOCATED( m_liq_eb_av ) )  THEN
                      ALLOCATE( m_liq_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   m_liq_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =      &
                                  tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'm_soil' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d2(:,:,:)
                   m_soil(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                          tmp_3d2(nzb_soil:nzt_soil,nysf-nbgp:nynf             &
                          +nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'm_soil_av' )
                   IF ( .NOT. ALLOCATED( m_soil_av ) )  THEN
                      ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d2(:,:,:)
                   m_soil_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =      &
                                    tmp_3d2(nzb_soil:nzt_soil,nysf             &
                                    -nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws_eb_av' )
                   IF ( .NOT. ALLOCATED( qsws_eb_av ) )  THEN
                      ALLOCATE( qsws_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws_liq_eb_av' )
                   IF ( .NOT. ALLOCATED( qsws_liq_eb_av ) )  THEN
                      ALLOCATE( qsws_liq_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_liq_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                CASE ( 'qsws_soil_eb_av' )
                   IF ( .NOT. ALLOCATED( qsws_soil_eb_av ) )  THEN
                      ALLOCATE( qsws_soil_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_soil_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws_veg_eb_av' )
                   IF ( .NOT. ALLOCATED( qsws_veg_eb_av ) )  THEN
                      ALLOCATE( qsws_veg_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_veg_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =  &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_eb_av' )
                   IF ( .NOT. ALLOCATED( shf_eb_av ) )  THEN
                      ALLOCATE( shf_eb_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_eb_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 't_soil' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   t_soil(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,               &
                                                nxlf-nbgp:nxrf+nbgp)

                CASE ( 't_soil_av' )
                   IF ( .NOT. ALLOCATED( t_soil_av ) )  THEN
                      ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d2(:,:,:)
                   t_soil_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =      &
                                    tmp_3d(:,nysf-nbgp:nynf+nbgp,             &
                                    nxlf-nbgp:nxrf+nbgp)


               CASE DEFAULT
                  WRITE( message_string, * ) 'unknown variable named "',       &
                                        TRIM( field_char ), '" found in',      &
                                        '&data from prior run on PE ', myid
                  CALL message( 'lsm_read_restart_data', 'PA0441', 1, 2, 0, 6, &
                                 0 )

            END SELECT

         ENDDO

         READ ( 13 )  field_char

      ENDDO
   ENDIF

 END SUBROUTINE lsm_read_restart_data

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of roughness length for open water (lakes, ocean). The
!> parameterization follows Charnock (1955). Two different implementations
!> are available: as in ECMWF-IFS (Beljaars 1994) or as in FLake (Subin et al.
!> 2012)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_z0_water_surface

       USE control_parameters,                                                 &
           ONLY: g, kappa, molecular_viscosity

       IMPLICIT NONE

       INTEGER :: i  !< running index
       INTEGER :: j  !< running index

       REAL(wp), PARAMETER :: alpha_ch  = 0.018_wp !< Charnock constant (0.01-0.11). Use 0.01 for FLake and 0.018 for ECMWF
!       REAL(wp), PARAMETER :: pr_number = 0.71_wp !< molecular Prandtl number in the Charnock parameterization (differs from prandtl_number)
!       REAL(wp), PARAMETER :: sc_number = 0.66_wp !< molecular Schmidt number in the Charnock parameterization
!       REAL(wp) :: re_0 !< near-surface roughness Reynolds number


       DO  i = nxlg, nxrg    
          DO  j = nysg, nyng
             IF ( water_surface(j,i) )  THEN

!
!--             Disabled: FLake parameterization. Ideally, the Charnock 
!--             coefficient should depend on the water depth and the fetch
!--             length
!                re_0 = z0(j,i) * us(j,i) / molecular_viscosity
!        
!                z0(j,i) = MAX( 0.1_wp * molecular_viscosity / us(j,i),            &
!                              alpha_ch * us(j,i) / g )
!
!                z0h(j,i) = z0(j,i) * EXP( - kappa / pr_number * ( 4.0_wp * SQRT( re_0 ) - 3.2_wp ) )
!                z0q(j,i) = z0(j,i) * EXP( - kappa / pr_number * ( 4.0_wp * SQRT( re_0 ) - 4.2_wp ) )

!
!--              Set minimum roughness length for u* > 0.2
!                IF ( us(j,i) > 0.2_wp )  THEN
!                   z0h(j,i) = MAX( 1.0E-5_wp, z0h(j,i) )
!                   z0q(j,i) = MAX( 1.0E-5_wp, z0q(j,i) )
!                ENDIF

!
!--             ECMWF IFS model parameterization after Beljaars (1994). At low
!--             wind speed, the sea surface becomes aerodynamically smooth and
!--             the roughness scales with the viscosity. At high wind speed, the
!--             Charnock relation is used.
                z0(j,i) =   ( 0.11_wp * molecular_viscosity / us(j,i) )        &
                          + ( alpha_ch * us(j,i)**2 / g )

                z0h(j,i) = 0.40_wp * molecular_viscosity / us(j,i)
                z0q(j,i) = 0.62_wp * molecular_viscosity / us(j,i)

             ENDIF
          ENDDO
       ENDDO

    END SUBROUTINE calc_z0_water_surface


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of specific humidity of the skin layer (surface). It is assumend
!> that the skin is always saturated.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_q_surface

       IMPLICIT NONE

       INTEGER :: i              !< running index
       INTEGER :: j              !< running index
       INTEGER :: k              !< running index

       REAL(wp) :: resistance    !< aerodynamic and soil resistance term

       DO  i = nxlg, nxrg    
          DO  j = nysg, nyng
             k = nzb_s_inner(j,i)

!
!--          Calculate water vapour pressure at saturation
             e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp * ( t_surface_p(j,i)   &
                   - 273.16_wp ) / ( t_surface_p(j,i) - 35.86_wp ) )

!
!--          Calculate specific humidity at saturation
             q_s = 0.622_wp * e_s / surface_pressure

             resistance = r_a(j,i) / (r_a(j,i) + r_s(j,i))

!
!--          Calculate specific humidity at surface
             IF ( cloud_physics )  THEN
                q(k,j,i) = resistance * q_s + (1.0_wp - resistance)            &
                             * ( q(k+1,j,i) - ql(k+1,j,i) )
             ELSE
                q(k,j,i) = resistance * q_s + (1.0_wp - resistance)            &
                             * q(k+1,j,i)
             ENDIF

!
!--          Update virtual potential temperature
             vpt(k,j,i) = pt(k,j,i) * ( 1.0_wp + 0.61_wp * q(k,j,i) )

          ENDDO
       ENDDO

    END SUBROUTINE calc_q_surface


 END MODULE land_surface_model_mod
