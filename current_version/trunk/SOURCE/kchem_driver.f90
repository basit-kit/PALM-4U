!> @file kchem_driver.f90
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
! kk: Intial version (Klaus Ketelsen)
! bK: Changed initial values of chem_species
! RFo: Added tmp_fact
! FKa: Some formatting, todos added below, "_wp" added to REAL quantities' values,
!      the "if(var(1:3) == 'kc_')" is now included in check_parameters (as done
!      for "usm_" quantities
!
! Former revisions:
! -----------------
!
! Description:
! ------------
!> Interface PALM <-> kpp chemistry
!>
!> @todo Format this module according to PALM coding standard (see e.g. module
!>       template under http://palm.muk.uni-hannover.de/mosaik/downloads/8 or
!>       D3_coding_standard.pdf under https://palm.muk.uni-hannover.de/trac/downloads/16)
!> @todo Rename ssws, sswst to csws, cswst to avoid confusion with ssws/sswst 
!>       used together with passive scalar s
!> @todo Use same run index for all loops over NSPEC, e.g. nsp = 1, NSPEC
!>       (and: NSPEC-->nspec, i.e. no capital letters in variable names 
!>       [see formatting rules])
!> @todo "#ifdef KPP_CHEM" --> "#if defined( __chem)" (see other examples 
!>       in PALM)               #if defined( __parallel )
!> @todo „write(6...“ printing should eventually be handled in a separate 
!>       subroutine „kchem_header“ (see e.g. how it is done in routine 
!>       lsm_header in land_surface_model_mod.f90), which prints this 
!>       information to a MONITORING file named <jobname>_header (in this file, 
!>       all PALM setup information are included)
!> @todo Not all routines have meaningful names yet, e.g. "kchem_integrate" only
!>       refers to chemical reactions part, but the name itself could also refer
!>       to the other part of the chemistry integration (advection/diffusion)
!>       --> Please clarify
!> @todo Rename flag „use_kpp_chemistry“ to „chemistry“ (see e.g. flag 
!>       „plant_canopy“ used for the plant canopy model
!> @todo Please enter all of your TODOs here in this list
!> @todo Every subroutine needs a "Description" header 
!>       (as e.g. in urban_surface_mod.f90)
!> @todo Clear circular dependencies between check_parameters, kchem_initialize,
!>       and kchem_check_data_output

!------------------------------------------------------------------------------!

MODULE kchem_driver
   USE kinds,              ONLY: wp, iwp
   USE indices,            ONLY: nzb,nzt,nysg,nyng,nxlg,nxrg,nzb_s_inner,nys,nyn
   USE pegrid,             ONLY: myid, threads_per_task
   USE control_parameters, ONLY: dt_3d, ws_scheme_sca, initializing_actions  !ws_sch... added by bK
   USE arrays_3d,          ONLY: pt
   USE kchem_kpp,          ONLY: NSPEC, SPC_NAMES, NKPPCTRL, NMAXFIXSTEPS, t_steps, FILL_TEMP, kpp_integrate,     &
                                 NVAR, ATOL, RTOL
   USE cpulog,             ONLY: cpu_log, log_point_s

   IMPLICIT   NONE
   PRIVATE
   SAVE

!- Define chemical variables

   TYPE   species_def
      CHARACTER(LEN=8)                                   :: name
      CHARACTER(LEN=16)                                  :: unit
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc_p
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: tconc_m
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: rssws, rsswst           !bK changed ssws and sswst  to rssws and rsswst
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: flux_s, diss_s
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:)         :: flux_l, diss_l
   END TYPE species_def

   PUBLIC  species_def

   logical, PUBLIC                                               :: use_kpp_chemistry = .FALSE.

   TYPE(species_def),ALLOCATABLE,DIMENSION(:),TARGET, PUBLIC     :: chem_species

   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_1,spec_conc_2,spec_conc_3, spec_conc_av       !bK added spec_conc_av

   INTEGER,DIMENSION(NKPPCTRL)                           :: icntrl                            ! Fine tuning kpp
   REAL(kind=wp),DIMENSION(NKPPCTRL)                     :: rcntrl                            ! Fine tuning kpp

!------------------ For mean concentration over specified dt                         !bK add this block
   TYPE   species_av_def
      CHARACTER(LEN=8)                                   :: name
      CHARACTER(LEN=16)                                  :: unit
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc_av
   END TYPE species_av_def

   PUBLIC  species_av_def

   TYPE(species_av_def),ALLOCATABLE,DIMENSION(:),TARGET, PUBLIC     :: chem_species_av
!-------------------

   PUBLIC NSPEC
   PUBLIC NVAR       ! added NVAR for pe  bK, kd3
   PUBLIC SPC_NAMES  ! added for pe bk, kd4

!- Interface section
  INTERFACE kchem_boundary_conds
      MODULE PROCEDURE kchem_boundary_conds
  END INTERFACE kchem_boundary_conds

   INTERFACE kchem_initialize
      MODULE PROCEDURE kchem_initialize
   END INTERFACE kchem_initialize

   INTERFACE kchem_parin
      MODULE PROCEDURE kchem_parin
   END INTERFACE kchem_parin

   INTERFACE kchem_integrate
      MODULE PROCEDURE kchem_integrate_ij
   END INTERFACE kchem_integrate

   INTERFACE kchem_swap_timelevel
      MODULE PROCEDURE kchem_swap_timelevel
   END INTERFACE kchem_swap_timelevel

   INTERFACE kchem_define_netcdf_grid
      MODULE PROCEDURE kchem_define_netcdf_grid
   END INTERFACE kchem_define_netcdf_grid

   INTERFACE kchem_data_output_3d
      MODULE PROCEDURE kchem_data_output_3d
   END INTERFACE kchem_data_output_3d

   INTERFACE kchem_check_data_output
      MODULE PROCEDURE kchem_check_data_output
   END INTERFACE kchem_check_data_output

   INTERFACE kchem_check_data_output_pr
      MODULE PROCEDURE kchem_check_data_output_pr
   END INTERFACE kchem_check_data_output_pr

   INTERFACE kchem_3d_data_averaging
      MODULE PROCEDURE kchem_3d_data_averaging 
   END INTERFACE kchem_3d_data_averaging

   INTERFACE kchem_last_actions
      MODULE PROCEDURE kchem_last_actions 
   END INTERFACE kchem_last_actions

   INTERFACE kchem_read_restart_data
      MODULE PROCEDURE kchem_read_restart_data
   END INTERFACE kchem_read_restart_data



   PUBLIC kchem_check_data_output, kchem_data_output_3d,                       &
          kchem_define_netcdf_grid, kchem_initialize, kchem_integrate,         &
          kchem_parin, kchem_swap_timelevel, kchem_boundary_conds,             &
          kchem_3d_data_averaging, kchem_check_data_output_pr,                 &
          kchem_last_actions, kchem_read_restart_data
            


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize and set all boundary conditions for chemical species concentration
!------------------------------------------------------------------------------!

 SUBROUTINE kchem_boundary_conds( mode )                                        !Please note alphabetical order of listings, and also all other formatting
                                                                                !and continue in that exact way.
    USE control_parameters,                                                  &  !you might also have to add some other quantities here, since
        ONLY: bc_rs_b, bc_rs_t, bc_rs_t_val, chemistry, ibc_rs_b, ibc_rs_t,  &  ! bK added bc_rs_* and chemistry here. 
              outflow_l, outflow_n, outflow_r, outflow_s                        ! bK added ibc_rs_b and ibc_c_t, bc_rs_t_tal,reactive_scalar, 
    USE indices,                                                             &  !I didn't check for all that might be necessary to add to ONLY list
        ONLY:nxl, nxr,  nxlg, nxrg, nyng, nysg, nzb_s_inner, nzt                ! bK added nzt, nzb_c_inner, nxl, nxr


!    USE prognostic_equations_mod,                                            &

    USE arrays_3d, ONLY: dzu, rs_p, rs, trs_m                            !bK added dzu, rs_p, rs, trs_m
                                                             

    IMPLICIT NONE
    
    CHARACTER (len=*), INTENT(IN) :: mode
    INTEGER i,j                                                                 !bK added i, j

    SELECT CASE ( TRIM( mode ) )                                               !bK commented it


       CASE ( 'init' )                             ! bK commented it
!         This part below is from check_parameters.f90's SUBROUTINE set_bc_scalars
!         I already adjusted this part to work for chem species (see notes under CASE('set_bc_bottomtop'))
!
!--       Set Integer flags and check for possible erroneous settings for 
!--       bottom boundary condition  (the integer flags are then used in CASE('set_bc_bottomtop'))
!    if (myid == 0) print*,'entering kchem_boundary_conds ***** ', mode                                !bK      


          IF ( bc_rs_b == 'dirichlet' )  THEN
             ibc_rs_b = 0 
          ELSEIF ( bc_rs_b == 'neumann' )  THEN
             ibc_rs_b = 1 
          ELSE
              print*, 'unknown boundary condition: bc_rs_b ="'// TRIM( bc_rs_b ) // '"'            !bK added
!             message_string = 'unknown boundary condition: bc_rs_b ="' // TRIM( bc_rs_b ) // '"'  ! bK commented
             CALL message( 'kchem_boundary_conds', 'CHEM001', 1, 2, 0, 6, 0 )     !checks in kchem_driver have special error numbers --> "CHEM###",
          ENDIF                                                                   !so that we don't mix them up with PALM error numbers "PA0###"
!
!
!--       Set Integer flags and check for possible erroneous settings for top
!--       boundary condition. bK added *_rs_* here.
          IF ( bc_rs_t == 'dirichlet' )  THEN
             ibc_rs_t = 0 
          ELSEIF ( bc_rs_t == 'neumann' )  THEN
             ibc_rs_t = 1
          ELSEIF ( bc_rs_t == 'initial_gradient' )  THEN
             ibc_rs_t = 2
          ELSEIF ( bc_rs_t == 'nested' )  THEN
             ibc_rs_t = 3
          ELSE
             print*, 'unknown boundary condition: bc_rs_b ="'// TRIM( bc_rs_t ) // '"'            !bK added
!            message_string = 'unknown boundary condition: bc_rs_t ="' // TRIM( bc_rs_t ) // '"'  ! bK commented
             CALL message( 'check_parameters', 'CHEM002', 1, 2, 0, 6, 0 )
          ENDIF
  if(myid==0) print*,'fm kchem_boundary_conds ***** #3 ', mode                                !bK debug    

!
!--       We might later on need some parameter checkings as in check_parameters.f90's SUBROUTINE check_bc_scalars
!--       (these should be put under a new SUBROUTINE kchem_check_parameters, see e.g. urban_surface_mod.f90 also has SUBROUTINE usm_check_parameters)


       CASE ( 'set_bc_bottomtop' )                   !bK commented it


!--       Surface conditions for chemical species     !!!!this is the code for passive_scalar (copied from boundary_conds.f90) 
                                                      !!!!---> modify this to work for chemical species!!!
          IF ( ibc_rs_b == 0 ) THEN                   !use "ibc_rs_b" for chem species. you have to declare it in kchem_driver's global decleration part.
                                                      !and make it PUBLICly available. same holds for ibc_rs_t, bc_rs_b/bc_rs_t.
                                                      !don't confuse the integer flag "ibc_rs_b" with the &kpp_chem parameter "bc_rs_b" which needs to
                                                      !be defined as a parameter in kchem_parin. "bc_rs_b" can take the same values as "bc_s_b" (see PALM
                                                      !documentation)
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   rs_p(nzb_s_inner(j,i),j,i) = rs(nzb_s_inner(j,i),j,i)

                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
! if(myid == 0) print*,' *********  before ibc_rs_b  kchem_boundary_conds *****  ', ibc_rs_b                            !bK     
!                    print*,'the val of rs_p is : ', shape(rs_p), size(rs_p) 
                   rs_p(nzb_s_inner(j,i),j,i) = rs_p(nzb_s_inner(j,i)+1,j,i)
! if(myid == 0) print*,'********** after ibc_rs_bkchem_boundary_conds ********** ', ibc_rs_b                                !bK      

                ENDDO
             ENDDO
          ENDIF
!
!--       Top boundary conditions for chemical species
          IF ( ibc_rs_t == 0 )  THEN                    !use "ibc_rs_t" for chem species. you have to declare it in kchem_driver global decleration part
             rs_p(nzt+1,:,:) = rs(nzt+1,:,:)             !(see notes for "ibc_rs_b" above)
          ELSEIF ( ibc_rs_t == 1 )  THEN
             rs_p(nzt+1,:,:) = rs_p(nzt,:,:)
          ELSEIF ( ibc_rs_t == 2 )  THEN
             rs_p(nzt+1,:,:) = rs_p(nzt,:,:) + bc_rs_t_val * dzu(nzt+1)
          ENDIF


!
       CASE ( 'set_bc_lateral' )                       ! bK commented it
!       print*,'entering set bc lateral, kchem_boundary_conds', mode                                !bK      

!

!-- Lateral boundary conditions for chem species at outflow boundary (following code was copied from boundary_conds. It needs to be adjusted!)
          IF ( outflow_s )  THEN
             IF ( chemistry )  rs_p(:,nys-1,:) = rs_p(:,nys,:)    !bK added reactive_scalar in nested-if
          ELSEIF ( outflow_n )  THEN
             IF ( chemistry )  rs_p(:,nyn+1,:) = rs_p(:,nyn,:)
          ELSEIF ( outflow_l )  THEN
             IF ( chemistry )  rs_p(:,:,nxl-1) = rs_p(:,:,nxl)
          ELSEIF ( outflow_r )  THEN
             IF ( chemistry )  rs_p(:,:,nxr+1) = rs_p(:,:,nxr)
          ENDIF

! Lateral BC for inflow are automatically set when chem species are initialized (that's what Siggi told us in our discussion about BCs)
! So far, I think that's all we have to do.
!       print*,'leaving set bc lateral, kchem_boundary_conds', mode                                !bK      

       END SELECT

 END SUBROUTINE kchem_boundary_conds

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the kchem_driver
!------------------------------------------------------------------------------!
   SUBROUTINE kchem_initialize
      IMPLICIT   none
!--   local variables
      INTEGER                  :: i

!--   Allocate Memory for chemical species

      ALLOCATE(chem_species(NSPEC))
      ALLOCATE(spec_conc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
      ALLOCATE(spec_conc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
      ALLOCATE(spec_conc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))

      ALLOCATE(chem_species_av(NSPEC))          !bK added
      ALLOCATE(spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC)) 




      if(myid == 0)  then
         write(6,'(/,a,/)')  'kpp chemics >>>> List of species: '
      end if

      do i=1,NSPEC
         chem_species(i)%name = SPC_NAMES(i)
         chem_species_av(i)%name = SPC_NAMES(i)                !bK added

         if(myid == 0)  then
            write(6,'(a,i4,3x,a)')  '   Species: ',i,trim(SPC_NAMES(i))
         end if
         chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,i)
         chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,i)
         chem_species(i)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => spec_conc_3(:,:,:,i)
         chem_species_av(i)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => spec_conc_av(:,:,:,i)     !bK added

         ALLOCATE (chem_species(i)%rssws(nysg:nyng,nxlg:nxrg))                          !bK replaced ssws with rssws
         ALLOCATE (chem_species(i)%rsswst(nysg:nyng,nxlg:nxrg))                         !bK replaced sswst with rsswst 

         chem_species(i)%rssws  = 0.0_wp                                                !bK replaced ssws with rssws
         chem_species(i)%rsswst = 0.0_wp                                                !bK replaced sswst with rsswst                             


!          IF ( ws_scheme_sca )  THEN                                           !fk (needs to be revisited. due to
            ALLOCATE (chem_species(i)%flux_s(nzb+1:nzt,0:threads_per_task-1))   !circular dependencies between check_parameters,
            ALLOCATE (chem_species(i)%diss_s(nzb+1:nzt,0:threads_per_task-1))   !kchem_initialize, kchem_check_data_output,
            ALLOCATE (chem_species(i)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) !this IF clause creates problems,
            ALLOCATE (chem_species(i)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) !since ws_scheme_sca is first set "true"  
            chem_species(i)%flux_s = 0.0_wp                          !bK        !in check_parameters, which is so far called   
            chem_species(i)%flux_l = 0.0_wp                                     !after kchem_initialize
            chem_species(i)%diss_s = 0.0_wp
            chem_species(i)%diss_l = 0.0_wp
!          ENDIF         

      end do

      if(myid== 0) write(9,*) 'fm 100, chem_species_av%name is ', chem_species_av%name
      if(myid== 0) write(9,*) 'fm 100, shape of chem_species_av%conc_av is ', shape(chem_species_av(1)%conc_av)

!--   Set initial values

      IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN             !bK added if block
            CALL set_const_initial_values
      END IF

      if(myid == 0)  write(6,*) ' '

      RETURN

     CONTAINS
      SUBROUTINE set_const_initial_values
         IMPLICIT   none

!--      local variables
         INTEGER                  :: i

         if(myid == 0)  then
            write(6,'(/,a,/)')  'kpp chemics >>>> Set constant Initial Values: '
         end if

!        Default values are taken from smog.def from supplied kpp example
         do i=1,NSPEC
            if(trim(chem_species(i)%name) == 'NO')   then
!              chem_species(i)%conc = 8.725*1.0E+08
               chem_species(i)%conc = 0.5_wp                           !added by bK
!               chem_species(i)%conc = 0.01_wp                          !added by RFo
            else if(trim(chem_species(i)%name) == 'NO2') then
!              chem_species(i)%conc = 2.240*1.0E+08             
               chem_species(i)%conc = 0.1_wp                           !added by bK
!               chem_species(i)%conc = 0.05_wp                          !added by RFo
            else if(trim(chem_species(i)%name) == 'O3') then
               chem_species(i)%conc = 0.05_wp                          !added by bK
            else if(trim(chem_species(i)%name) == 'H2O') then
!              chem_species(i)%conc = 5.326*1.0E+11
               chem_species(i)%conc = 1.30*1.0E+4_wp                   !added by bK
            else if(trim(chem_species(i)%name) == 'O2') then
!              chem_species(i)%conc = 1.697*1.0E+16
               chem_species(i)%conc = 2.0*1.0E+5_wp                    !added by bK
            else if(trim(chem_species(i)%name) == 'RH') then
!               chem_species(i)%conc = 9.906*1.0E+01
                chem_species(i)%conc = 2.0_wp                          !added by bK
!                chem_species(i)%conc = 0.001_wp                        !added by RFo
            else if(trim(chem_species(i)%name) == 'RCHO') then
!              chem_species(i)%conc = 6.624*1.0E+08
              chem_species(i)%conc = 2.0_wp                           !added by bK
!               chem_species(i)%conc = 0.01_wp                          !added by RFo
!            else if(trim(chem_species(i)%name) == 'OH') then
!               chem_species(i)%conc = 1.0*1.0E-07_wp                   !added by bK
!            else if(trim(chem_species(i)%name) == 'HO2') then
!               chem_species(i)%conc = 1*1.0E-7_wp                      !added by bK
!            else if(trim(chem_species(i)%name) == 'RCCO2') then
!               chem_species(i)%conc = 1.0*1.0E-7_wp                    !added by bK
!            else if(trim(chem_species(i)%name) == 'RCOO2NO2') then
!               chem_species(i)%conc = 1.0*1.0E-7_wp                   !added by bK
            if(myid==0) print*, ' from kchem_initialize #1 '          !bK debug
           else
!   H2O = 2.0e+04;
               chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg) = 0.0_wp
            end if
!             chem_species(i)%conc_p  = 0.0_wp
            chem_species(i)%conc_p  = chem_species(i)%conc       ! added bK
            chem_species(i)%tconc_m = 0.0_wp

            if(myid == 0)  then
               write(6,'(a,3x,a,3x,a,e12.4)')  '   Species:     ',chem_species(i)%name(1:7),' Initial Value = ',chem_species(i)%conc(nzb,nysg,nxlg)
            end if
         end do

#if defined( __nopointer )
!kk      Hier mit message abbrechen
         if(myid == 0)  then
            write(6,*)  '   KPP does only run with POINTER Version'
         end if
         stop 'error'
#endif

         RETURN
      END SUBROUTINE set_const_initial_values


   END SUBROUTINE kchem_initialize

   SUBROUTINE kchem_parin
      IMPLICIT   none

      CHARACTER (LEN=80)                             ::  line                   ! < dummy string that contains the current line of the parameter file
      REAL(kind=wp), DIMENSION(NMAXFIXSTEPS)         ::  my_steps               ! List of fixed timesteps   my_step(1) = 0.0 automatic stepping

      NAMELIST /kpp_chem/ icntrl, rcntrl, my_steps

!--   Read kpp_chem namelist
      icntrl    = 0
      rcntrl    = 0.0_wp
      my_steps  = 0.0_wp
      line      = ' '
      icntrl(2) = 1                                   !Atol and Rtol are Scalar

      ATOL = 1.0_wp
      RTOL = 0.01_wp
!
!--   Try to find kpp chemistry package
      REWIND ( 11 )
      line = ' '
      DO   WHILE ( INDEX( line, '&kpp_chem' ) == 0 )
        READ ( 11, '(A)', END=10 )  line
      ENDDO
      BACKSPACE ( 11 )

!
!--   Read user-defined namelist
      READ ( 11, kpp_chem )

      use_kpp_chemistry = .TRUE.

 10   CONTINUE

      t_steps = my_steps

      if(myid <= 1)  then
         write(6,*) 'KPP Parin ',icntrl(2),icntrl(3),ATOL(1)
      end if

   write(9,*) 'fm kcd, leaving kc_parin #-1'
   flush(9)
      RETURN
   END SUBROUTINE kchem_parin

   SUBROUTINE kchem_integrate_ij (i, j)
      IMPLICIT   none
      INTEGER,INTENT(IN)       :: i,j

!--   local variables
      INTEGER                  :: k,m,istat
      INTEGER,dimension(20)    :: istatus
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt,NSPEC)                :: tmp_conc           
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt)                      :: tmp_temp
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt)                      :: tmp_fact   !! RFo

! ppm to molecules/cm**3
       tmp_fact = 10.e-6_wp*6.022e23_wp/22.414_wp/1000._wp     !*T/273.15*1013/p, 10.e-6_iwp replaced by 10.e-6_wp ..!bK
!      tmp_temp = 292
!      CALL fill_temp (istat, tmp_temp)                             ! Load constant temperature into kpp context
      CALL fill_temp (istat, pt(nzb_s_inner(j,i)+1:nzt,j,i))                             ! Load temperature into kpp context


      do m=1,NSPEC
         tmp_conc(:,m) = chem_species(m)%conc (nzb_s_inner(j,i)+1:nzt,j,i) * tmp_fact(:) ! RFo
      end do

      if(myid == 0 .and. i == 10 .and. j == 10)  then
         write(0,*) 'begin KPP step ',dt_3d
      end if
      if(myid == 0)  print*,'fm kc_driver, calling kpp_integrate in kchem_kpp #10.1'   !bK debug
      CALL cpu_log( log_point_s(80), 'kpp_integrate', 'start' )

      CALL kpp_integrate (dt_3d, tmp_conc, istatus=istatus)

      CALL cpu_log( log_point_s(80), 'kpp_integrate', 'stop' )

      do m=1,NSPEC
         chem_species(m)%conc (nzb_s_inner(j,i)+1:nzt,j,i) = tmp_conc(:,m) / tmp_fact(:)  ! RFo
      end do

      if(myid == 0 .and. i == 10 .and. j == 10)  then
         write(6,'(a,8i7)') ' KPP Status ',istatus(1:8)
      end if

      RETURN
   END SUBROUTINE kchem_integrate_ij

   SUBROUTINE kchem_swap_timelevel (level)
      IMPLICIT   none

      INTEGER,INTENT(IN)                  :: level

!--   local variables

      INTEGER               :: i

!        print*,' *** entering kchem_swap_timelevel ***) '
      if(level == 0)  then
         do i=1,NVAR                                        ! NSPEC replaced with NVAR bK kd1  
            chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,i)
            chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,i)
         end do
      else
         do i=1,NVAR                                        ! NSPEC replaced with NVAR bk  kd2
            chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_2(:,:,:,i)
            chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_1(:,:,:,i)
         end do
      end if

!        print*,' *** leaving kchem_swap_timelevel ***) '           !bK

      RETURN
   END SUBROUTINE kchem_swap_timelevel

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
   SUBROUTINE kchem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

      IMPLICIT NONE

      CHARACTER (LEN=*), INTENT(IN)  ::  var         !<
      LOGICAL, INTENT(OUT)           ::  found       !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

      found  = .TRUE.

      if(var(1:3) == 'kc_')   then                    ! always the same grid for chemistry variables
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'                             !kk Use same z axis as u variables. Has to be checked if OK
      else
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
      end if

      write(6,*) 'kchem_define_netcdf_grid ',TRIM(var),' ',trim(grid_x),' ',found

   END SUBROUTINE kchem_define_netcdf_grid

   SUBROUTINE kchem_check_data_output( var, unit, i, ilen, k )


      USE control_parameters,                                                 &
         ONLY: data_output, message_string

      IMPLICIT NONE

      CHARACTER (LEN=*) ::  unit     !<
      CHARACTER (LEN=*) ::  var      !<

      INTEGER(iwp) :: i
      INTEGER(iwp) :: ilen
      INTEGER(iwp) :: k

      INTEGER              :: n
      CHARACTER(len=16)    ::  spec_name

      unit = 'illegal'

      spec_name = TRIM(var(4:))                     !bK 1:3 is  'kc_', fm (4:), the array contains NO, NO2, RCHO, O3 fm namelist.

      if(myid == 0) print*,'fm #117@624 kc_check_data_output, var  = ', spec_name 

      do n=1,NSPEC
         if(TRIM(spec_name) == TRIM(chem_species(n)%name))   Then
            unit = 'ppm'
            if(myid == 0) print*,'fm #117@629 inside loop prnt spec_name  ',spec_name
         end if
     end do

      print*,'fm #117@633  kchem_check_data_output '       !bK debug
    

      RETURN
   END SUBROUTINE kchem_check_data_output

   SUBROUTINE kchem_check_data_output_pr( variable, var_count, unit, dopr_unit )


    USE arrays_3d,                                                             &
        ONLY: zu

    USE control_parameters,                                                    &
        ONLY: data_output_pr, message_string, chemistry

    USE indices

    USE profil_parameter

    USE statistics


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< 
    CHARACTER (LEN=*) ::  variable !< 
    CHARACTER (LEN=*) ::  dopr_unit
    CHARACTER(len=16) ::  spec_name
 
    INTEGER(iwp) ::  var_count, n  !<
    

    spec_name = TRIM(variable(4:))    
!             write(9,*) 'fm #32 .. chemistry ',chemistry

             IF (  .NOT.  chemistry )  THEN
                 message_string = 'data_output_pr = ' //                        &
                 TRIM( data_output_pr(var_count) ) // ' is not imp' // &
                          'lemented for chemistry = .FALSE.'
                 CALL message( 'check_parameters', 'PA0185', 1, 2, 0, 6, 0 )
            ELSE
                do n = 1, NSPEC
                    IF (TRIM( spec_name ) == TRIM( chem_species(n)%name ) ) THEN 
                        dopr_index(var_count) = 900 
                        dopr_unit  = 'ppm'
                        hom(:,2,900,:) = SPREAD( zu, 2, statistic_regions+1 )
                    ENDIF
                enddo
            ENDIF
!            if(myid ==0) write(9,*) 'fm #32 .. variable is ',variable,' ..and var_count is .. ', var_count

   END SUBROUTINE kchem_check_data_output_pr


   SUBROUTINE kchem_data_output_3d( av, variable, found, local_pf )


      USE indices

      USE kinds


      IMPLICIT NONE

      CHARACTER (LEN=*) ::  variable !<
      LOGICAL      ::  found !<
      INTEGER(iwp) ::  av    !<
      REAL(sp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !<

      !-- local variables

      INTEGER              ::  i, j, k,n
      CHARACTER(len=16)    ::  spec_name


      found = .FALSE.

      spec_name = TRIM(variable(4:))


! IF av

!av == 0

      do n=1,NSPEC
         if(TRIM(spec_name) == TRIM(chem_species(n)%name))   Then
            write(6,*) 'Output of species ',TRIM(variable),' ',TRIM(chem_species(n)%name)       !bK print output on screen
            
             IF (av == 0) THEN

             DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = chem_species(n)%conc(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
           
            ELSE

             DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = chem_species_av(n)%conc_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO


            END IF

            found = .TRUE.
         end if
      end do

! av == 1 


! local_pf = chem_species_av


!    if(myid == 0) print*, 'from kchem_data_output_3d #6'       !bK debug
!    flush(9)

      RETURN
   END SUBROUTINE kchem_data_output_3d

   SUBROUTINE kchem_3d_data_averaging ( mode, variable )

   USE control_parameters

   USE indices

   USE kinds
 
   IMPLICIT NONE
 
   CHARACTER (LEN=*) ::  mode    !< 
   CHARACTER (LEN=*) :: variable !< 
 
   INTEGER(iwp) ::  i !< 
   INTEGER(iwp) ::  j !< 
   INTEGER(iwp) ::  k !< 
   INTEGER(iwp) :: nn, kc_idx !<
  
   IF ( mode == 'allocate' )  THEN
 
        DO nn = 1, NSPEC
           IF (TRIM(variable(4:)) == TRIM(chem_species(nn)%name)) THEN
                kc_idx = nn 
!                ALLOCAqE(spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
                  chem_species_av(nn)%conc_av = 0.0_wp
                if(myid == 0) write(9,*) 'the kc_dix fm allocate is ', kc_idx,' ..and variable is .. ', variable

                flush(9)
            ENDIF
       ENDDO

   ELSEIF ( mode == 'sum' )  THEN
   
        DO nn = 1, NSPEC
           IF (TRIM(variable(4:)) == TRIM(chem_species(nn)%name)) THEN
                kc_idx = nn 
!                ALLOCATE(spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
        !         chem_species_av(nn)%conc_av = 0.0_wp
                if(myid == 0) write(9,*) 'the kc_dix fm sum is ', kc_idx,' ..and variable is .. ', variable
                flush(9)

              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                   chem_species_av(kc_idx)%conc_av(k,j,i) = chem_species_av(kc_idx)%conc_av(k,j,i) + chem_species(kc_idx)%conc(k,j,i)
                  ENDDO
               ENDDO
              ENDDO
              ! EXIT
            ENDIF
       ENDDO
 

   ELSEIF ( mode == 'average' )  THEN
        DO nn = 1, NSPEC
           IF (TRIM(variable(4:)) == TRIM(chem_species(nn)%name)) THEN
                kc_idx = nn 
!                ALLOCATE(spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
        !         chem_species_av(nn)%conc_av = 0.0_wp
                if(myid == 0) write(9,*) 'the kc_dix fm average is ', kc_idx,' ..and variable is .. ', variable

                flush(9)
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb, nzt+1
                       chem_species_av(kc_idx)%conc_av(k,j,i) = chem_species_av(kc_idx)%conc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                    ENDDO
                 ENDDO
              ENDDO
              ! EXIT
            ENDIF
       ENDDO
       
   ENDIF     


   END SUBROUTINE kchem_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Write restart data for land surface model
!------------------------------------------------------------------------------!
   SUBROUTINE kchem_last_actions
 

   USE control_parameters
    
   USE kinds

   IMPLICIT NONE 

   INTEGER(iwp) :: kc_n !< 
   CHARACTER(LEN=20)                                         :: chems_name
   REAL(kind=wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   :: chems_conc



   print*,'fm #110@862 write_binary =  ', write_binary
  
   IF ( write_binary(1:4) == 'true' )  THEN
   DO kc_n = 1, NSPEC
!
      chems_name = chem_species(kc_n)%name
!
      WRITE(14) chems_name; WRITE(14) chem_species(kc_n)%conc

      if(myid==0) print*,' fm #110@889 val of chem_species(n)%name is ', chems_name

   END DO
   
   WRITE ( 14 )  '*** end kchem ***   '

!       IF ( ALLOCATED( qsws_soil_eb_av ) )  THEN 
!          WRITE ( 14 )  'qsws_soil_eb_av     ';  WRITE ( 14 )  qsws_soil_eb_av
!       ENDIF 

   END IF

   END SUBROUTINE kchem_last_actions

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> read restart data for chemical species .. kchem_driver
!------------------------------------------------------------------------------!


   SUBROUTINE kchem_read_restart_data( i, nxlfa, nxl_on_file, nxrfa, nxr_on_file,   &
                                     nynfa, nyn_on_file, nysfa, nys_on_file,        &
                                     offset_xa, offset_ya, overlap_count, tmp_3d)    
                                     
 
    USE control_parameters
             
    USE indices
    
    USE kinds
    
    USE pegrid

    IMPLICIT NONE 

    CHARACTER (LEN=20) :: field_char   !<   

    INTEGER(iwp) ::  i, kc_n         !< 
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
!!
!!
    REAL(wp),                                                                  &
       DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp, nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_3d 


   IF ( initializing_actions == 'read_restart_data' )  THEN
      READ ( 13 )  field_char

      DO  WHILE ( TRIM( field_char ) /= '*** end kchem ***')

         DO  k = 1, overlap_count

            nxlf = nxlfa(i,k)
            nxlc = nxlfa(i,k) + offset_xa(i,k)
            nxrf = nxrfa(i,k)
            nxrc = nxrfa(i,k) + offset_xa(i,k)
            nysf = nysfa(i,k)
            nysc = nysfa(i,k) + offset_ya(i,k)
            nynf = nynfa(i,k)
            nync = nynfa(i,k) + offset_ya(i,k)



             DO kc_n = 1, NSPEC
!                if(myid==0) print*,'fm 110@974 kchem_restart_data ',field_char
                IF (TRIM( field_char ) == TRIM(chem_species(kc_n)%name) )  THEN

                   IF ( k == 1 )  READ ( 13 )  tmp_3d

                   chem_species(kc_n)%conc = tmp_3d
                if(myid==0) print*,'fm 110@982 kchem_restart_data ',field_char

                END IF

              END DO

         ENDDO

         READ ( 13 )  field_char
!
      ENDDO
   ENDIF
   if(myid==0) print*,'fm 110@994  before going back #41'
!
   END SUBROUTINE kchem_read_restart_data
!

END MODULE kchem_driver
