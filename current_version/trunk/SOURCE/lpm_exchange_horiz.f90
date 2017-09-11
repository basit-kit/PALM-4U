!> @file lpm_exchange_horiz.f90
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
! $Id: lpm_exchange_horiz.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1936 2016-06-13 13:37:44Z suehring
! Deallocation of unused memory
! 
! 1929 2016-06-09 16:25:25Z suehring
! Bugfixes: 
! - reallocation of new particles 
!   ( did not work for small number of min_nr_particle )
! - dynamical reallocation of north-south exchange arrays ( particles got lost )
! - north-south exchange ( nr_move_north, nr_move_south were overwritten by zero )
! - horizontal particle boundary conditions in serial mode
!
! Remove unused variables
! Descriptions in variable declaration blocks added
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. Unused variables removed.
!
! 1783 2016-03-06 18:36:17Z raasch
! new netcdf-module included
!
! 1691 2015-10-26 16:17:44Z maronga
! Formatting corrections.
!
! 1685 2015-10-08 07:32:13Z raasch
! bugfix concerning vertical index offset in case of ocean
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
! 
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf output queries
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 851 2012-03-15 14:32:58Z raasch
! Bugfix: resetting of particle_mask and tail mask moved from end of this
! routine to lpm
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
! Exchange of particles between the subdomains.
!------------------------------------------------------------------------------!
 MODULE lpm_exchange_horiz_mod
 

    USE control_parameters,                                                    &
        ONLY:  dz, message_string, simulated_time

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy, dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nyn, nys, nzb, nzt

    USE kinds

    USE lpm_pack_arrays_mod,                                                   &
        ONLY:  lpm_pack_arrays

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format

    USE particle_attributes,                                                   &
        ONLY:  alloc_factor, deleted_particles, grid_particles,                &
               ibc_par_lr, ibc_par_ns, min_nr_particle,                        &
               mpi_particle_type, number_of_particles,                         &
               offset_ocean_nzt, offset_ocean_nzt_m1, particles,               &
               particle_type, prt_count, trlp_count_sum,                       &
               trlp_count_recv_sum, trnp_count_sum, trnp_count_recv_sum,       &
               trrp_count_sum, trrp_count_recv_sum, trsp_count_sum,            &
               trsp_count_recv_sum, zero_particle

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp), PARAMETER ::  NR_2_direction_move = 10000 !<
    INTEGER(iwp)            ::  nr_move_north               !<
    INTEGER(iwp)            ::  nr_move_south               !<

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  move_also_north
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  move_also_south

    SAVE

    PRIVATE
    PUBLIC lpm_exchange_horiz, lpm_move_particle, realloc_particles_array,     &
           dealloc_particles_array

    INTERFACE lpm_exchange_horiz
       MODULE PROCEDURE lpm_exchange_horiz
    END INTERFACE lpm_exchange_horiz

    INTERFACE lpm_move_particle
       MODULE PROCEDURE lpm_move_particle
    END INTERFACE lpm_move_particle

    INTERFACE realloc_particles_array
       MODULE PROCEDURE realloc_particles_array
    END INTERFACE realloc_particles_array

    INTERFACE dealloc_particles_array
       MODULE PROCEDURE dealloc_particles_array
    END INTERFACE dealloc_particles_array
CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange between subdomains.
!> As soon as one particle has moved beyond the boundary of the domain, it
!> is included in the relevant transfer arrays and marked for subsequent
!> deletion on this PE.
!> First sweep for crossings in x direction. Find out first the number of
!> particles to be transferred and allocate temporary arrays needed to store
!> them.
!> For a one-dimensional decomposition along y, no transfer is necessary,
!> because the particle remains on the PE, but the particle coordinate has to
!> be adjusted.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_exchange_horiz

    IMPLICIT NONE

    INTEGER(iwp) ::  i                 !< grid index (x) of particle positition
    INTEGER(iwp) ::  ip                !< index variable along x
    INTEGER(iwp) ::  j                 !< grid index (y) of particle positition
    INTEGER(iwp) ::  jp                !< index variable along y
    INTEGER(iwp) ::  kp                !< index variable along z
    INTEGER(iwp) ::  n                 !< particle index variable 
    INTEGER(iwp) ::  trlp_count        !< number of particles send to left PE
    INTEGER(iwp) ::  trlp_count_recv   !< number of particles receive from right PE
    INTEGER(iwp) ::  trnp_count        !< number of particles send to north PE
    INTEGER(iwp) ::  trnp_count_recv   !< number of particles receive from south PE
    INTEGER(iwp) ::  trrp_count        !< number of particles send to right PE
    INTEGER(iwp) ::  trrp_count_recv   !< number of particles receive from left PE
    INTEGER(iwp) ::  trsp_count        !< number of particles send to south PE
    INTEGER(iwp) ::  trsp_count_recv   !< number of particles receive from north PE

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvlp  !< particles received from right PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvnp  !< particles received from south PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvrp  !< particles received from left PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvsp  !< particles received from north PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trlp  !< particles send to left PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trnp  !< particles send to north PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trrp  !< particles send to right PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trsp  !< particles send to south PE

    CALL cpu_log( log_point_s(23), 'lpm_exchange_horiz', 'start' )

#if defined( __parallel )

!
!-- Exchange between subdomains.
!-- As soon as one particle has moved beyond the boundary of the domain, it
!-- is included in the relevant transfer arrays and marked for subsequent
!-- deletion on this PE.
!-- First sweep for crossings in x direction. Find out first the number of
!-- particles to be transferred and allocate temporary arrays needed to store
!-- them.
!-- For a one-dimensional decomposition along y, no transfer is necessary,
!-- because the particle remains on the PE, but the particle coordinate has to
!-- be adjusted.
    trlp_count  = 0
    trrp_count  = 0

    trlp_count_recv   = 0
    trrp_count_recv   = 0

    IF ( pdims(1) /= 1 )  THEN
!
!--    First calculate the storage necessary for sending and receiving the data.
!--    Compute only first (nxl) and last (nxr) loop iterration.
       DO  ip = nxl, nxr, nxr - nxl
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt

                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask )  THEN
                      i = ( particles(n)%x + 0.5_wp * dx ) * ddx
!
!--                   Above calculation does not work for indices less than zero
                      IF ( particles(n)%x < -0.5_wp * dx )  i = -1

                      IF ( i < nxl )  THEN
                         trlp_count = trlp_count + 1
                      ELSEIF ( i > nxr )  THEN
                         trrp_count = trrp_count + 1
                      ENDIF
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO

       IF ( trlp_count  == 0 )  trlp_count  = 1
       IF ( trrp_count  == 0 )  trrp_count  = 1

       ALLOCATE( trlp(trlp_count), trrp(trrp_count) )

       trlp = zero_particle
       trrp = zero_particle

       trlp_count  = 0
       trrp_count  = 0

    ENDIF
!
!-- Compute only first (nxl) and last (nxr) loop iterration
    DO  ip = nxl, nxr, nxr-nxl
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 ) CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Only those particles that have not been marked as 'deleted' may
!--             be moved.
                IF ( particles(n)%particle_mask )  THEN

                   i = ( particles(n)%x + 0.5_wp * dx ) * ddx
!
!--                Above calculation does not work for indices less than zero
                   IF ( particles(n)%x < - 0.5_wp * dx )  i = -1

                   IF ( i <  nxl )  THEN
                      IF ( i < 0 )  THEN
!
!--                   Apply boundary condition along x
                         IF ( ibc_par_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               particles(n)%x        = ( nx + 1 ) * dx + particles(n)%x
                               particles(n)%origin_x = ( nx + 1 ) * dx + &
                               particles(n)%origin_x
                            ELSE
                               trlp_count = trlp_count + 1
                               trlp(trlp_count)   = particles(n)
                               trlp(trlp_count)%x = ( nx + 1 ) * dx + trlp(trlp_count)%x
                               trlp(trlp_count)%origin_x = trlp(trlp_count)%origin_x + &
                               ( nx + 1 ) * dx
                               particles(n)%particle_mask  = .FALSE.
                               deleted_particles = deleted_particles + 1

                               IF ( trlp(trlp_count)%x >= (nx + 0.5_wp)* dx - 1.0E-12_wp )  THEN
                                  trlp(trlp_count)%x = trlp(trlp_count)%x - 1.0E-10_wp
                                  !++ why is 1 subtracted in next statement???
                                  trlp(trlp_count)%origin_x = trlp(trlp_count)%origin_x - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%x       = -particles(n)%x
                            particles(n)%speed_x = -particles(n)%speed_x

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will be 
!--                      send to the neighbouring PE
                         trlp_count = trlp_count + 1
                         trlp(trlp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ELSEIF ( i > nxr )  THEN
                      IF ( i > nx )  THEN
!
!--                      Apply boundary condition along x
                         IF ( ibc_par_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               particles(n)%x = particles(n)%x - ( nx + 1 ) * dx
                               particles(n)%origin_x = particles(n)%origin_x - &
                               ( nx + 1 ) * dx
                            ELSE
                               trrp_count = trrp_count + 1
                               trrp(trrp_count) = particles(n)
                               trrp(trrp_count)%x = trrp(trrp_count)%x - ( nx + 1 ) * dx
                               trrp(trrp_count)%origin_x = trrp(trrp_count)%origin_x - &
                               ( nx + 1 ) * dx
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles = deleted_particles + 1

                            ENDIF

                         ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%x       = 2 * ( nx * dx ) - particles(n)%x
                            particles(n)%speed_x = -particles(n)%speed_x

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will be send
!--                      to the neighbouring PE
                         trrp_count = trrp_count + 1
                         trrp(trrp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Allocate arrays required for north-south exchange, as these
!-- are used directly after particles are exchange along x-direction.
    ALLOCATE( move_also_north(1:NR_2_direction_move) )
    ALLOCATE( move_also_south(1:NR_2_direction_move) )

    nr_move_north = 0
    nr_move_south = 0
!
!-- Send left boundary, receive right boundary (but first exchange how many
!-- and check, if particle storage must be extended)
    IF ( pdims(1) /= 1 )  THEN

       CALL MPI_SENDRECV( trlp_count,      1, MPI_INTEGER, pleft,  0, &
                          trrp_count_recv, 1, MPI_INTEGER, pright, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvrp(MAX(1,trrp_count_recv)))

       CALL MPI_SENDRECV( trlp(1)%radius, max(1,trlp_count), mpi_particle_type,&
                          pleft, 1, rvrp(1)%radius,                            &
                          max(1,trrp_count_recv), mpi_particle_type, pright, 1,&
                          comm2d, status, ierr )

       IF ( trrp_count_recv > 0 )  CALL Add_particles_to_gridcell(rvrp(1:trrp_count_recv))

       DEALLOCATE(rvrp)

!
!--    Send right boundary, receive left boundary
       CALL MPI_SENDRECV( trrp_count,      1, MPI_INTEGER, pright, 0, &
                          trlp_count_recv, 1, MPI_INTEGER, pleft,  0, &
                          comm2d, status, ierr )

       ALLOCATE(rvlp(MAX(1,trlp_count_recv)))

       CALL MPI_SENDRECV( trrp(1)%radius, max(1,trrp_count), mpi_particle_type,&
                          pright, 1, rvlp(1)%radius,                           &
                          max(1,trlp_count_recv), mpi_particle_type, pleft, 1, &
                          comm2d, status, ierr )

       IF ( trlp_count_recv > 0 )  CALL Add_particles_to_gridcell(rvlp(1:trlp_count_recv))

       DEALLOCATE( rvlp )
       DEALLOCATE( trlp, trrp )

    ENDIF

!
!-- Check whether particles have crossed the boundaries in y direction. Note 
!-- that this case can also apply to particles that have just been received
!-- from the adjacent right or left PE.
!-- Find out first the number of particles to be transferred and allocate
!-- temporary arrays needed to store them.
!-- For a one-dimensional decomposition along y, no transfer is necessary,
!-- because the particle remains on the PE.
    trsp_count  = nr_move_south
    trnp_count  = nr_move_north

    trsp_count_recv   = 0
    trnp_count_recv   = 0

    IF ( pdims(2) /= 1 )  THEN
!
!--    First calculate the storage necessary for sending and receiving the
!--    data
       DO  ip = nxl, nxr
          DO  jp = nys, nyn, nyn-nys    !compute only first (nys) and last (nyn) loop iterration
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask )  THEN
                      j = ( particles(n)%y + 0.5_wp * dy ) * ddy
!
!--                   Above calculation does not work for indices less than zero
                      IF ( particles(n)%y < -0.5_wp * dy )  j = -1

                      IF ( j < nys )  THEN
                         trsp_count = trsp_count + 1
                      ELSEIF ( j > nyn )  THEN
                         trnp_count = trnp_count + 1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       IF ( trsp_count  == 0 )  trsp_count  = 1
       IF ( trnp_count  == 0 )  trnp_count  = 1

       ALLOCATE( trsp(trsp_count), trnp(trnp_count) )

       trsp = zero_particle
       trnp = zero_particle

       trsp_count  = nr_move_south
       trnp_count  = nr_move_north
       
       trsp(1:nr_move_south) = move_also_south(1:nr_move_south)
       trnp(1:nr_move_north) = move_also_north(1:nr_move_north)

    ENDIF

    DO  ip = nxl, nxr
       DO  jp = nys, nyn, nyn-nys ! compute only first (nys) and last (nyn) loop iterration
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Only those particles that have not been marked as 'deleted' may
!--             be moved.
                IF ( particles(n)%particle_mask )  THEN

                   j = ( particles(n)%y + 0.5_wp * dy ) * ddy
!
!--                Above calculation does not work for indices less than zero
                   IF ( particles(n)%y < -0.5_wp * dy )  j = -1

                   IF ( j < nys )  THEN
                      IF ( j < 0 )  THEN
!
!--                      Apply boundary condition along y
                         IF ( ibc_par_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               particles(n)%y = ( ny + 1 ) * dy + particles(n)%y
                               particles(n)%origin_y = ( ny + 1 ) * dy + &
                                                     particles(n)%origin_y
                            ELSE
                               trsp_count         = trsp_count + 1
                               trsp(trsp_count)   = particles(n)
                               trsp(trsp_count)%y = ( ny + 1 ) * dy + &
                                                 trsp(trsp_count)%y
                               trsp(trsp_count)%origin_y = trsp(trsp_count)%origin_y &
                                                + ( ny + 1 ) * dy
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles = deleted_particles + 1

                               IF ( trsp(trsp_count)%y >= (ny+0.5_wp)* dy - 1.0E-12_wp )  THEN
                                  trsp(trsp_count)%y = trsp(trsp_count)%y - 1.0E-10_wp
                                  !++ why is 1 subtracted in next statement???
                                  trsp(trsp_count)%origin_y =                        &
                                                  trsp(trsp_count)%origin_y - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles          = deleted_particles + 1

                         ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%y       = -particles(n)%y
                            particles(n)%speed_y = -particles(n)%speed_y

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trsp_count = trsp_count + 1
                         trsp(trsp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ELSEIF ( j > nyn )  THEN
                      IF ( j > ny )  THEN
!
!--                       Apply boundary condition along y
                         IF ( ibc_par_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               particles(n)%y        = particles(n)%y - ( ny + 1 ) * dy
                               particles(n)%origin_y =                         &
                                          particles(n)%origin_y - ( ny + 1 ) * dy
                            ELSE
                               trnp_count         = trnp_count + 1
                               trnp(trnp_count)   = particles(n)
                               trnp(trnp_count)%y =                            &
                                          trnp(trnp_count)%y - ( ny + 1 ) * dy
                               trnp(trnp_count)%origin_y =                     &
                                         trnp(trnp_count)%origin_y - ( ny + 1 ) * dy
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles          = deleted_particles + 1
                            ENDIF

                         ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%y       = 2 * ( ny * dy ) - particles(n)%y
                            particles(n)%speed_y = -particles(n)%speed_y

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trnp_count = trnp_count + 1
                         trnp(trnp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Send front boundary, receive back boundary (but first exchange how many
!-- and check, if particle storage must be extended)
    IF ( pdims(2) /= 1 )  THEN

       CALL MPI_SENDRECV( trsp_count,      1, MPI_INTEGER, psouth, 0, &
                          trnp_count_recv, 1, MPI_INTEGER, pnorth, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvnp(MAX(1,trnp_count_recv)))
 
       CALL MPI_SENDRECV( trsp(1)%radius, trsp_count, mpi_particle_type,      &
                          psouth, 1, rvnp(1)%radius,                             &
                          trnp_count_recv, mpi_particle_type, pnorth, 1,   &
                          comm2d, status, ierr )

       IF ( trnp_count_recv  > 0 )  CALL Add_particles_to_gridcell(rvnp(1:trnp_count_recv))

       DEALLOCATE(rvnp)

!
!--    Send back boundary, receive front boundary
       CALL MPI_SENDRECV( trnp_count,      1, MPI_INTEGER, pnorth, 0, &
                          trsp_count_recv, 1, MPI_INTEGER, psouth, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvsp(MAX(1,trsp_count_recv)))

       CALL MPI_SENDRECV( trnp(1)%radius, trnp_count, mpi_particle_type,      &
                          pnorth, 1, rvsp(1)%radius,                          &
                          trsp_count_recv, mpi_particle_type, psouth, 1,   &
                          comm2d, status, ierr )

       IF ( trsp_count_recv > 0 )  CALL Add_particles_to_gridcell(rvsp(1:trsp_count_recv))

       DEALLOCATE(rvsp)

       number_of_particles = number_of_particles + trsp_count_recv

       DEALLOCATE( trsp, trnp )

    ENDIF

    DEALLOCATE( move_also_north )
    DEALLOCATE( move_also_south )

#else

    DO  ip = nxl, nxr, nxr-nxl
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Apply boundary conditions

                IF ( particles(n)%x < -0.5_wp * dx )  THEN

                   IF ( ibc_par_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%x = ( nx + 1 ) * dx + particles(n)%x

                   ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%x       = -dx - particles(n)%x
                      particles(n)%speed_x = -particles(n)%speed_x
                   ENDIF

                ELSEIF ( particles(n)%x >= ( nx + 0.5_wp ) * dx )  THEN

                   IF ( ibc_par_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%x = particles(n)%x - ( nx + 1 ) * dx

                   ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%x       = ( nx + 1 ) * dx - particles(n)%x
                      particles(n)%speed_x = -particles(n)%speed_x
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO  ip = nxl, nxr
       DO  jp = nys, nyn, nyn-nys
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                IF ( particles(n)%y < -0.5_wp * dy )  THEN

                   IF ( ibc_par_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%y = ( ny + 1 ) * dy + particles(n)%y

                   ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%y       = -dy - particles(n)%y
                      particles(n)%speed_y = -particles(n)%speed_y
                   ENDIF

                ELSEIF ( particles(n)%y >= ( ny + 0.5_wp ) * dy )  THEN

                   IF ( ibc_par_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%y = particles(n)%y - ( ny + 1 ) * dy

                   ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%y       = ( ny + 1 ) * dy - particles(n)%y
                      particles(n)%speed_y = -particles(n)%speed_y
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO
#endif

!
!-- Accumulate the number of particles transferred between the subdomains
#if defined( __parallel )
    trlp_count_sum      = trlp_count_sum      + trlp_count
    trlp_count_recv_sum = trlp_count_recv_sum + trlp_count_recv
    trrp_count_sum      = trrp_count_sum      + trrp_count
    trrp_count_recv_sum = trrp_count_recv_sum + trrp_count_recv
    trsp_count_sum      = trsp_count_sum      + trsp_count
    trsp_count_recv_sum = trsp_count_recv_sum + trsp_count_recv
    trnp_count_sum      = trnp_count_sum      + trnp_count
    trnp_count_recv_sum = trnp_count_recv_sum + trnp_count_recv
#endif

    CALL cpu_log( log_point_s(23), 'lpm_exchange_horiz', 'stop' )

 END SUBROUTINE lpm_exchange_horiz

!------------------------------------------------------------------------------!
! Description:
! ------------
!> If a particle moves from one processor to another, this subroutine moves 
!> the corresponding elements from the particle arrays of the old grid cells 
!> to the particle arrays of the new grid cells.
!------------------------------------------------------------------------------!
 SUBROUTINE Add_particles_to_gridcell (particle_array)

    IMPLICIT NONE

    INTEGER(iwp)        ::  ip        !< grid index (x) of particle
    INTEGER(iwp)        ::  jp        !< grid index (x) of particle
    INTEGER(iwp)        ::  kp        !< grid index (x) of particle
    INTEGER(iwp)        ::  n         !< index variable of particle
    INTEGER(iwp)        ::  pindex    !< dummy argument for new number of particles per grid box

    LOGICAL             ::  pack_done !<

    TYPE(particle_type), DIMENSION(:), INTENT(IN)  ::  particle_array !< new particles in a grid box
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  temp_ns        !< temporary particle array for reallocation

    pack_done     = .FALSE.

    DO n = 1, SIZE(particle_array)

       IF ( .NOT. particle_array(n)%particle_mask )  CYCLE

       ip = ( particle_array(n)%x + 0.5_wp * dx ) * ddx
       jp = ( particle_array(n)%y + 0.5_wp * dy ) * ddy
       kp =   particle_array(n)%z / dz + 1 + offset_ocean_nzt

       IF ( ip >= nxl  .AND.  ip <= nxr  .AND.  jp >= nys  .AND.  jp <= nyn    &
            .AND.  kp >= nzb+1  .AND.  kp <= nzt)  THEN ! particle stays on processor
          number_of_particles = prt_count(kp,jp,ip)
          particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

          pindex = prt_count(kp,jp,ip)+1
          IF( pindex > SIZE(grid_particles(kp,jp,ip)%particles) )  THEN
             IF ( pack_done )  THEN
                CALL realloc_particles_array (ip,jp,kp)
             ELSE
                CALL lpm_pack_arrays
                prt_count(kp,jp,ip) = number_of_particles
                pindex = prt_count(kp,jp,ip)+1
                IF ( pindex > SIZE(grid_particles(kp,jp,ip)%particles) )  THEN
                   CALL realloc_particles_array (ip,jp,kp)
                ENDIF
                pack_done = .TRUE.
             ENDIF
          ENDIF
          grid_particles(kp,jp,ip)%particles(pindex) = particle_array(n)
          prt_count(kp,jp,ip) = pindex
       ELSE
          IF ( jp <= nys - 1 )  THEN
             nr_move_south = nr_move_south+1
!
!--          Before particle information is swapped to exchange-array, check 
!--          if enough memory is allocated. If required, reallocate exchange
!--          array.
             IF ( nr_move_south > SIZE(move_also_south) )  THEN
!
!--             At first, allocate further temporary array to swap particle 
!--             information.
                ALLOCATE( temp_ns(SIZE(move_also_south)+NR_2_direction_move) )
                temp_ns(1:nr_move_south-1) = move_also_south(1:nr_move_south-1)
                DEALLOCATE( move_also_south )
                ALLOCATE( move_also_south(SIZE(temp_ns)) )
                move_also_south(1:nr_move_south-1) = temp_ns(1:nr_move_south-1)
                DEALLOCATE( temp_ns )

             ENDIF

             move_also_south(nr_move_south) = particle_array(n)

             IF ( jp == -1 )  THEN
                move_also_south(nr_move_south)%y =                             &
                   move_also_south(nr_move_south)%y + ( ny + 1 ) * dy
                move_also_south(nr_move_south)%origin_y =                      &
                   move_also_south(nr_move_south)%origin_y + ( ny + 1 ) * dy
             ENDIF
          ELSEIF ( jp >= nyn+1 )  THEN
             nr_move_north = nr_move_north+1
!
!--          Before particle information is swapped to exchange-array, check 
!--          if enough memory is allocated. If required, reallocate exchange
!--          array.
             IF ( nr_move_north > SIZE(move_also_north) )  THEN
!
!--             At first, allocate further temporary array to swap particle 
!--             information.
                ALLOCATE( temp_ns(SIZE(move_also_north)+NR_2_direction_move) )
                temp_ns(1:nr_move_north-1) = move_also_south(1:nr_move_north-1)
                DEALLOCATE( move_also_north )
                ALLOCATE( move_also_north(SIZE(temp_ns)) )
                move_also_north(1:nr_move_north-1) = temp_ns(1:nr_move_north-1)
                DEALLOCATE( temp_ns )

             ENDIF

             move_also_north(nr_move_north) = particle_array(n)
             IF ( jp == ny+1 )  THEN
                move_also_north(nr_move_north)%y =                             &
                   move_also_north(nr_move_north)%y - ( ny + 1 ) * dy
                move_also_north(nr_move_north)%origin_y =                      &
                   move_also_north(nr_move_north)%origin_y - ( ny + 1 ) * dy
             ENDIF
          ELSE
             WRITE(0,'(a,8i7)') 'particle out of range ',myid,ip,jp,kp,nxl,nxr,nys,nyn
          ENDIF
       ENDIF
    ENDDO

    RETURN

 END SUBROUTINE Add_particles_to_gridcell




!------------------------------------------------------------------------------!
! Description:
! ------------
!> If a particle moves from one grid cell to another (on the current 
!> processor!), this subroutine moves the corresponding element from the
!> particle array of the old grid cell to the particle array of the new grid 
!> cell.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_move_particle

    IMPLICIT NONE

    INTEGER(iwp)        ::  i           !< grid index (x) of particle position
    INTEGER(iwp)        ::  ip          !< index variable along x
    INTEGER(iwp)        ::  j           !< grid index (y) of particle position
    INTEGER(iwp)        ::  jp          !< index variable along y
    INTEGER(iwp)        ::  k           !< grid index (z) of particle position
    INTEGER(iwp)        ::  kp          !< index variable along z
    INTEGER(iwp)        ::  n           !< index variable for particle array
    INTEGER(iwp)        ::  np_old_cell !< number of particles per grid box before moving
    INTEGER(iwp)        ::  n_start     !< start index
    INTEGER(iwp)        ::  pindex      !< dummy argument for number of new particle per grid box

    LOGICAL             ::  pack_done   !<

    TYPE(particle_type), DIMENSION(:), POINTER  ::  particles_old_cell !< particles before moving

    CALL cpu_log( log_point_s(41), 'lpm_move_particle', 'start' )

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             np_old_cell = prt_count(kp,jp,ip)
             IF ( np_old_cell <= 0 )  CYCLE
             particles_old_cell => grid_particles(kp,jp,ip)%particles(1:np_old_cell)
             n_start = -1
             
             DO  n = 1, np_old_cell
                i = ( particles_old_cell(n)%x + 0.5_wp * dx ) * ddx
                j = ( particles_old_cell(n)%y + 0.5_wp * dy ) * ddy
                k = particles_old_cell(n)%z / dz + 1 + offset_ocean_nzt
!
!--             Check, if particle has moved to another grid cell.
                IF ( i /= ip  .OR.  j /= jp  .OR.  k /= kp )  THEN
!
!--                The particle has moved to another grid cell. Now check, if 
!--                particle stays on the same processor.
                   IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.      &
                        j <= nyn  .AND.  k >= nzb+1  .AND.  k <= nzt)  THEN
!
!--                   If the particle stays on the same processor, the particle
!--                   will be added to the particle array of the new processor.
                      number_of_particles = prt_count(k,j,i)
                      particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                      pindex = prt_count(k,j,i)+1
                      IF (  pindex > SIZE(grid_particles(k,j,i)%particles)  )  &
                      THEN
                         n_start = n
                         EXIT
                      ENDIF

                      grid_particles(k,j,i)%particles(pindex) = particles_old_cell(n)
                      prt_count(k,j,i) = pindex

                      particles_old_cell(n)%particle_mask = .FALSE.
                   ENDIF
                ENDIF
             ENDDO

             IF ( n_start >= 0 )  THEN
                pack_done = .FALSE.
                DO  n = n_start, np_old_cell
                   i = ( particles_old_cell(n)%x + 0.5_wp * dx ) * ddx
                   j = ( particles_old_cell(n)%y + 0.5_wp * dy ) * ddy
                   k = particles_old_cell(n)%z / dz + 1 + offset_ocean_nzt
                   IF ( i /= ip  .OR.  j /= jp  .OR.  k /= kp )  THEN
!
!--                   Particle is in different box
                      IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.   &
                           j <= nyn  .AND.  k >= nzb+1  .AND.  k <= nzt)  THEN
!
!--                      Particle stays on processor
                         number_of_particles = prt_count(k,j,i)
                         particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                         pindex = prt_count(k,j,i)+1
                         IF ( pindex > SIZE(grid_particles(k,j,i)%particles) ) &
                         THEN
                            IF ( pack_done )  THEN
                               CALL realloc_particles_array(i,j,k)
                            ELSE
                               CALL lpm_pack_arrays
                               prt_count(k,j,i) = number_of_particles
!
!--                            If number of particles in the new grid box 
!--                            exceeds its allocated memory, the particle array
!--                            will be reallocated
                               IF ( pindex > SIZE(grid_particles(k,j,i)%particles) )  THEN
                                  CALL realloc_particles_array(i,j,k)
                               ENDIF

                               pack_done = .TRUE.
                            ENDIF
                         ENDIF

                         grid_particles(k,j,i)%particles(pindex) = particles_old_cell(n)
                         prt_count(k,j,i) = pindex

                         particles_old_cell(n)%particle_mask = .FALSE.
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(41), 'lpm_move_particle', 'stop' )

    RETURN

 END SUBROUTINE lpm_move_particle

!------------------------------------------------------------------------------!
! Description:
! ------------
!> If the allocated memory for the particle array do not suffice to add arriving
!> particles from neighbour grid cells, this subrouting reallocates the 
!> particle array to assure enough memory is available. 
!------------------------------------------------------------------------------!
 SUBROUTINE realloc_particles_array (i,j,k,size_in)

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in)                       ::  i              !<
    INTEGER(iwp), INTENT(in)                       ::  j              !<
    INTEGER(iwp), INTENT(in)                       ::  k              !<
    INTEGER(iwp), INTENT(in), OPTIONAL             ::  size_in        !<

    INTEGER(iwp)                                   :: old_size        !<
    INTEGER(iwp)                                   :: new_size        !<
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE :: tmp_particles_d !<
    TYPE(particle_type), DIMENSION(500)            :: tmp_particles_s !<

    old_size = SIZE(grid_particles(k,j,i)%particles)

    IF ( PRESENT(size_in) )   THEN
       new_size = size_in
    ELSE
       new_size = old_size * ( 1.0_wp + alloc_factor / 100.0_wp )
    ENDIF

    new_size = MAX( new_size, min_nr_particle, old_size + 1 )

    IF ( old_size <= 500 )  THEN

       tmp_particles_s(1:old_size) = grid_particles(k,j,i)%particles(1:old_size)

       DEALLOCATE(grid_particles(k,j,i)%particles)
       ALLOCATE(grid_particles(k,j,i)%particles(new_size))

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_s(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

    ELSE

       ALLOCATE(tmp_particles_d(new_size))
       tmp_particles_d(1:old_size) = grid_particles(k,j,i)%particles

       DEALLOCATE(grid_particles(k,j,i)%particles)
       ALLOCATE(grid_particles(k,j,i)%particles(new_size))

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_d(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

       DEALLOCATE(tmp_particles_d)

    ENDIF
    particles => grid_particles(k,j,i)%particles(1:number_of_particles)

    RETURN
 END SUBROUTINE realloc_particles_array





 SUBROUTINE dealloc_particles_array

    IMPLICIT NONE

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j
    INTEGER(iwp) ::  k
    INTEGER(iwp) :: old_size        !<
    INTEGER(iwp) :: new_size        !<

    LOGICAL                                        :: dealloc  

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE :: tmp_particles_d !<
    TYPE(particle_type), DIMENSION(500)            :: tmp_particles_s !<

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Determine number of active particles
             number_of_particles = prt_count(k,j,i)
!
!--          Determine allocated memory size
             old_size = SIZE( grid_particles(k,j,i)%particles )
!
!--          Check for large unused memory
             dealloc = ( ( number_of_particles < min_nr_particle .AND.         &
                           old_size            > min_nr_particle )  .OR.       &
                         ( number_of_particles > min_nr_particle .AND.         &
                           old_size - number_of_particles *                    &
                              ( 1.0_wp + 0.01_wp * alloc_factor ) > 0.0_wp ) )
                         

             IF ( dealloc )  THEN
                IF ( number_of_particles < min_nr_particle )  THEN
                   new_size = min_nr_particle
                ELSE 
                   new_size = INT( number_of_particles * ( 1.0_wp + 0.01_wp * alloc_factor ) )
                ENDIF

                IF ( number_of_particles <= 500 )  THEN

                   tmp_particles_s(1:number_of_particles) = grid_particles(k,j,i)%particles(1:number_of_particles)

                   DEALLOCATE(grid_particles(k,j,i)%particles)
                   ALLOCATE(grid_particles(k,j,i)%particles(new_size))

                   grid_particles(k,j,i)%particles(1:number_of_particles)          = tmp_particles_s(1:number_of_particles)
                   grid_particles(k,j,i)%particles(number_of_particles+1:new_size) = zero_particle

                ELSE

                   ALLOCATE(tmp_particles_d(number_of_particles))
                   tmp_particles_d(1:number_of_particles) = grid_particles(k,j,i)%particles(1:number_of_particles)

                   DEALLOCATE(grid_particles(k,j,i)%particles)
                   ALLOCATE(grid_particles(k,j,i)%particles(new_size))

                   grid_particles(k,j,i)%particles(1:number_of_particles)          = tmp_particles_d(1:number_of_particles)
                   grid_particles(k,j,i)%particles(number_of_particles+1:new_size) = zero_particle

                   DEALLOCATE(tmp_particles_d)

                ENDIF

             ENDIF
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE dealloc_particles_array


END MODULE lpm_exchange_horiz_mod
