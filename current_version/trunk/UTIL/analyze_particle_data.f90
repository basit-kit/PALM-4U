 PROGRAM analyze_particle_data
!--------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation,
! either version 3 of the License, or (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2014  Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
! Former revisions:
! -----------------
! $Id: analyze_particle_data.f90 1310 2014-03-14 08:01:56Z raasch $
!
! 1046 2012-11-09 14:38:45Z maronga
! code put under GPL (PALM 3.9)
!
! Description:
! ------------
! This routine reads the particle data files generated by PALM
! and does some statistical analysis on these data.
!-------------------------------------------------------------------------------!

    IMPLICIT NONE

!
!-- Variable definitions
    CHARACTER (LEN=5)  ::  id_char
    CHARACTER (LEN=80) ::  run_description_header

    INTEGER, PARAMETER ::  spk = SELECTED_REAL_KIND( 6 )

    INTEGER ::  class, danz = 0, i, id, maximum_number_of_particles, n,   &
                number_of_intervalls, number_of_particles,                &
                number_of_vertical_levels, timelevel = 1, vertical_level, &
                total_number_of_particles

    INTEGER, DIMENSION(:), ALLOCATABLE ::  class_table
    
    LOGICAL ::  found

    REAL    ::  class_width, hdistance, km, km_x, km_y, particle_age, sigma, &
                sigma_local, sigma_local_x, sigma_local_y, sigma_x, sigma_y, &
                simulated_time, vertical_resolution
    REAL, DIMENSION(:,:), ALLOCATABLE ::  diffusivities

    TYPE particle_type
       SEQUENCE
       INTEGER ::  color, tailpoints
       REAL    ::  age, origin_x, origin_y, origin_z, size, speed_x, speed_y, &
                   speed_z, x, y, z
    END TYPE particle_type

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  particles


!
!-- Check, if file from PE0 exists and terminate program if it doesn't.
    WRITE (id_char,'(''_'',I4.4)')  danz
    INQUIRE ( FILE=id_char, EXIST=found )
!
!-- Find out the number of files (equal to the number of PEs which
!-- have been used in PALM) and open them
    DO  WHILE ( found )

       OPEN ( danz+110, FILE=id_char, FORM='UNFORMATTED' )
       danz = danz + 1
       WRITE (id_char,'(''_'',I4.4)')  danz
       INQUIRE ( FILE=id_char, EXIST=found )

    ENDDO
!
!-- Info-output
    PRINT*, ''
    PRINT*, '*** analyze_particle_data ***'
    IF ( danz /= 0 )  THEN
       PRINT*, '*** ', danz, ' file(s) found'
    ELSE
       PRINT*, '+++ file _0000 not found'
       PRINT*, '    program terminated'
       STOP
    ENDIF

!
!-- Loop over all timelevels of output 
    DO
       total_number_of_particles = 0
       sigma   = 0.0
       sigma_x = 0.0
       sigma_y = 0.0
!
!--    Loop over all files (reading data of the subdomains)
       DO  id = 0, danz-1
!
!--       Read file header
          IF ( timelevel == 1 )  THEN
             READ ( id+110 )  run_description_header
!
!--          Print header information
             IF ( id == 0 )  THEN
                PRINT*, '*** run: ', run_description_header
                PRINT*, ' '
                PRINT*, '--> enter class width in m:'
                READ*, class_width
                PRINT*, '--> enter number of class intervalls:'
                READ*, number_of_intervalls
                PRINT*, '--> enter vertical resolution in m:'
                READ*, vertical_resolution
                PRINT*, '--> enter number of vertical levels:'
                READ*, number_of_vertical_levels
!
!--             Allocate table space
                ALLOCATE( class_table( 0:number_of_intervalls ) )
                class_table = 0
                ALLOCATE( diffusivities(0:number_of_vertical_levels,5) )
                diffusivities = 0.0
             ENDIF
          ENDIF
!
!--       Read time information and indices
          READ ( id+110, END=10 )  simulated_time, maximum_number_of_particles,&
                                   number_of_particles
!
!--       Print timelevel and number of particles
          IF ( id == 0 )  THEN
             PRINT*, ' '
             PRINT*, '*** time: ', simulated_time
          ENDIF
          PRINT*, 'PE', id, ': ', number_of_particles, ' particles'
!
!--       Allocate array and read particle data
          ALLOCATE( particles(maximum_number_of_particles) )
          READ ( id+110 )  particles
!
!--       Analyze the particle data
          DO  n = 1, number_of_particles
!
!--          Calculate horizontal distance from of particle from its origin
             hdistance = SQRT( ( particles(n)%x - particles(n)%origin_x )**2 + &
                               ( particles(n)%y - particles(n)%origin_y )**2 )
             class = hdistance / class_width
             sigma_local   = hdistance**2
             sigma_local_x = ( particles(n)%x - particles(n)%origin_x )**2
             sigma_local_y = ( particles(n)%y - particles(n)%origin_y )**2

             vertical_level = particles(n)%origin_z / vertical_resolution
             IF ( vertical_level > number_of_vertical_levels )  THEN
                vertical_level = number_of_vertical_levels
             ENDIF

             IF ( class > number_of_intervalls )  THEN
                class = number_of_intervalls
!                PRINT*, 'x =',particles(n)%x,' y =',particles(n)%y
!                PRINT*, 'xo=',particles(n)%origin_x,' yo=',particles(n)%origin_y
             ENDIF

             class_table(class) = class_table(class) + 1

             diffusivities(vertical_level,1) = diffusivities(vertical_level,1) +&
                                               sigma_local
             diffusivities(vertical_level,2) = diffusivities(vertical_level,2) +&
                                               sigma_local_x
             diffusivities(vertical_level,3) = diffusivities(vertical_level,3) +&
                                               sigma_local_y
             diffusivities(vertical_level,4) = diffusivities(vertical_level,4) +&
                                               1.0

             vertical_level = particles(n)%z / vertical_resolution
             IF ( vertical_level > number_of_vertical_levels )  THEN
                vertical_level = number_of_vertical_levels
             ENDIF
             diffusivities(vertical_level,5) = diffusivities(vertical_level,5) +&
                                               1.0

!
!--          Summation for variances
             sigma = sigma + sigma_local
             sigma_x = sigma_x + sigma_local_x
             sigma_y = sigma_y + sigma_local_y
             total_number_of_particles = total_number_of_particles + 1

          ENDDO
!
!--       Store the particle age (it is provided that all particles have the
!--       same age)
          particle_age = particles(1)%age

!
!--       Deallocate particle array before data from next file are read
          DEALLOCATE( particles )

       ENDDO  ! next file
!
!--    Print statistics
       PRINT*, ' '
       PRINT*, '*** statistics for t = ', simulated_time
       DO  n = 0, number_of_intervalls-1
          WRITE ( *, 1 ) n*class_width, (n+1)*class_width, class_table(n)
   1      FORMAT (F6.1,' - ',F6.1, ' m   n = ',I7)
       ENDDO
       WRITE ( *, 2 )  (number_of_intervalls+1)*class_width, &
                       class_table(number_of_intervalls)
   2   FORMAT (6X,' > ',F6.1,' m   n = ',I7)

       sigma   = SQRT( sigma   / REAL( total_number_of_particles ) )
       km      = sigma**2 / ( 2.0 * particle_age )
       sigma_x = SQRT( sigma_x / REAL( total_number_of_particles ) )
       km_x    = sigma_x**2 / ( 2.0 * particle_age )
       sigma_y = SQRT( sigma_y / REAL( total_number_of_particles ) )
       km_y    = sigma_y**2 / ( 2.0 * particle_age )
       PRINT*, ' '
       WRITE ( *, 3 )  sigma, km, sigma_x, km_x, sigma_y, km_y
   3   FORMAT ('sigma   = ',F6.1,' m   Km   = ',F5.1,' m**2/s'/ &
               'sigma_x = ',F6.1,' m   Km_x = ',F5.1,' m**2/s'/ &
               'sigma_y = ',F6.1,' m   Km_y = ',F5.1,' m**2/s')

       PRINT*, ' '
       PRINT*, 'Height dependence of diffusivities:'
       DO  i = 0, number_of_vertical_levels-1
          IF ( diffusivities(i,4) == 0.0 )  diffusivities(i,4) = 1.0E-20
          WRITE ( *, 4 )  i*vertical_resolution, (i+1.0)*vertical_resolution,&
                          ( diffusivities(i,1) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          ( diffusivities(i,2) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          ( diffusivities(i,3) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          diffusivities(i,4), diffusivities(i,5)
   4      FORMAT (F6.1,'-',F6.1,' m  Km=',F5.1,' Km_x=',F5.1, &
                  ' Km_y=',F5.1,' n_o=',F7.0,' n=',F7.0)
       ENDDO
       IF ( diffusivities(number_of_vertical_levels,4) == 0.0 )  THEN
          diffusivities(number_of_vertical_levels,4) = 1.0E-20
       ENDIF
          i = number_of_vertical_levels
          WRITE ( *, 5 )  i*vertical_resolution,                        &
                          ( diffusivities(i,1) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          ( diffusivities(i,2) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          ( diffusivities(i,3) / diffusivities(i,4) ) / &
                                     ( 2.0 * particle_age ),            &
                          diffusivities(i,4), diffusivities(i,5)
   5   FORMAT (F6.1,'-...... m  Km=',F5.1,' Km_x=',F5.1, &
                  ' Km_y=',F5.1,' n_o=',F7.0,' n=',F7.0)

!
!--    Initialize class table for next timelevel
       class_table   =   0
       diffusivities = 0.0
       timelevel = timelevel + 1
       
    ENDDO  ! next timelevel

10  PRINT*, '*** EOF reached on file PARTICLE_DATA/_0000'

 END PROGRAM analyze_particle_data
