!> @file user_init.f90
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
! $Id: user_init.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1799 2016-04-05 08:35:55Z gronemeier
! Bugfix: added dots_num to variable list
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1320 2014-03-20 08:40:49Z raasch
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Execution of user-defined initializing actions
!------------------------------------------------------------------------------!
 SUBROUTINE user_init
 

    USE arrays_3d
    
    USE control_parameters
    
    USE grid_variables
    
    USE indices
    
    USE kinds
    
    USE netcdf_interface,                                                      &
        ONLY: dots_label, dots_unit, dots_num
    
    USE pegrid
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=20) ::  field_char   !< 
    CHARACTER(LEN=7)   ::  i_char = ''
    INTEGER(iwp), DIMENSION(1:2) ::  eso_ind      !< index values of emission_stripe_orig
    INTEGER(iwp)                 ::  esw_ind      !< index value of emission_stripe_width
    INTEGER(iwp)                 ::  i, j         !< running index

    
    
!
!-- Calculate indices for locations of additional timeseries output, used in
!-- user_statistics
    size_of_pos_array = SIZE(s_ts_pos_x)

    s_ts_pos_x_ind = NINT( s_ts_pos_x * ddx )
    s_ts_pos_y_ind = NINT( s_ts_pos_y * ddy )
    
!
!-- Additional timeseries output
    DO  i = 1, size_of_pos_array
       WRITE (i_char,'(''_'',I3.3)')  i
       dots_label(dots_num+i) = 's_pos' // TRIM(i_char)
       dots_unit(dots_num+i)  = 'kg'
    ENDDO

    dots_num_palm = dots_num
    dots_num = dots_num + size_of_pos_array
    
    
!
!-- Passive scalar shall only be released from a line source of width 
!-- emission_stripe_width, stripe origin is precribed by parameter 
!-- emission_stripe_orig under &userpar
    eso_ind = NINT( emission_stripe_orig_y * ddy )
    esw_ind = NINT( emission_stripe_width_y * ddy )
    WRITE(9,*) 'd', esw_ind
    FLUSH(9)

    DO  j = nysg, nyng
       IF ( ( eso_ind(1) <= j  .AND.  j < (eso_ind(1)+esw_ind) )  .OR.  ( eso_ind(2) <= j  .AND.  j < (eso_ind(2)+esw_ind) ) )  THEN
          ssws(j,:) = surface_scalarflux
       ELSE
          ssws(j,:) = 0.0_wp
       ENDIF
    ENDDO   

 END SUBROUTINE user_init

