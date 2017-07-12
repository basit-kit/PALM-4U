!> @file user_init_land_surface.f90
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
! $Id: user_init_land_surface.f90 2001 2016-08-20 18:41:22Z knoop $

! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1972 2016-07-26 07:52:02Z maronga
! Update of use statements
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Changed description text
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
! Description:
! ------------
!> Execution of user-defined actions to initiate the land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_land_surface
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE land_surface_model_mod

    USE netcdf_interface,                                                      &
        ONLY: dots_label, dots_unit, dots_num
    
    USE pegrid
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) :: i   !< running index
    INTEGER(iwp) :: j   !< running index

!
!-- Here the user-defined land surface initialization actions follow:



 END SUBROUTINE user_init_land_surface

