!> @file user_check_data_output_pr.f90
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
! $Id: user_check_data_output_pr.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + relatec changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
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
!> Set the unit of user defined profile output quantities. For those variables
!> not recognized by the user, the parameter unit is set to "illegal", which
!> tells the calling routine that the output variable is not defined and leads
!> to a program abort.
!------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output_pr( variable, var_count, unit )
 

    USE arrays_3d

    USE indices

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  dopr_unit

    USE profil_parameter

    USE statistics

    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  user_pr_index !< 
    INTEGER(iwp) ::  var_count     !< 

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    Add additional CASE statements depending on the number of quantities
!--    for which profiles are to be calculated. The respective calculations
!--    to be performed have to be added in routine user_statistics.
!--    The quantities are (internally) identified by a user-profile-number
!--    (see variable "user_pr_index" below). The first user-profile must be assigned
!--    the number "pr_palm+1", the second one "pr_palm+2", etc. The respective
!--    user-profile-numbers have also to be used in routine user_statistics!
       CASE ( 'kc_O3w' )                      ! quantity string as given in                 !bK uncommented/added this block
!                                            ! data_output_pr_user
          user_pr_index = pr_palm + 1  
          dopr_index(var_count)  = user_pr_index    ! quantities' user-profile-number       
          dopr_unit(var_count)   = 'ppmm/s'  ! quantity unit
          hom(:,2,user_pr_index,:)       = SPREAD( zu, 2, statistic_regions+1 )
!                                            ! grid on which the quantity is
!                                            ! defined (use zu or zw)                       !bK added/uncommented this block
       CASE ( 'kc_O3' )                      ! quantity string as given in
!                                            ! data_output_pr_user
          user_pr_index = pr_palm + 2                                                       
          dopr_index(var_count)  = user_pr_index    ! quantities' user-profile-number
          dopr_unit(var_count)   = 'ppm'  ! quantity unit
          hom(:,2,user_pr_index,:)       = SPREAD( zu, 2, statistic_regions+1 )
!                                            ! grid on which the quantity is
!                                            ! defined (use zu or zw)

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE user_check_data_output_pr

