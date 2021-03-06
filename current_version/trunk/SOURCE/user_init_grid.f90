!> @file user_init_grid.f90
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
! $Id: user_init_grid.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1968 2016-07-18 12:01:49Z suehring
! Change dimensions for nzb_local, which do not longer need to be set on 
! multigrid ghost points. 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 217 2008-12-09 18:00:48Z letzel
! +topography_grid_convention
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Execution of user-defined grid initializing actions
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_grid( nzb_local )
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  nzb_local   !< 

!
!-- Here the user-defined grid initializing actions follow:

!
!-- Set the index array nzb_local for non-flat topography.
!-- Here consistency checks concerning domain size and periodicity are necessary
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat', 'single_building', 'single_street_canyon' )
!
!--       Not allowed here since these are the standard cases used in init_grid.

       CASE ( 'user_defined_topography_1' )
!
!--       Here the user can define his own topography.
!--       After definition, please remove the following three lines!
          message_string = 'topography "' // topography // '" not available yet'
          CALL message( 'user_init_grid', 'UI0005', 1, 2, 0, 6, 0 )

       CASE DEFAULT
!
!--       The DEFAULT case is reached if the parameter topography contains a
!--       wrong character string that is neither recognized in init_grid nor 
!--       here in user_init_grid.
          message_string = 'unknown topography "' // topography // '"'
          CALL message( 'user_init_grid', 'UI0006', 1, 2, 0, 6, 0 )

    END SELECT


 END SUBROUTINE user_init_grid

