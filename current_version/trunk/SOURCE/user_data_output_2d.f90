!> @file user_data_output_2d.f90
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
! $Id: user_data_output_2d.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1551 2015-03-03 14:18:16Z maronga
! Replaced nzb and nzt+1 with the new array bounds nzb_do and nzt_do.
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
!> Resorts the user-defined output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) and sets the grid on which it is defined.
!> Allowed values for grid are "zu" and "zw".
!------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )
 

    USE indices

    USE kinds

    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av !< 
    INTEGER(iwp) ::  i  !< 
    INTEGER(iwp) ::  j  !< 
    INTEGER(iwp) ::  k  !< 
    INTEGER(iwp) ::  nzb_do !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av)
!--    have to be declared and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )
!          IF ( av == 0 )  THEN
!             DO  i = nxlg, nxrg
!                DO  j = nysg, nyng
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ELSE
!             DO  i = nxlg, nxrg
!                DO  j = nysg, nyng
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2_av(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ENDIF
!
!          grid = 'zu'


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE user_data_output_2d

